// Phobos primary impact Monte Carlo — production v3 angle diagnostics (2026-10): paired shielding ON/OFF
// Independent of step9_grav1.cpp and production v1; no overwrite of earlier programs.
//
// MODEL: stationary, isotropic, monoenergetic impactors at Mars-centric infinity.
// At a finite Mars-centric radius r, the Mars-only Liouville distribution has
//   v_loc(r)^2 = v_inf^2 + 2 GM_Mars/r,
//   n_loc(r)/n_inf = v_loc(r)/v_inf,
// and the velocity directions are isotropic (before removal of planet-blocked orbits).
//
// Particle injection from a moving, Phobos-centered local b-plane, perpendicular
// to its locally sampled velocity u.  The initial plane at t0=-L/v_loc(a) is
// r0 = r_Ph(0) + b - L*u.  b is uniform in area over [0, Bgate].
// The flux carried by each injected trajectory relative to the far-field number
// density is (v_loc(r0)/v_inf) * (v_loc(r0)-v_Ph(0) dot u).
// No 1/p_gate factor: sampling is DIRECT from the gate disk. Rate proxy per
// unit far-field number density is pi*Bgate^2*sum(weight)/NUM [m^3/s].
//
// This is a NEW physically motivated sampling model, not a validation of the
// gravity-on case. Before publishing, convergence in L, Bgate, dt, and the
// Mars-boundary treatment plus independent source validation are essential.
// No Mars atmosphere, noncircular Phobos orbit, or real size-speed distribution.
// Gravitational field of Phobos is a point mass even with OBJ surface model.
// Paired output: shield_ON excludes any impactor that crossed Mars (in the
// prior Mars-only Kepler orbit or during integration), shield_OFF includes it.
// Both are computed from IDENTICAL particle histories, not two simulations.
// The shield_OFF counterfactual uses uniform-density Mars interior gravity.
//
// Build: g++ -O3 -std=c++17 -fopenmp phobos_impact_production_v2.cpp -o phobos_impact_prod_v2
// Run: ./phobos_impact_prod_v2 --help

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>

constexpr double PI=3.1415926535897932384626433832795;
constexpr double MU_MARS=4.282837e13; // m^3/s^2
constexpr double MARS_RADIUS=3389500.0; // m
constexpr double MU_PHOBOS=7.11e5; // m^3/s^2
constexpr double V_INF_DEFAULT=7952.98; // m/s

struct Vec {double x=0,y=0,z=0;};
Vec operator+(Vec a,Vec b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
Vec operator-(Vec a,Vec b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
Vec operator*(Vec a,double q){return {a.x*q,a.y*q,a.z*q};}
Vec operator/(Vec a,double q){return a*(1.0/q);}
double dot(Vec a,Vec b){return a.x*b.x+a.y*b.y+a.z*b.z;}
Vec cross(Vec a,Vec b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double norm2(Vec a){return dot(a,a);}
double norm(Vec a){return std::sqrt(norm2(a));}
Vec unit(Vec a){return a/norm(a);}
Vec vmin(Vec a,Vec b){return {std::min(a.x,b.x),std::min(a.y,b.y),std::min(a.z,b.z)};}
Vec vmax(Vec a,Vec b){return {std::max(a.x,b.x),std::max(a.y,b.y),std::max(a.z,b.z)};}

struct Triangle {Vec a,b,c,lo,hi,centroid;};
struct Node {Vec lo,hi; int left=-1,right=-1,start=0,count=0;};
class Mesh {
public:
    std::vector<Triangle> tris;
    std::vector<int> order;
    std::vector<Node> nodes;
    double rbound=0;
    static bool segment_box(Vec a, Vec b, Vec lo, Vec hi){
        Vec d=b-a; double tlo=0.0,thi=1.0;
        for(int axis=0;axis<3;++axis){
            double x=(axis==0?a.x:axis==1?a.y:a.z);
            double dx=(axis==0?d.x:axis==1?d.y:d.z);
            double l=(axis==0?lo.x:axis==1?lo.y:lo.z);
            double h=(axis==0?hi.x:axis==1?hi.y:hi.z);
            if(std::fabs(dx)<1e-14){if(x<l||x>h)return false;}
            else {double t1=(l-x)/dx,t2=(h-x)/dx;if(t1>t2)std::swap(t1,t2);
                  tlo=std::max(tlo,t1);thi=std::min(thi,t2);if(tlo>thi)return false;}
        }
        return true;
    }
    static bool segment_tri(Vec p0,Vec p1,const Triangle &t,double &fraction){
        Vec dir=p1-p0,e1=t.b-t.a,e2=t.c-t.a;
        Vec pv=cross(dir,e2);double det=dot(e1,pv);
        if(std::fabs(det)<1e-14)return false;
        double inv=1/det;Vec tv=p0-t.a;
        double u=dot(tv,pv)*inv;if(u<-1e-12||u>1+1e-12)return false;
        Vec qv=cross(tv,e1);double v=dot(dir,qv)*inv;
        if(v<-1e-12||u+v>1+1e-12)return false;
        double s=dot(e2,qv)*inv;
        if(s<-1e-12||s>1+1e-12)return false;
        fraction=std::max(0.0,std::min(1.0,s));return true;
    }
    void add_triangle(Vec a,Vec b,Vec c){
        Triangle t; t.a=a;t.b=b;t.c=c;t.lo=vmin(vmin(a,b),c);t.hi=vmax(vmax(a,b),c);
        t.centroid=(a+b+c)/3.0;tris.push_back(t);
        rbound=std::max({rbound,norm(a),norm(b),norm(c)});
    }
    bool load(const std::string &path,double unit){
        std::ifstream in(path); if(!in){std::cerr<<"Could not open OBJ: "<<path<<"\n";return false;}
        std::vector<Vec> verts;
        std::string s;
        while(std::getline(in,s)){
            if(s.size()<2)continue;
            if(s[0]=='v'&&std::isspace((unsigned char)s[1])){
                std::istringstream ss(s.substr(1));Vec v;
                if(ss>>v.x>>v.y>>v.z) verts.push_back(v*unit);
            }else if(s[0]=='f'&&std::isspace((unsigned char)s[1])){
                std::istringstream ss(s.substr(1));std::vector<int> idx;std::string tok;
                while(ss>>tok){
                    // Accept positive OBJ vertex IDs, including f v/vt/vn.
                    int vi=std::atoi(tok.c_str());
                    if(vi<=0 || vi>(int)verts.size()){
                        std::cerr<<"Only positive, already-defined OBJ vertex IDs supported.\n";return false;
                    }
                    idx.push_back(vi-1);
                }
                for(size_t k=2;k<idx.size();++k)
                    add_triangle(verts[idx[0]],verts[idx[k-1]],verts[idx[k]]);
            }
        }
        if(tris.empty())return false;
        order.resize(tris.size());std::iota(order.begin(),order.end(),0);
        nodes.reserve(tris.size()*2); build(0,(int)order.size());return true;
    }
    int build(int start,int count){
        Node n; n.start=start;n.count=count;
        double inf=std::numeric_limits<double>::infinity();
        n.lo={inf,inf,inf};n.hi={-inf,-inf,-inf};
        for(int i=start;i<start+count;i++){
            const auto &t=tris[order[i]];n.lo=vmin(n.lo,t.lo);n.hi=vmax(n.hi,t.hi);
        }
        int id=(int)nodes.size();nodes.push_back(n);
        if(count<=8)return id;
        Vec ext=n.hi-n.lo;int axis=(ext.y>ext.x)?1:0;
        if(ext.z>(axis==0?ext.x:ext.y))axis=2;
        int mid=start+count/2;
        std::nth_element(order.begin()+start,order.begin()+mid,order.begin()+start+count,
            [&](int a,int b){
                Vec x=tris[a].centroid,y=tris[b].centroid;
                return (axis==0?x.x:axis==1?x.y:x.z)<(axis==0?y.x:axis==1?y.y:y.z);
            });
        int left=build(start,mid-start);int right=build(mid,start+count-mid);
        nodes[id].left=left;nodes[id].right=right;nodes[id].count=0;
        return id;
    }
    bool intersect(Vec a,Vec b,double &fraction,int &face_id)const{
        if(nodes.empty())return false;
        std::vector<int> stk;stk.reserve(64);stk.push_back(0);
        bool hit=false;fraction=1.0;face_id=-1;
        while(!stk.empty()){
            int ni=stk.back();stk.pop_back();const Node &n=nodes[ni];
            if(!segment_box(a,b,n.lo,n.hi))continue;
            if(n.left<0){
                for(int i=n.start;i<n.start+n.count;i++){
                    const auto &t=tris[order[i]];
                    if(!segment_box(a,b,t.lo,t.hi))continue;
                    double f;if(segment_tri(a,b,t,f)&&f<fraction){fraction=f;hit=true;face_id=order[i];}
                }
            }else{stk.push_back(n.left);stk.push_back(n.right);}
        }
        return hit;
    }
};

struct Cfg {
    int64_t n=1000000;
    double v_inf=V_INF_DEFAULT, orbit=9376000.0, radius=11100.0;
    double bgate=1500000.0, bmax=20000000.0, launch=2000000.0;
    double dt_scale=1.0,post=450.0;
    bool marsgrav=true,phgrav=true,gate_direct=true;
    bool mesh=false; std::string shape, out_on="",out_off="",out_through="",summary="";
    double obj_scale=1000.0;
    uint64_t seed=20261026;
};
void help(){
    std::cout<<"Usage: ./phobos_impact_prod_v2 [options]\n"
     <<" --num N --seed INTEGER --orbit-r METERS --radius METERS\n"
     <<" --v-inf MPS --bgate METERS --bmax METERS --launch-distance METERS\n"
     <<" --shape path.obj [--obj-unit km|m]  (omit for sphere)\n"
     <<" --mars-grav 0|1 --phobos-grav 0|1  (shield ON/OFF are always BOTH output)\n"
     <<" --sample gate|full   (gate=direct sampling within Bgate; full=from Bmax and skip)\n"
     <<" --dt-scale FLOAT --post-seconds FLOAT\n"
     <<" --out-on file.txt --out-off file.txt --out-through file.txt --summary file.txt\n";
}
bool parse(int argc,char**argv,Cfg &c){
    for(int i=1;i<argc;i++){
        std::string a=argv[i];if(a=="--help"){help();std::exit(0);}
        if(i+1>=argc){std::cerr<<"Missing value for "<<a<<"\n";return false;}
        std::string val=argv[++i];
        try{
            if(a=="--num")c.n=std::stoll(val);
            else if(a=="--seed")c.seed=std::stoull(val);
            else if(a=="--orbit-r")c.orbit=std::stod(val);
            else if(a=="--radius")c.radius=std::stod(val);
            else if(a=="--v-inf")c.v_inf=std::stod(val);
            else if(a=="--bgate")c.bgate=std::stod(val);
            else if(a=="--bmax")c.bmax=std::stod(val);
            else if(a=="--launch-distance")c.launch=std::stod(val);
            else if(a=="--post-seconds")c.post=std::stod(val);
            else if(a=="--dt-scale")c.dt_scale=std::stod(val);
            else if(a=="--mars-grav")c.marsgrav=(std::stoi(val)!=0);
            else if(a=="--phobos-grav")c.phgrav=(std::stoi(val)!=0);
            else if(a=="--shield"){std::cerr<<"--shield is replaced by paired ON/OFF output in v2; remove it.\n";return false;}
            else if(a=="--sample"){
                if(val!="full"&&val!="gate")return false;
                c.gate_direct=(val=="gate");
            }
            else if(a=="--shape"){c.shape=val;c.mesh=true;}
            else if(a=="--obj-unit"){
                if(val=="km")c.obj_scale=1000;
                else if(val=="m")c.obj_scale=1;
                else return false;
            }
            else if(a=="--out-on")c.out_on=val;
            else if(a=="--out-off")c.out_off=val;
            else if(a=="--out-through")c.out_through=val;
            else if(a=="--summary")c.summary=val;
            else{std::cerr<<"Unrecognized option "<<a<<"\n";return false;}
        }catch(const std::exception &e){std::cerr<<"Bad value: "<<a<<" "<<val<<"\n";return false;}
    }
    if(c.n<=0||c.v_inf<=0||c.orbit<MARS_RADIUS+c.radius||c.radius<=0||
       c.bgate<=0||c.bmax<c.bgate||c.launch<10*c.radius||c.dt_scale<=0||c.post<=0)return false;
    if(c.orbit-c.launch-c.bgate<=MARS_RADIUS){
        std::cerr<<"WARNING: injection disk may approach/intersect Mars; use a smaller Bgate/launch distance.\n";
    }
    if(c.out_on.empty()||c.out_off.empty()||c.out_through.empty()||c.summary.empty()){
        std::cerr<<"Provide --out-on, --out-off, --out-through, and --summary.\n";return false;
    }
    if(c.out_on==c.out_off||c.out_on==c.out_through||c.out_off==c.out_through){
        std::cerr<<"ON/OFF/through outputs must use different filenames.\n";return false;
    }
    return true;
}
struct Body {Vec r,v;};
Body phobos(double t,double phi,double orb,double omega){
    double th=phi+omega*t;
    return {{orb*std::cos(th),orb*std::sin(th),0},
            {-orb*omega*std::sin(th),orb*omega*std::cos(th),0}};
}
Vec bodypos(Vec rel,double t,double phi,double omega){
    double th=phi+omega*t,c=std::cos(th),s=std::sin(th);
    return {-rel.x*c-rel.y*s,rel.x*s-rel.y*c,rel.z};
}
struct State {Vec r,v;double t;};
Vec acceleration(const State &s,const Cfg &c,double phi,double omega){
    Vec a{};
    if(c.marsgrav){
        double rn=norm(s.r);
        if(rn>1.0){
            double d=std::max(rn,MARS_RADIUS);
            a=a+s.r*(-MU_MARS/(d*d*d)); // uniform-density interior, no singularity
        }
    }
    if(c.phgrav){
        Vec d=s.r-phobos(s.t,phi,c.orbit,omega).r;
        double r2=norm2(d),r=std::sqrt(r2);
        if(r>1.0)a=a+d*(-MU_PHOBOS/(r2*r));
    }
    return a;
}
State advance(State s,double dt,const Cfg &c,double phi,double omega){
    Vec k1r=s.v,k1v=acceleration(s,c,phi,omega);
    State a{s.r+k1r*(dt/2),s.v+k1v*(dt/2),s.t+dt/2};
    Vec k2r=a.v,k2v=acceleration(a,c,phi,omega);
    State b{s.r+k2r*(dt/2),s.v+k2v*(dt/2),s.t+dt/2};
    Vec k3r=b.v,k3v=acceleration(b,c,phi,omega);
    State d{s.r+k3r*dt,s.v+k3v*dt,s.t+dt};
    Vec k4r=d.v,k4v=acceleration(d,c,phi,omega);
    return {s.r+(k1r+(k2r+k3r)*2.0+k4r)*(dt/6.0),
            s.v+(k1v+(k2v+k3v)*2.0+k4v)*(dt/6.0),s.t+dt};
}
bool sphere_segment(Vec p0,Vec p1,double radius,double &frac){
    Vec d=p1-p0;double a=norm2(d),b=2*dot(p0,d),cc=norm2(p0)-radius*radius;
    if(cc<=0){frac=0;return true;} if(a<=0)return false;
    double disc=b*b-4*a*cc;if(disc<0)return false;
    double t=(-b-std::sqrt(disc))/(2*a);
    if(t<0||t>1)return false;
    frac=t;return true;
}
// If a particle at r0,v0 is OUTBOUND and its Mars-only hyperbolic orbit
// has periapsis < Mars radius, a physical planet must have blocked it earlier.
bool blocked_in_past(Vec r0,Vec v0,const Cfg &c){
    
    double r=norm(r0);if(r<=MARS_RADIUS)return true;
    if(!c.marsgrav){
        if(dot(r0,v0)<=0)return false;
        // Past straight path closest approach and positive past flight duration.
        double tau=dot(r0,v0)/norm2(v0);
        return norm(r0-v0*tau)<=MARS_RADIUS;
    }
    if(dot(r0,v0)<=0)return false;
    double h2=norm2(cross(r0,v0));
    double energy=0.5*norm2(v0)-MU_MARS/r;
    double ecc=std::sqrt(std::max(0.0,1.0+2.0*energy*h2/(MU_MARS*MU_MARS)));
    double peri=h2/(MU_MARS*(1+ecc));
    return peri<=MARS_RADIUS;
}
struct Result {
    bool hit=false,crossed_mars=false;
    double lat=0,lon=0,beta=0,vrel=0;
    // alpha: angle ABOVE local tangent plane (0=grazing, 90=normal)
    double alpha_deg=0,surface_speed=0,incidence_cos_raw=0;
    int face_id=-1; // zero-based OBJ triangle ID; -1 for sphere
    Vec velocity_body{},normal_body{}; // m/s and unit outward normal
};
Result trace(const Cfg &c,const Mesh *mesh,State s,double phi,double omega,double tend){
    const double rbound=mesh?mesh->rbound:c.radius;
    bool crossed_mars = norm(s.r)<=MARS_RADIUS;

    for(int j=0;j<500000 && s.t<tend;j++){
        Body ph=phobos(s.t,phi,c.orbit,omega);
        double dph=norm(s.r-ph.r);
        double rv=norm(s.v-ph.v);
        double dt=std::min(16.0,std::max(0.02,0.25*dph/std::max(1.0,rv)));
        if(dph<50*rbound)dt=std::min(dt,0.5);
        if(dph<8*rbound)dt=std::min(dt,0.07);
        {double mars_d=norm(s.r)-MARS_RADIUS;
            if(mars_d<400000)dt=std::min(dt,1.0);
        }
        dt=std::min(dt*c.dt_scale,tend-s.t);
        if(dt<=0)break;
        State next=advance(s,dt,c,phi,omega);
        if(!std::isfinite(next.r.x)||!std::isfinite(next.r.y)||!std::isfinite(next.r.z))break;
        Body ph1=phobos(next.t,phi,c.orbit,omega);
        Vec rel0=s.r-ph.r,rel1=next.r-ph1.r;
        double tf=1.0;bool collided=false;int face_id=-1;
        Vec p0=bodypos(rel0,s.t,phi,omega),p1=bodypos(rel1,next.t,phi,omega);
        if(mesh){
            // The chord approximation is resolved with dt-scale convergence runs.
            if(Mesh::segment_box(p0,p1,mesh->nodes[0].lo,mesh->nodes[0].hi))
                collided=mesh->intersect(p0,p1,tf,face_id);
        }else collided=sphere_segment(p0,p1,c.radius,tf);
        // One trajectory for BOTH shielding assumptions.  Mark Mars traversal
        // only when it occurred BEFORE the Phobos collision on this step.
        double mf=1.0;
        bool mars_hit=sphere_segment(s.r,next.r,MARS_RADIUS,mf);
        if(mars_hit && (!collided || mf<=tf))crossed_mars=true;
        if(collided){
            double hit_t=s.t+tf*dt;
            Vec p=p0+(p1-p0)*tf;
            double lat=std::atan2(p.z,std::hypot(p.x,p.y))*180/PI;
            double lon=std::atan2(p.y,p.x)*180/PI;
            // apex body-axis = -Y (lon -90 degrees)
            double beta=std::acos(std::max(-1.0,std::min(1.0,-p.y/norm(p))))*180/PI;
            State hitstate=advance(s,tf*dt,c,phi,omega);
            Body hitph=phobos(hit_t,phi,c.orbit,omega);
            // Preserve v2 vrel: inertial particle velocity relative to Phobos COM.
            // For crater incidence use the actual LOCAL co-rotating surface velocity:
            // v_surface = v_Ph_COM + omega cross (r_hit - r_Ph_COM).
            Vec vrel_com=hitstate.v-hitph.v;
            Vec rel_inertial=hitstate.r-hitph.r;
            Vec omega_vec{0,0,omega};
            Vec v_surface_inertial=vrel_com-cross(omega_vec,rel_inertial);
            Vec v_surface_body=bodypos(v_surface_inertial,hit_t,phi,omega);
            Vec n_body{};
            if(mesh){
                const Triangle &tri=mesh->tris.at((size_t)face_id);
                n_body=unit(cross(tri.b-tri.a,tri.c-tri.a));
            }else{
                n_body=unit(p); // radial outward normal for sphere
            }
            double vs=norm(v_surface_body);
            double mu= -dot(v_surface_body,n_body)/vs;
            double alpha=std::asin(std::clamp(mu,0.0,1.0))*180.0/PI;
            Result out;
            out.hit=true;out.crossed_mars=crossed_mars;
            out.lat=lat;out.lon=lon;out.beta=beta;
            out.vrel=norm(vrel_com); // EXACT v2 definition, unchanged
            out.alpha_deg=alpha;out.surface_speed=vs;
            out.incidence_cos_raw=mu;
            out.face_id=face_id;
            out.velocity_body=v_surface_body;out.normal_body=n_body;
            return out;
        }
        s=next;
        // A conservative safe escape: never terminate merely because distance grows.
    }
    Result out;out.crossed_mars=crossed_mars;return out;
}
struct Accum {
    uint64_t hits_off=0,hits_on=0,hits_through=0,prior=0,forward=0,skipped=0;
    double w_on=0,w2_on=0,w_off=0,w2_off=0,w_through=0,w2_through=0;
    std::string lines_on,lines_off,lines_through;
};
int main(int argc,char **argv){
    Cfg cfg;if(!parse(argc,argv,cfg)){help();return 2;}
    Mesh mesh;
    if(cfg.mesh){if(!mesh.load(cfg.shape,cfg.obj_scale)){return 3;}
        std::cerr<<"OBJ triangles="<<mesh.tris.size()<<" bound_km="<<mesh.rbound/1000<<"\n";
    }
    double v_orb=std::sqrt(MU_MARS/cfg.orbit),omega=v_orb/cfg.orbit;
    double v_ref=cfg.marsgrav?std::sqrt(cfg.v_inf*cfg.v_inf+2*MU_MARS/cfg.orbit):cfg.v_inf;
    double lead_time=cfg.launch/v_ref;
    double disk=cfg.gate_direct?cfg.bgate:cfg.bmax;
    // Each trajectory output weight is physical phase-space flux, not gate importance weight.
    // Rate proxy from full-disk proposals must instead use bmax^2; from gate use bgate^2.
    std::ofstream out_on(cfg.out_on),out_off(cfg.out_off),out_through(cfg.out_through);
    if(!out_on||!out_off||!out_through){
        std::cerr<<"Cannot open one of the paired output paths.\n";return 4;
    }
    const char *header="# lat_deg lon_deg beta_deg vrel_mps w_flux_mps b_initial_km u_x u_y u_z phase_rad alpha_deg v_surface_mps face_id v_surface_body_x_mps v_surface_body_y_mps v_surface_body_z_mps normal_body_x normal_body_y normal_body_z incidence_cos_raw\n";
    out_on<<header;out_off<<header;out_through<<header;
    const int nt=omp_get_max_threads();std::vector<Accum> a(nt);
    auto start=omp_get_wtime();
    #pragma omp parallel
    {
        int tid=omp_get_thread_num();auto &stats=a[tid];stats.lines_on.reserve(1<<17);stats.lines_off.reserve(1<<17);stats.lines_through.reserve(1<<17);
        // Deterministic reproducibility for a fixed thread count and OMP scheduling.
        std::mt19937_64 rng(cfg.seed+0x9e3779b97f4a7c15ULL*(uint64_t)(tid+1));
        std::uniform_real_distribution<double> U(0,1);
        #pragma omp for schedule(static)
        for(int64_t k=0;k<cfg.n;k++){
            double phi=2*PI*U(rng);
            double zz=2*U(rng)-1,theta=2*PI*U(rng),h=std::sqrt(std::max(0.0,1-zz*zz));
            Vec u{h*std::cos(theta),h*std::sin(theta),zz};
            Vec helper=(std::fabs(u.z)<0.9)?Vec{0,0,1}:Vec{1,0,0};
            Vec e1=unit(cross(helper,u)),e2=cross(u,e1);
            double bmag=disk*std::sqrt(U(rng)),psi=2*PI*U(rng);
            if(bmag>cfg.bgate){stats.skipped++;continue;}
            Vec bvec=e1*(bmag*std::cos(psi))+e2*(bmag*std::sin(psi));
            Body enc=phobos(0,phi,cfg.orbit,omega);
            Vec pos=enc.r+bvec-u*cfg.launch;
            double r=norm(pos);
            // Consistent Liouville speed in the counterfactual transparent,
            // uniform-density Mars potential, including injection points inside Mars.
            double grav_potential_twice=(r>=MARS_RADIUS)
              ? 2.0*MU_MARS/r
              : MU_MARS/MARS_RADIUS*(3.0-r*r/(MARS_RADIUS*MARS_RADIUS));
            double vloc=cfg.marsgrav?std::sqrt(cfg.v_inf*cfg.v_inf+grav_potential_twice):cfg.v_inf;
            Vec vel=u*vloc;
            bool prior=blocked_in_past(pos,vel,cfg);
            if(prior)stats.prior++;
            // Flux crossing the encounter-centered moving sampling surface.
            double w=(vloc/cfg.v_inf)*(vloc-dot(enc.v,u));
            if(w<=0)continue;
            Result res=trace(cfg,cfg.mesh?&mesh:nullptr,{pos,vel,-lead_time},phi,omega,cfg.post);
            if(res.crossed_mars)stats.forward++;
            if(res.hit){
                const bool through = prior || res.crossed_mars;
                stats.hits_off++;stats.w_off+=w;stats.w2_off+=w*w;
                char line[750];int l=std::snprintf(line,sizeof(line),
                     "%.8f %.8f %.8f %.8f %.8f %.8f %.9f %.9f %.9f %.9f %.8f %.8f %d %.8f %.8f %.8f %.9f %.9f %.9f %.9f\n",
                     res.lat,res.lon,res.beta,res.vrel,w,bmag/1000,u.x,u.y,u.z,phi,
                     res.alpha_deg,res.surface_speed,res.face_id,
                     res.velocity_body.x,res.velocity_body.y,res.velocity_body.z,
                     res.normal_body.x,res.normal_body.y,res.normal_body.z,
                     res.incidence_cos_raw);
                stats.lines_off.append(line,l);
                if(through){
                    stats.hits_through++;stats.w_through+=w;stats.w2_through+=w*w;
                    stats.lines_through.append(line,l);
                }else{
                    stats.hits_on++;stats.w_on+=w;stats.w2_on+=w*w;
                    stats.lines_on.append(line,l);
                }
            }
        }
    }
    uint64_t hits_on=0,hits_off=0,hits_through=0,prior=0,forward=0,skipped=0;
    double sw_on=0,sw2_on=0,sw_off=0,sw2_off=0,sw_through=0,sw2_through=0;
    for(auto &x:a){
        out_on<<x.lines_on;out_off<<x.lines_off;out_through<<x.lines_through;
        hits_on+=x.hits_on;hits_off+=x.hits_off;hits_through+=x.hits_through;
        prior+=x.prior;forward+=x.forward;skipped+=x.skipped;
        sw_on+=x.w_on;sw2_on+=x.w2_on;
        sw_off+=x.w_off;sw2_off+=x.w2_off;
        sw_through+=x.w_through;sw2_through+=x.w2_through;
    }
    out_on.close();out_off.close();out_through.close();
    const double area=PI*disk*disk;
    auto rate=[&](double sw){return area*sw/(double)cfg.n;};
    auto se=[&](double w2){return area*std::sqrt(w2)/(double)cfg.n;};
    std::ostringstream os;os<<std::setprecision(12);
    os<<"model=production_local_liouville_v3_angles_paired_shielding\n"
      <<"shape="<<(cfg.mesh?cfg.shape:"sphere")<<"\n"
      <<"surface_type="<<(cfg.mesh?"OBJ":"sphere")<<"\n"
      <<"surface_bound_km="<<(cfg.mesh?mesh.rbound/1000:cfg.radius/1000)<<"\n"
      <<"num="<<cfg.n<<"\nseed="<<cfg.seed<<"\nthreads="<<nt
      <<"\norbit_r_km="<<cfg.orbit/1000<<"\norbital_speed_mps="<<v_orb
      <<"\nphobos_radius_km="<<cfg.radius/1000
      <<"\nv_inf_mps="<<cfg.v_inf
      <<"\nBgate_km="<<cfg.bgate/1000<<"\nBmax_km="<<cfg.bmax/1000
      <<"\nsampling="<<(cfg.gate_direct?"gate":"full")
      <<"\nlaunch_distance_km="<<cfg.launch/1000
      <<"\nlaunch_time_seconds="<<lead_time
      <<"\npost_seconds="<<cfg.post<<"\ndt_scale="<<cfg.dt_scale
      <<"\nmars_gravity="<<cfg.marsgrav<<"\nphobos_gravity="<<cfg.phgrav
      <<"\nshielding_mode=paired_single_trajectory\n"
      <<"prior_mars_crossed_particles="<<prior
      <<"\nforward_mars_crossed_particles="<<forward
      <<"\nskipped_by_gate="<<skipped
      <<"\nhits_shield_on="<<hits_on
      <<"\nhits_shield_off="<<hits_off
      <<"\nhits_through_only="<<hits_through
      <<"\nweight_sum_on_mps="<<sw_on
      <<"\nweight_sum_off_mps="<<sw_off
      <<"\nweight_sum_through_mps="<<sw_through
      <<"\nrate_proxy_shield_on_m3_per_s_per_ninf="<<rate(sw_on)
      <<"\nrate_proxy_shield_off_m3_per_s_per_ninf="<<rate(sw_off)
      <<"\nrate_proxy_through_only_m3_per_s_per_ninf="<<rate(sw_through)
      <<"\nrate_1sigma_shield_on="<<se(sw2_on)
      <<"\nrate_1sigma_shield_off="<<se(sw2_off)
      <<"\nrate_1sigma_through_only="<<se(sw2_through)
      <<"\nimpact_angle_definition=from_local_tangent_plane_0_grazing_90_normal"
      <<"\nimpact_angle_surface_velocity=Phobos_COM_plus_corotation"
      <<"\nimpact_angle_normal=OBJ_face_winding_or_sphere_radial"
      <<"\nimpact_catalog_columns=20"
      <<"\noutput_shield_on="<<cfg.out_on
      <<"\noutput_shield_off="<<cfg.out_off
      <<"\noutput_through_only="<<cfg.out_through
      <<"\nelapsed_seconds="<<(omp_get_wtime()-start)<<"\n";
    std::cerr<<os.str();
    std::ofstream log(cfg.summary);
    if(!log){std::cerr<<"Cannot open summary file "<<cfg.summary<<"\n";return 5;}
    log<<os.str();log.close();
    if(hits_off!=hits_on+hits_through || std::fabs(rate(sw_off)-rate(sw_on)-rate(sw_through))>1e-8*std::max(1.0,rate(sw_off))){
        std::cerr<<"ERROR: Paired ON+through = OFF invariant failed.\n";return 6;
    }
    return 0;
}
