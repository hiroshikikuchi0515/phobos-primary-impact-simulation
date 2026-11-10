// Phobos primary-impact source-sphere production v5 (2026-11).
// Compiles against the unchanged phobos_impact_production_v3.cpp.
// v5: CLI safety/documentation/provenance only; dynamical and sampling algorithms match v4.
// Source sphere: uniform area + isotropic Mars-frame velocity direction,
// weighted by the inward flux of the local Mars-only Liouville distribution.
// IMPORTANT: the flux estimator assumes every collision trajectory crosses the
// source sphere inward exactly once before collision (single-inward-crossing).
// Phobos gravity, if enabled, is applied only during trace; the source
// distribution remains the Mars-only far-field stationary reservoir.
// Catalog: 20 columns aligned with production v3 except column 6, which is
// source_radius_km rather than b_initial_km; 6 extra source diagnostics.
// Compile in your production work directory using C++17 and OpenMP.

#define main phobos_unused_original_production_main
#include "phobos_impact_production_v3.cpp"
#undef main

#include <filesystem>
#include <cerrno>
#include <iomanip>
#include <memory>
#include <stdexcept>

namespace fs = std::filesystem;

struct CatalogMoments {
    uint64_t hits = 0;
    double sum = 0, sum2 = 0;
    double mars = 0, mars2 = 0, anti = 0, anti2 = 0;
    double lead = 0, lead2 = 0, trail = 0, trail2 = 0;

    void add(double w, double longitude_deg) {
        ++hits;
        sum += w;
        sum2 += w*w;
        const double longitude = longitude_deg * PI / 180.0;
        if (std::cos(longitude) > 0) {
            mars += w;
            mars2 += w*w;
        } else {
            anti += w;
            anti2 += w*w;
        }
        if (std::sin(longitude) < 0) {
            lead += w;
            lead2 += w*w;
        } else {
            trail += w;
            trail2 += w*w;
        }
    }

    void merge(const CatalogMoments& b) {
        hits += b.hits;
        sum += b.sum; sum2 += b.sum2;
        mars += b.mars; mars2 += b.mars2;
        anti += b.anti; anti2 += b.anti2;
        lead += b.lead; lead2 += b.lead2;
        trail += b.trail; trail2 += b.trail2;
    }
};

struct ThreadCatalogs {
    CatalogMoments on, off, through;
    uint64_t inward_proposals = 0, prior_inward = 0;
    uint64_t forward_inward = 0, recorded_lines = 0;
    std::ofstream fon, foff, fthrough;
};

struct RatioStat {
    double value = 0, se = 0;
};

static RatioStat weighted_ratio(double a, double a2, double b, double b2) {
    if (!(b > 0)) return {std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()};
    const double r = a / b;
    return {r, std::sqrt(a2 + r*r*b2)/b};
}

static double rate_sigma(const CatalogMoments& m, double area, double n) {
    if (n <= 1) return std::numeric_limits<double>::quiet_NaN();
    // Each proposal (including nonhits) is an independent random draw.
    // Unbiased sample variance of the mean, not Poisson-only approximation.
    const double centered = std::max(0.0, m.sum2 - m.sum*m.sum/n);
    return area*std::sqrt(centered/(n*(n-1)));
}

static void print_moments(std::ostream& os, const char* label,
                          const CatalogMoments& m, double area, double n) {
    const auto lt = weighted_ratio(m.lead,m.lead2,m.trail,m.trail2);
    const auto ma = weighted_ratio(m.mars,m.mars2,m.anti,m.anti2);
    os << label << " hits=" << m.hits << "\n"
       << "  Rate=" << std::setprecision(12) << area*m.sum/n
       << " +/- " << std::setprecision(7) << rate_sigma(m,area,n) << "\n"
       << "  L/T=" << std::setprecision(9) << lt.value << " +/- " << lt.se << "\n"
       << "  M/A=" << ma.value << " +/- " << ma.se << "\n";
}

static std::string part_name(const std::string& stem, int tid) {
    return stem + ".source_sphere_part_" + std::to_string(tid);
}

static bool merge_catalog(const std::string& dest, int nt,
                          const char* header) {
    std::ofstream out(dest, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    out << header << '\n';
    for (int i=0; i<nt; ++i) {
        const auto part = part_name(dest,i);
        std::ifstream in(part, std::ios::binary);
        if (!in) return false;
        // Inserting an empty streambuf sets failbit on some implementations.
        if (fs::file_size(part) > 0) out << in.rdbuf();
        if (!out) return false;
        in.close();
    }
    out.close();
    if (!out) return false;
    for (int i=0; i<nt; ++i) fs::remove(part_name(dest,i));
    return true;
}

int main(int argc, char** argv) {
    // Parse the one new flag ourselves; delegate every existing flag to v3.
    double source_radius_km = 30.0;
    std::string summary_path;
    std::vector<char*> args;
    args.reserve(argc);
    args.push_back(argv[0]);
    for (int i=1; i<argc; ++i) {
        const std::string opt(argv[i]);
        if (opt == "--help" || opt == "-h") {
            std::cout
              << "Source-sphere production v5 (Mars-centered Liouville injection)\n"
              << "Main options: --num N --seed S --orbit-r METERS --v-inf MPS\n"
              << "          --shape file.obj [--obj-unit km|m] (omit for sphere)\n"
              << "          --mars-grav 0|1 --phobos-grav 0|1\n"
              << "          --source-radius-km R --post-seconds T --dt-scale FLOAT\n"
              << "          --out-on PATH --out-off PATH --out-through PATH --summary PATH\n"
              << "Optional: --radius METERS (spherical target; default 11100)\n"
              << "Source-sphere proposals are always uniform in area; inward flux is weighted.\n"
              << "The Phobos gravitational field is a point mass and is switched ON in production.\n"
              << "Legacy gate flags (--sample, --bgate, --bmax, --launch-distance)\n"
              << "are prohibited here; they do not define the source-sphere sampler.\n";
            return 0;
        }
        if (opt == "--sample" || opt == "--bgate" || opt == "--bmax" ||
            opt == "--launch-distance") {
            std::cerr << "Gate-specific flag " << opt
                      << " is NOT applicable to source-sphere production.\n";
            return 2;
        }
        if (opt == "--mars-grav" || opt == "--phobos-grav") {
            if (i+1 >= argc || (std::string(argv[i+1]) != "0" &&
                                std::string(argv[i+1]) != "1")) {
                std::cerr << opt << " requires exactly 0 or 1\n";
                return 2;
            }
        }
        if (opt == "--source-radius-km") {
            if (++i >= argc) {
                std::cerr << "--source-radius-km needs a value\n";
                return 2;
            }
            char* end=nullptr;
            errno=0;
            source_radius_km=std::strtod(argv[i],&end);
            if (errno || !end || *end || !std::isfinite(source_radius_km)) {
                std::cerr << "Invalid --source-radius-km\n";
                return 2;
            }
        } else {
            if (opt == "--summary" && i+1 < argc) summary_path=argv[i+1];
            args.push_back(argv[i]);
        }
    }
    if (summary_path.empty()) {
        std::cerr << "Please supply --summary <new_summary_path>\n";
        return 2;
    }
    Cfg cfg;
    if (!parse(static_cast<int>(args.size()),args.data(),cfg)) return 2;

    Mesh mesh;
    const Mesh* shape = nullptr;
    if (cfg.mesh) {
        if (!mesh.load(cfg.shape,cfg.obj_scale)) return 3;
        shape = &mesh;
        std::cerr << "OBJ triangles=" << mesh.tris.size()
                  << " bound_km=" << mesh.rbound/1000.0 << '\n';
    }
    const double Rsrc=1000.0*source_radius_km;
    const double rbound=shape ? mesh.rbound : cfg.radius;
    // The OBJ header is descriptive metadata; independently calculate true facet area.
    double mesh_area_km2=0.0;
    if (shape) for (const auto& tri : mesh.tris)
        mesh_area_km2 += 0.5*norm(cross(tri.b-tri.a,tri.c-tri.a))/1e6;
    if (!(Rsrc > rbound) || !(cfg.orbit-Rsrc > MARS_RADIUS)) {
        std::cerr << "Invalid source radius: must enclose shape, avoid Mars\n";
        return 4;
    }
    if (cfg.n <= 1 || cfg.v_inf <= 0 || cfg.post <= 0) {
        std::cerr << "Invalid N, v-inf or post-seconds\n";
        return 4;
    }
    // Note: the surface density is Mars-only, even if Phobos gravity is ON.
    const double omega=std::sqrt(MU_MARS/(cfg.orbit*cfg.orbit*cfg.orbit));
    const double area=4.0*PI*Rsrc*Rsrc;
    const int nt=omp_get_max_threads();
    const std::string paths[4]={cfg.out_on,cfg.out_off,cfg.out_through,summary_path};
    for (const auto& path : paths) {
        if (path.empty() || fs::exists(path)) {
            std::cerr << "Output absent or already exists: " << path << '\n';
            return 5;
        }
        const fs::path parent=fs::path(path).parent_path();
        if (!parent.empty() && !fs::is_directory(parent)) {
            std::cerr << "Create output directory first: " << parent << '\n';
            return 5;
        }
    }
    if (paths[0]==paths[1] || paths[0]==paths[2] || paths[0]==paths[3] ||
        paths[1]==paths[2] || paths[1]==paths[3] || paths[2]==paths[3]) {
        std::cerr << "All output filenames must differ\n";
        return 5;
    }
    for (int tid=0;tid<nt;++tid) {
        for (int j=0;j<3;++j) {
            if (fs::exists(part_name(paths[j],tid))) {
                std::cerr << "Stale temporary file: " << part_name(paths[j],tid) << '\n';
                return 5;
            }
        }
    }
    std::vector<ThreadCatalogs> stats(nt);
    for (int tid=0; tid<nt; ++tid) {
        stats[tid].fon.open(part_name(cfg.out_on,tid),std::ios::binary);
        stats[tid].foff.open(part_name(cfg.out_off,tid),std::ios::binary);
        stats[tid].fthrough.open(part_name(cfg.out_through,tid),std::ios::binary);
        if (!stats[tid].fon || !stats[tid].foff || !stats[tid].fthrough) {
            std::cerr << "Cannot create temporary catalog files\n";
            return 5;
        }
    }

    const double start=omp_get_wtime();
    #pragma omp parallel
    {
        const int tid=omp_get_thread_num();
        auto& st=stats[tid];
        std::mt19937_64 rng(cfg.seed +
            0x9e3779b97f4a7c15ULL*(static_cast<uint64_t>(tid)+1));
        std::uniform_real_distribution<double> U(0,1);
        #pragma omp for schedule(static)
        for (int64_t i=0; i<cfg.n; ++i) {
            const double phi=2*PI*U(rng);
            const double zn=2*U(rng)-1,an=2*PI*U(rng);
            const double hn=std::sqrt(std::max(0.0,1-zn*zn));
            const Vec outward{hn*std::cos(an),hn*std::sin(an),zn};
            const double zu=2*U(rng)-1,au=2*PI*U(rng);
            const double hu=std::sqrt(std::max(0.0,1-zu*zu));
            const Vec u{hu*std::cos(au),hu*std::sin(au),zu};
            const Body ph=phobos(0.0,phi,cfg.orbit,omega);
            const Vec pos=ph.r+outward*Rsrc;
            const double r=norm(pos);
            const double vloc=cfg.marsgrav
                ? std::sqrt(cfg.v_inf*cfg.v_inf+2.0*MU_MARS/r)
                : cfg.v_inf;
            const Vec vel=u*vloc;
            const double inward=-dot(vel-ph.v,outward);
            if (!(inward > 0)) continue;
            ++st.inward_proposals;
            const double w=(vloc/cfg.v_inf)*inward;
            const bool prior=blocked_in_past(pos,vel,cfg);
            if (prior) ++st.prior_inward;
            const Result result=trace(cfg,shape,{pos,vel,0.0},phi,omega,cfg.post);
            if (result.crossed_mars) ++st.forward_inward;
            if (!result.hit) continue;
            const bool through=prior || result.crossed_mars;
            st.off.add(w,result.lon);
            if (through) st.through.add(w,result.lon);
            else st.on.add(w,result.lon);

            // Do not silently rename source radius to an impact parameter.
            // The first 20 columns match v3 apart from the sixth column.
            char line[1024];
            const int length=std::snprintf(line,sizeof(line),
                "%.8f %.8f %.8f %.8f %.8f %.8f "
                "%.9f %.9f %.9f %.9f %.8f %.8f %d "
                "%.8f %.8f %.8f %.9f %.9f %.9f %.9f "
                "%.9f %.9f %.9f %.8f %.8f %.8f\n",
                result.lat,result.lon,result.beta,result.vrel,w,source_radius_km,
                u.x,u.y,u.z,phi,result.alpha_deg,result.surface_speed,
                result.face_id,result.velocity_body.x,result.velocity_body.y,
                result.velocity_body.z,result.normal_body.x,result.normal_body.y,
                result.normal_body.z,result.incidence_cos_raw,
                outward.x,outward.y,outward.z,r/1000.0,vloc,inward);
            if (length < 0 || length >= static_cast<int>(sizeof(line))) continue;
            st.foff.write(line,length);
            if (through) st.fthrough.write(line,length);
            else st.fon.write(line,length);
            ++st.recorded_lines;
        }
    }
    const double elapsed=omp_get_wtime()-start;
    for (auto& st:stats) {
        st.fon.close(); st.foff.close(); st.fthrough.close();
        if (!st.fon || !st.foff || !st.fthrough) {
            std::cerr << "Catalog write failed; keep part files for recovery\n";
            return 6;
        }
    }
    CatalogMoments on,off,through;
    uint64_t inward=0, prior=0,forward=0,lines=0;
    for (const auto& st:stats) {
        on.merge(st.on); off.merge(st.off); through.merge(st.through);
        inward+=st.inward_proposals;
        prior+=st.prior_inward;
        forward+=st.forward_inward;
        lines+=st.recorded_lines;
    }
    if (lines != off.hits || on.hits + through.hits != off.hits) {
        std::cerr << "Internal catalog count mismatch\n";
        return 6;
    }
    const char* header =
        "# lat_deg lon_deg beta_deg vrel_mps w_flux_mps source_radius_km "
        "u_x u_y u_z phase_rad alpha_deg v_surface_mps face_id "
        "v_surface_body_x_mps v_surface_body_y_mps v_surface_body_z_mps "
        "normal_body_x normal_body_y normal_body_z incidence_cos_raw "
        "source_normal_inertial_x source_normal_inertial_y source_normal_inertial_z "
        "source_mars_distance_km source_vloc_mps inward_flux_mps";
    if (!merge_catalog(cfg.out_off,nt,header) ||
        !merge_catalog(cfg.out_on,nt,header) ||
        !merge_catalog(cfg.out_through,nt,header)) {
        std::cerr << "Failed to merge partial catalogs\n";
        return 6;
    }
    const double N=static_cast<double>(cfg.n);
    const double blocked_frac=off.sum>0 ? through.sum/off.sum :
        std::numeric_limits<double>::quiet_NaN();
    const double blocked_se=off.sum>0 ?
        std::sqrt((1-blocked_frac)*(1-blocked_frac)*through.sum2 +
                  blocked_frac*blocked_frac*on.sum2)/off.sum :
        std::numeric_limits<double>::quiet_NaN();
    std::cout << "Source sphere radius: " << source_radius_km << " km\n"
              << "N: " << cfg.n << "\n"
              << "v_inf: " << cfg.v_inf << " m/s\n"
              << "post_seconds: " << cfg.post << "\n";
    print_moments(std::cout,"SHIELD OFF",off,area,N);
    print_moments(std::cout,"SHIELD ON ",on,area,N);
    print_moments(std::cout,"SHIELD THROUGH",through,area,N);
    std::cout << "Shield fraction: " << std::setprecision(9)
              << 100*blocked_frac << " +/- " << 100*blocked_se << " %\n";
    std::ofstream summary(summary_path);
    if (!summary) {
        std::cerr << "Cannot open summary path\n";
        return 6;
    }
    summary << std::setprecision(15)
      << "model=source_sphere_production_v5\n"
      << "shape=" << (cfg.mesh?"obj":"sphere") << '\n'
      << "surface_bound_km=" << rbound/1000 << '\n'
      << "shape_file=" << (cfg.mesh ? cfg.shape : "none") << '\n'
      << "mesh_facets=" << (shape ? mesh.tris.size() : 0) << '\n'
      << "mesh_surface_area_km2=" << mesh_area_km2 << '\n'
      << "phobos_gravity_model=point_mass_COM\n"
      << "source_exterior_phobos_gravity=neglected_Mars_only_injection\n"
      << "obj_scale_m_per_obj_unit=" << cfg.obj_scale << '\n'
      << "num=" << cfg.n << '\n'
      << "seed=" << cfg.seed << '\n'
      << "threads=" << nt << '\n'
      << "orbit_r_km=" << cfg.orbit/1000 << '\n'
      << "phobos_radius_km=" << cfg.radius/1000 << '\n'
      << "v_inf_mps=" << cfg.v_inf << '\n'
      << "source_radius_km=" << source_radius_km << '\n'
      << "source_sphere_area_m2=" << area << '\n'
      << "post_seconds=" << cfg.post << '\n'
      << "dt_scale=" << cfg.dt_scale << '\n'
      << "mars_gravity=" << cfg.marsgrav << '\n'
      << "phobos_gravity=" << cfg.phgrav << '\n'
      << "shielding_mode=paired_single_trajectory\n"
      << "source_distribution=Mars_only_local_Liouville\n"
      << "source_crossing_assumption=single_inward_crossing\n"
      << "source_weights=local_number_density_ratio_times_inward_relative_speed\n"
      << "inward_proposals=" << inward << '\n'
      << "prior_mars_crossed_inward_proposals=" << prior << '\n'
      << "forward_mars_crossed_inward_proposals=" << forward << '\n'
      << "hits_shield_on=" << on.hits << '\n'
      << "hits_shield_off=" << off.hits << '\n'
      << "hits_through_only=" << through.hits << '\n'
      << "weight_sum_on_mps=" << on.sum << '\n'
      << "weight_sum_off_mps=" << off.sum << '\n'
      << "weight_sum_through_mps=" << through.sum << '\n'
      << "rate_proxy_shield_on_m3_per_s_per_ninf=" << area*on.sum/N << '\n'
      << "rate_proxy_shield_off_m3_per_s_per_ninf=" << area*off.sum/N << '\n'
      << "rate_proxy_through_only_m3_per_s_per_ninf=" << area*through.sum/N << '\n'
      << "rate_1sigma_shield_on=" << rate_sigma(on,area,N) << '\n'
      << "rate_1sigma_shield_off=" << rate_sigma(off,area,N) << '\n'
      << "rate_1sigma_through_only=" << rate_sigma(through,area,N) << '\n'
      << "LT_shield_on=" << weighted_ratio(on.lead,on.lead2,on.trail,on.trail2).value << '\n'
      << "LT_shield_on_1sigma=" << weighted_ratio(on.lead,on.lead2,on.trail,on.trail2).se << '\n'
      << "LT_shield_off=" << weighted_ratio(off.lead,off.lead2,off.trail,off.trail2).value << '\n'
      << "LT_shield_off_1sigma=" << weighted_ratio(off.lead,off.lead2,off.trail,off.trail2).se << '\n'
      << "MA_shield_on=" << weighted_ratio(on.mars,on.mars2,on.anti,on.anti2).value << '\n'
      << "MA_shield_on_1sigma=" << weighted_ratio(on.mars,on.mars2,on.anti,on.anti2).se << '\n'
      << "MA_shield_off=" << weighted_ratio(off.mars,off.mars2,off.anti,off.anti2).value << '\n'
      << "MA_shield_off_1sigma=" << weighted_ratio(off.mars,off.mars2,off.anti,off.anti2).se << '\n'
      << "shield_fraction_percent=" << 100*blocked_frac << '\n'
      << "shield_fraction_1sigma_percent=" << 100*blocked_se << '\n'
      << "impact_catalog_columns=26\n"
      << "impact_angle_definition=from_local_tangent_plane_0_grazing_90_normal\n"
      << "impact_angle_surface_velocity=Phobos_COM_plus_corotation\n"
      << "impact_angle_normal=OBJ_face_winding_or_sphere_radial\n"
      << "catalog_6th_column=source_radius_km_NOT_b_initial_km\n"
      << "flux_weight_unit=m_per_s\n"
      << "rate_normalization=4pi_source_radius_squared_over_num\n"
      << "output_shield_on=" << cfg.out_on << '\n'
      << "output_shield_off=" << cfg.out_off << '\n'
      << "output_through_only=" << cfg.out_through << '\n'
      << "elapsed_seconds=" << elapsed << '\n';
    summary.close();
    if (!summary) {
        std::cerr << "Failed to write summary\n";
        return 6;
    }
    std::cerr << "Completed in " << elapsed << " seconds\n";
    return 0;
}
