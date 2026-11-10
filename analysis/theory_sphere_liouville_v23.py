#!/usr/bin/env python3
"""Independent *local* Mars-only Liouville reference for a spherical Phobos.

This is NOT a finite-source-sphere propagation code. It integrates a sphere's
projected collision cross section over Mars-centered asymptotic isotropic
velocities, a forward geometric Mars-shadow criterion, and a uniform impact
disk. It is valid for a fixed circular Phobos orbit, Phobos gravity OFF,
R_body=11.1 km, and the adopted separable crater-production weight.

Body axes: +X toward Mars, -Y leading, +Z north. See README for assumptions.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.stats import qmc

GM_MARS=4.282837e13       # m^3/s^2
R_MARS=3389500.0         # m, absorbing sphere
R_PH=11100.0             # m, analytical target sphere
V_REF=8000.0             # m/s
CRATER_EXP=1.375         # q * mu = q * gamma


def shadow_solid_angle_fraction(a, v_inf):
    vloc=np.hypot(v_inf,np.sqrt(2*GM_MARS/a))
    bcrit=R_MARS*np.sqrt(v_inf**2+2*GM_MARS/R_MARS)
    s=np.clip(bcrit/(a*vloc),0,1)
    return float((1-np.sqrt(1-s*s))/2)


def calculate(a_km, vinf_kms, m=19, seed=12345):
    a=a_km*1000.0; v_inf=vinf_kms*1000.0
    v_orb=np.sqrt(GM_MARS/a)
    v_loc=np.sqrt(v_inf**2+2*GM_MARS/a)
    samples=qmc.Sobol(d=4,scramble=True,seed=seed).random_base2(m)
    mu=2*samples[:,0]-1
    phi=2*np.pi*samples[:,1]
    s=np.sqrt(np.maximum(0,1-mu**2))
    # Mars-centered velocity, expressed in Phobos body axes
    projectile_dir=np.column_stack((mu,s*np.cos(phi),s*np.sin(phi)))
    # Phobos velocity points toward -Y in body-fixed axes
    rel_vel=v_loc*projectile_dir
    rel_vel[:,1]+=v_orb
    speed=np.linalg.norm(rel_vel,axis=1)
    uhat=rel_vel/speed[:,None]
    # Make an orthonormal frame in the plane perpendicular to relative velocity
    helper=np.zeros_like(uhat)
    helper[:,2]=1.0
    helper[np.abs(uhat[:,2])>0.9]=[1.0,0.0,0.0]
    e1=np.cross(uhat,helper)
    e1/=np.linalg.norm(e1,axis=1)[:,None]
    e2=np.cross(uhat,e1)
    rho=np.sqrt(samples[:,2])
    ang=2*np.pi*samples[:,3]
    cosi=np.sqrt(np.maximum(0,1-rho*rho))
    # Disk-sampled impact normal for a sphere (Lambertian surface incidence)
    n=-uhat*cosi[:,None]+rho[:,None]*(np.cos(ang)[:,None]*e1+np.sin(ang)[:,None]*e2)
    # prior Mars-crossing condition: outbound after hyperbolic periapsis
    h=a*v_loc*np.sqrt(np.maximum(0,1-mu*mu))
    ecc=np.sqrt(1+(h*v_inf/GM_MARS)**2)
    rp=h*h/(GM_MARS*(1+ecc))
    through=(mu<0)&(rp<R_MARS)

    # All angles have equal probability; constant v_loc/v_inf cancels ratios.
    norm=(np.pi*R_PH**2)*(v_loc/v_inf)
    flux=speed
    crater=speed*(speed/V_REF)**CRATER_EXP * cosi**CRATER_EXP
    keep=(~through).astype(float)
    result={
        'orbit_km':a_km,'v_inf_km_s':vinf_kms,'sobol_m':m,
        'sobol_seed':seed,'n_samples':int(len(samples)),
        'shadow_solid_angle_exact_pct':100*shadow_solid_angle_fraction(a,v_inf),
        'shadow_solid_angle_numerical_pct':100*np.mean(through),
        'v_loc_mps':v_loc,'v_orb_mps':v_orb,
    }
    for name,wt in [('flux',flux),('crater',crater)]:
        off_rate=norm*wt.mean()
        on_rate=norm*np.mean(wt*keep)
        result[name+'_off_rate_per_ninf_m3s']=float(off_rate)
        result[name+'_on_rate_per_ninf_m3s']=float(on_rate)
        result[name+'_shield_pct']=float(100*(1-on_rate/off_rate))
        wt_on=wt*keep
        lt=wt_on[n[:,1]<0].sum()/wt_on[n[:,1]>=0].sum()
        ma=wt_on[n[:,0]>0].sum()/wt_on[n[:,0]<=0].sum()
        result[name+'_LT']=float(lt)
        result[name+'_MA']=float(ma)
        for nameq,sel in [('LM',(n[:,1]<0)&(n[:,0]>0)),
                          ('LA',(n[:,1]<0)&(n[:,0]<=0)),
                          ('TM',(n[:,1]>=0)&(n[:,0]>0)),
                          ('TA',(n[:,1]>=0)&(n[:,0]<=0))]:
            result[name+'_'+nameq+'_on_rate_per_ninf_m3s']=float(norm*wt_on[sel].sum()/len(samples))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--radii',nargs='+',type=float,default=[9376])
    p.add_argument('--speeds',nargs='+',type=float,default=[7.95298,40])
    p.add_argument('--m',type=int,default=19,help='Sobol points = 2^m per scramble')
    p.add_argument('--replicates',type=int,default=4)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for orbit in a.radii:
        for v in a.speeds:
            for rep in range(a.replicates):
                r=calculate(orbit,v,a.m,10101+rep*4999)
                r['replicate']=rep+1
                rows.append(r)
    with (a.out/'theory_sphere_individual_v23.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    keys=[k for k in rows[0] if k not in ('orbit_km','v_inf_km_s','sobol_m','sobol_seed','n_samples')]
    summary=[]
    for orbit in a.radii:
        for v in a.speeds:
            rr=[x for x in rows if x['orbit_km']==orbit and x['v_inf_km_s']==v]
            out={'orbit_km':orbit,'v_inf_km_s':v,'sobol_m':a.m,'replicates':len(rr)}
            for k in keys:
                vals=np.array([r[k] for r in rr],dtype=float)
                out[k+'_mean']=np.mean(vals)
                out[k+'_se']=np.std(vals,ddof=1)/np.sqrt(len(vals)) if len(vals)>1 else float('nan')
            summary.append(out)
    with (a.out/'theory_sphere_summary_v23.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    for r in summary:
        print(f"a={r['orbit_km']:.0f} km, vinf={r['v_inf_km_s']:.5g} km/s  "
              f"shadow(angular exact)={r['shadow_solid_angle_exact_pct_mean']:.5f}%  "
              f"flux shield={r['flux_shield_pct_mean']:.4f}%  "
              f"crater L/T={r['crater_LT_mean']:.5f}  "
              f"crater M/A={r['crater_MA_mean']:.5f}")
    print('Independent local reference saved to',a.out)

if __name__=='__main__': main()
