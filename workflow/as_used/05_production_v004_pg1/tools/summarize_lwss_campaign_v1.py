#!/usr/bin/env python3
"""Streaming four-sector flux and crater-proxy summary of completed LWSS-MC runs."""
from __future__ import annotations
from pathlib import Path
import argparse, math, csv, re, sys
import numpy as np

BASE=Path(__file__).resolve().parents[1]
LABELS=('LM','LA','TM','TA')

def areas_from_obj(path):
    vs=[];fs=[]
    with open(path) as f:
        for line in f:
            if line.startswith('v '): vs.append([float(x) for x in line.split()[1:4]])
            elif line.startswith('f '): fs.append([int(x.split('/')[0])-1 for x in line.split()[1:4]])
    v=np.asarray(vs,dtype=np.float64)
    f=np.asarray(fs,dtype=np.int64)
    tri=v[f]
    cent=tri.mean(axis=1)
    a=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    lm=(cent[:,1]<0)&(cent[:,0]>0)
    la=(cent[:,1]<0)&(cent[:,0]<=0)
    tm=(cent[:,1]>=0)&(cent[:,0]>0)
    ta=(cent[:,1]>=0)&(cent[:,0]<=0)
    areas={k:float(a[m].sum()) for k,m in zip(LABELS,(lm,la,tm,ta))}
    print('v004 physical areas [km²]:',areas,'global=',sum(areas.values()),flush=True)
    return areas

def parse_s(path):
    return dict(line.split('=',1) for line in path.read_text().splitlines() if '=' in line)

def ratio_mom(sums,squares,first,second,afirst,asecond):
    A=sum(sums[z] for z in first);B=sum(sums[z] for z in second)
    A2=sum(squares[z] for z in first);B2=sum(squares[z] for z in second)
    if A<=0 or B<=0: return (float('nan'),float('nan'))
    ratio=(A/afirst)/(B/asecond)
    se=ratio*math.sqrt(A2/A**2+B2/B**2)
    return (ratio,se)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--results-root',default='results/lwss_campaign_v1/v004pg1_n10000000_R20_dt1_post450')
    p.add_argument('--shape',default='shape_model/phobos_g_148m_spc_obj_0000n00000_v004.obj')
    p.add_argument('--out',default='plots/lwss_campaign_v1/sector_summary_v1.csv')
    p.add_argument('--q',type=float,default=2.5)
    p.add_argument('--mu',type=float,default=.55)
    p.add_argument('--gamma',type=float,default=.55)
    p.add_argument('--vref-mps',type=float,default=8000.)
    p.add_argument('--chunksize',type=int,default=100000)
    args=p.parse_args()
    try: import pandas as pd
    except ImportError:
        print('Requires pandas: install pandas in your mmx conda environment.',file=sys.stderr);raise
    if args.vref_mps<=0 or args.chunksize<1: p.error('invalid value')
    shape=(BASE/args.shape).resolve() if not Path(args.shape).is_absolute() else Path(args.shape).resolve()
    areas=areas_from_obj(shape)
    root=Path(args.results_root)
    if not root.is_absolute():root=BASE/root
    output=Path(args.out)
    if not output.is_absolute():output=BASE/output
    files=sorted(root.rglob('*.summary.txt'))
    if not files: raise FileNotFoundError(f'No complete summary files in {root}')
    rows=[]
    for summary in files:
        s=parse_s(summary)
        if s.get('model')!='source_sphere_production_v5' or s.get('phobos_gravity')!='1' or s.get('shape')!='obj':
            print('SKIP different condition:',summary);continue
        prefix=summary.name.removesuffix('.summary.txt')
        on=summary.with_name(prefix+'_on.txt')
        if not on.exists():
            print('SKIP no ON catalog:',on);continue
        sums={k:0. for k in LABELS};squares={k:0. for k in LABELS}
        fsums={k:0. for k in LABELS};fsquares={k:0. for k in LABELS}
        counts={k:0 for k in LABELS}
        for block in pd.read_csv(on,sep=r'\s+',comment='#',header=None,chunksize=args.chunksize):
            d=block.to_numpy(dtype=np.float64)
            if d.shape[1]!=26: raise ValueError(f'Expected 26 columns in {on}; got {d.shape[1]}')
            latlon=d[:,1]*math.pi/180.
            is_leading=np.sin(latlon)<0
            is_mars=np.cos(latlon)>0
            weights=d[:,4]
            crater_w=weights*(d[:,11]/args.vref_mps)**(args.mu*args.q)*np.maximum(0,d[:,19])**(args.gamma*args.q)
            sects=(is_leading & is_mars, is_leading & ~is_mars,~is_leading & is_mars,~is_leading & ~is_mars)
            for k,m in zip(LABELS,sects):
                x=crater_w[m];y=weights[m]
                sums[k]+=float(x.sum());squares[k]+=float((x*x).sum());fsums[k]+=float(y.sum());fsquares[k]+=float((y*y).sum());counts[k]+=int(m.sum())
        l_area=areas['LM']+areas['LA'];t_area=areas['TM']+areas['TA']
        m_area=areas['LM']+areas['TM'];a_area=areas['LA']+areas['TA']
        proxy_lt,proxy_lt_se=ratio_mom(sums,squares,['LM','LA'],['TM','TA'],l_area,t_area)
        proxy_ma,proxy_ma_se=ratio_mom(sums,squares,['LM','TM'],['LA','TA'],m_area,a_area)
        flux_lt,flux_lt_se=ratio_mom(fsums,fsquares,['LM','LA'],['TM','TA'],l_area,t_area)
        N=float(s['num']);Asrc=float(s['source_sphere_area_m2'])
        row={'orbit_km':s['orbit_r_km'],'v_inf_km_s':float(s['v_inf_mps'])/1000,'N':int(N),'seed':s['seed'],'Rsrc_km':s['source_radius_km'],'shield_percent':s['shield_fraction_percent'],'L_T_crater_proxy_area_normalized':proxy_lt,'L_T_crater_proxy_SE':proxy_lt_se,'M_A_crater_proxy_area_normalized':proxy_ma,'M_A_crater_proxy_SE':proxy_ma_se,'L_T_impact_flux_area_normalized':flux_lt,'L_T_impact_flux_SE':flux_lt_se,'hit_count':sum(counts.values()),'model_surface_area_km2':sum(areas.values())}
        for k in LABELS:
            row[f'{k}_area_km2']=areas[k]
            row[f'{k}_hits']=counts[k]
            row[f'{k}_crater_proxy_rate_per_area']=Asrc/N*sums[k]/(areas[k]*1e6)
        rows.append(row)
        print('DONE',prefix,'LT',f'{proxy_lt:.5f} ± {proxy_lt_se:.5f}',flush=True)
    if not rows:raise RuntimeError('No suitable v5 PG1 catalog in '+str(root))
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print('Written:',output,'cases:',len(rows))
if __name__=='__main__':main()
