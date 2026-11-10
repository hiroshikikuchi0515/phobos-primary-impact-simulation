#!/usr/bin/env python3
"""Phobos observational-frame rotation sensitivity for fixed v004 source fields.

Only rotate catalog crater centers; leave the model field, quadrant planes and
Stickney25 reference cap fixed in the adopted shape-body frame. The axes are
+X toward Mars, -Y leading, +Z north. This is a sensitivity stress test, NOT an
established PCK10/PCK11 frame transform. Requires original crater coordinates.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2

EXPECTED={
    'S':{'global':np.array([107,150,127,126]),'masked':np.array([93,150,127,126])},
    'K':{'global':np.array([168,216,177,164]),'masked':np.array([151,216,177,164])},
}
CATN={'S':'Salamuniccar','K':'Karachevtseva'}
QUADS=['LM','LA','TM','TA']
ANGLES=[-5,-2,-1,0,1,2,5]
CANDIDATE_COLUMNS={
    'lon':['lon_e_deg','lon_east_deg','longitude_east','longitude_e','lon','longitude','Longitude','lon_deg'],
    'lat':['lat_deg','latitude','Latitude','lat'],
    'diam':['diameter_km','diam_km','diameter','diam','D_km','D','diameter_m'],
}

def colpick(df,specified,name):
    if specified:
        if specified not in df.columns:raise ValueError(f'No column {specified!r} in input columns: {list(df.columns)}')
        return specified
    possible=[x for x in CANDIDATE_COLUMNS[name] if x in df.columns]
    if len(possible)==1:return possible[0]
    raise ValueError(f'Please specify --*-{name}: found {possible}; columns={list(df.columns)}')


def load_catalog(path,loncol,latcol,dcol,diam_unit='km',lon_positive='east'):
    df=pd.read_csv(path)
    cols=[colpick(df,loncol,'lon'),colpick(df,latcol,'lat'),colpick(df,dcol,'diam')]
    vals=df[cols].apply(pd.to_numeric,errors='coerce')
    if vals.isna().any().any():raise ValueError(f'Non-numeric coordinate or diameter in {path}')
    lon=np.asarray(vals.iloc[:,0],float)
    lat=np.asarray(vals.iloc[:,1],float)
    diam=np.asarray(vals.iloc[:,2],float)
    if lon_positive=='west':lon=-lon
    if diam_unit=='m':diam/=1000
    if np.any(abs(lat)>90):raise ValueError('Latitude outside [-90,90] deg')
    if np.any(diam<=0):raise ValueError('Nonpositive diameters')
    keep=diam>=.3
    lon=np.deg2rad(lon[keep]);lat=np.deg2rad(lat[keep])
    xyz=np.column_stack((np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)))
    return xyz,int(len(df)),int(keep.sum()),cols


def rotate(points,axis,ang_deg):
    if ang_deg==0:return points.copy()
    a=np.deg2rad(ang_deg);c=np.cos(a);s=np.sin(a)
    if axis=='X':R=np.array([[1,0,0],[0,c,-s],[0,s,c]])
    elif axis=='Y':R=np.array([[c,0,s],[0,1,0],[-s,0,c]])
    elif axis=='Z':R=np.array([[c,-s,0],[s,c,0],[0,0,1]])
    else:raise ValueError(axis)
    return points@R.T


def categorize(v):
    # The model definition puts LM: (+X,-Y) and TA: (-X,+Y)
    sec=np.select([
        (v[:,0]>0)&(v[:,1]<0),
        (v[:,0]<=0)&(v[:,1]<0),
        (v[:,0]>0)&(v[:,1]>=0),
    ],[0,1,2],default=3)
    st_center=np.deg2rad(-50.0);st_lat=np.deg2rad(1.3)
    center=np.array([np.cos(st_lat)*np.cos(st_center),np.cos(st_lat)*np.sin(st_center),np.sin(st_lat)])
    in_st=(v@center)>=np.cos(np.deg2rad(25))
    # 25deg Stickney cap must be entirely in LM for this adopted geometry
    if np.any(in_st & (sec!=0)):
        raise ValueError('Stickney25 center or coordinate axes inconsistent with LM convention')
    return sec,in_st


def g_stat(obs,p):
    if not np.isclose(sum(p),1):raise ValueError('bad model fractions')
    if np.any(np.asarray(p)<=0):raise ValueError('zero/negative model probabilities')
    O=np.array(obs,dtype=float);N=O.sum()
    if N<=0:return float('nan'),float('nan')
    g=2*np.sum(np.where(O>0,O*np.log(np.maximum(O,1e-300)/(N*p)),0))
    return float(g),float(chi2.sf(g,3))


def model_from_288(path):
    df=pd.read_csv(path)
    req={'orbit_km','v_inf_km_s','mask','sector','sum_crater_w'}
    if not req.issubset(df.columns):raise ValueError(f'Model CSV missing columns {req-set(df.columns)}')
    result={}
    for (a,v,mask),q in df.groupby(['orbit_km','v_inf_km_s','mask']):
        arr=q.set_index('sector').loc[QUADS,'sum_crater_w'].to_numpy(dtype=float)
        result[(float(a),float(v),'global' if mask=='Global' else 'masked')]=arr/arr.sum()
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--salam',type=Path,required=True,help='Individual Salamuniccar crater records CSV')
    p.add_argument('--kara',type=Path,required=True,help='Individual Karachevtseva crater records CSV')
    p.add_argument('--model',type=Path,required=True,help='288-row v004 production CSV')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--s-lon');p.add_argument('--s-lat');p.add_argument('--s-diam')
    p.add_argument('--k-lon');p.add_argument('--k-lat');p.add_argument('--k-diam')
    p.add_argument('--s-lon-positive',choices=['east','west'],default='east')
    p.add_argument('--k-lon-positive',choices=['east','west'],default='east')
    p.add_argument('--s-diam-unit',choices=['km','m'],default='km')
    p.add_argument('--k-diam-unit',choices=['km','m'],default='km')
    p.add_argument('--allow-baseline-mismatch',action='store_true')
    p.add_argument('--radii',nargs='+',type=float,default=[9376.,20000.])
    p.add_argument('--speeds',nargs='+',type=float,default=[7.95298,40.])
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    cats={};reports=[]
    for key,path,L,B,D,LP,DU in [
        ('S',args.salam,args.s_lon,args.s_lat,args.s_diam,args.s_lon_positive,args.s_diam_unit),
        ('K',args.kara,args.k_lon,args.k_lat,args.k_diam,args.k_lon_positive,args.k_diam_unit)]:
        xyz,nraw,nsel,cols=load_catalog(path,L,B,D,DU,LP)
        sec,st=categorize(xyz)
        global_counts=np.bincount(sec,minlength=4)
        masked_counts=np.bincount(sec[~st],minlength=4)
        baseline={'catalog':key,'input':str(path),'records_all':nraw,'records_Dge0p3':nsel,
                  'columns':str(cols),'lon_positive':LP,'diam_unit':DU,
                  'N_stickney25':int(sum(st)),'global_counts':str(global_counts.tolist()),
                  'masked_counts':str(masked_counts.tolist()),
                  'baseline_exact_match':bool(np.array_equal(global_counts,EXPECTED[key]['global']) and
                                            np.array_equal(masked_counts,EXPECTED[key]['masked']))}
        reports.append(baseline);cats[key]=(xyz,sec,st)
    pd.DataFrame(reports).to_csv(args.out/'catalog_baseline_audit_v23.csv',index=False)
    for r in reports:
        print(f"{r['catalog']} D>=0.3 N={r['records_Dge0p3']} global {r['global_counts']}; "
              f"masked {r['masked_counts']}, St25={r['N_stickney25']}, "
              f"baseline={'PASS' if r['baseline_exact_match'] else 'FAIL'}")
    if not all(x['baseline_exact_match'] for x in reports) and not args.allow_baseline_mismatch:
        raise SystemExit('ABORT: baseline does not reproduce published v22 quadrant and Stickney25 counts. '
                         'Fix coordinate conventions/CSV selections; --allow-baseline-mismatch for diagnostics only.')
    mod=model_from_288(args.model)
    outputs=[]
    for k,(xyz,sec_base,st_base) in cats.items():
        for axis,angle in [('none',0)]+[(ax,an) for ax in 'XYZ' for an in ANGLES if an!=0]:
            v=rotate(xyz,axis,angle) if axis!='none' else xyz
            sec,st=categorize(v)
            count_global=np.bincount(sec,minlength=4)
            count_masked=np.bincount(sec[~st],minlength=4)
            for mask,obs in [('global',count_global),('masked',count_masked)]:
                for a in args.radii:
                    for speed in args.speeds:
                        mk=(a,speed,mask)
                        if mk not in mod:raise ValueError(f'No model: {mk}')
                        g,pvalue=g_stat(obs,mod[mk])
                        outputs.append({
                            'catalog':k,'axis':axis,'angle_deg':angle,'mask':mask,
                            'orbit_km':a,'v_inf_km_s':speed,
                            'N':int(sum(obs)),'N_stickney25':int(sum(st)),
                            'num_quadrant_switches':int(np.sum(sec!=sec_base)),
                            'num_Stickney_switches':int(np.sum(st!=st_base)),
                            **{f'obs_{quad}':int(obs[i]) for i,quad in enumerate(QUADS)},
                            'G':g,'p_nominal':pvalue,'non_rejected_nominal':bool(pvalue>=.05)
                        })
    tbl=pd.DataFrame(outputs)
    tbl.to_csv(args.out/'rotation_gof_full_v23.csv',index=False)
    baseline=tbl[tbl.axis=='none'][['catalog','mask','orbit_km','v_inf_km_s','G']]
    baseline=baseline.rename(columns={'G':'G_baseline'})
    tbl=tbl.merge(baseline,on=['catalog','mask','orbit_km','v_inf_km_s'],validate='many_to_one')
    tbl['delta_G']=tbl.G-tbl.G_baseline
    tbl.to_csv(args.out/'rotation_gof_delta_v23.csv',index=False)
    by=tbl[(tbl['mask']=='masked')&(tbl.orbit_km==9376)]
    print('Range of masked-p at current orbit (includes ±1, ±2, ±5 deg about X,Y,Z):')
    for (k,speed),df in by.groupby(['catalog','v_inf_km_s']):
        print(f' {k} vinf={speed:.5g} km/s:  p_min={df.p_nominal.min():.6g}, '
              f'p_max={df.p_nominal.max():.6g}, max abs deltaG={df.delta_G.abs().max():.3f}')
    print('Saved',args.out)

if __name__=='__main__': main()
