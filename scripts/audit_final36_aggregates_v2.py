#!/usr/bin/env python3
"""Check FINAL 36-case exported summary consistency; no raw event files needed.
Requires pandas; does not prove physical correctness or full reproduction.
"""
from __future__ import annotations
from pathlib import Path
import argparse, hashlib, math, sys
import pandas as pd

EXPECTED_ORBITS=(9376,12000,15000,20000)
EXPECTED_SPEEDS=(0.5,1,2,3,4,7.95298,20,40,60)
EXPECTED_SHAPE_SHA='1ee0d9e331c165ac890e7823fd8fa80f019e386017dc6d4286993268518ddc16'
SECTORS=('LM','LA','TM','TA')

def k(orbit,speed):
    # Decimal serialization of exactly the same physical v∞ differs across CSVs.
    return int(round(float(orbit))),round(float(speed),5)

def frame_lookup(df):
    dd={}
    for _,row in df.iterrows():
        key=k(row.orbit_km,row.v_inf_km_s)
        if key in dd: raise AssertionError(f'Duplicate case: {key}')
        dd[key]=row
    return dd

def check(root):
    base=root/'data'/'derived'
    files={
        'manifest':base/'manifest_36cases_portable_v2.csv',
        'sector_summary':base/'sector_summary_36cases_v2.csv',
        'model_global':base/'model_global_four_sector_36cases_v2.csv',
        'model_288':base/'model_288_reference.csv',
    }
    for name,path in files.items():
        if not path.is_file():raise FileNotFoundError(f'{name}: {path}')
    m=pd.read_csv(files['manifest']); ss=pd.read_csv(files['sector_summary'])
    fg=pd.read_csv(files['model_global']); f288=pd.read_csv(files['model_288'])
    assert (len(m),len(ss),len(fg),len(f288))==(36,36,144,288),'Unexpected row counts'
    valid={(i,v) for i in EXPECTED_ORBITS for v in EXPECTED_SPEEDS}
    mk=frame_lookup(m);sk=frame_lookup(ss)
    assert set(mk)==set(sk)==valid,'Case matrix incomplete or duplicated'
    assert m.status.eq('complete').all(),'Incomplete cases'
    assert (m.num==10_000_000).all() and (m.source_radius_km==20).all(),'Wrong proposal count or source radius'
    assert (m.threads==8).all(),'Unexpected thread count'
    assert m.seed.nunique()==36,'Seeds reused across cases'
    assert m.model_file_sha256.eq(EXPECTED_SHAPE_SHA).all(),'Wrong shape hash'
    for col in ['model_file','on_file','off_file','through_file','summary_file']:
        assert all(Path(str(x)).name==str(x) for x in m[col]), f'Absolute/local path detected in {col}'
    for z in ('mask','sector'):
        assert z in f288, z
    assert set(f288['mask'])=={'Global','Outside Stickney25'},'Unexpected masking'
    assert set(fg.sector)==set(f288.sector)==set(SECTORS),'Unexpected sectors'
    a=f288[f288['mask']=='Global'].copy()
    assert len(a)==144
    sdict={(k(r.orbit_km,r.v_inf_km_s),str(r.sector)):r for _,r in a.iterrows()}
    gdict={(k(r.orbit_km,r.v_inf_km_s),str(r.sector)):r for _,r in fg.iterrows()}
    assert set(sdict)==set(gdict)=={(kk,s) for kk in valid for s in SECTORS},'Global sector keys differ'
    for key in sdict:
        x=sdict[key];y=gdict[key]
        assert int(x.hits)==int(y.raw_hits),f'hit count mismatch {key}'
        assert math.isclose(float(x.production_fraction),float(y.crater_proxy_share),rel_tol=1e-11,abs_tol=1e-11),f'fraction mismatch {key}'
        assert math.isclose(float(x.area_km2),float(y.area_km2),rel_tol=1e-11,abs_tol=1e-11),f'area mismatch {key}'
    for kk in valid:
        xr={s:sdict[(kk,s)] for s in SECTORS}
        local=sk[kk]
        if int(round(float(local.hit_count)))!=sum(int(x.hits) for x in xr.values()):raise AssertionError(f'Hit sum mismatch: {kk}')
        def ratio(s1,s2):
            v1=sum(float(xr[s].crater_proxy_rate_per_area)*float(xr[s].area_km2) for s in s1)/sum(float(xr[s].area_km2) for s in s1)
            v2=sum(float(xr[s].crater_proxy_rate_per_area)*float(xr[s].area_km2) for s in s2)/sum(float(xr[s].area_km2) for s in s2)
            return v1/v2
        lt=ratio(('LM','LA'),('TM','TA'))
        ma=ratio(('LM','TM'),('LA','TA'))
        assert math.isclose(lt,float(local.L_T_crater_proxy_area_normalized),rel_tol=1e-10),f'L/T mismatch {kk}'
        assert math.isclose(ma,float(local.M_A_crater_proxy_area_normalized),rel_tol=1e-10),f'M/A mismatch {kk}'
        assert math.isclose(sum(float(x.area_km2) for x in xr.values()),1632.15946129769,abs_tol=1e-5)
        assert math.isclose(float(local.N),10_000_000)
        assert int(local.seed)==int(mk[kk].seed),f'Seed mismatch {kk}'
        assert math.isclose(float(local.Rsrc_km),20.,abs_tol=1e-9)
        assert not pd.isna(local.shield_percent)
    print('PASS: 36/36 production cases, 144 global rows, 288 masked rows')
    print('PASS: unique cases/seeds, v004 shape SHA, Rsrc=20 km, N=10M, threads=8')
    print('PASS: exact global sector count/fraction consistency; derived L/T and M/A values')
    print('PASS: portable manifest contains no absolute paths')
    print('CAUTION: No raw-event or shape-model reproduction is tested here.')

if __name__=='__main__':
    pa=argparse.ArgumentParser(description=__doc__)
    pa.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[1])
    args=pa.parse_args()
    check(args.repo_root.expanduser().resolve())
