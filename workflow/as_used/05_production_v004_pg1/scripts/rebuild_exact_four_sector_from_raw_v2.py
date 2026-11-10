#!/usr/bin/env python3
"""v004 + Phobos point-gravity ON: matched model/observation masks and diameter bins.

Run in user's WSL after supplying an OBSERVATION catalog normalized to:
 catalog,lon_deg,lat_deg,diameter_km
where lon_deg is body-fixed EAST-positive signed, leading ~ -90, trailing ~ +90.

Unlike the aggregate-only outputs, this command uses actual 26-column ON collision
catalogs and individual observed crater center coordinates for exact quadrant/mask tests.

This does not assume the two overlapping crater catalogs are independent.
"""
from __future__ import annotations
import argparse, csv, math, json, sys, hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2, beta

SECTORS=('LM','LA','TM','TA')
DIAM_BINS=(('0.3–<0.6',.3,.6),('0.6–<1.0',.6,1.),('≥1.0',1.,np.inf),('≥0.3',.3,np.inf),('≥0.6',.6,np.inf))
S_LON,S_LAT=-50.,1.3
STICK_R=25.


def parse_summary(path):
 return dict(line.split('=',1) for line in Path(path).read_text().splitlines() if '=' in line)


def src_regions(obj):
 vs=[]; fs=[]
 with Path(obj).open() as f:
  for line in f:
   if line.startswith('v '):vs.append([float(x) for x in line.split()[1:4]])
   elif line.startswith('f '):fs.append([int(x.split('/')[0])-1 for x in line.split()[1:4]])
 v=np.array(vs,dtype=np.float64)
 f=np.array(fs,dtype=np.int64)
 t=v[f]
 c=t.mean(axis=1)
 area=.5*np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)
 lon=np.arctan2(c[:,1],c[:,0]);lat=np.arctan2(c[:,2],np.hypot(c[:,0],c[:,1]))
 sec=sectors(lon)
 st=inside_stickney(lon,lat)
 areas={}
 for mask,m in [('Global',np.ones_like(st,bool)),('Outside Stickney25',~st),('Stickney25',st)]:
  areas[mask]={s:float(area[(sec==s)&m].sum()) for s in SECTORS}
 if sum(areas['Stickney25'][s] for s in ('LA','TM','TA'))>1e-8:raise ValueError('Stickney25 no longer wholly LM: audit coordinate frame')
 if not np.isclose(sum(areas['Global'].values()),1632.1594613,atol=.02):raise ValueError('Model area inconsistent with v004 148m reference')
 return areas


def sectors(lon_rad):
 lead=np.sin(lon_rad)<0; mars=np.cos(lon_rad)>0
 return np.select([lead & mars,lead & ~mars,~lead & mars],['LM','LA','TM'],default='TA')


def inside_stickney(lon,lat):
 lon0=math.radians(S_LON);lat0=math.radians(S_LAT)
 return np.sin(lat)*math.sin(lat0)+np.cos(lat)*math.cos(lat0)*np.cos(lon-lon0)>=math.cos(math.radians(STICK_R))


def parse_observations(path):
 p=Path(path)
 if p.suffix.lower() in ('.xlsx','.xls'):
  df=pd.read_excel(p)
 else:
  df=pd.read_csv(p)
 needed={'catalog','lon_deg','lat_deg','diameter_km'}
 if not needed.issubset(df.columns):
  raise ValueError(f'Expected columns {sorted(needed)}; got {df.columns.tolist()}. Export from the original catalogs to a canonical coordinate CSV first.')
 df=df.copy()
 df['catalog']=df.catalog.astype(str).str.strip().replace({'Salamuniccar':'Salamunićcar','Salamunićcar':'Salamunićcar','K':'Karachevtseva','S':'Salamunićcar'})
 df=df[df.catalog.isin(['Salamunićcar','Karachevtseva'])].copy()
 for col in ('lon_deg','lat_deg','diameter_km'):
  df[col]=pd.to_numeric(df[col],errors='raise')
 if not np.isfinite(df[['lon_deg','lat_deg','diameter_km']].to_numpy()).all():raise ValueError('Nonfinite observation coordinate/diameter')
 if df.lat_deg.abs().max()>90.0001:raise ValueError('Latitude outside [-90,90]')
 if df.diameter_km.min()<=0:raise ValueError('Nonpositive crater diameter; units incorrect?')
 # DO NOT reverse a west-positive source here. Require already harmonized signed EAST.
 df['lon_deg']=(df.lon_deg+180)%360-180
 df['sector']=sectors(np.deg2rad(df.lon_deg.to_numpy()))
 df['inside_Stickney25']=inside_stickney(np.deg2rad(df.lon_deg.to_numpy()),np.deg2rad(df.lat_deg.to_numpy()))
 counts=[]
 for cat in ('Salamunićcar','Karachevtseva'):
  g=df[df.catalog==cat]
  for dname,dlo,dhi in DIAM_BINS[:3]:
   t=g[(g.diameter_km>=dlo)&(g.diameter_km<dhi)]
   counts.append((cat,dname,len(t),int(t.inside_Stickney25.sum())))
 expected={('Salamunićcar','0.3–<0.6'):(378,10),('Salamunićcar','0.6–<1.0'):(76,2),('Salamunićcar','≥1.0'):(56,2),
           ('Karachevtseva','0.3–<0.6'):(552,12),('Karachevtseva','0.6–<1.0'):(101,3),('Karachevtseva','≥1.0'):(72,2)}
 bad=[(cat,di,n,i,expected[(cat,di)]) for cat,di,n,i in counts if (n,i)!=expected[(cat,di)]]
 if bad:
  raise ValueError('Raw observation normalization does not reproduce the archived exact global/Stickney25 counts. '
                   'STOP: check longitude convention, diameter conversion, crater selection, catalog edition, and mask center. '+str(bad))
 return df


def obs_tables(obs,areas):
 counts=[]
 for cat in ('Salamunićcar','Karachevtseva'):
  d=obs[obs.catalog==cat]
  for mask,keep in [('Global',np.ones(len(d),bool)),('Outside Stickney25',~d.inside_Stickney25.to_numpy())]:
   sub=d.iloc[np.flatnonzero(keep)]
   for dbin,lo,hi in DIAM_BINS:
    z=sub[(sub.diameter_km>=lo)&(sub.diameter_km<hi)]
    record={'catalog':cat,'mask':mask,'diameter_bin':dbin,'N':len(z)}
    for sec in SECTORS:record[sec]=int((z.sector==sec).sum())
    for area_label,ss in [('leading',('LM','LA')),('trailing',('TM','TA')),('mars',('LM','TM')),('anti',('LA','TA'))]:
     record[area_label]=sum(record[s] for s in ss)
    aa=areas[mask]
    def ratio_num(n,d,an,ad):return (n/an)/(d/ad) if d else math.nan
    record['LT_density']=ratio_num(record['leading'],record['trailing'],aa['LM']+aa['LA'],aa['TM']+aa['TA'])
    record['MA_density']=ratio_num(record['mars'],record['anti'],aa['LM']+aa['TM'],aa['LA']+aa['TA'])
    counts.append(record)
 return pd.DataFrame(counts)


def accumulate_model(on,areas,stats,N,Rsrc_km,q,mu,gamma,vref,chunksize):
 # Accumulate observed-source-flux moment sums directly; no sector factorization.
 sums={m:{s:0.0 for s in SECTORS} for m in ('Global','Outside Stickney25')}
 w2={m:{s:0.0 for s in SECTORS} for m in sums}
 hits={m:{s:0 for s in SECTORS} for m in sums}
 for block in pd.read_csv(on,sep=r'\s+',comment='#',header=None,chunksize=chunksize):
  d=block.to_numpy(dtype=float)
  if d.shape[1]!=26:raise ValueError(f'{on}: invalid 26-column catalog; got {d.shape[1]}')
  if not np.isfinite(d).all():raise ValueError(f'{on}: non-finite values')
  lon=np.deg2rad(d[:,1]);lat=np.deg2rad(d[:,0]);sector=sectors(lon)
  inside=inside_stickney(lon,lat)
  if np.any(d[:,4]<0) or np.any(d[:,11]<=0) or np.any(d[:,19]<-1e-7) or np.any(d[:,19]>1.000001):
   raise ValueError(f'{on}: non-physical flux/speed/incidence weights')
  w=d[:,4]*(d[:,11]/vref)**(mu*q)*np.clip(d[:,19],0,1)**(gamma*q)
  for mask,keep in [('Global',np.ones(len(d),bool)),('Outside Stickney25',~inside)]:
   for s in SECTORS:
    selector=keep&(sector==s)
    x=w[selector]
    sums[mask][s]+=float(np.sum(x))
    w2[mask][s]+=float(np.sum(x*x))
    hits[mask][s]+=int(np.sum(selector))
 scale=4*math.pi*(Rsrc_km*1000)**2/N
 rows=[]
 for mask in sums:
  tot=sum(sums[mask].values());area=areas[mask]
  if tot<=0:raise ValueError(f'No collisions for {mask}')
  for s in SECTORS:
   rows.append(dict(mask=mask,sector=s,hits=hits[mask][s],sum_crater_w=sums[mask][s],sum_crater_w2=w2[mask][s],
     area_km2=area[s],crater_proxy_rate_per_area=scale*sums[mask][s]/(area[s]*1e6),
     production_fraction=sums[mask][s]/tot))
 return rows


def compare_model_obs(md,ob):
 arr=np.array([float(md.loc[s,'production_fraction']) for s in SECTORS])
 obs=np.array([int(ob[s]) for s in SECTORS])
 n=obs.sum();expected=n*arr
 if np.any(expected<=0):raise ValueError('Zero model expected count')
 pearson=float(np.sum((obs-expected)**2/expected))
 g=float(2*np.sum(np.where(obs>0,obs*np.log(np.maximum(obs,1e-300)/expected),0)))
 return dict(pearson_chisq=pearson,pearson_p_df3=float(chi2.sf(pearson,3)),G_deviance=g,G_p_df3=float(chi2.sf(g,3)),min_expected=float(expected.min()))


def main():
 pa=argparse.ArgumentParser(description=__doc__)
 pa.add_argument('--observations',required=True,help='Canonical CSV with catalog,lon_deg,lat_deg,diameter_km (east-positive signed)')
 pa.add_argument('--model-root',required=True,help='Directory holding manifest_v1.csv and ON catalogs')
 pa.add_argument('--shape',required=True,help='Full path to phobos_g_148m_spc_obj_0000n00000_v004.obj')
 pa.add_argument('--out',required=True,help='New output directory; this script refuses to overwrite existing outputs')
 pa.add_argument('--q',type=float,default=2.5)
 pa.add_argument('--mu',type=float,default=.55)
 pa.add_argument('--gamma',type=float,default=.55)
 pa.add_argument('--vref',type=float,default=8000)
 pa.add_argument('--chunksize',type=int,default=100000)
 pa.add_argument('--max-cases',type=int,default=0,help='For smoke test only; 0 means all 36')
 a=pa.parse_args()
 if not a.q>0 or not a.mu>=0 or not a.gamma>=0 or not a.vref>0:pa.error('Invalid scaling parameters')
 out=Path(a.out).expanduser()
 for n in ['v004_area.csv','observed_four_sectors.csv','model_four_sectors_masked.csv','fixed_model_four_sector_tests.csv','run_info.json']:
  if (out/n).exists():raise FileExistsError(f'Will not overwrite: {out/n}')
 areas=src_regions(a.shape)
 print('Geometry:',areas,flush=True)
 obs=parse_observations(a.observations)
 print('Historical crater total and mask counts verified.',flush=True)
 observed=obs_tables(obs,areas)
 root=Path(a.model_root)
 manifest=pd.read_csv(root/'manifest_v1.csv')
 if len(manifest)!=36:raise ValueError('Expected 36 production cases')
 if not manifest.status.eq('complete').all():raise ValueError('Incomplete manifest')
 shape_hash=hashlib.sha256(Path(a.shape).read_bytes()).hexdigest()
 if not (manifest.model_file_sha256==shape_hash).all():raise ValueError('Shape hash mismatch')
 out.mkdir(parents=True,exist_ok=True)
 pd.DataFrame([dict(mask=m,sector=s,area_km2=areas[m][s]) for m in areas for s in SECTORS]).to_csv(out/'v004_area.csv',index=False)
 observed.to_csv(out/'observed_four_sectors.csv',index=False)
 allrows=[];tests=[]
 for i,entry in manifest.iterrows():
  if a.max_cases>0 and i>=a.max_cases:break
  path=Path(entry.on_file)
  if not path.exists():path=root/path.name
  if not path.exists():raise FileNotFoundError(f'Missing ON catalog for {entry.identifier}: {path}')
  sp=Path(entry.summary_file)
  if not sp.exists():sp=root/sp.name
  summary=parse_summary(sp)
  if summary.get('phobos_gravity')!='1' or summary.get('shape')!='obj' or summary.get('mars_gravity')!='1':
   raise ValueError(f'{entry.identifier}: changed physical parameters')
  sect_rows=accumulate_model(path,areas,summary,int(entry.num),float(entry.source_radius_km),a.q,a.mu,a.gamma,a.vref,a.chunksize)
  for r in sect_rows:r.update(orbit_km=entry.orbit_km,v_inf_km_s=entry.v_inf_km_s,seed=entry.seed,source_file=str(path))
  allrows.extend(sect_rows)
  frame=pd.DataFrame(sect_rows)
  for mask in ('Global','Outside Stickney25'):
   model=frame[frame['mask']==mask].set_index('sector')
   if mask=='Global':
    ref=Path('plots/lwss_campaign_v1/sector_summary_v1.csv')
    # compare internally with summary rates from same run if available in user's original folder
    if int(model.hits.sum())!=int(summary['hits_shield_on']):raise ValueError('Global hits do not match summary')
   for _,r in observed[observed['mask']==mask].iterrows():
    test=compare_model_obs(model,r)
    test.update(orbit_km=entry.orbit_km,v_inf_km_s=entry.v_inf_km_s,seed=entry.seed,
        mask=mask,catalog=r.catalog,diameter_bin=r.diameter_bin,N=int(r.N),
        model_LM_fraction=float(model.loc['LM','production_fraction']),model_LA_fraction=float(model.loc['LA','production_fraction']),
        model_TM_fraction=float(model.loc['TM','production_fraction']),model_TA_fraction=float(model.loc['TA','production_fraction']))
    tests.append(test)
  print('DONE',entry.identifier,flush=True)
 pd.DataFrame(allrows).to_csv(out/'model_four_sectors_masked.csv',index=False)
 pd.DataFrame(tests).to_csv(out/'fixed_model_four_sector_tests.csv',index=False)
 (out/'run_info.json').write_text(json.dumps(dict(shape=str(a.shape),shape_sha256=shape_hash,observations=str(a.observations),
    model_root=str(root),completed_cases=min(len(manifest),a.max_cases) if a.max_cases else len(manifest),
    q=a.q,mu=a.mu,gamma=a.gamma,vref=a.vref,mask_center_lon=S_LON,mask_center_lat=S_LAT,mask_radius_deg=STICK_R,
    caveat='Nominal 3-df tests: no crater detection/secondary correction, model Monte Carlo uncertainty neglected, catalog overlap; model remains isotropic and stationary at infinity.'),indent=2,ensure_ascii=False))
 print('WROTE',out,flush=True)

if __name__=='__main__':main()
