#!/usr/bin/env python3
"""Direct v004-PG1 Stickney25 masked four-sector crater-production summary from existing 26-col ON catalogs.
Unlike the full observation script, this works without access to observed crater coordinates.
"""
from pathlib import Path
import argparse,hashlib,math,json
import pandas as pd
from rebuild_exact_four_sector_from_raw_v2 import src_regions,parse_summary,accumulate_model,SECTORS

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--model-root',required=True)
 p.add_argument('--shape',required=True)
 p.add_argument('--out',required=True)
 p.add_argument('--max-cases',type=int,default=0)
 p.add_argument('--chunksize',type=int,default=100000)
 p.add_argument('--q',type=float,default=2.5)
 p.add_argument('--mu',type=float,default=.55)
 p.add_argument('--gamma',type=float,default=.55)
 p.add_argument('--vref',type=float,default=8000)
 a=p.parse_args()
 root=Path(a.model_root); out=Path(a.out)
 if out.exists():raise FileExistsError('Refusing to overwrite output: '+str(out))
 areas=src_regions(a.shape)
 manifest=pd.read_csv(root/'manifest_v1.csv')
 if len(manifest)!=36 or not manifest.status.eq('complete').all():raise ValueError('Expected 36 complete runs')
 hash=hashlib.sha256(Path(a.shape).read_bytes()).hexdigest()
 if not manifest.model_file_sha256.eq(hash).all():raise ValueError('Model geometry hash mismatch')
 rows=[]
 for i,m in manifest.iterrows():
  if a.max_cases and i>=a.max_cases:break
  on=Path(m.on_file)
  if not on.exists():on=root/on.name
  if not on.exists():raise FileNotFoundError(str(on))
  ss=Path(m.summary_file)
  if not ss.exists():ss=root/ss.name
  s=parse_summary(ss)
  if s.get('phobos_gravity')!='1' or s.get('mars_gravity')!='1' or s.get('shape')!='obj':
   raise ValueError('Unexpected gravity/shape in '+str(ss))
  r=accumulate_model(on,areas,s,int(m.num),float(m.source_radius_km),a.q,a.mu,a.gamma,a.vref,a.chunksize)
  assert sum(int(z['hits']) for z in r if z['mask']=='Global')==int(s['hits_shield_on'])
  for x in r:x.update(orbit_km=m.orbit_km,v_inf_km_s=m.v_inf_km_s,seed=m.seed)
  rows.extend(r)
  print('DONE',m.identifier,flush=True)
 out.parent.mkdir(parents=True,exist_ok=True)
 pd.DataFrame(rows).to_csv(out,index=False)
 print('WROTE',out,'rows',len(rows),flush=True)

if __name__=='__main__':main()
