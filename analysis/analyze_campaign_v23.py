#!/usr/bin/env python3
"""Compare externally measured MC run results, source-radius convergence,
and an independently computed local Liouville spherical reference.

Input metrics must be from ACTUAL C++ runs, not synthetic or reference data.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

MEAS=['flux_off_rate','flux_on_rate','crater_off_rate','crater_on_rate','lt_crater','ma_crater']
BASE_KEYS=['shape','phobos_gravity','v_inf_km_s','r_src_km']

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--plan',type=Path,required=True)
 p.add_argument('--metrics',type=Path,required=True)
 p.add_argument('--theory',type=Path,required=True)
 p.add_argument('--out',type=Path,required=True)
 p.add_argument('--allow-partial',action='store_true',help='Analyze only rows with complete measured metrics')
 args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
 plan=pd.read_csv(args.plan)
 metrics=pd.read_csv(args.metrics)
 if metrics.run_id.duplicated().any():raise SystemExit('Duplicate run_id in metrics')
 data=plan.merge(metrics,on='run_id',how='inner',validate='one_to_one')
 if data.empty:raise SystemExit('No matching run_id between plan/metrics')
 for col in MEAS:
  if col not in data.columns:raise SystemExit(f'Missing measured column: {col}')
  data[col]=pd.to_numeric(data[col],errors='coerce')
 if data[MEAS].isna().any().any():
  miss=data.loc[data[MEAS].isna().any(axis=1),'run_id'].tolist()
  if not args.allow_partial:
   raise SystemExit('ABORT: missing measured rates or L/T; do not use partly filled table: '+str(miss[:8])+'; add --allow-partial to skip')
  print(f'NOTICE: skipping {len(miss)} incomplete measurement rows')
  data=data.dropna(subset=MEAS).copy()
  if data.empty:raise SystemExit('No complete measured rows available')
 if (data[['flux_off_rate','crater_off_rate','flux_on_rate','crater_on_rate']]<=0).any().any():
  raise SystemExit('Nonpositive flux/crater rates')
 if (data.flux_on_rate > data.flux_off_rate*(1+1e-12)).any():raise SystemExit('Flux ON > OFF')
 if (data.crater_on_rate > data.crater_off_rate*(1+1e-12)).any():raise SystemExit('Crater ON > OFF')
 data['flux_shield_pct']=100*(1-data.flux_on_rate/data.flux_off_rate)
 data['crater_shield_pct']=100*(1-data.crater_on_rate/data.crater_off_rate)
 data.to_csv(args.out/'verified_measured_runs_v23.csv',index=False)
 group=data.groupby(BASE_KEYS,sort=True,dropna=False)
 agg=group[['lt_crater','ma_crater','flux_shield_pct','crater_shield_pct']].agg(['mean','std','count'])
 agg.columns=['_'.join(x) for x in agg.columns]
 agg.reset_index().to_csv(args.out/'run_group_summary_v23.csv',index=False)
 # Matched-seed radius comparisons, both metrics and speeds
 paired=[]
 for (shape,pg,v),sub in data.groupby(['shape','phobos_gravity','v_inf_km_s']):
  if str(shape).lower()!='obj' or int(pg)!=1:continue
  for radius in [15,30]:
   small=sub[sub.r_src_km==radius].set_index('seed')
   base=sub[sub.r_src_km==20].set_index('seed')
   shared=small.index.intersection(base.index)
   if len(shared)==0:continue
   for metric in ['lt_crater','ma_crater','flux_shield_pct','crater_shield_pct']:
    delta=small.loc[shared,metric]-base.loc[shared,metric]
    relative=100*delta/base.loc[shared,metric].values
    paired.append({'shape':shape,'pg':pg,'v_inf_km_s':v,'comparison':f'R{radius}-R20',
                   'metric':metric,'n_matched_seeds':len(shared),
                   'mean_absolute_difference':delta.mean(),
                   'se_absolute_difference':delta.std(ddof=1)/np.sqrt(len(shared)) if len(shared)>1 else np.nan,
                   'mean_relative_difference_pct':relative.mean(),
                   'min_absolute_difference':delta.min(),'max_absolute_difference':delta.max()})
 if paired:pd.DataFrame(paired).to_csv(args.out/'paired_radius_differences_v23.csv',index=False)
 # Independent sphere benchmark vs locally integrated theory, do not label this exact finite-R solution
 theory=pd.read_csv(args.theory)
 bench=[]
 for (shape,pg,v,R),sub in group:
  if str(shape).lower()!='sphere' or int(pg)!=0 or float(R)!=20:continue
  q=theory[(np.isclose(theory.orbit_km,9376))&(np.isclose(theory.v_inf_km_s,v))]
  if len(q)!=1:raise SystemExit(f'Missing/ambiguous theoretical reference vinf={v}')
  q=q.iloc[0]
  for mc_key,thkey in [('lt_crater','crater_LT'),('ma_crater','crater_MA'),
                       ('flux_shield_pct','flux_shield_pct'),('crater_shield_pct','crater_shield_pct')]:
   x=sub[mc_key].to_numpy();theory_mean=float(q[thkey+'_mean'])
   mcmean=x.mean();semc=x.std(ddof=1)/np.sqrt(len(x)) if len(x)>1 else np.nan
   setheory=float(q[thkey+'_se'])
   bench.append({'v_inf_km_s':v,'r_src_km':R,'n_seeds':len(x),'metric':mc_key,
                 'measured_MC_mean':mcmean,'MC_seed_SE':semc,
                 'independent_local_theory':theory_mean,'theory_QMC_SE':setheory,
                 'MC_minus_theory':mcmean-theory_mean,
                 'relative_residual_pct':100*(mcmean-theory_mean)/theory_mean,
                 'note':'Local Mars-only analytical sphere limit, not exact finite-20-km LWSS benchmark'})
 if bench:pd.DataFrame(bench).to_csv(args.out/'sphere_theory_comparison_v23.csv',index=False)
 print('Checked',len(data),'supplied complete measurement rows (provenance must be independently verified).')
 print('Saved comparison outputs to',args.out)

if __name__=='__main__':main()
