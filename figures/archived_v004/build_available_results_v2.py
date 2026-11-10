#!/usr/bin/env python3
"""Reconstruct exactly supported v004 area/marginal-observation products.
Does NOT invent joint quadrant counts or Stickney-masked impact production from aggregate data.
"""
from pathlib import Path
import sys, math
import numpy as np
import pandas as pd
from scipy.stats import beta, binomtest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
MESH=ROOT/'input/shape_v004_148m.obj'
SOURCE=ROOT/'input/sector_summary_v1.csv'
MARG=ROOT/'input/Phobos_crater_diameter_LT_bins_v1.csv'
OUT=ROOT/'data';FIG=ROOT/'figures'
for p in (OUT,FIG):p.mkdir(parents=True,exist_ok=True)

def mesh_regions(path):
 v=[];f=[]
 with path.open() as h:
  for line in h:
   if line.startswith('v '):v.append([float(x) for x in line.split()[1:4]])
   if line.startswith('f '):f.append([int(x.split('/')[0])-1 for x in line.split()[1:4]])
 v=np.asarray(v);f=np.asarray(f,dtype=int)
 tri=v[f]
 cent=tri.mean(axis=1)
 area=0.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
 lon=np.arctan2(cent[:,1],cent[:,0])
 lat=np.arctan2(cent[:,2],np.hypot(cent[:,0],cent[:,1]))
 sector=np.select([(np.sin(lon)<0)&(np.cos(lon)>0), (np.sin(lon)<0)&(np.cos(lon)<=0), (np.sin(lon)>=0)&(np.cos(lon)>0)],['LM','LA','TM'],default='TA')
 slon,slat=np.deg2rad(-50),np.deg2rad(1.3)
 stick=(np.sin(lat)*np.sin(slat)+np.cos(lat)*np.cos(slat)*np.cos(lon-slon)>=np.cos(np.deg2rad(25)))
 d={}
 for mask,k in [(np.ones_like(stick,bool),'Global'),(~stick,'Outside Stickney25'),(stick,'Stickney25')]:
  d[k]={s:float(area[(sector==s)&mask].sum()) for s in ('LM','LA','TM','TA')}
 return d

regions=mesh_regions(MESH)
area_rows=[]
for m,d in regions.items():
 for s,a in d.items():area_rows.append(dict(mask=m,sector=s,area_km2=a))
area=pd.DataFrame(area_rows)
area.to_csv(OUT/'v004_area_Stickney25_four_sectors_v2.csv',index=False)
full=regions['Global'];outside=regions['Outside Stickney25'];
assert all(regions['Stickney25'][s]==0 for s in ('LA','TM','TA'))
assert abs(sum(full.values())-1632.1594612976864)<1e-8

marg=pd.read_csv(MARG)
# input N and counts are independently archival and were transcribed from published-product analysis
assert set(marg['mask'])=={'Global','Outside Stickney25'}
assert all(marg.N==marg.leading_count+marg.trailing_count)
assert sum(marg[marg['mask']=='Global'].N)==1235

rows=[]
for _,r in marg.iterrows():
 aa=regions[r['mask']]
 A_L=aa['LM']+aa['LA']; A_T=aa['TM']+aa['TA']
 k,n=int(r.leading_count),int(r.N)
 l=n-k
 lo=0 if k==0 else beta.ppf(.025,k,n-k+1)
 hi=1 if k==n else beta.ppf(.975,k+1,n-k)
 fn=lambda p: (p/(1-p))*(A_T/A_L) if p<1 else float('inf')
 rows.append(dict(catalog=r.catalog,mask=r['mask'],diameter_interval=r.diameter_bin_km,
   leading=k,trailing=l,total=n,A_leading_km2=A_L,A_trailing_km2=A_T,
   area_normalized_LT=(k/l)*(A_T/A_L) if l else float('nan'),
   CI95_low=fn(lo),CI95_high=fn(hi),
   source='archived disjoint diameter marginal counts; NOT individual-catalog remeasurement'))
df=pd.DataFrame(rows)
df.to_csv(OUT/'observed_disjoint_diameter_LT_v004_area_v2.csv',index=False)

allrows=[]
for catalog in df.catalog.unique():
 for mask in ['Global','Outside Stickney25']:
  d=df[(df.catalog==catalog)&(df['mask']==mask)]
  aa=regions[mask];AL=aa['LM']+aa['LA'];AT=aa['TM']+aa['TA']
  for s,sel in [('D>=0.3',d),('D>=0.6',d[d.diameter_interval!='0.3–0.6']),('D>=1.0',d[d.diameter_interval=='≥1.0'])]:
   nL=int(sel.leading.sum());nT=int(sel.trailing.sum());n=nL+nT
   lo=beta.ppf(.025,nL,nT+1) if nL else 0
   hi=beta.ppf(.975,nL+1,nT) if nT else 1
   conv=lambda p:p/(1-p)*AT/AL if p<1 else float('inf')
   allrows.append(dict(catalog=catalog,mask=mask,diameter_threshold=s,leading=nL,trailing=nT,N=n,
       area_normalized_LT=nL/nT*AT/AL,CI95_low=conv(lo),CI95_high=conv(hi)))
cum=pd.DataFrame(allrows);cum.to_csv(OUT/'observed_cumulative_diameter_LT_v004_area_v2.csv',index=False)

counts=[]
for cat in marg.catalog.unique():
 g=marg[(marg.catalog==cat)&(marg['mask']=='Global')].set_index('diameter_bin_km')
 o=marg[(marg.catalog==cat)&(marg['mask']=='Outside Stickney25')].set_index('diameter_bin_km')
 for d in g.index:
  z=g.loc[d];zz=o.loc[d]
  inside=int(z.N-zz.N)
  assert inside==int(z.leading_count-zz.leading_count)
  assert z.trailing_count==zz.trailing_count
  counts.append(dict(catalog=cat,diameter_interval=d,global_N=int(z.N),inside_Stickney25=inside,
     outside_Stickney25=int(zz.N),fraction_inside=inside/int(z.N),
     uniform_area_fraction=sum(regions['Stickney25'].values())/sum(full.values())))
st=pd.DataFrame(counts)
st.to_csv(OUT/'observed_Stickney25_disjoint_counts_v2.csv',index=False)

sec=pd.read_csv(SOURCE)
assert len(sec)==36 and sec.N.eq(10000000).all()
secs=('LM','LA','TM','TA')
new=[]
for _,r in sec.iterrows():
 d={s:float(r[f'{s}_crater_proxy_rate_per_area'])*float(r[f'{s}_area_km2']) for s in secs}
 tot=sum(d.values())
 af=sum(full.values())
 for s in secs:
  assert abs(float(r[f'{s}_area_km2'])-full[s])<1e-6
  new.append(dict(orbit_km=int(r.orbit_km),v_inf_km_s=float(r.v_inf_km_s),
   sector=s,area_km2=full[s],crater_proxy_share=d[s]/tot,model_relative_density= (d[s]/full[s])/(tot/af),
   raw_hits=int(r[f'{s}_hits']),relative_production_rate_per_area=float(r[f'{s}_crater_proxy_rate_per_area'])))
new=pd.DataFrame(new)
new.to_csv(OUT/'model_v004_pg1_36cases_four_sector_global_v2.csv',index=False)

# Global (not masked!) binomial comparisons for each model/size/catalog. One-dimensional only.
tests=[]
for _,m in sec.iterrows():
 AL=full['LM']+full['LA'];AT=full['TM']+full['TA']
 prodR=float(m.L_T_crater_proxy_area_normalized)
 p=prodR*AL/(prodR*AL+AT)
 for _,obs in cum[cum['mask']=='Global'].iterrows():
  k=int(obs.leading);n=int(obs.N)
  tests.append(dict(catalog=obs.catalog,diameter_threshold=obs.diameter_threshold,N=n,observed_leading=k,
   orbit_km=int(m.orbit_km),v_inf_km_s=float(m.v_inf_km_s),
   observed_LT_density=obs.area_normalized_LT,model_LT_density=prodR,
   expected_leading_fraction=p,exact_p_conditional_binomial=binomtest(k,n,p).pvalue,
   interpretation='global L/T marginal test ONLY; selected hypotheses, no multiplicity correction'))
pd.DataFrame(tests).to_csv(OUT/'global_model_observation_LT_binomial_36cases_v2.csv',index=False)

plt.rcParams.update({'font.size':10,'axes.labelsize':11,'xtick.labelsize':9,'ytick.labelsize':9,
                     'legend.fontsize':9,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':240})
col=['#28608d','#af5f25']
fig,axs=plt.subplots(1,2,figsize=(10.8,4.25),sharey=True)
for i,cat in enumerate(['Salamunićcar','Karachevtseva']):
 ax=axs[i]
 for j,(mask,marker) in enumerate([('Global','o'),('Outside Stickney25','s')]):
  d=df[(df.catalog==cat)&(df['mask']==mask)]
  bmap={'0.3–0.6':0,'0.6–1.0':1,'≥1.0':2}
  d=d.assign(x=d.diameter_interval.map(bmap)).sort_values('x')
  dx=(-.08 if j==0 else .08)
  ax.errorbar(d.x.to_numpy()+dx,d.area_normalized_LT,yerr=[d.area_normalized_LT-d.CI95_low,d.CI95_high-d.area_normalized_LT],
              color=col[j],marker=marker,ls='none',capsize=3,elinewidth=1.4,markersize=6,
              label='Global' if j==0 else 'Outside Stickney25')
 ax.set_xticks(range(3),['0.3–0.6','0.6–1.0','≥1.0'])
 ax.axhline(1,color='.5',ls=':',lw=1)
 ax.grid(alpha=.16,axis='y');ax.set_xlim(-.32,2.32);ax.set_ylim(.37,2.9)
 ax.set_xlabel('Crater diameter D (km)')
 ax.text(.01,1.03,f'({chr(97+i)})',transform=ax.transAxes,fontweight='bold')
axs[0].set_ylabel('Leading / trailing crater density')
fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1.02),ncol=2,frameon=False)
fig.subplots_adjust(top=.83,bottom=.2,left=.09,right=.98,wspace=.12)
fig.savefig(FIG/'Fig_obs_diameter_and_Stickney25_v004_v2.png',bbox_inches='tight')
fig.savefig(FIG/'Fig_obs_diameter_and_Stickney25_v004_v2.pdf',bbox_inches='tight');plt.close(fig)

# New 36-case PG1 OBJ direct four-sector production shares and spatial relative rates, 7.953 reference speed
sref=new[np.isclose(new.v_inf_km_s,7.95298)].copy()
fig,axs=plt.subplots(1,2,figsize=(10.8,4.15))
colors={'LM':'#33679e','LA':'#0e8c7b','TM':'#a66734','TA':'#993c72'}
for s in secs:
 z=sref[sref.sector==s].sort_values('orbit_km')
 axs[0].plot(z.orbit_km/1000,z.crater_proxy_share,marker='o',lw=1.75,color=colors[s],label=s)
 axs[1].plot(z.orbit_km/1000,z.model_relative_density,marker='o',lw=1.75,color=colors[s])
for i,a in enumerate(axs):
 a.set_xticks([9.376,12,15,20],['9.376','12','15','20']);a.set_xlabel('Phobos orbital radius (10³ km)')
 a.grid(alpha=.17);a.text(.01,1.03,f'({chr(97+i)})',transform=a.transAxes,fontweight='bold')
axs[0].set_ylabel('Fraction of global crater-production proxy')
axs[1].set_ylabel('Local crater-production density / global mean')
fig.legend(*axs[0].get_legend_handles_labels(),ncol=4,loc='upper center',bbox_to_anchor=(.5,1.02),frameon=False)
fig.subplots_adjust(top=.82,bottom=.2,left=.09,right=.98,wspace=.22)
fig.savefig(FIG/'Fig_model_v004_four_sectors_7p953_v2.png',bbox_inches='tight')
fig.savefig(FIG/'Fig_model_v004_four_sectors_7p953_v2.pdf',bbox_inches='tight');plt.close(fig)

# Descriptive Stickney25: observed centers inside vs physical-area uniform rate
fig,axs=plt.subplots(1,2,figsize=(10.8,4.0),sharey=True)
for i,cat in enumerate(['Salamunićcar','Karachevtseva']):
 d=st[st.catalog==cat].copy().reset_index(drop=True)
 d['idx']=d.diameter_interval.map({'0.3–0.6':0,'0.6–1.0':1,'≥1.0':2})
 d=d.sort_values('idx')
 x=d.idx.to_numpy()
 p=d.fraction_inside.to_numpy()
 n=d.global_N.to_numpy()
 k=d.inside_Stickney25.to_numpy()
 ci_lo=beta.ppf(.025,np.maximum(k,1),n-k+1);ci_lo[k==0]=0
 ci_hi=beta.ppf(.975,k+1,np.maximum(n-k,1));ci_hi[k==n]=1
 axs[i].errorbar(x,p,yerr=[p-ci_lo,ci_hi-p],marker='o',capsize=3,ls='none',label='Observed fraction',color='#315f8b')
 axs[i].axhline(d.uniform_area_fraction.iloc[0],ls='--',lw=1.5,color='#b36a35',label='Uniform physical-area fraction')
 axs[i].set_xticks([0,1,2],['0.3–0.6','0.6–1.0','≥1.0']);axs[i].set_xlim(-.3,2.3)
 axs[i].set_xlabel('Crater diameter D (km)');axs[i].grid(alpha=.14,axis='y')
 axs[i].text(.01,1.03,f'({chr(97+i)})',transform=axs[i].transAxes,fontweight='bold')
axs[0].set_ylabel('Fraction of catalog craters inside Stickney25')
axs[0].set_ylim(0,.15)
fig.legend(*axs[0].get_legend_handles_labels(),ncol=2,loc='upper center',bbox_to_anchor=(.5,1.02),frameon=False)
fig.subplots_adjust(top=.81,bottom=.2,left=.1,right=.98,wspace=.16)
fig.savefig(FIG/'Fig_Stickney25_observed_fraction_v004_v2.png',bbox_inches='tight')
fig.savefig(FIG/'Fig_Stickney25_observed_fraction_v004_v2.pdf',bbox_inches='tight');plt.close(fig)

print('v004 geometry',full,'Stickney',regions['Stickney25'])
print('S and K global >=1:')
print(cum[(cum['mask']=='Global')&(cum.diameter_threshold=='D>=1.0')].to_string(index=False))
print('model v7.953 current and 20k sectors:')
print(sref[sref.orbit_km.isin([9376,20000])][['orbit_km','sector','crater_proxy_share','model_relative_density']].to_string(index=False))
