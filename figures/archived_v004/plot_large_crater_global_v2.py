#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1];data=root/'data';out=root/'figures'
model=pd.read_csv(root/'input/sector_summary_v1.csv')
obs=pd.read_csv(data/'observed_cumulative_diameter_LT_v004_area_v2.csv')
plt.rcParams.update({'font.size':10,'axes.labelsize':11,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':8,'pdf.fonttype':42,'savefig.dpi':240})
fig,ax=plt.subplots(figsize=(7.6,4.8))
for speed,color,mark in [(7.95298,'#24598b','o'),(40,'#ac6935','s')]:
 t=model[np.isclose(model.v_inf_km_s,speed)].sort_values('orbit_km')
 ax.errorbar(t.orbit_km/1000,t.L_T_crater_proxy_area_normalized,yerr=t.L_T_crater_proxy_SE,
   color=color,marker=mark,lw=1.9,markersize=5,capsize=1.7,label=f'LWSS-MC: $v_\\infty$ = {speed:g} km s$^{{-1}}$')
for cat,col in [('Salamunićcar','#238c6d'),('Karachevtseva','#aa5c92')]:
 row=obs[(obs.catalog==cat)&(obs['mask']=='Global')&(obs.diameter_threshold=='D>=1.0')].iloc[0]
 ax.axhspan(row.CI95_low,row.CI95_high,alpha=.12,color=col)
 ax.axhline(row.area_normalized_LT,color=col,lw=1.4,ls='--',label=f'{cat}: D ≥ 1 km')
ax.set(xlabel='Phobos orbital radius (10³ km)',ylabel='Leading / trailing density ratio',xlim=(9.05,20.5),ylim=(.42,2.58))
ax.set_xticks([9.376,12,15,20],['9.376','12','15','20'])
ax.grid(alpha=.12)
fig.legend(*ax.get_legend_handles_labels(),ncol=2,loc='upper center',bbox_to_anchor=(.5,1.03),frameon=False)
fig.subplots_adjust(top=.74,bottom=.18,left=.13,right=.98)
fig.savefig(out/'Fig_large_craters_vs_orbit_global_v004_v2.png',bbox_inches='tight')
fig.savefig(out/'Fig_large_craters_vs_orbit_global_v004_v2.pdf',bbox_inches='tight')
