#!/usr/bin/env python3
"""Independent reproduction of 144-cell rotation robustness from rotation_gof_full_v25.csv."""
import argparse
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

ap=argparse.ArgumentParser()
ap.add_argument('--input',type=Path,required=True)
ap.add_argument('--outdir',type=Path,default=Path('rotation_check_v30'))
a=ap.parse_args();a.outdir.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(a.input)
key=['catalog','mask','orbit_km','v_inf_km_s']
assert len(df)==2736 and len(df.groupby(key))==144
assert len(df[['axis','angle_deg']].drop_duplicates())==19
baseline=df[(df.axis=='none')&(df.angle_deg==0)][key+['non_rejected_nominal']].rename(columns={'non_rejected_nominal':'baseline_nonreject'})
rows=[]
for identifiers,g in df.groupby(key,sort=True):
    base=g[(g.axis=='none')&(g.angle_deg==0)].iloc[0]
    arr=g.non_rejected_nominal.astype(bool).to_numpy()
    state='Always non-rejected' if arr.all() else ('Always rejected' if not arr.any() else 'Threshold-changing')
    rows.append(dict(zip(key,identifiers))|dict(state=state,p_min=g.p_nominal.min(),p_max=g.p_nominal.max(),n_flips=int((g.non_rejected_nominal.astype(bool)!=bool(base.non_rejected_nominal)).sum())))
res=pd.DataFrame(rows)
res.to_csv(a.outdir/'rotation_144case_robustness_v30.csv',index=False)
print(res.groupby(['catalog','mask','state']).size().to_string())
fig,axs=plt.subplots(2,2,figsize=(10.2,9.1),sharex=True,sharey=True)
fig.subplots_adjust(left=.09,right=.985,top=.955,bottom=.16,wspace=.18,hspace=.28)
vs=sorted(res.v_inf_km_s.unique()); aa=sorted(res.orbit_km.unique());color=['#c9d4de','#e5a761','#337f73']
mapping={'Always rejected':0,'Threshold-changing':1,'Always non-rejected':2}
for ax,(cat,mask,lbl) in zip(axs.flat,[('S','global','(a) Salamunićcar — global'),('S','masked','(b) Salamunićcar — outside Stickney25'),('K','global','(c) Karachevtseva — global'),('K','masked','(d) Karachevtseva — outside Stickney25')]):
    s=res[(res.catalog==cat)&(res['mask']==mask)]
    mat=np.array([[mapping[s[(s.orbit_km==o)&np.isclose(s.v_inf_km_s,v)].state.iloc[0]] for o in aa] for v in vs])
    ax.imshow(mat,origin='lower',interpolation='nearest',aspect='auto',vmin=-.5,vmax=2.5,cmap=ListedColormap(color))
    ax.set_title(lbl,loc='left',fontsize=10,pad=8);ax.set_xticks(range(4),['9,376','12,000','15,000','20,000']);ax.set_yticks(range(9),[('%g'%v) for v in vs]);ax.set_xlabel('Orbital radius (km)');ax.set_ylabel(r'$v_\infty$ (km s$^{-1}$)')
    ax.set_xticks(np.arange(-.5,4,1),minor=True);ax.set_yticks(np.arange(-.5,9,1),minor=True);ax.grid(which='minor',color='white',lw=1.6);ax.tick_params(which='minor',length=0)
    for yy in range(9):
        for xx in range(4):
            if mat[yy,xx]==1:ax.text(xx,yy,'△',ha='center',va='center',fontsize=13,color='#382b20')
fig.legend(handles=[Patch(facecolor=color[0],label='Always rejected'),Patch(facecolor=color[1],label='Threshold-changing'),Patch(facecolor=color[2],label='Always non-rejected')],loc='lower center',bbox_to_anchor=(.5,.02),ncol=3,frameon=False)
fig.savefig(a.outdir/'rotation_robustness_v30.png',dpi=260,bbox_inches='tight');plt.close(fig)
