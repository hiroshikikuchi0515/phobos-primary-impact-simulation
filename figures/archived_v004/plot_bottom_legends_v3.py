#!/usr/bin/env python3
"""Publication figure layout v3: legends below all axes. Numeric data unchanged."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'; INPUT=ROOT/'input'; FIG=ROOT/'figures'
FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.labelsize':10.5,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':250})

def save(fig,stem):
    for ext in ('png','pdf'):
        fig.savefig(FIG/(stem+'.'+ext),bbox_inches='tight',facecolor='white')
    plt.close(fig)

def panels(axs):
    for i,ax in enumerate(axs):
        ax.text(-.08,1.04,'('+chr(97+i)+')',transform=ax.transAxes,fontsize=11,
                fontweight='bold',ha='left',va='bottom',clip_on=False)
        ax.grid(alpha=.16)

def legend_bottom(fig,handles,labels,ncol=2,fontsize=9,title=None):
    fig.legend(handles,labels,ncol=ncol,loc='lower center',bbox_to_anchor=(.5,.012),
               frameon=False,fontsize=fontsize,title=title,handlelength=2.1,columnspacing=1.7)

# 1. 36 velocity/orbit models, direct shape production vs raw impact density
sector=pd.read_csv(INPUT/'sector_summary_v1.csv')
speeds=[.5,1,2,3,4,7.95298,20,40,60]
colors=plt.cm.viridis(np.linspace(.04,.91,len(speeds)))
fig,axs=plt.subplots(1,2,figsize=(12.0,5.15),sharex=True)
handles=[];labels=[]
for v,c in zip(speeds,colors):
    z=sector[np.isclose(sector.v_inf_km_s,v)].sort_values('orbit_km')
    label='7.953' if abs(v-7.95298)<.001 else f'{v:g}'
    hh=axs[0].errorbar(z.orbit_km/1000,z.L_T_crater_proxy_area_normalized,yerr=z.L_T_crater_proxy_SE,
                       color=c,marker='o',markersize=4.4,lw=1.6,capsize=1.6,label=label)
    axs[1].errorbar(z.orbit_km/1000,z.L_T_impact_flux_area_normalized,yerr=z.L_T_impact_flux_SE,
                    color=c,marker='o',markersize=4.4,lw=1.6,capsize=1.6)
    handles.append(hh);labels.append(label)
for ax in axs:
    ax.set_xlabel('Phobos orbital radius (10³ km)')
    ax.set_xticks([9.376,12,15,20],['9.376','12','15','20'])
    ax.axhline(1,color='.6',lw=.9,ls=':')
axs[0].set_ylabel('Leading / trailing crater-production proxy')
axs[1].set_ylabel('Leading / trailing impact-flux density')
panels(axs)
fig.subplots_adjust(left=.085,right=.985,top=.91,bottom=.31,wspace=.19)
legend_bottom(fig,handles,labels,ncol=5,fontsize=9,title=r'$v_\infty$ (km s$^{-1}$)')
save(fig,'Fig_orbit_LT_36cases_v004_pg1_v3')

# 2. D >= 1km vs orbit
obs=pd.read_csv(DATA/'observed_cumulative_diameter_LT_v004_area_v2.csv')
fig,ax=plt.subplots(figsize=(8.0,5.0))
for speed,color,marker in [(7.95298,'#24598b','o'),(40,'#ac6935','s')]:
    z=sector[np.isclose(sector.v_inf_km_s,speed)].sort_values('orbit_km')
    ax.errorbar(z.orbit_km/1000,z.L_T_crater_proxy_area_normalized,yerr=z.L_T_crater_proxy_SE,
                color=color,marker=marker,lw=1.9,markersize=5,capsize=1.7,
                label=rf'Model: $v_\infty$ = {speed:g} km s$^{{-1}}$')
for cat,col in [('Salamunićcar','#238c6d'),('Karachevtseva','#aa5c92')]:
    r=obs[(obs.catalog==cat)&(obs['mask']=='Global')&(obs.diameter_threshold=='D>=1.0')].iloc[0]
    ax.axhspan(r.CI95_low,r.CI95_high,alpha=.12,color=col)
    ax.axhline(r.area_normalized_LT,color=col,lw=1.5,ls='--',label=f'{cat}: D ≥ 1 km')
ax.set(xlabel='Phobos orbital radius (10³ km)',ylabel='Leading / trailing density ratio',xlim=(9.05,20.5),ylim=(.42,2.58))
ax.set_xticks([9.376,12,15,20],['9.376','12','15','20'])
ax.grid(alpha=.15)
h,l=ax.get_legend_handles_labels()
fig.subplots_adjust(left=.14,right=.98,top=.97,bottom=.32)
legend_bottom(fig,h,l,ncol=2,fontsize=9)
save(fig,'Fig_large_craters_vs_orbit_global_v004_v3')

# 3. diameter group obs L/T for global and outside Stickney
obs_dis=pd.read_csv(DATA/'observed_disjoint_diameter_LT_v004_area_v2.csv')
fig,axs=plt.subplots(1,2,figsize=(11.1,4.85),sharey=True)
for i,cat in enumerate(['Salamunićcar','Karachevtseva']):
    ax=axs[i]
    for j,(mask,marker,col) in enumerate([('Global','o','#28608d'),('Outside Stickney25','s','#af5f25')]):
        d=obs_dis[(obs_dis.catalog==cat)&(obs_dis['mask']==mask)].copy()
        d['x']=d.diameter_interval.map({'0.3–0.6':0,'0.6–1.0':1,'≥1.0':2})
        d=d.sort_values('x')
        xx=d.x.to_numpy()+(-.08 if j==0 else .08)
        ax.errorbar(xx,d.area_normalized_LT,yerr=[d.area_normalized_LT-d.CI95_low,d.CI95_high-d.area_normalized_LT],
                    marker=marker,color=col,ls='none',capsize=3,markersize=6,elinewidth=1.3,label=mask)
    ax.set_xticks(range(3),['0.3–0.6','0.6–1.0','≥1.0'])
    ax.set_xlim(-.32,2.32);ax.set_ylim(.37,2.90)
    ax.axhline(1,color='.5',ls=':',lw=1)
    ax.set_xlabel('Crater diameter D (km)')
axs[0].set_ylabel('Leading / trailing crater density')
panels(axs)
fig.subplots_adjust(left=.095,right=.985,top=.91,bottom=.26,wspace=.12)
legend_bottom(fig,*axs[0].get_legend_handles_labels(),ncol=2)
save(fig,'Fig_obs_diameter_and_Stickney25_v004_v3')

# 4. 4 sectors unmasked model
sec4=pd.read_csv(DATA/'model_v004_pg1_36cases_four_sector_global_v2.csv')
ref=sec4[np.isclose(sec4.v_inf_km_s,7.95298)]
fig,axs=plt.subplots(1,2,figsize=(11.3,4.85))
cc={'LM':'#33679e','LA':'#0e8c7b','TM':'#a66734','TA':'#993c72'}
for s in ('LM','LA','TM','TA'):
    z=ref[ref.sector==s].sort_values('orbit_km')
    axs[0].plot(z.orbit_km/1000,z.crater_proxy_share,marker='o',lw=1.8,color=cc[s],label=s)
    axs[1].plot(z.orbit_km/1000,z.model_relative_density,marker='o',lw=1.8,color=cc[s])
for ax in axs:
    ax.set_xticks([9.376,12,15,20],['9.376','12','15','20'])
    ax.set_xlabel('Phobos orbital radius (10³ km)')
axs[0].set_ylabel('Fraction of global crater-production proxy')
axs[1].set_ylabel('Local crater-production density / global mean')
panels(axs)
fig.subplots_adjust(left=.092,right=.985,top=.91,bottom=.26,wspace=.23)
legend_bottom(fig,*axs[0].get_legend_handles_labels(),ncol=4)
save(fig,'Fig_model_v004_four_sectors_7p953_v3')

# 5. Observed Stickney fractions and area-uniform reference
st=pd.read_csv(DATA/'observed_Stickney25_disjoint_counts_v2.csv')
fig,axs=plt.subplots(1,2,figsize=(11.1,4.8),sharey=True)
for i,cat in enumerate(['Salamunićcar','Karachevtseva']):
    d=st[st.catalog==cat].copy()
    d['idx']=d.diameter_interval.map({'0.3–0.6':0,'0.6–1.0':1,'≥1.0':2})
    d=d.sort_values('idx')
    x=d.idx.to_numpy();p=d.fraction_inside.to_numpy();n=d.global_N.to_numpy();k=d.inside_Stickney25.to_numpy()
    lo=beta.ppf(.025,np.maximum(k,1),n-k+1);lo[k==0]=0
    hi=beta.ppf(.975,k+1,np.maximum(n-k,1));hi[k==n]=1
    ax=axs[i]
    ax.errorbar(x,p,yerr=[p-lo,hi-p],marker='o',capsize=3,ls='none',color='#315f8b',
                label='Observed fraction')
    ax.axhline(d.uniform_area_fraction.iloc[0],ls='--',lw=1.5,color='#b36a35',label='Uniform physical-area fraction')
    ax.set_xticks([0,1,2],['0.3–0.6','0.6–1.0','≥1.0'])
    ax.set_xlim(-.3,2.3);ax.set_ylim(0,.15)
    ax.set_xlabel('Crater diameter D (km)')
axs[0].set_ylabel('Fraction of catalog craters inside Stickney25')
panels(axs)
fig.subplots_adjust(left=.11,right=.985,top=.91,bottom=.28,wspace=.17)
legend_bottom(fig,*axs[0].get_legend_handles_labels(),ncol=2)
save(fig,'Fig_Stickney25_observed_fraction_v004_v3')
print('Generated 5 pairs PNG/PDF (v3), bottom legends outside plot areas.')
