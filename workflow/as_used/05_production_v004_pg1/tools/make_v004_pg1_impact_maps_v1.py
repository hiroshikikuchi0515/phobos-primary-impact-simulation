#!/usr/bin/env python3
"""Generate paper maps from VERIFIED PG1/v004 LWSS 26-column event catalogs.

By default reads only two existing 10M-production cases (9376 and 20000 km,
v_inf=7.95298 km/s). NO simulation rerun. Always uses true v004 triangle area.
Creates PNG/PDF figures and tiny 15-degree gridded CSVs (not 26-column copies).

Modes: ON, OFF and THROUGH are paired event catalogs; shielding map is
100 * THROUGH_crater_weight / OFF_crater_weight, with low-statistic bins masked.
Density maps are divided by the actual mesh area per longitude-latitude bin and
normalized to each case's GLOBAL SHIELDING-OFF mean, for matched visual scales.
This is crater-production *proxy*, NOT a measured crater-density map.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

FIELDS='lat_deg lon_deg beta_deg vrel_mps w_flux_mps source_radius_km u_x u_y u_z phase_rad alpha_deg v_surface_mps face_id v_surface_body_x_mps v_surface_body_y_mps v_surface_body_z_mps normal_body_x normal_body_y normal_body_z incidence_cos_raw source_normal_inertial_x source_normal_inertial_y source_normal_inertial_z source_mars_distance_km source_vloc_mps inward_flux_mps'.split()
EXPECTED_SHA='1ee0d9e331c165ac890e7823fd8fa80f019e386017dc6d4286993268518ddc16'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):h.update(part)
    return h.hexdigest()

def geometry(shape,lon_edges,lat_edges):
    vertices=[];faces=[]
    with shape.open('r') as f:
        for line in f:
            if line.startswith('v '):vertices.append(tuple(map(float,line.split()[1:4])))
            elif line.startswith('f '):faces.append([int(x.split('/')[0])-1 for x in line.split()[1:4]])
    v=np.asarray(vertices,dtype=np.float64);ind=np.asarray(faces,dtype=np.int64)
    if len(ind)<100000:raise ValueError('Expected full 148-m v004 OBJ mesh; wrong shape?')
    tris=v[ind];vec=np.cross(tris[:,1]-tris[:,0],tris[:,2]-tris[:,0]); areas=np.linalg.norm(vec,axis=1)/2
    center=tris.mean(axis=1)
    lat=np.degrees(np.arctan2(center[:,2],np.hypot(center[:,0],center[:,1])))
    lon=np.degrees(np.arctan2(center[:,1],center[:,0]))
    area2,_,_=np.histogram2d(lat,lon,bins=[lat_edges,lon_edges],weights=areas)
    total=area2.sum()
    if abs(total-1632.15946129769)>0.02:raise ValueError(f'Unexpected mesh area: {total:.6f} km2 (expected 1632.15946)')
    return area2,{'mesh_vertices':len(v),'mesh_faces':len(ind),'mesh_area_km2':total,
                  'area_method':'full v004 OBJ facet areas, each assigned to 15-degree bin by triangle centroid'}


def input_path(x,manifest):
    p=Path(str(x))
    if p.is_file():return p
    # Robust if archives moved but filenames retained
    p=manifest.parent/p.name
    if p.is_file():return p
    raise FileNotFoundError(f'Event catalog missing: {x} ; also tried {p}')


def read_catalog(path, lat_edges,lon_edges,scatter_max,rng,chunk):
    with path.open() as f:header=f.readline().strip().split()
    if not header or header[0]!='#' or header[1:]!=FIELDS:
        raise ValueError(f'Unexpected 26-column header in {path}. Refuse to guess columns.')
    shp=(len(lat_edges)-1,len(lon_edges)-1)
    acc={key:np.zeros(shp,dtype=np.float64) for key in ('flux','crater','flux2','crater2','count')}
    # Exact uniform-priority sample: the scatter plot shows at most scatter_max
    # uniformly selected *event rows*, not a rate-weighted sample.
    chosen_keys=np.empty((0,),dtype=np.float64);chosen_points=np.empty((0,2),dtype=np.float64)
    n=0
    for frame in pd.read_csv(path,sep=r'\s+',comment='#',header=None,chunksize=chunk):
        a=frame.to_numpy(dtype=np.float64,copy=False)
        if a.shape[1]!=26:raise ValueError(f'Expected 26 columns, got {a.shape[1]} in {path}')
        if not np.isfinite(a).all():raise ValueError(f'NaN/Inf found in {path}')
        lat=a[:,0];lon=a[:,1];w=a[:,4]
        if np.any((lat<-90)|(lat>90)|(lon<-180)|(lon>180)|(w<0)):raise ValueError(f'Invalid location/weight {path}')
        if np.any(a[:,11]<=0):raise ValueError(f'Invalid surface speed {path}')
        w_cr=w*np.power(a[:,11]/8000,1.375)*np.power(np.clip(a[:,19],0,1),1.375)
        n+=len(a)
        for name,weights in (('flux',w),('crater',w_cr),('flux2',w*w),('crater2',w_cr*w_cr),('count',None)):
            hist,_,_=np.histogram2d(lat,lon,bins=(lat_edges,lon_edges),weights=weights)
            acc[name]+=hist
        keys=rng.random(len(a))
        pts=a[:,[1,0]]
        if len(chosen_keys)==scatter_max:
            selector=keys<chosen_keys.max()
            if not selector.any():continue
            keys=keys[selector];pts=pts[selector]
        newkeys=np.concatenate((chosen_keys,keys));newpts=np.concatenate((chosen_points,pts),axis=0)
        if len(newkeys)>scatter_max:
            ii=np.argpartition(newkeys,scatter_max-1)[:scatter_max]
            chosen_keys=newkeys[ii];chosen_points=newpts[ii]
        else:
            chosen_keys=newkeys;chosen_points=newpts
    if n==0: raise ValueError(f'Empty catalog {path}')
    acc['n']=n;acc['sample']=chosen_points
    return acc


def accumulate_by_case(row,manifest,lat_edges,lon_edges,args):
    dat={}
    case_id=str(row.identifier)
    for md in ('on','off','through'):
        cat=input_path(row[f'{md}_file'],manifest)
        print(f'[{case_id}] READING {md} ({cat.stat().st_size/1e6:.1f} MB): {cat}',flush=True)
        dat[md]=read_catalog(cat,lat_edges,lon_edges,args.scatter_max,np.random.default_rng(20261120+int(row.seed)%100000+{'on':1,'off':2,'through':3}[md]),args.chunksize)
    for k in ('flux','crater','flux2','crater2','count'):
        if not np.allclose(dat['off'][k],dat['on'][k]+dat['through'][k],rtol=1e-7,atol=0.01):
            raise ValueError(f'Paired identity OFF=ON+THROUGH fails for {case_id}, {k}')
    return dat


def rel_map(dat,area,kind,mode):
    avg=np.sum(dat['off'][kind])/np.sum(area)
    result=np.divide(dat[mode][kind],area,out=np.full_like(area,np.nan),where=area>0)/avg
    return result


def blocked_map(dat,area,neff_min=40):
    off=dat['off']['crater'];on=dat['through']['crater'];w2=dat['off']['crater2']
    eff=np.divide(off**2,w2,out=np.zeros_like(off),where=w2>0)
    frac=100*np.divide(on,off,out=np.full_like(off,np.nan),where=(off>0)&(area>0))
    frac[eff<neff_min]=np.nan
    return frac,eff


def cap_coords():
    ang=np.linspace(0,2*np.pi,360)
    lon0=np.deg2rad(-50);lat0=np.deg2rad(1.3);r=np.deg2rad(25)
    lat=np.arcsin(np.sin(lat0)*np.cos(r)+np.cos(lat0)*np.sin(r)*np.cos(ang))
    lon=lon0+np.arctan2(np.sin(ang)*np.sin(r)*np.cos(lat0),np.cos(r)-np.sin(lat0)*np.sin(lat))
    return np.degrees(lon),np.degrees(lat)

def axes(ax,lab=None):
    ax.set_xlim(-180,180);ax.set_ylim(-90,90)
    ax.set_xticks(np.arange(-180,181,90));ax.set_yticks(np.arange(-90,91,45))
    ax.set_xlabel('Body-fixed longitude (deg)');ax.set_ylabel('Latitude (deg)')
    ax.grid(alpha=.12,lw=.55)
    if lab:ax.text(-.09,1.05,lab,transform=ax.transAxes,fontsize=11,fontweight='bold',ha='left',va='bottom',clip_on=False)
    cx,cy=cap_coords();ax.plot(cx,cy,color='0.2',lw=.9,linestyle='--',alpha=.8)
    ax.axvline(-90,color='gray',lw=.5,ls=':',alpha=.5)
    ax.axvline(+90,color='gray',lw=.5,ls=':',alpha=.5)


def save(fig,outstem,figures):
    png=outstem.with_suffix('.png');pdf=outstem.with_suffix('.pdf')
    if png.exists() or pdf.exists():raise FileExistsError(f'{png} already exists. Supply a NEW --out directory (v2), do not overwrite.')
    fig.savefig(png,dpi=270,facecolor='white')
    fig.savefig(pdf,facecolor='white')
    plt.close(fig)
    figures.extend((str(png),str(pdf)))
    print('FIGURE',png,flush=True)


def draw_main(case,area,edges,out,figures,scale):
    dat=case['data'];lon_edges,lat_edges=edges
    fig,axs=plt.subplots(1,3,figsize=(13.7,4.85))
    fig.subplots_adjust(left=.075,right=.988,top=.85,bottom=.31,wspace=.30)
    ax=axs[0];pt=dat['on']['sample'];ax.scatter(pt[:,0],pt[:,1],marker='x',s=4,alpha=.38,linewidths=.37,c='black',rasterized=True)
    axes(ax,'(a)')
    fl=rel_map(dat,area,'flux','on');cr=rel_map(dat,area,'crater','on')
    ims=[None]
    im=axs[1].pcolormesh(lon_edges,lat_edges,fl,cmap='viridis',vmin=0,vmax=scale['flux'],rasterized=True,shading='flat');axes(axs[1],'(b)');ims.append(im)
    im=axs[2].pcolormesh(lon_edges,lat_edges,cr,cmap='viridis',vmin=0,vmax=scale['crater'],rasterized=True,shading='flat');axes(axs[2],'(c)');ims.append(im)
    fig.text(.17,.165,f'Unweighted subsample: {len(pt):,} of {dat["on"]["n"]:,} events',ha='center',fontsize=8.2)
    for i,label in [(1,'Impact-flux density / global OFF mean'),(2,'Crater-production proxy / global OFF mean')]:
        pos=axs[i].get_position(); cbax=fig.add_axes([pos.x0,pos.y0-.12,pos.width,.028]);cb=fig.colorbar(ims[i],cax=cbax,orientation='horizontal');cb.set_label(label,fontsize=8)
    save(fig,out/f'FigMain_impact_sites_flux_crater_{case["identifier"]}_v1',figures)


def draw_shield(case,area,edges,out,figures):
    lon_edges,lat_edges=edges;frac,neff=blocked_map(case['data'],area)
    fig,ax=plt.subplots(figsize=(10.6,5.1))
    fig.subplots_adjust(left=.095,right=.98,top=.86,bottom=.29)
    valid=frac[np.isfinite(frac)]
    vmax=max(12,float(np.percentile(valid,97))) if len(valid) else 20.0
    im=ax.pcolormesh(lon_edges,lat_edges,frac,cmap='magma',vmin=0,vmax=vmax,shading='flat',rasterized=True)
    if not len(valid):ax.text(.5,.5,'All bins below effective-count threshold',transform=ax.transAxes,ha='center',va='center',fontsize=9,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
    axes(ax,'(a)')
    pos=ax.get_position(); cax=fig.add_axes([pos.x0,pos.y0-.12,pos.width,.03]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_label('Locally Mars-blocked crater-production proxy: THROUGH / OFF (%)',fontsize=9)
    save(fig,out/f'FigSupp_shielding_crater_{case["identifier"]}_v1',figures)


def draw_pair_comparison(cases,area,edges,out,figures,scale):
    if len(cases)<2:return
    c1,c2=cases[:2];lon_edges,lat_edges=edges
    fig,axs=plt.subplots(2,2,figsize=(10.3,7.3),sharex=True,sharey=True)
    fig.subplots_adjust(left=.10,right=.95,top=.92,bottom=.16,hspace=.28,wspace=.19)
    lab=['(a)','(b)','(c)','(d)']
    for row,ca in enumerate((c1,c2)):
        for col,kind in enumerate(('flux','crater')):
            ax=axs[row,col]
            vals=rel_map(ca['data'],area,kind,'on')
            ax.pcolormesh(lon_edges,lat_edges,vals,cmap='viridis',vmin=0,vmax=scale[kind],shading='flat',rasterized=True)
            axes(ax,lab[row*2+col])
            ax.text(.02,.06,f'a = {ca["orbit_km"]:g} km',ha='left',va='bottom',transform=ax.transAxes,fontsize=8,
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':2})
    for col,kind in enumerate(('flux','crater')):
        pos=axs[-1,col].get_position();cax=fig.add_axes([pos.x0,pos.y0-.10,pos.width,.023])
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,scale[kind]),cmap='viridis'),cax=cax,orientation='horizontal')
        cb.set_label(('Impact-flux' if kind=='flux' else 'Crater-production proxy')+' density / respective global OFF mean',fontsize=8)
    vtag=f'{cases[0]["v_inf_km_s"]:g}'.replace('.','p')
    save(fig,out/f'FigSupp_orbit_compare_v{vtag}_density_v1',figures)


def draw_raw_pair(cases,edges,out,figures):
    if len(cases)<2:return
    c1,c2=cases[:2]
    fig,axs=plt.subplots(2,2,figsize=(10.1,7.25),sharex=True,sharey=True)
    fig.subplots_adjust(left=.10,right=.96,top=.92,bottom=.11,hspace=.25,wspace=.18)
    for ri,c in enumerate((c1,c2)):
        for ci,md in enumerate(('off','on')):
            ax=axs[ri,ci];pts=c['data'][md]['sample']
            ax.scatter(pts[:,0],pts[:,1],marker='x',s=3.2,alpha=.36,linewidths=.34,color='black',rasterized=True)
            axes(ax,f'({chr(97+ri*2+ci)})')
            ax.text(.02,.055,f'{md.upper()} | a={c["orbit_km"]:g} km | {len(pts):,} shown',transform=ax.transAxes,
                    ha='left',va='bottom',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':1})
    save(fig,out/'FigSupp_raw_impact_x_ON_OFF_two_orbits_v1',figures)


def grid_export(case,area,lat_edges,lon_edges,out):
    dat=case['data'];frac,eff=blocked_map(dat,area)
    p=out/f'grid15_{case["identifier"]}_v1.csv'
    if p.exists():raise FileExistsError(p)
    with p.open('w',newline='') as f:
        names=['lat_min','lat_max','lon_min','lon_max','area_km2','ON_hits','OFF_hits','THROUGH_hits','ON_flux_weight_sum','OFF_flux_weight_sum','THROUGH_flux_weight_sum','ON_crater_weight_sum','OFF_crater_weight_sum','THROUGH_crater_weight_sum','ON_flux_relative_global_OFF_mean','OFF_flux_relative_global_OFF_mean','ON_crater_relative_global_OFF_mean','OFF_crater_relative_global_OFF_mean','Mars_blocked_crater_percent_masked_neff_lt40','OFF_crater_effective_count']
        w=csv.writer(f);w.writerow(names)
        fl_on=rel_map(dat,area,'flux','on');fl_off=rel_map(dat,area,'flux','off')
        cr_on=rel_map(dat,area,'crater','on');cr_off=rel_map(dat,area,'crater','off')
        for i in range(len(lat_edges)-1):
            for j in range(len(lon_edges)-1):
                w.writerow([lat_edges[i],lat_edges[i+1],lon_edges[j],lon_edges[j+1],area[i,j],
                            *[int(round(dat[m]['count'][i,j])) for m in ('on','off','through')],
                            *[dat[m]['flux'][i,j] for m in ('on','off','through')],
                            *[dat[m]['crater'][i,j] for m in ('on','off','through')],
                            fl_on[i,j],fl_off[i,j],cr_on[i,j],cr_off[i,j],frac[i,j],eff[i,j]])
    return p


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',type=Path,default=Path('results/lwss_campaign_v1/v004pg1_n10000000_R20_dt1_post450/manifest_v1.csv'))
    p.add_argument('--shape',type=Path,default=Path('shape_model/phobos_g_148m_spc_obj_0000n00000_v004.obj'))
    p.add_argument('--out',type=Path,default=Path('figures_for_paper_v004_pg1_v1/02_pg1_maps'))
    p.add_argument('--orbits',type=float,nargs='+',default=[9376,20000])
    p.add_argument('--speed',type=float,default=7.95298)
    p.add_argument('--bin-deg',type=int,default=15)
    p.add_argument('--scatter-max',type=int,default=14000)
    p.add_argument('--chunksize',type=int,default=80000)
    p.add_argument('--dry-run',action='store_true')
    args=p.parse_args()
    if 360%args.bin_deg or 180%args.bin_deg:raise ValueError('bin-deg must divide 180 and 360')
    manifest=args.manifest.resolve();shape=args.shape.resolve();out=args.out.resolve()
    if not manifest.is_file() or not shape.is_file():raise FileNotFoundError(f'Could not find manifest={manifest} or shape={shape}')
    m=pd.read_csv(manifest)
    rows=m[np.isclose(m.v_inf_km_s,args.speed)&np.isin(m.orbit_km,args.orbits)].sort_values('orbit_km')
    if len(rows)!=len(args.orbits):raise ValueError(f'Matched {len(rows)} cases but requested {len(args.orbits)}: {rows.identifier.tolist()}')
    print('INPUTS: v004 + PG ON; reading existing 26-col catalogs, no trajectory propagation.')
    for _,r in rows.iterrows():
        print('CASE',r.identifier,'speed',r.v_inf_km_s,'orbit',r.orbit_km)
        if str(r.status)!='complete' or int(r.num)!=10000000 or str(r.model_file_sha256)!=EXPECTED_SHA:
            raise ValueError(f'Unexpected provenance for {r.identifier}')
        for mode in ('on','off','through'):
            path=input_path(r[f'{mode}_file'],manifest);print(' ',mode,path.name,round(path.stat().st_size/1e6,1),'MB')
    if args.dry_run:return
    if out.exists() and any(out.iterdir()):raise FileExistsError(f'Output directory {out} has contents. Change to a versioned --out ending in v2.')
    out.mkdir(parents=True,exist_ok=True)
    sh=sha(shape)
    if sh!=EXPECTED_SHA:raise ValueError(f'Wrong v004 OBJ SHA256 {sh}')
    lat_edges=np.arange(-90,90+args.bin_deg,args.bin_deg,dtype=np.float64)
    lon_edges=np.arange(-180,180+args.bin_deg,args.bin_deg,dtype=np.float64)
    area,geo=geometry(shape,lon_edges,lat_edges)
    print('PHYSICAL AREA',geo)
    cases=[]
    for _,r in rows.iterrows():
        dat=accumulate_by_case(r,manifest,lat_edges,lon_edges,args)
        cases.append({'identifier':str(r.identifier),'orbit_km':float(r.orbit_km),'v_inf_km_s':float(r.v_inf_km_s),'N':int(r.num),'seed':int(r.seed),'data':dat})
    scales={}
    for kind in ('flux','crater'):
        arr=np.concatenate([rel_map(x['data'],area,kind,'on')[np.isfinite(rel_map(x['data'],area,kind,'on'))].ravel() for x in cases])
        scales[kind]=float(max(1.0,np.percentile(arr,99)))
    print('COLOR SCALES',scales)
    plt.rcParams.update({'font.size':9,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'savefig.dpi':280,'pdf.fonttype':42,'ps.fonttype':42})
    figs=[];small=[]
    for c in cases:
        small.append(str(grid_export(c,area,lat_edges,lon_edges,out)))
        draw_main(c,area,(lon_edges,lat_edges),out,figs,scales)
        draw_shield(c,area,(lon_edges,lat_edges),out,figs)
    draw_raw_pair(cases,(lon_edges,lat_edges),out,figs)
    draw_pair_comparison(cases,area,(lon_edges,lat_edges),out,figs,scales)
    meta={'description':'Maps generated from archived 26-column event catalog (PG1/Ernst v004), not rerun.',
      'shape_SHA256':sh,'shape':str(shape),'manifest':str(manifest),'formula':'w_crater=w_flux*(v_surface/8000 m/s)^1.375*clip(incidence_cos_raw,0,1)^1.375',
      'normalization':'Area per 15-degree cell from v004 triangles, local cell density / global OFF mean separately for each flux/crater proxy and each case',
      'shield_fraction_formula':'100*sum(w_crater THROUGH)/sum(w_crater OFF); masked if OFF weighted effective count <40',
      'scatter':'uniform unweighted event-row reservoir sample, max '+str(args.scatter_max),
      'bin_deg':args.bin_deg,'mesh':geo,'scales':scales,'cases':[{'identifier':x['identifier'],'orbit_km':x['orbit_km'],'v_inf_km_s':x['v_inf_km_s'],'N':x['N'],'seed':x['seed'],'hit_count_by_mode':{md:x['data'][md]['n'] for md in ('on','off','through')}} for x in cases],
      'figures':figs,'grids':small,'warning':'Do not reuse raw old PG OFF maps as results of the new PG ON calculation. The maps depict model impact/production, not observed crater densities.'}
    (out/'figure_provenance_v1.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False)+'\n')
    captions=[
        'MAIN: Individual primary-impact sites (unweighted x-marks), area-normalized impact-flux density, and illustrative crater-production proxy at 7.95298 km/s. PG ON; Ernst SPC v004; shielding ON; 15-degree binning with true mesh surface area; all densities normalized to the global shielding-OFF mean. Dashed contour is the 25-degree Stickney cap. Scatter is a uniform subset of event rows for visualization, not a density estimator. No spatial smoothing.',
        'SUPP: Local production lost to Mars shielding computed from paired THROUGH and OFF crater-production sums, cell by cell. Cells with fewer than 40 OFF weighted-effective impacts are masked. Geometric 25-degree Stickney contour is indicated; this map is not evidence of crater preservation.',
        'SUPP: Raw primary-impact site maps at two orbital radii and for shielding ON/OFF. Uniform unweighted event subsamples; marker counts do not measure absolute flux. Both target mesh and force model are PG ON, Ernst SPC v004.',
        'SUPP: Two orbital radii compared for incident impact-flux and crater-production proxy per physical area, each relative to its own global shielding-OFF mean. Same display scale across orbital radii; 15-degree bins, no smoothing. PG ON, v004.'
    ]
    (out/'CAPTIONS_FOR_MANUSCRIPT_v1.txt').write_text('\n\n'.join(captions)+'\n')
    print('DONE',len(figs)//2,'figures (PNG+PDF),',len(small),'reduced gridded CSVs')
    print('IMPORTANT: upload the PNGs and grid15_*.csv to ChatGPT; NOT the huge event files.')
if __name__=='__main__':main()
