#!/usr/bin/env python3
"""Reproducible, no-overwrite LWSS-MC campaign for v004 OBJ and Phobos point-gravity ON."""
from __future__ import annotations
import argparse, subprocess, sys, os, csv, hashlib, json, re
from pathlib import Path

BASE=Path(__file__).resolve().parents[1]
DEFAULT_SPEEDS='0.5,1,2,3,4,7.95298,20,40,60'
DEFAULT_ORBITS='9376,12000,15000,20000'

def comma_floats(value):
    try:
        out=[float(x.strip()) for x in value.split(',') if x.strip()]
    except ValueError as e:
        raise argparse.ArgumentTypeError(str(e)) from e
    if not out or any(x <= 0 for x in out) or len(set(out)) != len(out):
        raise argparse.ArgumentTypeError('Speeds/orbits must be distinct positive numbers')
    return out

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()

def parse_summary(path):
    return dict(line.split('=',1) for line in path.read_text().splitlines() if '=' in line)

def tag_number(value):
    return ('%g' % value).replace('.','p')

def main():
    p=argparse.ArgumentParser(description='LWSS-MC v004+PG1: default is safe dry run; --execute performs runs.')
    p.add_argument('--execute',action='store_true',help='Actually run the campaign. Without flag only print plan.')
    p.add_argument('--num',type=int,default=10_000_000,help='Proposals per case, default 10000000')
    p.add_argument('--threads',type=int,default=8)
    p.add_argument('--seed-base',type=int,default=2026110600)
    p.add_argument('--source-radius-km',type=float,default=20.)
    p.add_argument('--dt-scale',type=float,default=1.)
    p.add_argument('--post-seconds',type=float,default=450.)
    p.add_argument('--speeds-km-s',type=comma_floats,default=comma_floats(DEFAULT_SPEEDS))
    p.add_argument('--orbits-km',type=comma_floats,default=comma_floats(DEFAULT_ORBITS))
    p.add_argument('--max-cases',type=int,default=0,help='Stop after this many NEW runs; 0 means no limit.')
    p.add_argument('--shape',default='shape_model/phobos_g_148m_spc_obj_0000n00000_v004.obj')
    a=p.parse_args()
    if a.num<=1 or a.threads<=0 or a.source_radius_km<=0 or a.dt_scale<=0 or a.post_seconds<=0 or a.max_cases<0:
        p.error('Invalid numeric options')
    shape=(BASE/a.shape).resolve() if not Path(a.shape).is_absolute() else Path(a.shape).resolve()
    src=BASE/'src/phobos_impact_source_sphere_production_v5.cpp'
    dep=BASE/'src/phobos_impact_production_v3.cpp'
    binary=BASE/'bin/phobos_source_sphere_v5_lwss'
    for f in [shape,src,dep]:
        if not f.is_file(): p.error('Missing required file: '+str(f))
    if 'v004' not in shape.name:
        p.error('Refusing other than v004 model: '+str(shape))
    if a.source_radius_km <= 14.: p.error('Source sphere must enclose v004 shape (rbound ~13.843 km) with margin')
    if min(a.orbits_km)*1000 - a.source_radius_km*1000 <= 3389500:
        p.error('Source sphere intersects Mars')
    campaign=f'v004pg1_n{a.num}_R{tag_number(a.source_radius_km)}_dt{tag_number(a.dt_scale)}_post{tag_number(a.post_seconds)}'
    out=BASE/'results/lwss_campaign_v1'/campaign
    logs=BASE/'logs/lwss_campaign_v1'/campaign
    table=out/'manifest_v1.csv'
    cases=[]
    for oi,orbit in enumerate(a.orbits_km):
        for vi,speed in enumerate(a.speeds_km_s):
            seed=a.seed_base+oi*100+vi+1
            ident=f'a{tag_number(orbit)}_v{tag_number(speed)}_seed{seed}'
            files={suffix:out/f'{ident}_{suffix}.txt' for suffix in ('on','off','through')}
            files['summary']=out/f'{ident}.summary.txt'
            complete=all(path.exists() for path in files.values())
            anyfile=any(path.exists() for path in files.values())
            if anyfile and not complete:
                raise RuntimeError('Partial output; manually inspect before resuming: '+ident)
            if complete:
                data=parse_summary(files['summary'])
                if any([
                  data.get('model')!='source_sphere_production_v5',
                  int(data.get('num',-1))!=a.num,
                  abs(float(data.get('v_inf_mps','nan'))-1000*speed)>1e-6,
                  abs(float(data.get('orbit_r_km','nan'))-orbit)>1e-6,
                  data.get('phobos_gravity')!='1',
                  data.get('mars_gravity')!='1',
                  int(data.get('threads',-1))!=a.threads,
                  data.get('shape_file')!=str(shape),
                  abs(float(data.get('dt_scale','nan'))-a.dt_scale)>1e-9,
                  abs(float(data.get('post_seconds','nan'))-a.post_seconds)>1e-9,
                  abs(float(data.get('source_radius_km','nan'))-a.source_radius_km)>1e-6
                ]): raise RuntimeError('Existing summary does not match requested case: '+ident)
            cases.append(dict(identifier=ident,orbit_km=orbit,v_inf_km_s=speed,seed=seed,model_file=str(shape),model_file_sha256='',num=a.num,source_radius_km=a.source_radius_km,threads=a.threads,status='complete' if complete else 'pending',**{k+'_file':str(v) for k,v in files.items()}))
    print('Campaign:',campaign,'shape=',shape,'cases=',len(cases),'complete=',sum(x['status']=='complete' for x in cases))
    print('Speed v_inf km/s:',a.speeds_km_s,'orbits km:',a.orbits_km)
    if not a.execute:
        for x in cases: print(' ',x['status'],x['identifier'])
        print('Dry-run only. Add --execute to start. Use --max-cases for a controlled pilot.')
        return
    out.mkdir(parents=True,exist_ok=True)
    logs.mkdir(parents=True,exist_ok=True)
    binary.parent.mkdir(parents=True,exist_ok=True)
    source_shas={str(k):sha256(k) for k in [src,dep,shape]}
    if not binary.exists() or max(src.stat().st_mtime,dep.stat().st_mtime)>binary.stat().st_mtime:
        cmd=['g++','-O3','-std=c++17','-fopenmp',str(src),'-o',str(binary)]
        print('Build:', ' '.join(cmd),flush=True)
        subprocess.run(cmd,check=True)
    print('Binary sha256:',sha256(binary),flush=True)
    changed=0
    try:
        for x in cases:
            if x['status']=='complete':
                print('SKIP:',x['identifier'],flush=True);continue
            if a.max_cases and changed>=a.max_cases: break
            ident=x['identifier'];cmd=[str(binary),'--num',str(a.num),'--seed',str(x['seed']),'--orbit-r',str(int(round(x['orbit_km']*1000))),'--radius','11100','--v-inf',str(x['v_inf_km_s']*1000),'--mars-grav','1','--phobos-grav','1','--shape',str(shape),'--obj-unit','km','--source-radius-km',str(a.source_radius_km),'--dt-scale',str(a.dt_scale),'--post-seconds',str(a.post_seconds),'--summary',x['summary_file'],'--out-on',x['on_file'],'--out-off',x['off_file'],'--out-through',x['through_file']]
            print('RUN:',ident,flush=True)
            with (logs/f'{ident}.log').open('w') as logfile:
                logfile.write('CMD: '+' '.join(cmd)+'\n')
                logfile.write('PROVENANCE: '+json.dumps(source_shas,ensure_ascii=False)+'\n')
                logfile.flush()
                try: subprocess.run(cmd,stdout=logfile,stderr=subprocess.STDOUT,env={**os.environ,'OMP_NUM_THREADS':str(a.threads)},check=True)
                except subprocess.CalledProcessError:
                    print('ERROR. Inspect:',logs/f'{ident}.log','and temporary .source_sphere_part files',file=sys.stderr)
                    raise
            x['status']='complete';changed+=1
            print('COMPLETE:',ident,flush=True)
    finally:
        for x in cases: x['model_file_sha256']=source_shas[str(shape)]
        with table.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(cases[0].keys()));writer.writeheader();writer.writerows(cases)
        (out/'provenance_v1.json').write_text(json.dumps({'source_files':source_shas,'binary_sha256':sha256(binary),'options':vars(a)},ensure_ascii=False,indent=2))
        print('Manifest:',table,flush=True)
    print('New runs:',changed,'completed:',sum(x['status']=='complete' for x in cases),'/',len(cases),flush=True)

if __name__=='__main__': main()
