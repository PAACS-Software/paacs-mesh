#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Regenerate archive-backed ABC, sizing and cylinder figures from native tools.

No configure/build step or mesher change. Run from any directory; explicit
executable and evidence paths make the generating inputs reviewable.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')


def invoke(command, log, runs):
    result=subprocess.run([str(x) for x in command],capture_output=True,text=True)
    log.write_text(result.stdout,encoding='utf-8')
    log.with_suffix('.stderr.log').write_text(result.stderr,encoding='utf-8')
    runs.append({'command':[str(x) for x in command],'returncode':result.returncode,
                 'stdout_path':str(log),'stdout_sha256':digest(log)})
    if result.returncode:
        raise RuntimeError(f'command failed ({result.returncode}): {command[0]}: {result.stderr[-1000:]}')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive-verification',type=Path,required=True)
    parser.add_argument('--spatial-executable',type=Path,required=True)
    parser.add_argument('--primitive-executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--skip-generation',action='store_true',help='use existing example inputs in output/inputs')
    args=parser.parse_args()
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=True)
    inputs=output/'inputs'; inputs.mkdir(exist_ok=True)
    runs=[]
    records=json.loads(args.archive_verification.read_text(encoding='utf-8-sig'))
    selected=[r for r in records if r['repeat']==0 and r['verification']['verified']]
    abc={'title':'PAACS / Five supported ABC planar bodies',
         'subtitle':'Real STEP imports · Native audited archives · Cyan exterior / gold exact interior section',
         'columns':3,'cases':[]}
    for r in selected:
        archive=Path(r['command'][2]).resolve(); source=Path(r['command'][1]).resolve()
        if digest(archive) != r['archive_sha256']:
            raise ValueError('ABC archive hash differs from native verification receipt')
        abc['cases'].append({'name':'ABC '+str(int(Path(r['path']).parts[0])),
                             'mesh':str(archive),'source':str(source),'units':'mm','show_bounds':True,
                             'verification':r['verification'],
                             'metadata':{'native_receipt_path':str(args.archive_verification.resolve()),
                                         'native_receipt_sha256':digest(args.archive_verification)},
                             'section':{'axis':'y','fraction':.5}})
    write(output/'abc-manifest.json',abc)
    spatial=inputs/'spatial'; spatial.mkdir(exist_ok=True)
    executable=args.spatial_executable.resolve()
    if not args.skip_generation:
        invoke([executable,'--gallery',spatial],inputs/'spatial-native.log',runs)
    spatial_cases=[('point','Point influence'),('segment','Segment influence'),
                   ('triangle','Triangle patch'),('box','Box region'),('ball','Ball region'),
                   ('minimum','Minimum: either region'),('maximum','Maximum: overlap'),
                   ('smooth','Smooth transition + clamp')]
    sizing={'title':'PAACS / Composable spatial size controls',
            'subtitle':'Reproduced native tests · Actual unit-cube Tet4 meshes · Exact y = 0.5 sections',
            'gallery_mode':'sections','columns':4,'cases':[]}
    for key,name in spatial_cases:
        sizing['cases'].append({'name':name,'mesh':str(spatial/(key+'.paacsmesh')),
            'field_csv':str(spatial/(key+'-field.csv')),'report':False,'units':'model units',
            'section':{'axis':'y','value':.5,'expected_section_area':1},
            'metadata':{'generator_path':str(executable),'generator_sha256':digest(executable),
                        'native_test_scope':'19 configurations/23 refusals, size-field audit, independent volume/boundary/interior-size/thread identity checks'}})
    write(output/'spatial-manifest.json',sizing)
    cylinder=inputs/'cylinder'; cylinder.mkdir(exist_ok=True)
    primitive=args.primitive_executable.resolve()
    comparison={'title':'PAACS / Cylinder volume uniformity',
                'subtitle':'Radius 1 · Height 2 · Requested size 0.3 · Same chord/volume tolerances · Straight-sided geometry',
                'columns':3,'distributions':True,'cases':[]}
    for key,family,strategy,name in [('swept-tet4','tet4','swept','Swept Tet4'),
           ('uniform-tet4','tet4','uniform','Uniform Tet4'),('uniform-tet10','tet10','uniform','Uniform Tet10')]:
        mesh=cylinder/(key+'.paacsmesh'); log=cylinder/(key+'.json')
        command=[primitive,'cylinder','.3',mesh,'--elements',family,'--cylinder-strategy',strategy,
                 '--threads','1','--chord-error','.01','--volume-error','.03']
        if not args.skip_generation:
            invoke(command,log,runs)
        native=json.loads(log.read_text(encoding='utf-8-sig'))
        if not native['audit_passed']:
            raise ValueError('cylinder native audit failed')
        comparison['cases'].append({'name':name,'mesh':str(mesh),'units':'model units',
                'show_midsides':family=='tet10','section':{'axis':'y','value':0},
                'verification':{'verified':True,'mesh_volume':native['mesh_volume'],
                                'minimum_quality':native['minimum_quality']},
                'metadata':{'native_result':native,'command':[str(x) for x in command],
                            'generator_sha256':digest(primitive)}})
    write(output/'cylinder-manifest.json',comparison)
    renderer=Path(__file__).with_name('mesh_visual_report.py')
    for name in ('abc','spatial','cylinder'):
        invoke([sys.executable,renderer,'--manifest',output/(name+'-manifest.json'),
                '--output',output/name],output/(name+'-render.log'),runs)
    write(output/'reproduction.json',{'schema':'paacs-mesh-visual-reproduction-v1',
                                    'driver_sha256':digest(__file__),'runs':runs})
    print(output/'reproduction.json')


if __name__=='__main__':
    main()
