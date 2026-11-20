#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Reproducible, measured mesh figures; see mesh_visual_report.md for inputs.

Requires numpy and matplotlib. Reads PAACS native v1 archives without a lossy
text export, or PAACS legacy ASCII VTK. Scientific plotting only: no ImageGen.
"""
import argparse
import json
import math
from pathlib import Path
import re
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm, Normalize, to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from mesh_visual_data import (CORNERS, NAMES, clip_polygon, read_geometry,
                              read_mesh, section, self_test, sha256, volume_data)

BG, INK, SUB = '#f3f6f9', '#18394d', '#526b7c'
CYAN, GOLD, EDGE = '#70bec9', '#e7bf7b', '#294d5c'
plt.rcParams.update({'font.family':'DejaVu Sans', 'text.color':INK,
                     'axes.labelcolor':SUB, 'xtick.color':SUB, 'ytick.color':SUB})


def slug(name):
    return re.sub(r'[^a-z0-9]+', '-', name.lower()).strip('-') or 'mesh'


def resolve(base, path):
    value = Path(path)
    return value.resolve() if value.is_absolute() else (base/value).resolve()


def simple_axes(ax):
    ax.set_facecolor(BG)
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_color('#a6b6c0')
    ax.tick_params(labelsize=9)


def histogram_bins(values, count=45):
    low, high = float(values.min()), float(values.max())
    # A genuinely uniform volume can still span a handful of floating-point
    # ulps. Fixed-width histogram bins then need a visible nonzero interval.
    if high-low < max(abs(low),abs(high),1)*1e-10:
        padding = max(abs((high+low)*.5),1)*.025
        low -= padding; high += padding
    return np.linspace(low,high,count+1)


def frame3d(ax, mesh, elevation=27, azimuth=-58):
    lo,hi = mesh.bounds
    extent = hi-lo
    margin = max(float(np.max(extent))*.03, np.finfo(float).tiny)
    ax.set(xlim=(lo[0]-margin,hi[0]+margin), ylim=(lo[1]-margin,hi[1]+margin),
           zlim=(lo[2]-margin,hi[2]+margin))
    ax.set_box_aspect(np.maximum(extent, np.max(extent)*.015))
    ax.view_init(elev=elevation, azim=azimuth)
    ax.set_proj_type('ortho')
    ax.set_axis_off()


def cutaway(ax, item):
    mesh, polygons, boundary = item['mesh'], item['polygons'], item['boundary']
    axis, value, keep = item['axis'], item['value'], item['keep']
    el,az = item['case'].get('camera', {}).get('elevation',27), \
            item['case'].get('camera', {}).get('azimuth',-58)
    view = np.array([math.cos(math.radians(el))*math.cos(math.radians(az)),
                     math.cos(math.radians(el))*math.sin(math.radians(az)),
                     math.sin(math.radians(el))])
    light = np.array([.5,-.8,1.]); light /= np.linalg.norm(light)
    faces, colors, visible_boundary = [], [], set()
    # Mesh bounds scan every node. The loaded mesh is fixed throughout this
    # drawing pass; computing scale per face made dense galleries quadratic.
    minimum_normal_length = mesh.scale**2*1e-24
    for ids in boundary:
        p = clip_polygon(mesh.points[list(ids)], axis, value, keep)
        if len(p) < 3:
            continue
        normal = np.cross(p[1]-p[0],p[2]-p[0])
        length = np.linalg.norm(normal)
        if length <= minimum_normal_length:
            continue
        normal /= length
        if normal@view <= 0:
            continue
        visible_boundary.update(ids)
        faces.append(p)
        colors.append(np.asarray(to_rgb(CYAN))*(.75+.25*max(0,float(normal@light))))
    # The cut polygons are the precise original-cell intersections, including
    # quadrilateral cuts through tetrahedra. They are never triangulated anew.
    if -keep*view[axis] > 0:
        for p in polygons:
            faces.append(p); colors.append(to_rgb(GOLD))
    ax.add_collection3d(Poly3DCollection(faces, facecolors=colors,
                         edgecolors=EDGE, linewidths=.27, zsort='average'))
    if item['case'].get('show_midsides', False):
        mids = set()
        for cell in mesh.cells:
            if cell.kind == 6:
                for e,(a,b) in enumerate(((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))):
                    if cell.ids[a] in visible_boundary and cell.ids[b] in visible_boundary:
                        mids.add(cell.ids[4+e])
        if mids:
            p = mesh.points[sorted(mids)]
            p = p[keep*(p[:,axis]-value) >= 0]
            ax.scatter(*p.T, s=2.2, color=INK, depthshade=False)
    frame3d(ax, mesh, el, az)


def section_axes(ax, item, norm=None, cmap='YlGnBu'):
    axes = [i for i in range(3) if i != item['axis']]
    sizes = np.cbrt(item['volumes'])[item['section_cells']]
    if norm is None:
        low, high = float(sizes.min()), float(sizes.max())
        norm = LogNorm(low*.95, high*1.05) if high-low < high*1e-10 else LogNorm(low,high)
    colors = plt.get_cmap(cmap)(norm(sizes))
    projected = [p[:,axes] for p in item['polygons']]
    ax.add_collection(PolyCollection(projected, facecolors=colors,
                         edgecolors=EDGE, linewidths=.25))
    lo,hi = item['mesh'].bounds
    ax.set(xlim=(lo[axes[0]],hi[axes[0]]), ylim=(lo[axes[1]],hi[axes[1]]), aspect='equal')
    # Keep narrow physical sections at the panel's left edge so their title
    # remains inside the figure after the equal-aspect axes shrink.
    ax.set_anchor('W')
    unit = item['case'].get('units','model units')
    ax.set_xlabel(f'{"xyz"[axes[0]]} ({unit})', fontsize=9)
    ax.set_ylabel(f'{"xyz"[axes[1]]} ({unit})', fontsize=9)
    simple_axes(ax)
    return norm


def load_case(case, base):
    path = resolve(base, case['mesh'])
    mesh = read_mesh(path)
    v,q,ids,boundary,metrics = volume_data(mesh)
    plane = case.get('section', {})
    axis = 'xyz'.index(plane.get('axis','y'))
    lo,hi = mesh.bounds
    value = plane.get('value', lo[axis]+plane.get('fraction',.5)*(hi[axis]-lo[axis]))
    polygons, cell_numbers, sections = section(mesh, ids, v, axis, value)
    if not polygons:
        raise ValueError(f'{case["name"]}: selected plane intersects no cell interiors')
    keep = plane.get('keep', -1 if axis == 2 else 1)
    if keep not in (-1,1):
        raise ValueError('section.keep must be -1 or +1')
    record = {'name':case['name'], 'mesh_path':str(path), 'mesh_sha256':sha256(path),
              'mesh_bytes':path.stat().st_size, 'units':case.get('units','model units'),
              'metrics':metrics, 'section':sections, 'case_metadata':case.get('metadata',{})}
    warning_target = case.get('quality_warning_target')
    if warning_target is not None:
        if not isinstance(warning_target, (int, float)) or not math.isfinite(warning_target) or not 0 <= warning_target <= 1:
            raise ValueError('quality_warning_target must be finite and between zero and one')
        record['quality_warning'] = {'target':warning_target,
                                    'volume_cells_below_target':int(np.count_nonzero(q < warning_target)),
                                    'cell_fraction_below_target':float(np.mean(q < warning_target))}
    native = case.get('verification')
    if native:
        record['native_verification'] = native
        if native.get('verified') is False:
            raise ValueError('manifest marks this native mesh unverified')
        expected = native.get('divergence_volume', native.get('mesh_volume'))
        if expected is not None:
            error = abs(metrics['volume_sum']-expected)/abs(expected)
            record['independent_volume_relative_delta_from_native'] = error
            if error > 1e-9:
                raise ValueError('independent corner volume disagrees with native receipt')
        expected = native.get('minimum_quality')
        if expected is not None:
            delta = abs(metrics['subtet_mean_ratio_min']-expected)
            record['independent_quality_delta_from_native'] = delta
            if delta > 1e-9:
                raise ValueError('independent native-defined quality disagrees with native receipt')
    if 'expected_section_area' in plane:
        expected = plane['expected_section_area']
        delta = abs(sections['section_area']-expected)/abs(expected)
        record['analytic_section_relative_error'] = delta
        if delta > 1e-9:
            raise ValueError('actual section area disagrees with analytic reference')
    if case.get('source'):
        source = resolve(base,case['source'])
        record.update(source_path=str(source), source_sha256=sha256(source))
    if case.get('field_csv'):
        field = resolve(base,case['field_csv'])
        record.update(field_csv_path=str(field), field_csv_sha256=sha256(field))
    return dict(case=case,mesh=mesh,volumes=v,quality=q,volume_cells=ids,
                boundary=boundary,polygons=polygons,section_cells=cell_numbers,
                axis=axis,value=value,keep=keep,record=record)


def stats_text(item):
    m = item['record']['metrics']
    minimum = m['tet_mean_ratio_min'] if m['tet_mean_ratio_samples'] == m['volume_cells'] else m['subtet_mean_ratio_min']
    measured = f'{minimum:.2e}' if 0 < minimum < .001 else f'{minimum:.3f}'
    quality = (f'Min. tet quality {measured}' if m['tet_mean_ratio_samples'] == m['volume_cells']
               else f'Min. SubTet quality {measured}')
    text = (f"{m['volume_cells']:,} volume cells · {m['nodes']:,} nodes\n"
            f"{quality} · Volume CV {100*m['volume_cv']:.1f}%")
    if item['case'].get('show_bounds',False):
        lo,hi = item['mesh'].bounds
        text += '\nBounds ' + ' × '.join(f'{x:.4g}' for x in hi-lo) + ' ' + item['case'].get('units','')
    if item['case'].get('stats_note'):
        text += '\n' + str(item['case']['stats_note'])
    return text


def size_bar(fig, scalar, **kwargs):
    bar = fig.colorbar(scalar, **kwargs)
    ticks = np.geomspace(scalar.norm.vmin,scalar.norm.vmax,5)
    bar.set_ticks(ticks); bar.set_ticklabels([f'{v:.3g}' for v in ticks])
    bar.minorticks_off()
    return bar


def save(fig, path, figures):
    fig.savefig(path, dpi=160, facecolor=BG)
    plt.close(fig)
    figures.append({'path':str(path.resolve()),'sha256':sha256(path)})


def individual(item, output, figures):
    fig = plt.figure(figsize=(13.4,10.6), facecolor=BG)
    fig.text(.055,.962,item['case']['name'],fontsize=22,weight='bold')
    fig.text(.055,.923,'Actual mesh · Cyan exterior / gold interior · Exact cell–plane section',fontsize=11,color=SUB)
    grid = fig.add_gridspec(2,2,left=.065,right=.90,top=.865,bottom=.16,
                           hspace=.44,wspace=.3,height_ratios=[1.65,1])
    ax = fig.add_subplot(grid[0,0],projection='3d',facecolor=BG)
    cutaway(ax,item)
    ax.set_title('Mesh cutaway',loc='left',weight='bold',fontsize=12)
    ax.text2D(0,-.03,stats_text(item),transform=ax.transAxes,color=SUB,fontsize=10,linespacing=1.5)
    ax = fig.add_subplot(grid[0,1],facecolor=BG)
    norm = section_axes(ax,item)
    ax.set_title(f'Actual section: {"xyz"[item["axis"]]} = {item["value"]:.6g}',
                 loc='left',fontsize=12,weight='bold')
    bar = size_bar(fig,plt.cm.ScalarMappable(norm=norm,cmap='YlGnBu'),ax=ax,fraction=.045,pad=.05)
    bar.set_label('Cell size: ∛volume',fontsize=9)
    ax = fig.add_subplot(grid[1,0]); simple_axes(ax)
    ratio = item['volumes']/item['volumes'].mean()
    ax.hist(ratio,bins=histogram_bins(ratio),weights=np.full(len(ratio),100/len(ratio)),
            color=CYAN,edgecolor=EDGE,linewidth=.3)
    ax.axvline(1,color='#b68442',lw=1.3,ls='--')
    ax.set(xlabel='Cell volume / mean cell volume',ylabel='Cells (%)')
    ax.set_title('Volume distribution',loc='left',fontsize=12,weight='bold')
    ax = fig.add_subplot(grid[1,1]); simple_axes(ax)
    q = item['quality']
    ax.hist(q,bins=np.linspace(0,1,41),weights=np.full(len(q),100/len(q)),
            color=GOLD,edgecolor=EDGE,linewidth=.3)
    mixed = item['record']['metrics']['tet_mean_ratio_samples'] != len(q)
    label = 'Native fixed SubTetMeanRatio' if mixed else 'Tetrahedron mean ratio (ideal = 1)'
    ax.set(xlabel=label,ylabel='Volume cells (%)',xlim=(0,1))
    ax.set_title('Quality distribution',loc='left',fontsize=12,weight='bold')
    warning = item['record'].get('quality_warning')
    if warning:
        ax.axvline(warning['target'],color='#b68442',lw=1.3,ls='--')
        ax.set_title(f"Quality · {100*warning['cell_fraction_below_target']:.1f}% below {warning['target']:g}",
                     loc='left',fontsize=12,weight='bold')
    short_hash = item['record']['mesh_sha256'][:16]
    fig.text(.055,.084,f"Source archive SHA256 {short_hash}… · {len(item['polygons']):,} section polygons",fontsize=9,color=SUB)
    fig.text(.055,.056,'CV = population standard deviation / mean. Low CV is uniformity; it does not prove correct geometry or requested sizing.',fontsize=9,color=SUB)
    fig.text(.055,.031,'Mixed quality uses the native fixed subtet diagnostic. Straight Tet10 uses corners. Native geometry/audit receipts are separate.',fontsize=8.6,color=SUB)
    save(fig,output/(slug(item['case']['name'])+'.png'),figures)


def gallery(items, spec, output, figures):
    mode = spec.get('gallery_mode','cutaways')
    cols = int(spec.get('columns',min(3,len(items))))
    rows = math.ceil(len(items)/cols)
    height = 4.8*rows+1.7
    fig = plt.figure(figsize=(5.2*cols,height),facecolor=BG)
    fig.text(.047,1-.46/height,spec.get('title','PAACS / Actual mesh validation'),fontsize=22,weight='bold')
    subtitle = spec.get('subtitle', 'Actual mesh cutaways · Cyan exterior / gold interior · Measured cells and quality')
    fig.text(.047,1-.84/height,subtitle,fontsize=11,color=SUB)
    grid = fig.add_gridspec(rows,cols,left=.047,right=.925 if mode == 'sections' else .975,
                           bottom=.30 if mode == 'sections' else .2,
                           top=1-1.43/height,hspace=.42,wspace=.23)
    norm = None
    if mode == 'sections':
        values = np.concatenate([np.cbrt(i['volumes'])[i['section_cells']] for i in items])
        norm = LogNorm(float(values.min())*.999,float(values.max())*1.001)
    for index,item in enumerate(items):
        ax = fig.add_subplot(grid[index//cols,index%cols],
                             projection='3d' if mode != 'sections' else None,facecolor=BG)
        if mode == 'sections':
            section_axes(ax,item,norm)
        else:
            cutaway(ax,item)
        ax.set_title(item['case']['name'],loc='left',fontsize=12,weight='bold',pad=12)
        ax.text2D(0,-.04,stats_text(item),transform=ax.transAxes,fontsize=9.4,color=SUB,
                  linespacing=1.5) if mode != 'sections' else \
            ax.text(0,-.16,stats_text(item),transform=ax.transAxes,fontsize=9,color=SUB,
                    linespacing=1.5,va='top')
    if mode == 'sections':
        cax = fig.add_axes([.948,.24,.012,.49])
        bar = size_bar(fig,plt.cm.ScalarMappable(norm=norm,cmap='YlGnBu'),cax=cax)
        bar.set_label('Actual cell size: cube root of volume',fontsize=10)
    fig.text(.047,.071,'Sections intersect original cells; four-sided polygons can be tetrahedron cuts. No replacement mesh or invented elements.',fontsize=9.5,color=SUB)
    footer = ('Sizing controls intentionally vary cell sizes. Mixed quality uses native SubTetMeanRatio; straight Tet10 remains faceted.'
              if mode == 'sections' else
              'True coordinate aspect is retained. Lower CV means more uniform volumes; mixed quality uses native SubTetMeanRatio.')
    fig.text(.047,.04,spec.get('footer',footer),fontsize=9.5,color=SUB)
    save(fig,output/'gallery.png',figures)


def comparison(items, output, figures):
    fig,axes = plt.subplots(1,3,figsize=(15.5,4.9),facecolor=BG)
    fig.subplots_adjust(left=.055,right=.97,bottom=.2,top=.76,wspace=.34)
    fig.text(.055,.92,'Measured element distributions',fontsize=19,weight='bold')
    colors = ['#cc9455','#268d9e','#264f7a','#789049','#8856a3']
    values = [i['volumes']/i['volumes'].mean() for i in items]
    high = max(float(v.max()) for v in values)
    bins = np.linspace(0,high*1.015,56)
    for index,item in enumerate(items):
        v = values[index]
        color = colors[index%len(colors)]
        weights = np.full(len(v),100/len(v))
        axes[0].hist(v,bins=bins,weights=weights,histtype='step',lw=1.8,label=item['case']['name'],color=color)
        sizes = np.cbrt(item['volumes']); sizes /= sizes.mean()
        axes[1].hist(sizes,bins=histogram_bins(sizes),weights=weights,histtype='step',lw=1.8,color=color)
        q = item['quality'][np.isfinite(item['quality'])]
        if len(q):
            axes[2].hist(q,bins=np.linspace(0,1,41),weights=np.full(len(q),100/len(q)),
                         histtype='step',lw=1.8,color=color)
    for ax,title in zip(axes,['Volume / mean volume','∛volume / mean ∛volume','Tet mean ratio (ideal = 1)']):
        simple_axes(ax); ax.set(xlabel=title,ylabel='Cells (%)')
    axes[0].legend(frameon=False,fontsize=8)
    axes[2].set_xlim(0,1)
    fig.text(.055,.066,'Cell-weighted distributions. Same requested controls must be established by the manifest; equal pictures alone do not establish equivalence.',fontsize=9,color=SUB)
    save(fig,output/'distributions.png',figures)


def geometry_gallery(cases, base, output, figures):
    cols,rows = min(3,len(cases)), math.ceil(len(cases)/3)
    fig = plt.figure(figsize=(cols*5.2,rows*4.8+1.65),facecolor=BG)
    fig.text(.047,.965,'CAD cases awaiting supported volume meshing',fontsize=22,weight='bold')
    fig.text(.047,.922,'Actual OCCT tessellation · These images show imported geometry, without claiming a PAACS volume mesh',fontsize=11,color=SUB)
    grid = fig.add_gridspec(rows,cols,left=.047,right=.97,top=.855,bottom=.145,hspace=.42,wspace=.18)
    records=[]
    for i,case in enumerate(cases):
        path=resolve(base,case['geometry']); mesh=read_geometry(path)
        ax=fig.add_subplot(grid[i//cols,i%cols],projection='3d',facecolor=BG)
        polys=[mesh.points[list(c.ids[:CORNERS[c.kind]])] for c in mesh.cells if c.kind in (2,3,4)]
        ax.add_collection3d(Poly3DCollection(polys,facecolors=CYAN,edgecolors=EDGE,linewidths=.17))
        frame3d(ax,mesh)
        ax.set_title(case['name'],loc='left',fontsize=12,weight='bold')
        reason='\n'.join(textwrap.wrap(case.get('reason','unsupported / unqualified'),width=44))
        ax.text2D(0,-.035,f"{case.get('status','Refused')} · {len(polys):,} tessellation faces\n{reason}",
                  transform=ax.transAxes,fontsize=9,color=SUB,linespacing=1.5)
        record={'name':case['name'],'geometry_path':str(path),'geometry_sha256':sha256(path),
                'tessellation_faces':len(polys),'status':case.get('status'),
                'reason':case.get('reason'),'metadata':case.get('metadata',{})}
        if case.get('source'):
            source=resolve(base,case['source']); record.update(source_path=str(source),source_sha256=sha256(source))
        records.append(record)
    fig.text(.047,.061,'OCCT display tessellation is an approximation to CAD surfaces, not a PAACS mesh or a geometric error certificate.',fontsize=9.5,color=SUB)
    save(fig,output/'geometry-overview.png',figures)
    return records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(),indent=2)); return
    if not args.manifest or not args.output:
        parser.error('--manifest and --output are required unless --self-test is set')
    manifest=args.manifest.resolve(); output=args.output.resolve(); output.mkdir(parents=True,exist_ok=True)
    spec=json.loads(manifest.read_text(encoding='utf-8-sig'))
    items=[load_case(c,manifest.parent) for c in spec.get('cases',[])]
    figures=[]
    for item in items:
        if item['case'].get('report',True):
            individual(item,output,figures)
    if items:
        gallery(items,spec,output,figures)
        if spec.get('distributions',False):
            comparison(items,output,figures)
    geometry=geometry_gallery(spec['geometry_cases'],manifest.parent,output,figures) if spec.get('geometry_cases') else []
    receipt={'schema':'paacs-mesh-visual-report-v1','manifest_path':str(manifest),
             'manifest_sha256':sha256(manifest),'python':sys.version,
             'numpy':np.__version__,'matplotlib':matplotlib.__version__,
             'renderer_sha256':sha256(Path(__file__)),
             'data_module_sha256':sha256(Path(__file__).with_name('mesh_visual_data.py')),
             'self_tests':self_test(),'meshes':[i['record'] for i in items],
             'geometry_cases':geometry,'figures':figures,
             'scope':'Rendered actual retained cells and their convex planar intersections; Tet10 verified straight-sided. Native audits/geometry verification are separate receipts.'}
    (output/'provenance.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'meshes':len(items),'geometry_cases':len(geometry),'figures':len(figures),
                      'provenance':str(output/'provenance.json')},indent=2))


if __name__=='__main__':
    try:
        main()
    except (ValueError,KeyError,IndexError,OSError) as error:
        print(f'mesh visual report: {error}',file=sys.stderr)
        sys.exit(1)
