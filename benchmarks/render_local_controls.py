#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Geometric cuts through exported Tet4/Prism6/Hex8 cells, never invented meshes."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm, to_rgb
from matplotlib.patches import Rectangle, Circle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from render_cylinder_comparison import read, statistics

BG = '#f4f7fa'
INK = '#19374b'
SUB = '#516879'
EDGES = {5: [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],
         7: [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)],
         8: [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]}

def section(points, cells, plane=.5):
    polygons, sizes, heights, seen = [], [], [], set()
    for kind in EDGES:
        ids = np.array([v for t,d,_,v in cells if d == 3 and t == kind], dtype=np.int64)
        if not len(ids):
            continue
        p = points[ids]
        selected = np.flatnonzero((p[:,:,1].min(axis=1) <= plane+1e-12) &
                                  (p[:,:,1].max(axis=1) >= plane-1e-12))
        for i in selected:
            vertices = []
            for a,b in EDGES[kind]:
                pa,pb = p[i,a],p[i,b]
                da,db = pa[1]-plane,pb[1]-plane
                if abs(da) < 1e-12:
                    vertices.append(pa[[0,2]])
                if da*db < 0:
                    vertices.append((pa+(pb-pa)*(-da)/(db-da))[[0,2]])
                if abs(db) < 1e-12:
                    vertices.append(pb[[0,2]])
            if len(vertices)<3:
                continue
            vertices=np.unique(np.round(vertices,13),axis=0)
            if len(vertices)<3:
                continue
            center=vertices.mean(axis=0)
            vertices=vertices[np.argsort(np.arctan2(vertices[:,1]-center[1],vertices[:,0]-center[0]))]
            area=.5*abs(np.sum(vertices[:,0]*np.roll(vertices[:,1],-1)-vertices[:,1]*np.roll(vertices[:,0],-1)))
            if area<1e-14:
                continue
            key=tuple(sorted(map(tuple,vertices)))
            if key in seen:
                continue
            seen.add(key)
            polygons.append(vertices)
            heights.append(float(p[i,:,2].mean()))
            if kind==5:
                vol=np.linalg.det(p[i,1:]-p[i,0])/6
                assert vol>0
                sizes.append(float(np.cbrt(vol)))
            else:
                sizes.append(float(p[i,:,2].max()-p[i,:,2].min()))
    return polygons,np.array(sizes),np.array(heights)

def clip(poly,plane=.5):
    out=[]
    for i,p in enumerate(poly):
        q=poly[(i+1)%len(poly)]
        inside=p[1]>=plane-1e-12
        if inside:
            out.append(p)
        if inside != (q[1]>=plane-1e-12):
            out.append(p+(q-p)*((plane-p[1])/(q[1]-p[1])))
    return np.array(out)

def cutaway(ax,points,cells,polygons):
    az,el=-57,26
    view=np.array([np.cos(np.deg2rad(el))*np.cos(np.deg2rad(az)),
                   np.cos(np.deg2rad(el))*np.sin(np.deg2rad(az)),np.sin(np.deg2rad(el))])
    light=np.array([.6,-.7,1.]);light/=np.linalg.norm(light)
    faces,colors=[],[]
    for _,dim,_,ids in cells:
        if dim!=2:
            continue
        poly=clip(points[list(ids)])
        if len(poly)<3:
            continue
        normal=np.cross(poly[1]-poly[0],poly[2]-poly[0])
        length=np.linalg.norm(normal)
        if length<1e-14:
            continue
        normal/=length
        if normal@view<=0:
            continue
        faces.append(poly)
        colors.append(np.array(to_rgb('#78c4cf'))*(.76+.24*max(0,normal@light)))
    for p in polygons:
        faces.append(np.column_stack((p[:,0],np.full(len(p),.5),p[:,1])))
        colors.append(to_rgb('#edc985'))
    ax.add_collection3d(Poly3DCollection(faces,facecolors=colors,edgecolors='#344d5e',linewidths=.22))
    ax.set(xlim=(0,1),ylim=(0,1),zlim=(0,1))
    ax.set_box_aspect((1,1,1));ax.view_init(el,az);ax.set_axis_off()

def annotate_region(ax,kind):
    if kind=='cube':
        ax.add_patch(Rectangle((.375,.375),.25,.25,fill=False,edgecolor='#ae3344',ls='--',lw=1.6))
    elif kind=='sphere':
        ax.add_patch(Circle((.5,.5),.16,fill=False,edgecolor='#ae3344',ls='--',lw=1.6))
    else:
        ax.plot([0,0],[0,1],color='#ae3344',lw=3,clip_on=False)

def axis2d(ax):
    ax.set(xlim=(0,1),ylim=(0,1),aspect='equal',xticks=[0,.5,1],yticks=[0,.5,1])
    ax.tick_params(labelsize=8,colors=SUB)
    for s in ax.spines.values():
        s.set_color('#91a2af')

def sizing_gallery(folder,output):
    fig=plt.figure(figsize=(16,15),facecolor=BG)
    fig.suptitle('PAACS rc.8  /  Local size controls through the volume',x=.055,y=.974,ha='left',fontsize=23,weight='bold',color=INK)
    fig.text(.055,.94,'Unit cube · Background edge cap 0.25 · Local edge cap 0.0625 · Actual exported Tet4 meshes',fontsize=12,color=SUB)
    fig.text(.055,.916,'Top: cutaways at field growth 1.3. Below: exact y = 0.5 sections at two growth settings.',fontsize=11,color=SUB)
    grid=fig.add_gridspec(3,3,left=.055,right=.91,bottom=.13,top=.875,hspace=.3,wspace=.18,height_ratios=[1,1,1])
    norm=LogNorm(vmin=.0625/np.cbrt(24),vmax=.25/np.cbrt(24))
    cmap=plt.get_cmap('YlOrBr')
    records=[]
    for col,(kind,title) in enumerate([('cube','Finer cube inside'),('sphere','Spherical influence region'),('face','Refine the x-min face')]):
        for grow,row in [(130,1),(170,2)]:
            path=folder/f'{kind}-{grow}.paacsmesh';points,cells=read(path)
            polygons,sizes,_=section(points,cells)
            stats,_=statistics(points,cells)
            count=sum(d==3 for _,d,_,_ in cells)
            records.append(dict(case=kind,growth=grow/100,volume_cells=count,nodes=len(points),
                                archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),metrics=stats))
            if row==1:
                ax=fig.add_subplot(grid[0,col],projection='3d',facecolor=BG)
                cutaway(ax,points,cells,polygons)
                ax.set_title(title,loc='left',fontsize=13,weight='bold',color=INK,pad=5)
                ax.text2D(0,-.02,f'{count:,} tetrahedra · min. quality {stats["quality_min"]:.3f}',transform=ax.transAxes,fontsize=10,color=SUB)
            ax=fig.add_subplot(grid[row,col],facecolor=BG)
            ax.add_collection(PolyCollection(polygons,facecolors=cmap(norm(sizes)),edgecolors='#624f38',linewidths=.22))
            annotate_region(ax,kind);axis2d(ax)
            ax.set_title(f'Field growth {grow/100:g}  ·  {count:,} cells',fontsize=11,loc='left',color=INK)
            if col==0:
                ax.set_ylabel('z',color=SUB)
            if row==2:
                ax.set_xlabel('x',color=SUB)
    cax=fig.add_axes([.928,.19,.012,.4])
    bar=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
    bar.set_label('Actual cell size: cube root of volume',fontsize=10,color=SUB)
    bar.set_ticks([.0625/np.cbrt(24),.125/np.cbrt(24),.25/np.cbrt(24)])
    bar.set_ticklabels(['0.0217','0.0433','0.0867']);bar.minorticks_off();bar.ax.tick_params(labelsize=9)
    fig.text(.055,.078,'Red outline = sizing influence region / selected face, not a material interface. No holes or hanging faces at refinement transitions.',fontsize=10.5,color=SUB)
    fig.text(.055,.052,'Smaller field growth spreads refinement farther. Discrete cell sizes are quantized; measured neighbor cbrt(volume) ratio is about 1.260 here.',fontsize=10.5,color=SUB)
    fig.text(.055,.027,'Quality is tetrahedron mean ratio (ideal 1). Cyan = exterior; gold = exposed interior. Cuts intersect the exported elements geometrically.',fontsize=10,color=SUB)
    fig.savefig(output,dpi=150,facecolor=BG);plt.close(fig)
    output.with_suffix('.json').write_text(json.dumps(records,indent=2)+'\n')

def layer_gallery(folder,output):
    fig=plt.figure(figsize=(16,10.6),facecolor=BG)
    fig.suptitle('PAACS rc.8  /  Wall layers into a coarser interior',x=.055,y=.967,ha='left',fontsize=23,weight='bold',color=INK)
    fig.text(.055,.925,'First thickness 0.01 · 5 prescribed wall layers · Growth 1.2 · Total depth 1 · Core thickness cap 0.15',fontsize=11.5,color=SUB)
    grid=fig.add_gridspec(2,3,left=.055,right=.97,bottom=.16,top=.86,hspace=.22,wspace=.15)
    records=[]
    for col,(kind,title) in enumerate([('structured','Structured triangles → wedges'),('unstructured','Unstructured triangles → wedges'),('mixed','Mixed triangles + quads → wedges + hexes')]):
        path=folder/f'{kind}.paacsmesh';points,cells=read(path);polygons,sizes,heights=section(points,cells)
        counts={str(t):sum(kind==t and d==3 for kind,d,_,_ in cells) for t in [7,8]}
        ax=fig.add_subplot(grid[0,col],projection='3d',facecolor=BG);cutaway(ax,points,cells,polygons)
        ax.set_title(title,fontsize=11.5,weight='bold',loc='left',color=INK)
        ax.text2D(0,-.035,f'{counts["7"]:,} wedges  ·  {counts["8"]:,} hexes  ·  {len(points):,} nodes',transform=ax.transAxes,fontsize=10,color=SUB)
        ax=fig.add_subplot(grid[1,col],facecolor=BG)
        colors=['#78c4cf' if h<.074416-1e-12 else '#edc985' for h in heights]
        ax.add_collection(PolyCollection(polygons,facecolors=colors,edgecolors='#344d5e',linewidths=.45));axis2d(ax)
        ax.axhline(.074416,color='#ae3344',ls='--',lw=1.2)
        ax.set_title('Actual section at y = 0.5',loc='left',fontsize=11,color=INK);ax.set_xlabel('x',color=SUB)
        if col==0:
            ax.set_ylabel('z / depth from wall',color=SUB)
        records.append(dict(case=kind,wedges=counts['7'],hexes=counts['8'],nodes=len(points),archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),section_polygons=len(polygons)))
    fig.text(.055,.096,'Section colors: cyan = 5 wall layers; gold = 12 transition/core layers. Dashed line marks the prescribed wall stack height 0.074416.',fontsize=10.5,color=SUB)
    fig.text(.055,.065,'Prism6 and Wedge6 name the same six-node element. An unstructured source still produces an ordered stack in the extrusion direction.',fontsize=10.5,color=SUB)
    fig.text(.055,.034,'This interior continues the wedge/hex columns. Curved wall collisions and connection to a separate tetrahedral core remain future work.',fontsize=10.5,color=SUB)
    fig.savefig(output,dpi=150,facecolor=BG);plt.close(fig)
    output.with_suffix('.json').write_text(json.dumps(records,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sizing_archives',type=Path);p.add_argument('layer_archives',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    sizing_gallery(a.sizing_archives,a.output/'local_sizing_rc8.png')
    layer_gallery(a.layer_archives,a.output/'prism_layers_rc8.png')
    print('Wrote sizing and layer figures from exported mesh archives.')
