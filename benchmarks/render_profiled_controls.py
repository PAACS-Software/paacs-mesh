#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Render exact cuts and boundary faces from qualified PAACS mesh archives."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize, to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from render_cylinder_comparison import read
from render_local_controls import section, clip

BG, INK, SUB = '#f4f7fa', '#18384b', '#526e80'
CMAP = plt.get_cmap('YlGnBu_r')
NORM = Normalize(.12, .36)
PLANE = .013

def target(z, mode):
    z=np.asarray(z)
    h=np.full_like(z,.36,dtype=float)
    intervals={'start':[(0,0)], 'end':[(4.8,4.8)], 'band':[(2,2.8)],
               'both':[(0,0),(4.8,4.8)]}.get(mode,[])
    for lo,hi in intervals:
        h=np.minimum(h,.12+(2*(1.25-1)/(1.25+1))*np.maximum.reduce([lo-z,z-hi,np.zeros_like(z)]))
    return h

def counts(cells):
    return sum(d==3 for _,d,_,_ in cells)

def meshfile(folder, mode, family, strategy):
    return folder/f'cylinder-{mode}-{family}-{strategy}.paacsmesh'

def record(path, points, cells):
    return dict(archive=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                nodes=len(points),volume_cells=counts(cells))

def sizing(folder,output):
    fig=plt.figure(figsize=(16,9),facecolor=BG)
    fig.text(.045,.955,'PAACS rc.9  /  One reusable size profile',fontsize=23,color=INK,weight='bold')
    fig.text(.045,.91,'Select a cap or an axial band · Local 0.12 → growth 1.25 → global cap 0.36',fontsize=12,color=SUB)
    grid=fig.add_gridspec(3,3,left=.045,right=.965,bottom=.17,top=.84,wspace=.22,hspace=.72,width_ratios=[1.75,1.75,1])
    records=[]
    for row,(mode,title) in enumerate([('start','Refine the start cap'),('end','Refine the end cap'),('band','Refine a central band')]):
        stations=None
        for col,(family,label) in enumerate([('tet4','Tet4'),('hex-dominant','Hex8 from tetrahedra')]):
            path=meshfile(folder,mode,family,'volume');points,cells=read(path)
            polygons,_,heights=section(points,cells,PLANE)
            ax=fig.add_subplot(grid[row,col],facecolor='white')
            ax.add_collection(PolyCollection([p[:,::-1] for p in polygons],facecolors=CMAP(NORM(target(heights,mode))),edgecolors='#3b5867',linewidths=.28))
            ax.set(xlim=(0,4.8),ylim=(-.62,.62),aspect='equal',xticks=[0,1.2,2.4,3.6,4.8],yticks=[-.6,0,.6])
            ax.tick_params(labelsize=8,colors=SUB)
            ax.set_title(f'{title}  ·  {label}\n{counts(cells):,} volume cells',loc='left',fontsize=10.5,color=INK,pad=8)
            ax.set_xlabel('Axial position z',fontsize=9,color=SUB)
            if col==0:ax.set_ylabel('Radial x',fontsize=9,color=SUB)
            if mode=='start':ax.plot([0,0],[-.6,.6],color='#bc4151',lw=3,clip_on=False)
            if mode=='end':ax.plot([4.8,4.8],[-.6,.6],color='#bc4151',lw=3,clip_on=False)
            if mode=='band':
                for x in [2,2.8]:ax.axvline(x,color='#bc4151',ls='--',lw=1)
            records.append(record(path,points,cells))
            if col==0:stations=np.unique(points[:,2])
        ax=fig.add_subplot(grid[row,2],facecolor='white')
        z=np.linspace(0,4.8,600)
        ax.plot(z,target(z,mode),color='#21698b',lw=2,label='Requested size')
        ax.step((stations[:-1]+stations[1:])/2,np.diff(stations),where='mid',color='#c0782e',lw=1.7,label='Actual slab width')
        ax.set(xlim=(0,4.8),ylim=(0,.39),xticks=[0,2.4,4.8],yticks=[.12,.24,.36])
        ax.tick_params(labelsize=8,colors=SUB);ax.grid(alpha=.16)
        ax.set_title('Target and axial spacing',fontsize=10.5,color=INK,loc='left')
        ax.set_xlabel('Axial position z',fontsize=9,color=SUB)
        if row==0:ax.legend(fontsize=7.5,loc='lower right',frameon=False)
    cax=fig.add_axes([.16,.085,.44,.014])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cax,orientation='horizontal')
    cb.set_ticks([.12,.24,.36]);cb.ax.tick_params(labelsize=9)
    fig.text(.63,.082,'Color = requested size, not element quality',fontsize=10,color=SUB)
    fig.text(.045,.042,'Exact y = 0.013 cuts through exported meshes. Cut polygons may have more than four sides, including in a valid Hex8 mesh.',fontsize=10,color=SUB)
    fig.text(.045,.018,'Growth is guaranteed for axial slab widths. It is not a bound on every 3D edge or on neighboring equivalent-volume sizes.',fontsize=10,color=SUB)
    fig.savefig(output,dpi=160,facecolor=BG);plt.close(fig)
    return records

def cutaway(ax,points,cells,mode):
    polygons,_,heights=section(points,cells,PLANE)
    faces,colors=[],[]
    for kind,dim,_,ids in cells:
        if dim!=2:continue
        q=clip(points[list(ids)],PLANE)
        if len(q)<3:continue
        faces.append(np.column_stack([q[:,2],q[:,0],-q[:,1]]))
        colors.append(to_rgb('#87cbd2'))
    for polygon,z in zip(polygons,heights):
        faces.append(np.column_stack([polygon[:,1],polygon[:,0],np.full(len(polygon),-PLANE)]))
        colors.append(CMAP(NORM(target(z,mode))))
    ax.add_collection3d(Poly3DCollection(faces,facecolors=colors,edgecolors='#375567',linewidths=.23))
    ax.set(xlim=(0,4.8),ylim=(-.6,.6),zlim=(-.6,.6))
    ax.set_box_aspect((4.8,1.2,1.2),zoom=1.15);ax.view_init(24,-65);ax.set_proj_type('ortho');ax.set_axis_off()

def choices(folder,output):
    fig=plt.figure(figsize=(16,8.8),facecolor=BG)
    fig.text(.045,.951,'PAACS rc.9  /  Choose the mesh strategy',fontsize=23,weight='bold',color=INK)
    fig.text(.045,.907,'The same fine start-cap control drives each mesh. Exterior faces and interior cuts come from exported geometry.',fontsize=11.5,color=SUB)
    grid=fig.add_gridspec(2,3,left=.035,right=.975,bottom=.16,top=.86,hspace=.12,height_ratios=[1,1])
    records=[]
    for col,(family,strategy,label) in enumerate([('hex-dominant','sweep','Hex8 · fixed-section sweep'),('tet4','volume','Tet4 · adaptive sections'),('hex-dominant','volume','Hex8 · converted adaptive sections')]):
        path=meshfile(folder,'start',family,strategy);points,cells=read(path);records.append(record(path,points,cells))
        ax=fig.add_subplot(grid[0,col],projection='3d',facecolor=BG);cutaway(ax,points,cells,'start')
        ax.set_title(label,fontsize=12,color=INK,loc='left',weight='bold')
        ax.text2D(.06,.04,f'{counts(cells):,} volume cells  ·  {len(points):,} nodes',transform=ax.transAxes,fontsize=10,color=SUB)
        caps=grid[1,col].subgridspec(1,2,wspace=.09)
        for side,tag in enumerate([2,3]):
            selected=[points[list(ids),:2] for _,dim,t,ids in cells if dim==2 and t==tag]
            ax=fig.add_subplot(caps[0,side],facecolor=BG)
            ax.add_collection(PolyCollection(selected,facecolors=CMAP(NORM(.12 if side==0 else .36)),edgecolors='#3b5867',linewidths=.32))
            ax.set(xlim=(-.65,.65),ylim=(-.65,.65),aspect='equal');ax.set_axis_off()
            ax.set_title(('Fine start' if side==0 else 'Remote end')+f'\n{len(selected):,} boundary faces',fontsize=10,color=INK,pad=6)
    fig.text(.045,.1,'Sweep retains its fine cross-section at both ends. Volume mode coarsens the interior section while preserving a common outer polygon.',fontsize=10.5,color=SUB)
    fig.text(.045,.065,'Converted Hex8 gives four hexes per tetrahedron. It is a flexible conforming route; it can cost more memory and have poorer quality than the sweep.',fontsize=10.5,color=SUB)
    fig.text(.045,.03,'Quality warnings remain visible and can be changed to strict rejection. These are analytic-cylinder examples, not a general CAD qualification.',fontsize=10.5,color=SUB)
    fig.savefig(output,dpi=160,facecolor=BG);plt.close(fig)
    return records

def main():
    parser=argparse.ArgumentParser();parser.add_argument('archives',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    records=sizing(args.archives,args.output/'paacs-rc9-size-controls.png')
    records+=choices(args.archives,args.output/'paacs-rc9-cylinder-choices.png')
    (args.output/'profiled-gallery-manifest.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
