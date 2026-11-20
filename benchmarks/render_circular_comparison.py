#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Actual archive cutaways; no synthetic geometry or illustrated elements."""
import argparse,json,subprocess,tempfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from render_cylinder_comparison import read,statistics

def surface(ax,pts,cells,quadratic,shape):
    boundary={tuple(sorted(ids[:3])) for kind,dim,tag,ids in cells if dim==2}
    exposed={};faces=[[1,2,3],[0,3,2],[0,1,3],[0,2,1]]
    for kind,dim,tag,ids in cells:
        if dim!=3 or pts[list(ids[:4]),2].mean()>1e-12:continue
        for f in faces:
            conn=tuple(ids[j] for j in f);key=tuple(sorted(conn))
            if key in exposed:del exposed[key]
            else:exposed[key]=(conn,key in boundary,ids)
    az=-58 if shape=='torus' else -48;el=44
    view=np.array([np.cos(np.deg2rad(el))*np.cos(np.deg2rad(az)),np.cos(np.deg2rad(el))*np.sin(np.deg2rad(az)),np.sin(np.deg2rad(el))])
    light=np.array([.5,-.7,1]);light/=np.linalg.norm(light)
    polygons=[];colors=[];mids=set()
    for conn,bound,ids in exposed.values():
        poly=pts[list(conn)];normal=np.cross(poly[1]-poly[0],poly[2]-poly[0]);normal/=np.linalg.norm(normal)
        if normal@view<=0:continue
        polygons.append(poly)
        colors.append(np.array(to_rgb('#70c5cc' if bound else '#eccb87'))*(.73+.27*max(0,normal@light)))
        if quadratic:
            for e,(a,b) in enumerate([(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]):
                if ids[a] in conn and ids[b] in conn:mids.add(ids[4+e])
    ax.add_collection3d(Poly3DCollection(polygons,facecolors=colors,edgecolors='#2d4f60',linewidths=.18))
    if mids:
        p=pts[sorted(mids)];ax.scatter(p[:,0],p[:,1],p[:,2],s=.7,c='#234454',depthshade=False)
    if shape=='torus':ax.set_xlim(-3,3);ax.set_ylim(-3,3);ax.set_box_aspect((6,6,2))
    else:ax.set_xlim(0,3);ax.set_ylim(0,3);ax.set_box_aspect((3,3,2))
    ax.set_zlim(-1,1);ax.view_init(el,az);ax.set_axis_off()

def main():
    p=argparse.ArgumentParser();p.add_argument('executable',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--size',type=float,default=.2);args=p.parse_args()
    fig=plt.figure(figsize=(17,11.5),facecolor='#f3f6f9')
    fig.suptitle('PAACS rc.7  /  Balanced tetrahedra around curves',x=.045,y=.97,ha='left',fontsize=24,weight='bold',color='#18334b')
    fig.text(.045,.924,f'Solid torus and 90° elbow · Minor radius 1 · Major radius 2 · Requested size {args.size:g}',fontsize=13,color='#526478')
    fig.text(.045,.893,'Actual mesh cutaways · Cyan: exterior · Gold: exposed interior · Tet10 dots: shared midside nodes',fontsize=11,color='#526478')
    grid=fig.add_gridspec(2,3,left=.04,right=.97,bottom=.18,top=.86,hspace=.42,wspace=.12)
    records=[]
    with tempfile.TemporaryDirectory() as temp:
        for row,shape in enumerate(['torus','elbow']):
            for col,(family,strategy,title) in enumerate([('tet4','swept','Previous swept Tet4'),('tet4','uniform','Balanced Tet4'),('tet10','uniform','Balanced Tet10')]):
                archive=Path(temp)/(shape+family+strategy+'.paacsmesh')
                run=subprocess.run([str(args.executable.resolve()),shape,str(args.size),str(archive),'--elements',family,'--circular-strategy',strategy,'--detect-features','on'],capture_output=True,text=True,check=True)
                result=json.loads(run.stdout);pts,cells=read(archive);stats,_=statistics(pts,cells)
                records.append(dict(shape=shape,family=family,strategy=strategy,size=args.size,result=result,metrics=stats))
                ax=fig.add_subplot(grid[row,col],projection='3d',facecolor='#f3f6f9')
                surface(ax,pts,cells,family=='tet10',shape)
                ax.set_title(f'{shape.title()}  /  {title}',loc='left',fontsize=13,weight='bold',color='#18334b')
                ax.text2D(0,-.02,f'{result[family+"s"]:,} cells · {result["nodes"]:,} nodes',transform=ax.transAxes,fontsize=11,color='#526478')
                ax.text2D(0,-.105,f'Min. quality {stats["quality_min"]:.3f} · Volume CV {100*stats["volume_cv"]:.1f}%',transform=ax.transAxes,fontsize=11,color='#526478')
                print(shape,family,strategy,stats['quality_min'],stats['volume_cv'],flush=True)
    fig.text(.045,.080,'Lower volume CV means more consistent cell sizes. Higher mean-ratio quality means less distorted tetrahedra (ideal = 1).',fontsize=11,color='#526478')
    fig.text(.045,.050,'Equal requested size and tolerances; cell counts differ. Bends still compress inner cells. Tet10 geometry remains faceted.',fontsize=11,color='#526478')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=150,facecolor=fig.get_facecolor());plt.close(fig)
    args.output.with_suffix('.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
