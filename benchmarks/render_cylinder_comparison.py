#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Render actual mesh archives and independently measure their corner tets."""
import argparse,json,struct,subprocess,tempfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def read(path):
    data=path.read_bytes();assert data[:8]==b'PAACSM8\0'
    n,c,_=struct.unpack_from('<QQQ',data,16)
    pts=np.array([struct.unpack_from('<ddd',data,40+32*i) for i in range(n)])
    cells=[];at=40+32*n;arities=[2,3,3,6,4,4,10,6,8,5]
    for i in range(c):
        kind,dim,tag=struct.unpack_from('<B3xii',data,at)
        ids=struct.unpack_from('<'+'q'*arities[kind],data,at+12)
        cells.append((kind,dim,tag,ids));at+=12+8*arities[kind]
    assert at==len(data)
    return pts,cells

def statistics(pts,cells):
    ids=np.array([ids[:4] for kind,dim,tag,ids in cells if dim==3])
    p=pts[ids];v=np.linalg.det(p[:,1:]-p[:,0,None])/6
    assert np.all(v>0)
    edge2=np.stack([np.sum((p[:,i]-p[:,j])**2,axis=1) for i in range(4) for j in range(i)],axis=1)
    quality=12*np.power(3*v,2/3)/np.sum(edge2,axis=1)
    equivalent=np.cbrt(v)
    metrics=dict(volume_cv=float(v.std()/v.mean()),volume_max_min=float(v.max()/v.min()),
                 size_cv=float(equivalent.std()/equivalent.mean()),
                 quality_min=float(quality.min()),quality_p10=float(np.quantile(quality,.1)),
                 quality_median=float(np.median(quality)),
                 edge_length_cell_weighted_quantiles=np.quantile(np.sqrt(edge2),[.1,.5,.9]).tolist(),
                 volume_quantiles=np.quantile(v,[.1,.5,.9]).tolist())
    return metrics,v/v.mean()

def surface(ax,pts,cells,high):
    boundary={tuple(sorted(ids[:3])) for kind,dim,tag,ids in cells if dim==2}
    exposed={};faces=[[1,2,3],[0,3,2],[0,1,3],[0,2,1]]
    for kind,dim,tag,ids in cells:
        if dim!=3 or pts[list(ids[:4]),0].mean()>1e-12:continue
        for f in faces:
            conn=tuple(ids[j] for j in f);key=tuple(sorted(conn))
            if key in exposed:del exposed[key]
            else:exposed[key]=(conn,key in boundary,ids)
    az=-55;el=30
    view=np.array([np.cos(np.deg2rad(el))*np.cos(np.deg2rad(az)),np.cos(np.deg2rad(el))*np.sin(np.deg2rad(az)),np.sin(np.deg2rad(el))])
    light=np.array([.5,-.7,1]);light/=np.linalg.norm(light)
    polys=[];colors=[];mids=set()
    for conn,bound,ids in exposed.values():
        poly=pts[list(conn)];normal=np.cross(poly[1]-poly[0],poly[2]-poly[0]);normal/=np.linalg.norm(normal)
        if normal@view<=0:continue
        polys.append(poly);colors.append(np.array(to_rgb('#64bcc6' if bound else '#e2ae67'))*(.7+.3*max(0,normal@light)))
        if high:
            for e,(a,b) in enumerate([(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]):
                if ids[a] in conn and ids[b] in conn:mids.add(ids[4+e])
    ax.add_collection3d(Poly3DCollection(polys,facecolors=colors,edgecolors='#25485b',linewidths=.25))
    if mids:
        p=pts[sorted(mids)];ax.scatter(p[:,0],p[:,1],p[:,2],s=2,c='#17394d',depthshade=False)
    ax.set_xlim(-1,1);ax.set_ylim(-1,1);ax.set_zlim(0,2);ax.set_box_aspect((2,2,2))
    ax.view_init(el,az);ax.set_axis_off()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('executable',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--size',type=float,default=.2);a=parser.parse_args()
    fig=plt.figure(figsize=(15,8.8),facecolor='#f3f6f9');records=[];volumes=[]
    fig.suptitle('PAACS rc.6  /  More uniform cylinder tetrahedra',x=.045,y=.97,ha='left',fontsize=22,weight='bold',color='#18334b')
    fig.text(.045,.919,f'Radius 1 · Height 2 · Requested size {a.size:g} · Same geometry tolerances · Actual mesh cutaways',fontsize=12,color='#526478')
    grid=fig.add_gridspec(2,3,height_ratios=[3.6,1.2],left=.045,right=.97,bottom=.13,top=.86,hspace=.53,wspace=.12)
    with tempfile.TemporaryDirectory() as temp:
        for i,(family,strategy,title) in enumerate([('tet4','swept','Previous swept Tet4'),('tet4','uniform','Uniform Tet4'),('tet10','uniform','Uniform Tet10')]):
            path=Path(temp)/(str(i)+'.paacsmesh')
            run=subprocess.run([str(a.executable.resolve()),'cylinder',str(a.size),str(path),'--elements',family,'--cylinder-strategy',strategy],capture_output=True,text=True,check=True)
            result=json.loads(run.stdout);pts,cells=read(path);metrics,volume=statistics(pts,cells)
            records.append(dict(family=family,strategy=strategy,size=a.size,result=result,metrics=metrics));volumes.append(volume)
            ax=fig.add_subplot(grid[0,i],projection='3d',facecolor='#f3f6f9');surface(ax,pts,cells,family=='tet10')
            ax.set_title(title,loc='left',fontsize=15,weight='bold',color='#18334b')
            ax.text2D(0,-.04,f'{result[family+"s"]:,} cells · {result["nodes"]:,} nodes',transform=ax.transAxes,fontsize=11,color='#526478')
            ax.text2D(0,-.12,f'Min. quality {metrics["quality_min"]:.3f} · Volume CV {100*metrics["volume_cv"]:.1f}%',transform=ax.transAxes,fontsize=11,color='#526478')
        ax=fig.add_subplot(grid[1,:2],facecolor='#f3f6f9')
        bins=np.linspace(0,max(v.max() for v in volumes[:2])*1.04,55)
        for v,label,color in zip(volumes[:2],['Previous swept Tet4','Uniform Tet4'],['#db8b56','#2a8e9c']):
            ax.hist(v,bins=bins,weights=np.full(len(v),100/len(v)),histtype='step',linewidth=2,label=label,color=color)
        ax.set_xlabel('Cell volume / mean cell volume (closer to 1 = more uniform)',fontsize=10)
        ax.set_ylabel('Cells (%)');ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False,fontsize=9)
        ax=fig.add_subplot(grid[1,2]);ax.set_axis_off()
        ax.text(0,.95,'What this proves',weight='bold',fontsize=12,color='#18334b')
        ax.text(0,.68,'More consistent cell sizes in this cylinder.\nPositive cells and a complete boundary audit.\nTet10 shares midsides; geometry stays faceted.\nThis is a structured sweep, not general CAD.',va='top',linespacing=1.7,fontsize=10,color='#526478')
    fig.text(.045,.035,'Gold = exposed interior · Cyan = exterior · Tet10 dots = midside nodes · CV = standard deviation / mean volume. Lower CV is more uniform.',fontsize=9.5,color='#526478')
    a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=170,facecolor=fig.get_facecolor())
    a.output.with_suffix('.json').write_text(json.dumps(records,indent=2)+'\n');print(a.output)
if __name__=='__main__':main()
