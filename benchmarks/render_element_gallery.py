#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Actual archive geometry, element choices and cone quality. No synthetic meshes."""
import argparse,json,struct,subprocess,tempfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
p=argparse.ArgumentParser();p.add_argument('executable',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
arities=[2,3,3,6,4,4,10,6,8,5]
faces={5:[[1,2,3],[0,3,2],[0,1,3],[0,2,1]],6:[[1,2,3],[0,3,2],[0,1,3],[0,2,1]],
8:[[0,3,2,1],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]],
9:[[0,3,2,1],[0,1,4],[1,2,4],[2,3,4],[3,0,4]]}
def read(path):
    data=path.read_bytes();assert data[:8]==b'PAACSM8\0';n,c,_=struct.unpack_from('<QQQ',data,16)
    pts=np.array([struct.unpack_from('<ddd',data,40+32*i) for i in range(n)]);cells=[];at=40+32*n
    for i in range(c):
        kind,dim,tag=struct.unpack_from('<B3xii',data,at);ids=struct.unpack_from('<'+'q'*arities[kind],data,at+12)
        cells.append((kind,dim,tag,ids));at+=12+8*arities[kind]
    assert at==len(data);return pts,cells
cases=[('cylinder','Cylinder · Hex8','.3',[]),('cylinder','Cylinder · Tet4','.3',['--elements','tet4']),
       ('cylinder','Cylinder · Tet10','.3',['--elements','tet10']),
       ('cone','Cone · swept Tet4','.15',['--elements','tet4','--cone-strategy','swept']),
       ('cone','Cone · graded Tet4','.15',['--elements','tet4']),
       ('cone','Sharp cone · quality warnings','.3',['--radius','.3','--height','3','--elements','tet4'])]
fig=plt.figure(figsize=(15,10),facecolor='#f3f6f9');records=[]
fig.suptitle('PAACS 1.0.0-rc.5   /   Element choice & sharp tips',x=.04,y=.975,ha='left',fontsize=23,weight='bold',color='#18334b')
fig.text(.04,.93,'Actual volume cells • complete boundary audit • Tet10 Jacobian certificates • C++26 / GCC 16',fontsize=12,color='#526478')
with tempfile.TemporaryDirectory() as temp:
 for i,(shape,title,size,extra) in enumerate(cases):
    path=Path(temp)/(str(i)+'.paacsmesh');run=subprocess.run([str(a.executable.resolve()),shape,size,str(path),*extra],capture_output=True,text=True,check=True)
    result=json.loads(run.stdout);records.append(dict(title=title,size=size,options=extra,result=result))
    pts,cells=read(path);boundary={}
    for kind,dim,tag,ids in cells:
        if dim==2:boundary[tuple(sorted(ids[:3] if kind==3 else ids))]=tag
    exposed={}
    for ordinal,(kind,dim,tag,ids) in enumerate(cells):
        if dim!=3:continue
        corners=ids[:4] if kind in [5,6] else ids
        if pts[list(corners),0].mean()>1e-12:continue
        quality=1
        if kind in [5,6]:
            q=pts[list(corners)];vol=np.linalg.det((q[1:]-q[0]).T)/6
            quality=12*(3*vol)**(2/3)/sum(np.sum((q[x]-q[y])**2) for x in range(4) for y in range(x))
        for face in faces[kind]:
            conn=tuple(ids[j] for j in face);key=tuple(sorted(conn))
            if key in exposed:del exposed[key]
            else:exposed[key]=(conn,boundary.get(key,99),quality,kind,ids)
    polygons=[];colors=[];midpoints=[];az=-55;el=27
    view=np.array([np.cos(np.deg2rad(el))*np.cos(np.deg2rad(az)),np.cos(np.deg2rad(el))*np.sin(np.deg2rad(az)),np.sin(np.deg2rad(el))])
    light=np.array([.5,-.7,1]);light/=np.linalg.norm(light)
    for conn,tag,q,kind,ids in exposed.values():
        poly=pts[list(conn)];normal=np.cross(poly[1]-poly[0],poly[2]-poly[0]);normal/=np.linalg.norm(normal)
        if normal@view<=0:continue
        color='#db6961' if q<.1 else '#e0ad65' if tag==99 else '#6dbfc9'
        polygons.append(poly);colors.append(np.array(to_rgb(color))*(.68+.32*max(0,normal@light)))
        if kind==6:
            edges=[(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]
            for e,(u,v) in enumerate(edges):
                if ids[u] in conn and ids[v] in conn:midpoints.append(ids[4+e])
    ax=fig.add_subplot(2,3,i+1,projection='3d',facecolor='#f3f6f9')
    ax.add_collection3d(Poly3DCollection(polygons,facecolors=colors,edgecolors='#224659',linewidths=.22,antialiased=True,zsort='average'))
    if midpoints:
        mid=pts[sorted(set(midpoints))];ax.scatter(mid[:,0],mid[:,1],mid[:,2],s=2,c='#17394d',depthshade=False)
    low=pts.min(axis=0);high=pts.max(axis=0);ax.set_xlim(low[0],high[0]);ax.set_ylim(low[1],high[1]);ax.set_zlim(low[2],high[2]);ax.set_box_aspect(high-low)
    ax.view_init(elev=el,azim=az);ax.set_axis_off();ax.set_title(title,loc='left',fontsize=14,color='#18334b',weight='bold',pad=3)
    total=sum(result[k] for k in ['hexes','pyramids','tet4s','tet10s'])
    ax.text2D(0,-.025,f'{total:,} cells  ·  {result["nodes"]:,} nodes  ·  min. quality {result["minimum_quality"]:.3f}',transform=ax.transAxes,fontsize=9.5,color='#526478')
    ax.text2D(0,-.092,f'{result["low_quality_cells"]:,} cells below 0.1  |  volume error {100*result["relative_volume_error"]:.3f}%',transform=ax.transAxes,fontsize=9,color='#526478')
fig.text(.04,.022,'Gold = exposed interior. Red = quality < 0.1. Tet10 dots show shared midside nodes; geometry is straight-sided. Cone comparisons use equal requested sizes.',fontsize=9.5,color='#526478')
fig.subplots_adjust(left=.04,right=.97,bottom=.115,top=.875,wspace=.12,hspace=.31)
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=180,facecolor=fig.get_facecolor());a.output.with_suffix('.json').write_text(json.dumps(records,indent=2)+'\n');print(a.output)
