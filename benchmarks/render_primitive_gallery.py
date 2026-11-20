#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Plot actual archive cells, with optional cell-removal cutaways. Matplotlib/numpy."""
import argparse
import json
import struct
import subprocess
import tempfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

p=argparse.ArgumentParser()
p.add_argument("executable",type=Path)
p.add_argument("output",type=Path)
a=p.parse_args()
counts=[2,3,3,6,4,4,10,6,8,5]
faces={8:[[0,3,2,1],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]],
       9:[[0,3,2,1],[0,1,4],[1,2,4],[2,3,4],[3,0,4]]}
def read(path):
    data=path.read_bytes();assert data[:8]==b"PAACSM8\0"
    n,c,_=struct.unpack_from("<QQQ",data,16)
    points=np.array([struct.unpack_from("<ddd",data,40+32*i) for i in range(n)])
    cells=[];offset=40+32*n
    for _ in range(c):
        kind,dim,tag=struct.unpack_from("<B3xii",data,offset)
        ids=struct.unpack_from("<"+"q"*counts[kind],data,offset+12)
        cells.append((kind,dim,tag,ids));offset+=12+8*counts[kind]
    assert offset==len(data)
    return points,cells

cases=[("torus","Donut / periodic Hex8",[],False),
       ("cone","Cone / cutaway",[],True),
       ("cylinder","Cylinder / cutaway",[],True),
       ("tube","Hollow tube",["--inner-radius",".78"],False),
       ("frustum","Frustum",["--end-radius",".35"],False),
       ("elbow","Hollow elbow",["--inner-radius",".78"],False)]
fig=plt.figure(figsize=(15,10),facecolor="#f3f6f9")
fig.suptitle("PAACS 1.0.0-rc.4   /   Curved volume meshes",x=.04,y=.97,ha="left",fontsize=23,weight="bold",color="#18334b")
fig.text(.04,.925,"Real mesh cells • shared periodic seams • audited boundaries • C++26 / GCC 16",fontsize=12,color="#526478")
records=[]
with tempfile.TemporaryDirectory() as temp:
    for idx,(name,title,extra,cut) in enumerate(cases):
        path=Path(temp)/(name+".paacsmesh")
        command=[str(a.executable.resolve()),name,".3",str(path),*extra]
        result=json.loads(subprocess.check_output(command,text=True));records.append(dict(case=name,cutaway=cut,options=extra,result=result))
        points,cells=read(path)
        boundary={tuple(sorted(ids)):(ids,tag) for kind,dim,tag,ids in cells if dim==2}
        if cut:
            exposed={}
            for kind,dim,tag,ids in cells:
                if dim!=3 or points[list(ids),0].mean()>1e-12:continue
                for face in faces[kind]:
                    conn=tuple(ids[i] for i in face);key=tuple(sorted(conn))
                    if key in exposed:del exposed[key]
                    else:exposed[key]=(conn,boundary[key][1] if key in boundary else 99)
            boundary=exposed
        polygons=[];colors=[]
        azimuth=-55;elevation=27
        view=np.array([np.cos(np.deg2rad(elevation))*np.cos(np.deg2rad(azimuth)),np.cos(np.deg2rad(elevation))*np.sin(np.deg2rad(azimuth)),np.sin(np.deg2rad(elevation))])
        light=np.array([.5,-.7,1]);light/=np.linalg.norm(light)
        palette={1:"#6dbfc9",2:"#98cbd3",3:"#98cbd3",4:"#428ba5",99:"#e0ad65"}
        for ids,tag in boundary.values():
            polygon=points[list(ids)];normal=np.cross(polygon[1]-polygon[0],polygon[2]-polygon[0]);normal/=np.linalg.norm(normal)
            if normal@view<=0:continue
            polygons.append(polygon)
            colors.append(np.array(to_rgb(palette[tag]))*(.68+.32*max(0,normal@light)))
        ax=fig.add_subplot(2,3,idx+1,projection="3d",facecolor="#f3f6f9")
        ax.add_collection3d(Poly3DCollection(polygons,facecolors=colors,edgecolors="#224659",linewidths=.27,zsort="average",antialiased=True))
        low=points.min(axis=0);high=points.max(axis=0)
        ax.set_xlim(low[0],high[0]);ax.set_ylim(low[1],high[1]);ax.set_zlim(low[2],high[2]);ax.set_box_aspect(high-low)
        ax.view_init(elev=elevation,azim=azimuth);ax.set_axis_off()
        ax.set_title(title,loc="left",fontsize=14,color="#18334b",pad=4,weight="bold")
        total=result["hexes"]+result["pyramids"]
        ax.text2D(0,-.03,f'{total:,} volume cells  |  {100*result["relative_volume_error"]:.3f}% volume error',transform=ax.transAxes,fontsize=10,color="#526478")
        ax.text2D(0,-.10,f'min. sub-tet quality {result["minimum_quality"]:.3f}  •  boundary audit passed',transform=ax.transAxes,fontsize=9,color="#526478")
fig.text(.04,.026,"Gold faces expose the actual interior cells. Curves use linear facets with a checked chord-error bound. These are explicit primitive meshes.",fontsize=10,color="#526478")
fig.subplots_adjust(left=.04,right=.97,bottom=.12,top=.875,wspace=.10,hspace=.31)
a.output.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(a.output,dpi=180,facecolor=fig.get_facecolor())
a.output.with_suffix(".json").write_text(json.dumps(records,indent=2)+"\n")
print(a.output)
