#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Exact faces and geometric sections from the complex-geometry test archives."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,to_rgb
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from render_cylinder_comparison import read,statistics
from render_local_controls import section,clip
BG,INK,SUB='#f4f7fa','#19374b','#516879'
CMAP=plt.get_cmap('YlGnBu_r');NORM=Normalize(.12,.45)
SWEEP={'i_beam':(0,0,4),'multi_bore':(.25,-.1,2.8),'star_duct':(0,0,3.2),'cooling_channels':(0,0,2.4),'notched_frame':(.2,.15,2.5),'thin_bracket':(0,0,2.2)}

def field(z,length,mode):
 z=np.asarray(z);a,b=(.38*length,.6*length) if mode=='band' else ((length,length) if mode=='end' else (0,0))
 d=np.maximum.reduce([a-z,z-b,np.zeros_like(z)])
 return np.minimum(.45,.12+2*(1.25-1)/(1.25+1)*d)

def get(folder,name):
 path=folder/(name+'.paacsmesh');pts,cells=read(path)
 count=sum(d==3 for _,d,_,_ in cells)
 return pts,cells,dict(archive=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),nodes=len(pts),volume_cells=count)

def cut(ax,pts,cells,plane,length,height,mode):
 polygons,_,z=section(pts,cells,plane)
 faces,colors=[],[]
 light=np.array([.5,-.7,1.]);light/=np.linalg.norm(light)
 az,el=-58,25;view=np.array([np.cos(np.deg2rad(el))*np.cos(np.deg2rad(az)),np.cos(np.deg2rad(el))*np.sin(np.deg2rad(az)),np.sin(np.deg2rad(el))])
 for _,dim,tag,ids in cells:
  if dim!=2:continue
  poly=clip(pts[list(ids)],plane)
  if len(poly)<3:continue
  normal=np.cross(poly[1]-poly[0],poly[2]-poly[0]);norm=np.linalg.norm(normal)
  if norm<1e-14:continue
  normal/=norm
  if normal@view<=0:continue
  faces.append(poly);colors.append(np.array(to_rgb('#88cbd1'))*(.73+.27*max(0,normal@light)))
 for poly,h in zip(polygons,z):
  faces.append(np.column_stack([poly[:,0],np.full(len(poly),plane),poly[:,1]]));colors.append(CMAP(NORM(field(h/height*length,length,mode))))
 ax.add_collection3d(Poly3DCollection(faces,facecolors=colors,edgecolors='#3b5564',linewidths=.22))
 lo=pts.min(axis=0);hi=pts.max(axis=0);ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2]));ax.set_box_aspect(hi-lo,zoom=1.05);ax.view_init(el,az);ax.set_proj_type('ortho');ax.set_axis_off()

def gallery(folder,out):
 entries=[('i_beam','start','native','I-beam · Prism6'),('multi_bore','start','tet4','Six polygonal bores · Tet4'),('star_duct','band','tet4','Star duct · Tet4'),('cooling_channels','start','tet4','Twelve cooling channels · Tet4'),('notched_frame','band','native','Notched frame · Prism6'),('thin_bracket','start','tet4','Thin bracket · Tet4')]
 fig=plt.figure(figsize=(16,11),facecolor=BG);fig.text(.04,.957,'PAACS rc.10  /  Controls on complex polygonal geometry',fontsize=23,weight='bold',color=INK)
 fig.text(.04,.919,'Actual exported meshes · Selected-cap / internal-band sizing · Preserved face tags and through-holes',fontsize=12,color=SUB)
 grid=fig.add_gridspec(2,3,left=.035,right=.975,top=.86,bottom=.14,hspace=.18,wspace=.04);records=[]
 for i,(name,mode,family,label) in enumerate(entries):
  pts,cells,record=get(folder,f'{name}-{mode}-{family}');records.append(record)
  sweep=np.array(SWEEP[name]);length=np.linalg.norm(sweep)
  plane={'i_beam':1.3,'multi_bore':.8,'star_duct':.04,'cooling_channels':1.5,'notched_frame':1.05,'thin_bracket':.057}[name]
  ax=fig.add_subplot(grid[i//3,i%3],projection='3d',facecolor=BG);cut(ax,pts,cells,plane,length,sweep[2],mode)
  ax.set_title(label,loc='left',fontsize=13,color=INK,weight='bold',pad=3)
  ax.text2D(.03,-.018,f'{record["volume_cells"]:,} cells · {record["nodes"]:,} nodes · '+('fine band' if mode=='band' else 'fine start cap'),transform=ax.transAxes,fontsize=10,color=SUB)
 cax=fig.add_axes([.17,.083,.4,.013]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cax,orientation='horizontal');cb.set_ticks([.12,.3,.45]);cb.ax.tick_params(labelsize=9)
 fig.text(.60,.08,'Cut color = requested axial size',fontsize=10,color=SUB)
 fig.text(.04,.035,'Cuts intersect the exported cells exactly. Cyan = exterior. Cross-section sizing stays fixed; the selected control changes spacing along the sweep.',fontsize=10,color=SUB)
 fig.text(.04,.013,'Polygonal holes approximate bores; these examples do not establish arbitrary curved-CAD support. Native triangular sweeps produce Prism6 (Wedge6).',fontsize=10,color=SUB)
 fig.savefig(out,dpi=160,facecolor=BG);plt.close(fig);return records

def controls(folder,out):
 name='multi_bore';sweep=np.array(SWEEP[name]);length=np.linalg.norm(sweep);height=sweep[2]
 fig=plt.figure(figsize=(16,10),facecolor=BG);fig.text(.04,.959,'PAACS rc.10  /  The same controls on a part with six holes',fontsize=22,weight='bold',color=INK)
 fig.text(.04,.917,'Selected start cap, opposite cap, or internal band · Axial local size 0.12 → growth 1.25 → global cap 0.45',fontsize=11.5,color=SUB)
 grid=fig.add_gridspec(3,3,left=.05,right=.97,top=.846,bottom=.15,wspace=.3,hspace=.65,width_ratios=[1.25,1.25,.85]);records=[]
 for row,mode in enumerate(['start','end','band']):
  stations=None
  for col,(family,label) in enumerate([('native','Prism6'),('tet4','Tet4')]):
   pts,cells,record=get(folder,f'{name}-{mode}-{family}');records.append(record);polys,_,heights=section(pts,cells,.8)
   ax=fig.add_subplot(grid[row,col],facecolor='white')
   polys=[np.column_stack([poly[:,1]/height*length,poly[:,0]]) for poly in polys]
   ax.add_collection(PolyCollection(polys,facecolors=CMAP(NORM(field(heights/height*length,length,mode))),edgecolors='#3b5564',linewidths=.32))
   ax.set(xlim=(0,length),ylim=(0,4.3),xticks=[0,1,2,length],yticks=[0,2,4]);ax.tick_params(labelsize=8,colors=SUB);ax.set_xlabel('Distance along sweep',fontsize=9,color=SUB)
   if col==0:ax.set_ylabel('x · exact y = 0.8 cut',fontsize=9,color=SUB)
   ax.set_title(f'{mode.capitalize()} control · {label} · {record["volume_cells"]:,} cells',loc='left',fontsize=10.5,color=INK,pad=8)
   if col==0:stations=np.unique(pts[:,2])/height*length
  ax=fig.add_subplot(grid[row,2],facecolor='white');z=np.linspace(0,length,500)
  ax.plot(z,field(z,length,mode),c='#21698b',lw=2,label='Requested size');ax.step((stations[:-1]+stations[1:])/2,np.diff(stations),where='mid',c='#c0782e',lw=1.6,label='Actual layer width')
  ax.set(xlim=(0,length),ylim=(0,.49),xticks=[0,length/2,length],yticks=[.12,.3,.45]);ax.tick_params(labelsize=8,colors=SUB);ax.grid(alpha=.18)
  ax.set_title('Axial size and measured spacing',loc='left',fontsize=10.5,color=INK);ax.set_xlabel('Distance along sweep',fontsize=9,color=SUB)
  if row==0:ax.legend(frameon=False,fontsize=8,loc='lower right')
 cax=fig.add_axes([.14,.076,.44,.013]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cax,orientation='horizontal');cb.set_ticks([.12,.3,.45]);cb.ax.tick_params(labelsize=9)
 fig.text(.63,.073,'Color = requested size, not quality',fontsize=10,color=SUB)
 fig.text(.04,.033,'White channels are geometric holes, not missing cells. Prism-to-tet conversion adds no nodes and preserves the same layers and exterior geometry.',fontsize=10,color=SUB)
 fig.text(.04,.012,'The guarantee applies to axial layer widths and adjacent growth. Tangential resolution is a separate control; a short plateau can remain below the cap.',fontsize=10,color=SUB)
 fig.savefig(out,dpi=160,facecolor=BG);plt.close(fig);return records

def main():
 p=argparse.ArgumentParser();p.add_argument('archives',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 records=gallery(a.archives,a.output/'paacs-rc10-complex-shapes.png')
 records+=controls(a.archives,a.output/'paacs-rc10-complex-controls.png')
 (a.output/'complex-gallery-manifest.json').write_text(json.dumps(records,indent=2)+'\n')
 from PIL import Image
 for path in a.output.glob('paacs-rc10-*.png'):
  with Image.open(path) as im:im.verify()
  print(path,path.stat().st_size,flush=True)
if __name__=='__main__':main()
