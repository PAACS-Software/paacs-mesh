#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Plot the retained isolated-process measurements, without extrapolation."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

p=argparse.ArgumentParser();p.add_argument('measurements',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
data=json.loads(a.measurements.read_text());rows=data['summary']
fig,axes=plt.subplots(2,2,figsize=(13.5,9.2),facecolor='#f4f7fa')
fig.suptitle('PAACS rc.8  /  Measured scaling of local controls and layers',x=.065,y=.965,ha='left',fontsize=20,weight='bold',color='#19374b')
fig.text(.065,.915,'Three fresh processes per point · Median full mesher call · Maximum peak RSS · Linux / GCC 16 / C++26',fontsize=10.5,color='#516879')
for row,(cases,title) in enumerate([(['uniform-fine','cube','sphere','face'],'Tet4 local sizing'),(['structured','unstructured','mixed'],'Layered patches')]):
    for kind,color in zip(cases,['#687385','#bc6232','#258aa0','#735db4']):
        for threads,style in [(1,'--'),(4,'-')]:
            values=sorted([r for r in rows if r['case']==kind and r['threads']==threads],key=lambda r:r['volumeCells'])
            x=[r['volumeCells'] for r in values]
            for col,key in enumerate(['wallSeconds','peak_rss_kib']):
                y=[r[key]/(1024 if col==1 else 1) for r in values]
                axes[row,col].plot(x,y,style,marker='o',markersize=4,lw=1.7,color=color,label=kind if threads==4 else None,alpha=1 if threads==4 else .65)
    for col in range(2):
        ax=axes[row,col];ax.set_facecolor('#f4f7fa');ax.set_xscale('log');ax.set_yscale('log')
        ax.set_title(title+(' / time' if col==0 else ' / process memory'),loc='left',color='#19374b',fontsize=12,weight='bold')
        ax.set_xlabel('Generated volume cells');ax.set_ylabel('Mesher seconds' if col==0 else 'Peak RSS (MiB)')
        ax.spines[['top','right']].set_visible(False);ax.grid(True,which='major',alpha=.18);ax.minorticks_off()
        ax.legend(frameon=False,fontsize=9)
fig.legend(handles=[Line2D([],[],color='#516879',ls='--',label='1 thread'),Line2D([],[],color='#516879',ls='-',label='4 threads')],loc='lower center',bbox_to_anchor=(.5,.072),ncol=2,frameon=False)
fig.text(.065,.045,'Uniform-fine refines the whole box to the local target; adaptive cases use different size fields and cell counts. No equal-solver-error comparison.',fontsize=9.2,color='#516879')
fig.text(.065,.019,'Time excludes field compilation, source triangulation, hashing and file I/O. RSS includes source generation and startup; small runs have a memory floor.',fontsize=9.2,color='#516879')
fig.subplots_adjust(left=.08,right=.965,top=.865,bottom=.18,hspace=.4,wspace=.24)
fig.savefig(a.output,dpi=160,facecolor=fig.get_facecolor());plt.close(fig)
