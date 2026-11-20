#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Plot measured cylinder scaling; no fitted/extrapolated performance."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
rows=json.loads(a.input.read_text())['summary']
fig,axes=plt.subplots(1,2,figsize=(12.5,5.5),facecolor='#f3f6f9')
fig.suptitle('PAACS rc.6  /  Cylinder cost at equal requested sizes',x=.065,ha='left',y=.97,fontsize=19,weight='bold',color='#18334b')
for family,color in [('tet4','#2a8e9c'),('tet10','#b77939')]:
    for strategy,style in [('swept','--'),('uniform','-')]:
        r=sorted([v for v in rows if v['element_family']==family and v['strategy']==strategy and v['threads']==4],key=lambda v:v['size'],reverse=True)
        x=[v['size'] for v in r]
        for ax,key,scale in [(axes[0],'mesher_seconds',1),(axes[1],'peak_rss_kib',1/1024)]:
            ax.plot(x,[v[key]*scale for v in r],style,marker='o',color=color,label=f'{family.upper()} · {strategy}',linewidth=2,markersize=5)
for ax in axes:
    ax.set_xscale('log',base=2);ax.set_yscale('log');ax.invert_xaxis();ax.set_facecolor('#f3f6f9')
    ax.set_xticks([.3,.15,.075,.0375],['0.3','0.15','0.075','0.0375']);ax.set_xlabel('Requested size h → finer mesh',fontsize=11)
    ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.18,which='major');ax.tick_params(labelsize=10)
axes[0].set_ylabel('Median mesher time (s)');axes[0].set_title('Four CPU threads · 3 repeats',loc='left',fontsize=13,color='#18334b')
axes[1].set_ylabel('Peak test-process RSS (MiB)');axes[1].set_title('Maximum across repeats',loc='left',fontsize=13,color='#18334b')
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.11),ncol=4,frameon=False,fontsize=10)
fig.text(.065,.035,'96 fresh-process runs, including 1-thread checks. Different cell counts at the same geometry tolerances. RSS includes independent test oracles.\nCoarse uniform meshes can cost more because curvature constrains their global spacing.',fontsize=9,color='#526478')
fig.subplots_adjust(left=.08,right=.965,bottom=.27,top=.82,wspace=.26)
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=175,facecolor=fig.get_facecolor());print(a.output)
