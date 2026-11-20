#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Standalone performance figure from the retained isolated-process samples."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

def main():
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();data=json.loads(a.input.read_text())['summary']
    fig,axes=plt.subplots(2,2,figsize=(13,9),facecolor='#f3f6f9')
    fig.subplots_adjust(left=.085,right=.97,top=.80,bottom=.16,hspace=.40,wspace=.23)
    fig.suptitle('PAACS rc.7  /  Circular mesh size sweep',x=.06,y=.967,ha='left',fontsize=22,weight='bold',color='#18334b')
    fig.text(.06,.918,'4 threads · 3 fresh-process samples per point · Same nominal size and geometry limits',fontsize=12,color='#526478')
    for row,shape in enumerate(['torus','elbow']):
        for col,(field,divisor,label) in enumerate([('mesher_seconds',1,'Mesher time (seconds)'),('peak_rss_kib',1024,'Peak process RSS (MiB)')]):
            ax=axes[row,col];ax.set_facecolor('#f3f6f9')
            for family,style in [('tet4','-'),('tet10','--')]:
                for strategy,color in [('swept','#c57642'),('uniform','#208e9e')]:
                    group=sorted([r for r in data if r['case']==shape and r['element_family']==family and r['strategy']==strategy and r['threads']==4],key=lambda r:r['size'],reverse=True)
                    assert len(group)==3
                    ax.plot([r['size'] for r in group],[r[field]/divisor for r in group],style,marker='o',color=color,lw=2,label=f'{family.upper()} / {"Previous" if strategy=="swept" else "Balanced"}')
            ax.set_xscale('log');ax.set_yscale('log');ax.invert_xaxis()
            ax.set_xticks([.3,.15,.075]);ax.set_xticklabels(['0.3','0.15','0.075'])
            ax.set_xlabel('Requested size (finer →)');ax.set_ylabel(label)
            ax.set_title(f'{shape.title()} · {"Runtime" if col==0 else "Memory"}',loc='left',weight='bold',color='#18334b')
            ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.16,which='major')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.54,.88),ncol=4,frameon=False,fontsize=10)
    fig.text(.06,.078,'Time: median generation + validation + audit, excluding test oracles and file I/O. RSS: maximum of the 3 processes, including test oracles.',fontsize=9,color='#526478')
    fig.text(.06,.051,'Strategies produce different cell counts. These are two analytic primitives, not general CAD or equal-solver-error comparisons.',fontsize=9,color='#526478')
    a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=160,facecolor=fig.get_facecolor());plt.close(fig)
if __name__=='__main__':main()
