#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('results',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
rows=json.loads(a.results.read_text())['summary']
fig,axes=plt.subplots(1,3,figsize=(15,4.9),facecolor='#f3f6f9')
colors={('cylinder','tet4','automatic'):'#268699',('cylinder','tet10','automatic'):'#254e84',('cone','tet4','automatic'):'#dc9855',('cone','tet10','automatic'):'#b95e43',('cone','tet4','swept'):'#927dab'}
for key,color in colors.items():
    selected=sorted([r for r in rows if (r['case'],r['element_family'],r['strategy'])==key and r['threads']==4],key=lambda r:r['volume_cells'])
    x=[r['volume_cells'] for r in selected];label=f'{key[0]} {key[1]}'+(' · swept' if key[2]=='swept' else '')
    for ax,values in zip(axes,[[r['mesher_seconds'] for r in selected],[r['peak_rss_kib']/1024 for r in selected],[r['minimum_quality'] for r in selected]]):
        ax.plot(x,values,'o-',lw=1.8,ms=4,label=label,color=color);ax.set_xscale('log');ax.set_yscale('log')
for ax,title,y in zip(axes,['Meshing time · 4 threads','Whole-process peak memory','Minimum quality under refinement'],['Seconds · median of 3','MiB · maximum of 3','Minimum corner/sub-tet mean ratio']):
    ax.set_title(title,loc='left',fontsize=12,weight='bold',color='#18334b');ax.set_xlabel('Volume cells');ax.set_ylabel(y);ax.grid(alpha=.15);ax.set_facecolor('#f3f6f9');ax.spines[['top','right']].set_visible(False)
axes[2].axhline(.1,color='#b95e43',lw=.8,ls='--');axes[0].legend(fontsize=8,loc='upper left')
fig.suptitle('PAACS rc.5   /   Measured Tet4 & Tet10 scaling',x=.04,y=.97,ha='left',fontsize=19,weight='bold',color='#18334b')
fig.text(.04,.025,'GCC 16 / C++26 · Linux · fresh processes · CPU only · RSS includes geometric test oracles · hashes match across 1 and 4 threads',fontsize=9,color='#526478')
fig.subplots_adjust(left=.06,right=.98,bottom=.18,top=.81,wspace=.3);fig.savefig(a.output,dpi=180,facecolor=fig.get_facecolor());print(a.output)
