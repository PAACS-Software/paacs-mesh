#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Measured runtime, RSS and mesh-quality curves; no invented data."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument("summary",type=Path);p.add_argument("output",type=Path);a=p.parse_args()
rows=json.loads(a.summary.read_text());names=list(dict.fromkeys(r["case"] for r in rows))
fig,axes=plt.subplots(1,3,figsize=(15,4.8),facecolor="#f3f6f9")
colors=["#167b96","#c57531","#6c64ad","#427c5b"]
for name,color in zip(names,colors):
    selected=sorted([r for r in rows if r["case"]==name and r["threads"]==1],key=lambda r:r["volume_cells"])
    x=[r["volume_cells"] for r in selected]
    for ax,key,mult in [(axes[0],"mesher_seconds",1),(axes[1],"process_peak_rss_kib",1/1024),(axes[2],"minimum_quality",1)]:
        ax.plot(x,[r[key]*mult for r in selected],"o-",label=name,color=color,linewidth=2,markersize=5)
for ax,title,ylabel in zip(axes,["Mesher time / 1 thread","Peak process memory","Quality while refining"],["Seconds, median of 3 runs","MiB, maximum of 3 runs","Minimum sub-tet mean ratio"]):
    ax.set_facecolor("#f3f6f9");ax.set_title(title,loc="left",fontsize=14,weight="bold");ax.set_xscale("log");ax.set_yscale("log");ax.set_xlabel("Volume cells");ax.set_ylabel(ylabel);ax.grid(alpha=.18);ax.spines[["top","right"]].set_visible(False)
axes[0].legend(frameon=False,fontsize=9)
fig.suptitle("PAACS rc.4  /  Measured scaling",x=.055,ha="left",fontsize=19,weight="bold",y=.98)
fig.text(.055,.018,"Generation + cell validation + boundary audit. RSS also includes independent qualification. Cone-tip quality remains a limitation.",fontsize=10,color="#526478")
fig.subplots_adjust(left=.06,right=.98,bottom=.20,top=.80,wspace=.29)
fig.savefig(a.output,dpi=170,facecolor=fig.get_facecolor())
