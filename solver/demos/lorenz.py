#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

os.makedirs("example6",exist_ok=True)

make = subprocess.Popen(["make","lorenz"])
make.wait()

N = 10

args = ["./lorenz"]
process = subprocess.Popen(args)
process.wait()

data = np.loadtxt("example6/lorenz.dat")
T,x,y,z = data[:,0],data[:,1],data[:,2],data[:,3]

fig = plt.figure(constrained_layout=True,figsize=(8,8))
gs = fig.add_gridspec(6,1)
axlist = []

ax = fig.add_subplot(gs[:3,0],projection='3d')
ax.set_title("Lorenz attractor")
ax.plot(x,y,z,label="Lorenz attractor",linewidth=0.375)
ax.set_xlim3d([-20,20])
ax.set_ylim3d([-20,20])
ax.set_zlim3d(bottom=0,top=50)
axlist.append(ax)

for i,(lab,data,ylim) in enumerate(zip(("x","y","z"),(x,y,z),([-20,20],[-20,20],[0,50]))):
  ax = fig.add_subplot(gs[3+i,0])
  ax.set_ylabel("${}$".format(lab))
  if lab == "z" :
    ax.set_xlabel("$t$")
  ax.set_xlim((T[0],T[-1]))
  ax.set_ylim(ylim)
  ax.plot( T , data , label="${}$".format(lab) , color=cm.tab10(i+1) )
  axlist.append(ax)

handles, labels = map(lambda x:sum(x,start=[]),zip(*map( lambda ax:ax.get_legend_handles_labels() , axlist )))
#fig.legend(handles, labels,loc="center",bbox_to_anchor=(0.26,0.55),ncol=len(labels))
plt.show()
