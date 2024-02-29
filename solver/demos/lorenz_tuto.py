# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

name = "lorenz_tuto"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()


fig = plt.figure(constrained_layout=True, figsize=(8, 8))
gs = fig.add_gridspec(6, 1)
axlist = []

ax1 = fig.add_subplot(gs[:, 0], projection='3d')
ax1.set_title("Lorenz attractor")
ax1.set_xlim3d([-20, 20])
ax1.set_ylim3d([-20, 20])
ax1.set_zlim3d(bottom=0, top=50)

for j, (meth, label) in enumerate([("rk44", "RK(4,4)"), ("lrk44", "LRK(4,4)"), ("lie", "Lie"), ("strang", "Strang")]):
    filename = os.path.join(data_dir, meth+".dat")
    data = np.loadtxt(filename)
    T, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    ax1.plot(x, y, z, label=label, color=f"C{j}", linewidth=0.375)

ax1.legend()
axlist.append(ax1)

handles, labels = map(lambda x: sum(x, start=[]), zip(
    *map(lambda ax: ax.get_legend_handles_labels(), axlist)))
# fig.legend(handles, labels,loc="center",bbox_to_anchor=(0.26,0.55),ncol=len(labels))
plt.show()
