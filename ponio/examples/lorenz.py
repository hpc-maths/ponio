# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

name = "lorenz"
data_dir = name+"_data"

make = subprocess.Popen(["make", name])
make.wait()

args = [os.path.join(".", name)]
process = subprocess.Popen(args)
process.wait()

data = np.loadtxt(os.path.join(data_dir, "lorenz.dat"))
t, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# 3d plot
fig = plt.figure(constrained_layout=True, figsize=(8, 8))
gs = fig.add_gridspec(6, 1)
axlist = []

ax = fig.add_subplot(gs[:, 0], projection='3d')
ax.set_title("Lorenz attractor")
ax.set_xlim3d([-20, 20])
ax.set_ylim3d([-20, 20])
ax.set_zlim3d(bottom=0, top=50)

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')
# ax.grid(False)

ax.plot(x, y, z, linewidth=0.5)

axlist.append(ax)

plt.savefig(os.path.join(data_dir, "01.png"))

ax.set_title("Lorenz Attractor")
plt.show()

# 2d plot
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3)

ax1.plot(t, x, linewidth=0.5)
ax1.set_ylabel("$x$")
# print(ax1.get_xticklabels()[0])
xticklabels = ax1.get_xticklabels()
for xlab in xticklabels:
    xlab.set_text("")

ax1.set_xticks(ax1.get_xticks(), xticklabels)
ax1.set_xlim(t[0], t[-1])

ax2.plot(t, y, linewidth=0.5)
ax2.set_ylabel("$y$")
ax2.set_xticks(ax2.get_xticks(), xticklabels)
ax2.set_xlim(t[0], t[-1])

ax3.plot(t, z, linewidth=0.5)
ax3.set_ylabel("$z$")

ax3.set_xlabel("$t$")
ax3.set_xlim(t[0], t[-1])

plt.savefig(os.path.join(data_dir, "02.png"))

fig.suptitle("Composantes of Lorenz Attractor")
plt.show()
