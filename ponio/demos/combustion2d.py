# Copyright 2022 PONIO TEAM. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

#!/usr/bin/env python

import subprocess
import os
import glob

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

name = "combustion2d"
data_dir = name+"_data"

# make = subprocess.Popen(["make", name])
# make.wait()

# args = [os.path.join(".", name)]
# process = subprocess.Popen(args)
# process.wait()


def load_data(filename):
    i = int(os.path.basename(filename)[2:-4])
    data = np.loadtxt(filename)

    return (i, data)


print("load data...")
data = []
for f in glob.glob(os.path.join(data_dir, "u_*.dat")):
    i = int(os.path.basename(f)[2:-4])
    if i % 100 == 0:
        data.append(load_data(f))
print(f"sort {len(data)} data...")
data = sorted(data, key=lambda x: x[0])

print("prepare visu...")
fig = plt.figure()
ax = plt.axes()
txt_frame = ax.text(.5, 1.05, '', transform=ax.transAxes, va='center')

im = plt.imshow(data[0][1], vmin=1-1e-2, vmax=1+1e-1)


def animate(i):
    txt_frame.set_text(f"frame: {i}")
    im.set_array(data[i][1])
    return [im]


anim = animation.FuncAnimation(
    fig, animate, frames=len(data), interval=1000 / 30)

plt.show()
