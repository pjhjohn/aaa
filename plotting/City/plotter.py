x, y, index = [], [], []
num = 0

with open('cycle.in.101') as text :
	num = (int)(text.readline())
	for i in range(num) :
		coord = text.readline().split(' ')
		x.append(float(coord[0]))
		y.append(float(coord[1]))

postfix = ''
with open('cycle.out.101' + postfix) as txt :
	index = txt.readline().split()
	for i in range(len(index)) :
		index[i] = int(index[i]) - 1

zipped = zip(x, y)
verts = []
for i in range(num) :
	verts.append(zipped[index[i]])
verts.append(zipped[index[0]])

import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
codes = [Path.MOVETO] + [Path.LINETO]*(num-1) + [Path.CLOSEPOLY]

path = Path(verts, codes)
fig = plt.figure()
ax = fig.add_subplot(111)
patch = patches.PathPatch(path, facecolor='white', lw=2)
ax.add_patch(patch)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
plt.savefig('plot_cycle_out_101%s.png' % postfix)