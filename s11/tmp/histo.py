import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import matplotlib.path as path

fig, ax = plt.subplots()

# charge data
# previous do: awk {'print $2'} velocidades_al_final.dat > vy_al_final.dat

data = open('vx_al_final.dat','r')
# Convert from a list of strings to a list of floats
data = map(float, data) 

# compare values from the list and array
ls = max(data)
li = min(data)
print ls, li
# we charge nat (as n) and number of bins
n, bins = np.histogram(data, 100)

left = np.array(bins[:-1])
print left
right = np.array(bins[1:])
print right

bottom = np.zeros(len(left))
top = bottom + n

XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T

# Rectangles with path and patch
barpath = path.Path.make_compound_path_from_polys(XY)

patch = patches.PathPatch(barpath, facecolor='green', edgecolor='gray', alpha=0.6)
ax.add_patch(patch)

# update the view limits
ax.set_xlim(left[0], right[-1])
ax.set_ylim(bottom.min(), top.max())

plt.show()
