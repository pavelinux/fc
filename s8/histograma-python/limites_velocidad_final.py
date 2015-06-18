#!/usr/bin/env python
# To be executed in ipython
import subprocess
vx = subprocess.call(["cat", "vx_al_final.dat"])
vx = map(float, vx)
limite_s = max(vx)
limite_i = min(vx)

print limite_s, limite_i
