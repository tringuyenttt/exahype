#!/usr/bin/python

import numpy as np

# picture size in pixel
picsize = np.array([1252, 165])

# patchsize or order+1 subcell grid size
basissize = 10

# This script gives you a width, height and maximum mesh
# size to resolve picsize pixel by pixel on a mesh with basissize
# subpoints per patch. It takes into account the partitioning
# of peano, i.e.

partitioning = 3.

# and tries not to waste too much points around the picture.
# It tells you how much it has to waste. Wasting is at least better
# than over-sampling the picture, thus creating crazily big data.

def find_nearest(array,value):
    A = array-value
    A[A < 0] = 1e10 # penalize negative values to find positives
    idx = A.argmin()
    return np.unravel_index(A.argmin(), A.shape)

basis = np.arange(10)
level = np.arange(4)

BASIS,LEVEL = np.meshgrid(basis,level)
resolutions = (partitioning * BASIS)**LEVEL

rx = find_nearest(resolutions, picsize[0])
ry = find_nearest(resolutions, picsize[1])

print "Optimal simulation domain to cover a picture of extends"
print picsize
print "is in x direction: domainsize=",resolutions[rx]
print "                   patchsize=",BASIS[rx], "level=",LEVEL[rx]
print "is in y direction: domainsize=",resolutions[ry]
print "                   patchsize=",BASIS[ry], "level=",LEVEL[ry]

rmax = tuple(map(max, zip(rx,ry)))

print "so the best we can do is patchsize=",BASIS[rmax],"level=",LEVEL[rmax]
print "  giving us a domain with extends ", resolutions[rmax]


