#!/bin/env python

from spatial_index import SpatialIndex

spt = SpatialIndex()


print('------')

retval = spt.cone_search(34., 45., 0.4)

print(retval['index_constraint'])
print(retval['geom_constraint'])

print('------')

retval = spt.cone_search(34., 45., 0.4, mode=spt.HPX, level=14)

print(retval['index_constraint'])
print(retval['geom_constraint'])

print('------')

retval = spt.cone_search(129.4, 43.7, 0.4, level=7)

print(retval['index_constraint'])
print(retval['geom_constraint'])

print('------')

retval = spt.polygon_search(4, [159.903, 159.903, 159.902, 159.902], [43.103, 43.102, 43.102, 43.103], level=7)

print(retval['index_constraint'])
print(retval['geom_constraint'])

print('------')
