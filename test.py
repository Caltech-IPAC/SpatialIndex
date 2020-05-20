from spatial_index import SpatialIndex

spt = SpatialIndex()


retval = spt.cone_search(34., 45., 0.4)

print(retval['index_constraint'])
print(retval['geom_constraint'])


retval = spt.cone_search(34., 45., 0.4, mode=spt.HPX, level=14)

print(retval['index_constraint'])
print(retval['geom_constraint'])


retval = spt.cone_search(129.4, 43.7, 0.4, level=7)

print(retval['index_constraint'])
print(retval['geom_constraint'])

