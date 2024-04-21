from bte_projection import BTE_PROJECTION

geo_coord = [40.74843814459844, -73.98566440289457]
xy_coord = [-8525873.069135161, -6026164.9710848285]

print(BTE_PROJECTION.fromGeo(geo_coord[0], geo_coord[1]))
print(BTE_PROJECTION.toGeo(xy_coord[0], xy_coord[1]))