### OIB_track_query.py ###
""" read in OIB-AK-tracks.gpkg to query
all tracks within geographic area or date range

env: ragu

Brandon S. Tober
01MAR2021
"""
### imports ###
import os,sys
import fiona
from shapely.geometry import Point, LineString, Polygon
from shapely import wkt

# inputs
pfix = "/mnt/g/MARS/targ/supl/UAF/radar/2021/"
fpath = pfix + "OIB_2021_tracks.gpkg"
# x = [-141.5, -139.7]        # x lower, upper bounds
# y = [59.6, 60.5]            # y lower, upper bounds
out_fname = "kenneth_2021_tracks.txt"

# initialize output file list
out = []
fout = open("/mnt/c/Users/btobers/Documents/data/radar/kennicott/" + out_fname, "w")

# # create area polygon
# points = [[x[0],y[0]],
#         [x[1],y[0]],
#         [x[1],y[1]],
#         [x[0],y[1]],
#         [x[0],y[0]]]
# poly = Polygon(points)
# poly = wkt.loads("POLYGON ((-140.237214236105 60.2730619913547,\
#                             -139.774371622576 59.8105500244558,\
#                             -140.625452260143 59.6734001530028,\
#                             -141.337733991742 59.9006519264868,\
#                             -140.86909574261 60.2714177201349,\
#                             -140.237214236105 60.2730619913547))")

poly = wkt.loads("MULTIPOLYGON (((-143.28780902582 61.7795298326555,-142.771338593303 61.4678247363248,-142.942685254444 61.3645306498214,-143.448826278311 61.6981097879999,-143.28780902582 61.7795298326555)))")

if not os.path.isfile(fpath):
    sys.exit(1)

layers = fiona.listlayers(fpath)

for layer in layers:
    print(layer)
    with fiona.open(fpath, layer=layer) as layer:

        for feature in layer:
            ls_3d = LineString(feature["geometry"]["coordinates"])
            ls_2d = LineString([xy[0:2] for xy in list(ls_3d.coords)]) 

            if ls_2d.intersects(poly):
                fn = feature["properties"]["fname"]
                year = fn[:4]
                if os.path.isfile(pfix + "/hdf5/" + fn):
                    fout.write(pfix + "/hdf5/" + fn + "\n")
                elif os.path.isfile(pfix + + "/hdf5/" + fn):
                    fout.write(pfix + "/hdf5/" + fn + "\n")
                else:
                    print(fn)
                    pass
        
fout.close()