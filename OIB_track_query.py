### OIB_track_query.py ###
""" read in OIB-AK-tracks.gpkg to query
all tracks within geographic area or date range
and export a list of the queried tracks


Brandon S. Tober
20210301
updated: 20220223
"""
### imports ###
import os,sys,argparse
import fiona
from shapely.geometry import Point, LineString, Polygon
from shapely import wkt

def query(gpkg,l,r,b,t,out):
    """
    open gpkg file and use lat, lon inputs to create bounding poly from which to query tracks and export file names
    this is currently based on the hierarchy data organization with NSIDC
    """
    # initialize output file list
    fout = open(out, "w")

    # create area polygon using input lat lon as wkt inputs
    # poly = wkt.loads("MULTIPOLYGON (((-144.2 60, -144.2 60.72, -139.8 60.7,-139.9 59.5, -144.2 60)))")
    poly = wkt.loads(f"MULTIPOLYGON ((({l} {b}, {l} {t}, {r} {t}, {r} {b}, {l} {b})))")

    # get parent folder 
    path = os.path.dirname(gpkg)

    layers = fiona.listlayers(gpkg)

    print("-------------Query Results:-------------")
    for layer in layers:
        # print(layer)
        with fiona.open(gpkg, layer=layer) as layer:

            for feature in layer:
                ls_3d = LineString(feature["geometry"]["coordinates"])
                ls_2d = LineString([xy[0:2] for xy in list(ls_3d.coords)]) 

                if ls_2d.intersects(poly):
                    fn = feature["properties"]["fname"]
                    year = fn[:4]
                    if os.path.isfile(path + "/hdf5/" + fn):
                        print(path + "/hdf5/" + fn)
                        fout.write(path + "/hdf5/" + fn + "\n")

                    else:
                        # print(fn)
                        pass
            
    fout.close()

    return



def main():

    # Set up CLI
    parser = argparse.ArgumentParser(
    description="Program for querying NASA OIB radar tracks based on latitude/longitude bounds"
    )
    parser.add_argument("gpkg", help="path to geopackage containing radar profiles", nargs="+")
    parser.add_argument("lon0", help="western longitude bound", nargs="+")
    parser.add_argument("lon1", help="eastern longitude bound", nargs="+")
    parser.add_argument("lat0", help="southern latitude bound", nargs="+")
    parser.add_argument("lat1", help="northern latitude bound", nargs="+")
    parser.add_argument("out", help="path of output query results text file", nargs="+")
    args = parser.parse_args()
    out = args.out[0]
    if not out.endswith(".txt"):
        out += ".txt"

    # check if path exists
    if not os.path.isfile(args.gpkg[0]):
        print(f"Geopackage not found: {args.gpkg[0]}")
        exit(1)

    # call the query function to go through the geopackage and query tracks based on bounding box and save file with track paths
    query(args.gpkg[0], args.lon0[0], args.lon1[0], args.lat0[0], args.lat1[0], out)


    print(f"Query results exported to:\t {os.path.abspath(out)}")


# execute if run as a script
if __name__ == "__main__":
    main()