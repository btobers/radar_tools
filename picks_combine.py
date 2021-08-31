import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import glob
import sys, os

# script designed to merge the picks from one directory in to a single csv file
# BST, 20191127
# updated 20210317

# gpkg is a funciton for saving picks to a geopackage/shapefile
def gpkg(fpath, df, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"):
    # fpath is the path for where the exported csv pick file should be saved [str]
    # df pick output dataframe
    # crs is the coordinate reference system for the shapefile output
    df_copy = df.copy()
    if df_copy["lon"].isnull().all() or df_copy["lat"].isnull().all():
        print("no geopackage was exported due to missing gps data")
        return
    # convert lon, lat to shapely points
    geometry = [Point(xy) for xy in zip(df_copy["lon"], df_copy["lat"])]
    df_copy.drop(["lon", "lat"], axis=1)

    # create geopandas df and export
    gdf = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
    gdf.to_file(fpath, driver="GPKG")

    print("geopackage exported successfully:\t" + fpath)


in_dir = "/mnt/c/Users/btobers/Documents/data/radar/kennicott/pick/"
# in_dir = "/media/btober/chupacabra/MARS/targ/supl/UAF/radar/IRUAFHF1B/pick/"
out_fname = "kennicott_wrst_v2"

track = np.array([]).astype(str)
trace = np.array([]).astype(int)
lon = np.array([]).astype(np.float)
lat = np.array([]).astype(np.float)
elev = np.array([]).astype(np.float)
srf_sample = np.array([]).astype(np.float)
srf_twtt = np.array([]).astype(np.float)
srf_elev = np.array([]).astype(np.float)
srf_amp = np.array([]).astype(np.float)
bed_sample = np.array([]).astype(np.float)
bed_twtt = np.array([]).astype(np.float)
bed_elev = np.array([]).astype(np.float)
bed_amp = np.array([]).astype(np.float)
srf_bed_thick = np.array([]).astype(np.float)

for file in glob.glob(in_dir + "*pk*.csv"):
    fname = file.split('/')[-1].rstrip("_pk_bst.csv")
    print(fname)

    dat = np.genfromtxt(file, delimiter=",", dtype = None, names = True)
    track = np.append(track, np.repeat(fname,dat.shape[0]))
    trace = np.append(trace,dat["trace"].astype(int))
    lon = np.append(lon,dat["lon"])
    lat = np.append(lat,dat["lat"])
    elev = np.append(elev,dat["elev"])
    srf_sample = np.append(srf_sample,dat["srf_sample"])
    srf_twtt = np.append(srf_twtt, dat["srf_twtt"])
    srf_elev = np.append(srf_elev,dat["srf_elev"])
    srf_amp = np.append(srf_amp,dat["srf_amp"])
    bed_sample = np.append(bed_sample,dat["bed_sample"])
    bed_twtt = np.append(bed_twtt,dat["bed_twtt"])
    bed_elev = np.append(bed_elev,dat["bed_elev"])
    bed_amp = np.append(bed_amp,dat["bed_amp"])
    srf_bed_thick = np.append(srf_bed_thick,dat["srf_bed_thick"])


out = pd.DataFrame({"track": track, "trace": trace, "lon": lon, "lat": lat, "elev": elev, "srf_sample": srf_sample,
                    "srf_twtt": srf_twtt, "srf_elev": srf_elev, "srf_amp": srf_amp, 
                    "bed_sample": bed_sample, "bed_twtt": bed_twtt, 
                    "bed_elev": bed_elev, "bed_amp": bed_amp, "srf_bed_thick": srf_bed_thick})

out.to_csv(in_dir + out_fname + ".csv", index=False)
if os.path.isfile(in_dir + out_fname + ".gpkg"):
    os.remove(in_dir + out_fname + ".gpkg")
gpkg(in_dir + out_fname + ".gpkg", out)

print("Combined picks exported to:\t" + in_dir + out_fname)