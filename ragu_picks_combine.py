import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import sys, os, glob, argparse, fnmatch

# script designed to merge RAGU picks from one directory in to a single csv file
# BST, 20191127
# updated 20220222

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


def merge(dir):
    """
    merge() goes through a file directory and reads in all RAGU .csv pick files
    then merges them into a single set of numpy arrays to be saved
    """
    # we'll instatiate a dictionary to hold all our merged arrays
    merged_dict = {}
    # get list of pick files
    flist = glob.glob(dir + "*pk*.csv")
    # first we'll need to open a file and get all the unique fields that we'll need to make arrays for
    df = pd.read_csv(flist[0])
    cols = list(df.columns)
    dtypes = df.dtypes

    # instatiate an array for each column
    merged_dict["track"] = np.array([]).astype(str)
    for i, col in enumerate(cols):
        merged_dict[col] = np.array([]).astype(dtypes[i])

    for f in flist:
        # parse track name - this is a little kludgy but it does the trick
        fname = f.split("/")[-1].split('\\')[-1].rsplit("_pk_",1)[0]
        print(fname)

        # read file to dataframe
        df = pd.read_csv(f)

        # repeat track name for each trace
        merged_dict["track"] = np.append(merged_dict["track"], np.repeat(fname,df.shape[0]))
        # loop through all fields and append to array in merged_dict
        for i, col in enumerate(cols):
            merged_dict[col] = np.append(merged_dict[col], df[col])

    return merged_dict


def main():

    # Set up CLI
    parser = argparse.ArgumentParser(
    description="Program for merging a directory containing multiple RAGU pick files into a single merged pick file"
    )
    parser.add_argument("path", help="Path to RAGU pickfiles", nargs="+")
    parser.add_argument("fname", help="Output merged file name", nargs="+")
    args = parser.parse_args()

    # check if path exists
    if not os.path.isdir(args.path[0]):
        print(f"Path not found: {args.path[0]}")
        exit(1)

    # call the merge function to go through all files and merge all fields in pick files contained in path
    data = merge(args.path[0])

    # create pandas dataframe
    out = pd.DataFrame(data)

    # save output as text file and as geopackage
    out.to_csv(args.path[0] + args.fname[0] + ".csv", index=False)
    if os.path.isfile(args.path[0] + args.fname[0] + ".gpkg"):
        os.remove(args.path[0] + args.fname[0] + ".gpkg")
    gpkg(args.path[0] + args.fname[0] + ".gpkg", out)

    print(f"Combined picks exported to:\t {args.path[0] + args.fname[0]}.csv, {args.path[0] + args.fname[0]}.gpkg")


# execute if run as a script
if __name__ == "__main__":
    main()