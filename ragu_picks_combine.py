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

    return

def merge(dir, merged=False):
    """
    merge() goes through a file directory and reads in all RAGU .csv pick files
    then merges them into a single set of numpy arrays to be saved
    """
    # we'll instatiate a dictionary to hold all our merged arrays
    merged_dict = {}
    # get list of pick files
    flist = glob.glob(dir + "*.csv")

    if not merged:
        # instatiate an array to hold track name for each trace in all files - this won't be encountered in the individual files so we must do it ahead of time
        merged_dict["fname"] = np.array([]).astype(str)

    for f in flist:
        # skip any note file
        if "note" in f:
            continue
        # skip combined pick files
        if "combine" in f:
            continue
    
        print(f)
        # read file to dataframe
        df = pd.read_csv(f)
        cols = list(df.columns)
        dtypes = df.dtypes

        if not merged:
            # parse track name - this is a little kludgy but it does the trick
            fname = f.split("/")[-1].split('\\')[-1].replace('.csv','')
            # repeat track name for each trace
            merged_dict["fname"] = np.append(merged_dict["fname"], np.repeat(fname, df.shape[0]))
            # hold total length of merged dictionary thus far - we'll need this later in the loop
            l0 = merged_dict["fname"].shape[0]

        # loop through all fields and append to array in merged_dict
        for i, col in enumerate(cols):

            # instatiate column 
            if col not in merged_dict.keys():
                merged_dict[col] = np.array([]).astype(dtypes[i])

            # append series from df to merged dictionary
            merged_dict[col] = np.append(merged_dict[col], df[col])

            if not merged:
                # since some ragu exports may contain more horizons than others, we'll need to account for opening a file later that has additional fields - this would cause length inconsistencies in the merged dictionary
                # account for length issues if a new column was added to the merged dictionary after the first file in the loop - add nan's to have equal length
                l = merged_dict[col].shape[0]
                if l != l0:
                    merged_dict[col] = np.append(np.repeat(np.nan, l0-l), merged_dict[col])

    return merged_dict


def main():
    # Set up CLI
    parser = argparse.ArgumentParser(
    description='''Program for merging a directory containing multiple RAGU pick files into a single merged pick file\n\n$python ragu_picks_combine.py F:/MARS/orig/xtra/OIB-AK/radar/2015/pk_bst/ 2015_pk_bst''', 
    formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("path", help="Path to RAGU pickfiles", nargs="+")
    parser.add_argument("fname", help="Output merged file name", nargs="+")
    parser.add_argument('-merged', help='Flag: input files are merged from multiple tracks (not individual track pick files)', default=False, action='store_true')
    args = parser.parse_args()

    # check if path exists
    if not os.path.isdir(args.path[0]):
        print(f"Path not found: {args.path[0]}")
        exit(1)

    # call the merge function to go through all files and merge all fields in pick files contained in path
    data = merge(args.path[0], merged=args.merged)

    # create pandas dataframe
    out = pd.DataFrame(data)

    # get number of tracks
    print(f"{len(out['fname'].unique())} pick files combined.")

    # save output as text file and as geopackage
    out.to_csv(args.path[0] + args.fname[0] + ".csv", index=False)
    if os.path.isfile(args.path[0] + args.fname[0] + ".gpkg"):
        os.remove(args.path[0] + args.fname[0] + ".gpkg")
    gpkg(args.path[0] + args.fname[0] + ".gpkg", out)

    print(f"Combined picks exported to:\t {args.path[0] + args.fname[0]}.csv, {args.path[0] + args.fname[0]}.gpkg")


# execute if run as a script
if __name__ == "__main__":
    main()