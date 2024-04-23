import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import sys, os, argparse, fnmatch, configparser
from ragu import ingest

# script designed to build geopackage from multiple radar datafiles
# BST, 20240423

# gpkg is a funciton for saving picks to a geopackage/shapefile
def gpkg(df, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"):
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
    return gdf

def nav_merge(datpath, navinfo):

    # initialize list to hold all nav data of interest
    lons = []
    lats = []
    elevs = []
    fns = []
    
    # loop through files in datpath, try to use ragu to open
    for f in os.listdir(datpath):
        try:
            rdata = ingest.ingest(datpath + f).read(simpath=None, navcrs=navinfo["crs"], body=navinfo["body"])
            navdf = rdata.navdf
        except:
            print(f'Unable to load file: {datpath + f}')
            continue
            
        # append values
        lons.extend(navdf.lon.values.flatten())
        lats.extend(navdf.lat.values.flatten())
        elevs.extend(navdf.elev.values.flatten())
        fns.extend(np.repeat(rdata.fn, rdata.tnum).flatten())

    # instatiate a dictionary to hold all our merged arrays
    merged_nav = {}
    merged_nav['fn'] = fns
    merged_nav['lon'] = lons
    merged_nav['lat'] = lats
    merged_nav['elev'] = elevs

    return merged_nav


def main():
    # Set up CLI
    parser = argparse.ArgumentParser(
    description='''Program for merging a directory containing multiple RAGU pick files into a single merged pick file\n\n$python ragu_picks_combine.py --pickfiles /home/user/data/flist.txt --outdir /tmp --outname merged_picks\nalternatively:\n$python ragu_picks_combine.py --pickfiles /home/user/data/pickfile1.csv /home/user/data/pickfile2.csv --outdir /tmp --outname merged_picks''', 
    formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-datpath", help="Radar data file path", required=True)
    parser.add_argument("-outpath", help="Output navigation data file path", required=True)
    args = parser.parse_args()

    assert os.path.isdir(args.datpath), f'radar data directory does not exist: {args.datpath}'
    assert os.path.isdir(os.path.split(args.outpath)[0]), f'output directory does not exist: {os.path.split(args.outpath)[0]}'

    # get ragu base config
    basedir = os.path.join(os.path.expanduser('~'),'RAGU')
    configPath = basedir + '/config.ini'
    assert os.path.exists(configPath), f'ragu configuration file does not exist at the expected location - have you run ragu previously? {configPath}'
    conf = configparser.ConfigParser()
    conf.read(configPath)

    nav = nav_merge(datpath = args.datpath, navinfo = conf['nav'])
    nav = pd.DataFrame(nav)
    gdf = gpkg(df = nav, crs = conf['nav']['crs'])
    # save output as text file and as geopackage
    nav.to_csv(args.outpath + '/merged_nav.csv',index=False)
    gdf.to_file(args.outpath + '/merged_nav.gpkg', driver="GPKG")

    print(f"Combined nav exported to:\t {args.outpath}")


# execute if run as a script
if __name__ == "__main__":
    main()