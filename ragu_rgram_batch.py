""" 
Batch radargram script
 
This script loops through a directory containing radar datafiles and RAGU generated pick files and creates a radargram for each set

Created by: Brandon S. Tober
Date: 20220520

Dependencies:
- ragu (https://github.com/btobers/ragu)
- numpy
- h5py
- pandas
- matplotlib
"""

### impots ###
import sys, os, itertools, glob, argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from geopandas import GeoDataFrame
from shapely.geometry import Point, LineString
import rasterio as rio
import pyproj

# append ragu path to access required modules
sys.path.append("C:/Users/btober/OneDrive/Documents/code/radar/ragu/code")
from ingest import ingest
from tools import utils

# function to load radargram and any accompanying picks
def load(f, datPath, pickPath):
    # ingest radar data
    fPath = datPath + "/" + f.rstrip("_pk_bst.csv") + ".h5"
    if os.path.isfile(fPath):
        igst = ingest(fpath)
        rdata = igst.read("", navcrs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", body="earth")

    # initialize dict to hold any import horizons
    horizons = {}
    pkPath = pickPath + "/" + f
    if os.path.isfile(pkPath):
        dat = pd.read_csv(pkPath)
        if dat.shape[0] != rdata.tnum:
            raise ValueError("import_pick error:\t pick file size does not match radar data")
        else:
            # get name of any horizons in pick file
            keys = fnmatch.filter(dat.keys(), "*sample*")
            for horizon in keys:
                # pick files should have ["horizon_sample"] keys
                horizon = horizon.split("_")[0]
                # add to horizon dict
                horizons[horizon] = dat[horizon + "_sample"]

    return rdata, horizons

    # make a nice radargram
    def rgram(rdata, horizons, params):

        # if clutter sim exists make 3 panels, otherwise 2
        if rdata.flags.sim:
            nPanels = 3
            # sim color limints
            vsmin = np.floor(np.nanpercentile(rdata.sim,10))
            vsmax = np.nanmax(rdata.sim)
        else:
            nPanels = 2

        # data color limits
        vdmin = np.floor(np.nanpercentile(rdata.proc.get_curr_dB(),10))
        vdmax = np.nanmax(rdata.proc.get_curr_dB())
        extent = [0, rdata.tnum, rdata.snum, 0]



        fig, ax = plt.subplots(panels, figsize=(params["pnlWidth"], nPanels*params["pnlHgt"]))
        fig.suptitle(rdata.fn)


srf_samps = rdata.pick.horizons["srf"]
bed_samps = rdata.pick.horizons["bed_imported"]

# plotting time
# first, take a quick look at rgram to get an idea of zoomed extent preference to set in next cell

extent = [0, rdata.tnum, rdata.snum, 0]
fig, ax = plt.subplots(figsize=(20,10))
ax.imshow(rdata.proc.get_curr_dB(), aspect="auto", extent=extent, cmap=cmap, vmin=vdmin, vmax=vdmax)


# plotting parameters - set zoomed extent preferences based on figure above
tl = 2500  # left trace
tr = rdata.tnum   # right trace
sb = 1500   # bottom sample
st = 0   # top sample
y_depth_s0 = 150    # sample number to set as zero-depth
extent = [-tl, tr-tl, rdata.snum, 0]




def main():

    # Set up CLI
    parser = argparse.ArgumentParser(
    description="Program for creating radargrams for a list of datafiles with corresponding pick files"
    )
    parser.add_argument("files", help="List of files to create radargrams for", nargs="+")
    parser.add_argument("datPath", help="Path to radar datafiles", nargs="+")
    parser.add_argument("pickPath", help="Path to RAGU radar pick files", nargs="+")
    parser.add_argument("outPath", help="Path to to output generated radargrams", nargs="+")
    args = parser.parse_args()

    # plotting parameters
    params = {}
    params["cmap"] = "Greys_r"                            # matplotlib.pyplot.imshow color map
    params["pnlHgt"] = 1.5                                # panel height in inches for each panel in the generated radargram
    params["pnlWidth"] = 6.5                              # panel width in inches for each panel in the generated radargram
    params["smplCutFact"] = .5                            # factor by which to trim the bottom half of the radargram (.5 will preserve the upper half of the samples across the radargram)

    # check if paths exists
    if not os.path.isdir(args.datPath[0]):
        print(f"Data file path not found: {args.datPath[0]}")
        exit(1)

    elif not os.path.isdir(args.pickPath[0]):
        print(f"Pick file path not found: {args.pickPath[0]}")
        exit(1)

    elif not os.path.isdir(args.outPath[0]):
        print(f"Output path not found: {args.outPath[0]}")
        exit(1)
    
    else:
        pass