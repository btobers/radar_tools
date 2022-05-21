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
    fPath = datPath + "/" + f[:-len("__pk_bst.csv")] + ".h5"
    if os.path.isfile(fPath):
        igst = ingest(fPath)
        rdata = igst.read("", navcrs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", body="earth")
    else:
        print("Data file not found: {}".format(fPath))
        return None, None

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
def radargram(rdata, horizons, params, outPath):
    nPanels = 3
    if rdata.flags.sim:
        # sim color limints
        vsmin = np.floor(np.nanpercentile(rdata.sim,10))
        vsmax = np.nanmax(rdata.sim)
    else:
        nPanels -= 1

    if len(horizons) == 0:
        nPanels -=1


    # data color limits
    vdmin = np.floor(np.nanpercentile(rdata.proc.get_curr_dB(),10))
    vdmax = np.nanmax(rdata.proc.get_curr_dB())
    # radargram figure extent
    extent = [0, rdata.tnum, rdata.snum, 0]


    # initialize figure
    fig, ax = plt.subplots(nPanels, figsize=(params["pnlWidth"], nPanels*params["pnlHgt"]))
    fig.suptitle(rdata.fn)


    ax[0].imshow(rdata.proc.get_curr_dB(), aspect="auto", extent=extent, cmap=params["cmap"], vmin=vdmin, vmax=vdmax)
    # ax[0].set_xlim(tl, tr)
    ax[0].set_ylim(int(rdata.snum*params["smplCutFact"]), 0)
    ax[0].get_xaxis().set_visible(False)
    ax[0].get_yaxis().set_visible(False)
    # secaxy = ax[0].twinx()
    # secaxy.yaxis.set_ticks_position("left")
    # secaxy.yaxis.set_label_position("left")
    # secaxy.set_ylim(utils.twtt2depth((sb - y_depth_s0) * rdata.dt),
    #                 -1*utils.twtt2depth(y_depth_s0 * rdata.dt))
    # secaxy.set_ylabel("Depth (m)")
    # secaxy.set_yticks([0,500,1000])
    secaxx = ax[0].twiny()
    secaxx.xaxis.set_ticks_position("bottom")
    secaxx.xaxis.set_label_position("bottom")
    secaxx.set_xlim(rdata.navdf["dist"].iloc[0]*1e-3, rdata.navdf["dist"].iloc[-1]*1e-3)
    secaxx.get_xaxis().set_visible(False)

    if rdata.flags.sim:
        ax[1].imshow(rdata.sim,  aspect="auto", extent=extent, cmap=params["cmap"], vmin=vsmin, vmax=vsmax)
        # ax[1].set_xlim(tl, tr)
        ax[1].set_ylim(int(rdata.snum*params["smplCutFact"]), 0)
        ax[1].get_xaxis().set_visible(False)
        ax[1].get_yaxis().set_visible(False)
        # secaxy = ax[1].twinx()
        # secaxy.yaxis.set_ticks_position("left")
        # secaxy.yaxis.set_label_position("left")
        # secaxy.set_ylim(utils.twtt2depth((sb - y_depth_s0) * rdata.dt),
        #                 -1*utils.twtt2depth(y_depth_s0 * rdata.dt))
        # secaxy.set_ylabel("Depth (m)")
        # secaxy.set_yticks([0,500,1000])
        secaxx = ax[1].twiny()
        secaxx.xaxis.set_ticks_position("bottom")
        secaxx.xaxis.set_label_position("bottom")
        secaxx.set_xlim(rdata.navdf["dist"].iloc[-1]*1e-3, rdata.navdf["dist"].iloc[-1]*1e-3)
        secaxx.get_xaxis().set_visible(False)

    ax[-1].imshow(rdata.proc.get_curr_dB(),  aspect="auto", extent=extent, cmap=params["cmap"], vmin=vdmin, vmax=vdmax)
    # ax[-1].set_xlim(tl, tr)
    ax[-1].set_ylim(int(rdata.snum*params["smplCutFact"]), 0)
    for hz in horizons.items():
        ax[-1].plot(hz)
    ax[-1].get_xaxis().set_visible(False)
    ax[-1].get_yaxis().set_visible(False)
    # secaxy = ax[-1].twinx()
    # secaxy.yaxis.set_ticks_position("left")
    # secaxy.yaxis.set_label_position("left")
    # secaxy.set_ylim(utils.twtt2depth((sb - y_depth_s0) * rdata.dt),
    #                 -1*utils.twtt2depth(y_depth_s0 * rdata.dt))
    # secaxy.set_ylabel("Depth (m)")
    # secaxy.set_yticks([0,500,1000])
    secaxx = ax[-1].twiny()
    secaxx.xaxis.set_ticks_position("bottom")
    secaxx.xaxis.set_label_position("bottom")
    secaxx.set_xlim(rdata.navdf["dist"].iloc[0]*1e-3, rdata.navdf["dist"].iloc[-1]*1e-3)
    secaxx.set_xlabel("Along-Track Distance (km)")

    fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0.1)
    fig.savefig(outPath + "/" + rdata.fn + ".png", dpi=500, bbox_inches='tight', pad_inches=0.05, transparent=True)# facecolor = "#d9d9d9")



    # srf_samps = rdata.pick.horizons["srf"]
    # bed_samps = rdata.pick.horizons["bed_imported"]

    return


def main():

    # Set up CLI
    parser = argparse.ArgumentParser(
    description="Program for creating radargrams for a list of datafiles with corresponding pick files"
    )
    parser.add_argument("fileList", help="List of files for which to create radargrams", nargs="+")
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

    # loop through list of files, load radargram and any picks, generate radargram
    with open(args.fileList[0]) as files:
        for f in files:
            rdata, horizons = load(
                f,
                args.datPath[0],
                args.pickPath[0]
            )
            print(rdata, horizons)
            

            radargram(
                rdata,
                horizons,
                params,
                args.outPath[0],
            )

# execute if run as a script
if __name__ == "__main__":
    main()