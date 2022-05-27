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
import sys, os, itertools, glob, argparse, fnmatch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Calibri"
plt.rcParams['font.size'] = 10

# append ragu path to access required modules
sys.path.append("C:/Users/btober/OneDrive/Documents/code/radar/ragu/code")
from ingest import ingest
from tools import utils

# function to load radargram and any accompanying picks
def load(f, datPath, pickPath):
    # ingest radar data
    fPath = datPath + "/" + f
    if os.path.isfile(fPath):
        igst = ingest(fPath)
        rdata = igst.read("", navcrs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", body="earth")
    else:
        print("Data file not found: {}".format(fPath))
        return None, None

    # initialize dict to hold any import horizons
    horizons = {}
    if pickPath:
        pickPath = pickPath + "/" + f[:-3] + "_pk_bst.csv"
        dat = pd.read_csv(pickPath)
        if dat.shape[0] != rdata.tnum:
            raise ValueError("import_pick error:\t pick file size does not match radar data")
        else:
            # get name of any horizons in pick file
            keys = fnmatch.filter(dat.keys(), "*sample*")
            for horizon in keys:
                # pick files should have ["horizon_sample"] keys
                horizon = horizon.split("_")[0]
                # add to horizon dict
                horizons[horizon] = dat[horizon + "_sample"].to_numpy()

    return rdata, horizons

# make a nice radargram
def radargram(rdata, horizons, params, outPath):
    # get number of horizons, if any
    nHzs = len(horizons)
    # declare number of figure panels
    nPanels = 3
    if rdata.flags.sim:
        # sim color limints
        vsmin = np.floor(np.nanpercentile(rdata.sim,10))
        vsmax = np.nanmax(rdata.sim)
    # if no sim, remove panel
    else:
        nPanels -= 1
    # if no horizons, remove panel
    if nHzs == 0:
        nPanels -=1

    # data color limits
    vdmin = np.floor(np.nanpercentile(rdata.proc.get_curr_dB(),10))
    vdmax = np.nanmax(rdata.proc.get_curr_dB())

    # radargram figure extent        
    extent = [0, rdata.tnum, rdata.snum, 0]
    if params["xAxis"] == "distance":
        extent[1] = rdata.navdf["dist"].iloc[-1]*1e-3

    if params["yAxis"] == "time":
        extent[2] = rdata.snum*rdata.dt*1e6
        # convert horizons from sample to microseconds
        for key, array in horizons.items():
            horizons[key] = array*rdata.dt*1e6

    # initialize figure
    fig = plt.figure(figsize=(params["pnlWidth"], nPanels*params["pnlHgt"]))
    ax0 = fig.add_subplot(111)   # big subplot
    fig.suptitle(rdata.fn)

    # initialize subplots
    ax = []
    if nPanels == 1:
        ax.append(fig.add_subplot(111))
    elif nPanels == 2:
        ax.append(fig.add_subplot(211))
        ax.append(fig.add_subplot(212))
    elif nPanels == 3:
        ax.append(fig.add_subplot(311))
        ax.append(fig.add_subplot(312))
        ax.append(fig.add_subplot(313))

    # Top panel - uninterpreted radargram
    ax[0].imshow(rdata.proc.get_curr_dB(), aspect="auto", extent=extent, cmap=params["cmap"], vmin=vdmin, vmax=vdmax)
    ax[0].set_ylim(int(extent[2]*params["yCutFact"]), 0)

    # Middle panel - clutter simulation, if present
    if rdata.flags.sim:
        ax[1].imshow(rdata.sim,  aspect="auto", extent=extent, cmap=params["cmap"], vmin=vsmin, vmax=vsmax)
        ax[1].set_ylim(int(extent[2]*params["yCutFact"]), 0)

    # Bottom plot - radargram with interpretations
    if nHzs > 0:
        ax[-1].imshow(rdata.proc.get_curr_dB(),  aspect="auto", extent=extent, cmap=params["cmap"], vmin=vdmin, vmax=vdmax)
        ax[-1].set_ylim(int(extent[2]*params["yCutFact"]), 0)
        for hz in horizons.values():
            ax[-1].plot(np.linspace(0,extent[1],rdata.tnum), hz)

    ### labels ###
    # Turn off axis lines and ticks of the big subplot
    ax0.set_yticks([extent[2], 0])
    ax0.spines['top'].set_color('none')
    ax0.spines['bottom'].set_color('none')
    ax0.spines['left'].set_color('none')
    ax0.spines['right'].set_color('none')
    ax0.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

    # create along-track distance secondary y axis
    for i, axis in enumerate(ax):
        axis.yaxis.set_ticks_position('both')
        axis.xaxis.set_ticks_position('both')

    
        # turn off xtick labels for all but bottom panel
        if i < nPanels - 1:
            axis.set_xticklabels([])
        else:
            if params["xAxis"] == "distance":
                axis.set_xlabel("Along-Track Distance (km)")
            elif params["xAxis"] == "trace":
                axis.set_xlabel("Trace")

    if params["yAxis"] == "sample":
        ax0.set_ylabel("Sample")
    elif params["yAxis"] == "time":
        ax0.set_ylabel(r'Two-Way Travel Time ($\mu$s)')

    fig.tight_layout()
    plt.subplots_adjust(hspace=0.125)
    # save figure
    fig.savefig(outPath + "/" + rdata.fn + ".png", dpi=500, bbox_inches='tight', pad_inches=0.05, transparent=True)
    # clear the figure
    plt.clf()
    plt.close("all")
    return


def main():

    # Set up CLI
    parser = argparse.ArgumentParser(
    description="Program for creating radargrams for a list of datafiles with corresponding pick files"
    )
    parser.add_argument("-f", dest = "fileList", help="List of files for which to create radargrams, alternately a single data file name", nargs="+")
    parser.add_argument("-datpath", dest = "datPath", help="Path to radar datafiles")
    parser.add_argument("-pkpath", dest = "pickPath", help="Path to RAGU radar pick files", nargs="?")
    parser.add_argument("-outpath", dest = "outPath", help="Path to to output generated radargrams")
    args = parser.parse_args()
    pkPath = args.pickPath

    # plotting parameters
    params = {}
    params["cmap"] = "Greys_r"                              # matplotlib.pyplot.imshow color map
    params["pnlHgt"] = 1.5                                  # panel height in inches for each panel in the generated radargram
    params["pnlWidth"] = 6.5                                # panel width in inches for each panel in the generated radargram
    params["yCutFact"] = 2/6                                 # factor by which to trim the bottom half of the radargram (.5 will preserve the upper half of the samples across the radargram)
    params["yAxis"] = "time"                                # y axis label unit ("sample" or "time")
    params["xAxis"] = "distance"                            # x axis label unit ("trace" or "distance")

    # check if paths exists
    if not os.path.isdir(args.datPath):
        print(f"Data file path not found: {args.datPath}")
        exit(1)

    flist = args.fileList
    if len(flist) == 1 and os.path.isfile(flist[0]):
        with open(flist[0]) as f:
            flist = [line.strip() for line in f]

    if pkPath:
        if not os.path.isdir(pkPath):
            pkPath = None
            print(f"Pick file path not found: {args.pickPath}")
            print("Radargrams will be generated without picks.")

    if not os.path.isdir(args.outPath):
        print(f"Output path not found: {args.outPath}")
        exit(1)


    # loop through list of files, load radargram and any picks, generate radargram
    for f in flist:
        rdata, horizons = load(
            f,
            args.datPath,
            pkPath
        )            

        radargram(
            rdata,
            horizons,
            params,
            args.outPath,
        )

# execute if run as a script
if __name__ == "__main__":
    main()