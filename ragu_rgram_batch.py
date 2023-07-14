""" 
Batch radargram script
 
This script loops through a directory containing radar datafiles and RAGU generated pick files and creates a radargram for each set

Created by: Brandon S. Tober
Date: 22222522

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
import matplotlib.gridspec as gridspec
plt.rcParams["font.family"] = "Arial"
plt.rcParams['font.size'] = 8

# append ragu path to access required modules
sys.path.append("C:/Users/btober/OneDrive/Documents/code/radar/ragu/code")
from ingest import ingest
from tools import utils

# get_ax_size returns the size of a matplotlib axis
def get_ax_size(fig, ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height

# function to load radargram and any accompanying picks
def load(f, datPath, pickPath, elevFlag):
    # ingest radar data
    fPath = datPath + f
    if os.path.isfile(fPath):
        igst = ingest(fPath)
        rdata = igst.read("", navcrs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", body="earth")
    else:
        fn = f.replace('IRUAFHF2','IRUAFHF1B').replace('IRARES2','IRARES1B')
        fPath = datPath + fn + '.h5'
        if os.path.isfile(fPath):
            igst = ingest(fPath)
            rdata = igst.read("", navcrs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", body="earth")
        else:
            print("Data file not found: {}".format(fPath))
            return None, None, None

    # initialize dict to hold any import horizons
    horizon_dict = {}
    elev_dict = {}
    if pickPath:
        pickPath = pickPath + "/" + f + ".csv"
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
                horizon_dict[horizon] = dat[horizon + "_sample"].to_numpy()
                if elevFlag:
                    elev_dict[horizon] = dat[horizon + "_hgt"].to_numpy()
                    # if horizon=='bed':
                    #     elev_dict['bed_amp'] = dat['bed_amp']
                    #     elev_dict['ice_thickness'] = dat['srf_bed_thick']
    else:
        # still retain surface if in rdata data file
        horizon_dict["srf"] = rdata.pick.horizons[rdata.pick.get_srf()]

    return rdata, horizon_dict, elev_dict

# make a nice radargram
def radargram(rdata, horizon_dict, elev_dict, params, simFlag, outPath):
    if not simFlag:
        rdata.flags.sim = False

    # determine number of figure panels
    nPanels = 4
    if rdata.flags.sim:
        # sim color limints
        vsmin = np.floor(np.nanpercentile(rdata.sim,10))
        vsmax = np.nanmax(rdata.sim)
    # if no sim, remove panel
    else:
        nPanels -= 1

    # determine if interpreted radargram panel is being created
    srfTop = False          # flag to indicat that surface horizon will be on upper panel
    if len(horizon_dict) == 0:
        horizon_dict = None
        nPanels -=1
    # if only a surface sample array exists, we'll plot it atop the top panel rather than adding an interpretation panel
    elif (len(horizon_dict) == 1) and (list(horizon_dict.keys())[0] == "srf"):
        srfTop = True
        nPanels -=1

    # determine if elevation panel is being created
    if len(elev_dict) == 0:
        elev_dict = None
        nPanels -= 1

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
        if horizon_dict:
            for key, array in horizon_dict.items():
                horizon_dict[key] = array*rdata.dt*1e6

    # set subplot relative heights
    hgt_ratios = []
    for i in range(nPanels):
        if elev_dict is not None:
            if i == nPanels - 1:
                hgt_ratios.append(1)
            else:
                hgt_ratios.append(1)
        else:
            hgt_ratios.append(1)

    # initialize figure
    fig, ax = plt.subplots(
                                        figsize=(params["pnlWidth"], nPanels*params["pnlHgt"]), 
                                        nrows = nPanels,
                                        ncols = 1,
                                        gridspec_kw={'height_ratios': hgt_ratios}
                                        )
    # ax[0].set_title(rdata.fn)

    # uninterpreted radargram
    ax[0].imshow(rdata.proc.get_curr_dB(), aspect="auto", extent=extent, cmap=params["cmap"], vmin=vdmin, vmax=vdmax)
    ax[0].set_ylim(int(extent[2]*params["yCutFact"]), 0)
    ax[0].set_xlim(extent[1]*params["xCutFact"][0], extent[1]*params["xCutFact"][1])

    # clutter simulation, if present
    if rdata.flags.sim:
        ax[1].imshow(rdata.sim,  aspect="auto", extent=extent, cmap=params["cmap"], vmin=vsmin, vmax=vsmax)
        ax[1].set_ylim(int(extent[2]*params["yCutFact"]), 0)
        ax[1].set_xlim(extent[1]*params["xCutFact"][0], extent[1]*params["xCutFact"][1])

    # radargram with interpretations, if present
    if horizon_dict is not None:
        # get appropriate panel number
        if srfTop:
            pnl = 0
        if not srfTop:
            pnl = 1
            if rdata.flags.sim:
                pnl = 2 
        ax[pnl].imshow(rdata.proc.get_curr_dB(),  aspect="auto", extent=extent, cmap=params["cmap"], vmin=vdmin, vmax=vdmax)
        for (name, arr)  in horizon_dict.items():
            ax[pnl].plot(np.linspace(0,extent[1],rdata.tnum), arr, lw = 1, label=name)
        ax[pnl].set_ylim(int(extent[2]*params["yCutFact"]), 0)
        ax[pnl].set_xlim(extent[1]*params["xCutFact"][0], extent[1]*params["xCutFact"][1])
        # ax[pnl].legend(labels = ['Surface','Bed'], fancybox=False, frameon=False, borderaxespad=0, loc=params["lgd_pos"], edgecolor=None, handlelength=0.8, frameon=True, framealpha=.6, ncol=2, columnspacing=1) 
        # leg = ax[pnl].legend(labels = ['Surface','Bed'], fancybox=False, frameon=True, framealpha=.6, borderaxespad=0, loc='lower right', edgecolor=None, handlelength=0.8, ncol=2,columnspacing=1) 
        # leg.get_frame().set_edgecolor('k')

    # elevation profile, if flagged
    if elev_dict is not None:
        # if elevation drops below zero, add 0 horizontal line
        hline = False
        for (name, arr)  in elev_dict.items():
            # if name == 'bed_amp':
            #     a2 = ax[-1].twinx()
            #     arr = pd.Series(arr).rolling(100).mean()
            #     a2.plot(np.linspace(0,extent[1],rdata.tnum),20*np.log10(arr),lw=1, ls='-',color='r')
            #     a2.set_ylabel('Return Power (dB)')
            #     a2.yaxis.label.set_color('red')
            #     a2.tick_params(axis='y', colors='red')
            #     a2.spines['right'].set_color('red')
                # a2.set_ylim([-70,130])
                # a2.set_yticks([-50,0])
                # a2.yaxis.label.set_color('C1')
                # a2.tick_params(axis='y', colors='C1')
                # a2.spines['right'].set_color('C1')
            # elif name == 'ice_thickness':
            #     ax[-1].plot(np.linspace(0,extent[1],rdata.tnum), arr, lw = 1, c='k')
            # # else:
            ax[-1].plot(np.linspace(0,extent[1],rdata.tnum), arr, lw = 1, label=name)
            if (arr<0).sum() > 0:
                hline = True
                # hline=False
        if hline:
            ax[-1].axhline(0, ls='--',c='k', alpha=.2,label='_nolegend_')           # zero-m WGS84 elevation (roughly sea level) 
        ax[-1].legend(labels = ['Surface','Bed'], fancybox=False, frameon=False, borderaxespad=0, loc=params["lgd_pos"], edgecolor=None, handlelength=0.8)#, ncol=2,columnspacing=1) 
        ax[-1].set_xlim(extent[1]*params["xCutFact"][0], extent[1]*params["xCutFact"][1])

        # get vertical exag. for elev profile - only display if greater than 1
        if params["xAxis"] == "distance":
            dx = (1e3*(ax[-1].get_xlim()[1]-ax[-1].get_xlim()[0]))
            dy = (ax[-1].get_ylim()[1]-ax[-1].get_ylim()[0])
            xl,yl = get_ax_size(fig, ax[-1])
            ve = round((dx/xl)/(dy/yl))
            if ve > 1:
                ax[-1].annotate("VE = " + str(ve) + "x",xy=([.025,.04]), xycoords = "axes fraction")



    ### labels ###
    # make hidden axes to hold axis labels so we don't have duplicated ylabels - will be two if elevation profile is being created
    if elev_dict is not None:
        spec = fig.add_gridspec(
                                nrows = 2,
                                ncols = 1,
                                height_ratios = [(nPanels - 1)*(0.125 + hgt_ratios[0]),hgt_ratios[-1]]
                                )
        ax0 = fig.add_subplot(spec[0])
        ax1 = fig.add_subplot(spec[1])

        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.spines['top'].set_color('none')
        ax1.spines['bottom'].set_color('none')
        ax1.spines['left'].set_color('none')
        ax1.spines['right'].set_color('none')
        ax1.tick_params(labelcolor='1', top=False, bottom=False, left=False, right=False)
        ax1.set_ylabel("Elevation (m)", labelpad=22)

    else:
        ax0 = fig.add_subplot(111)

    # Turn off axis lines and ticks of the big subplot
    ax0.set_xticklabels([])
    ax0.set_yticklabels([])
    ax0.spines['top'].set_color('none')
    ax0.spines['bottom'].set_color('none')
    ax0.spines['left'].set_color('none')
    ax0.spines['right'].set_color('none')
    ax0.tick_params(labelcolor='1', top=False, bottom=False, left=False, right=False)

    # set axes for all subplots
    for i, axis in enumerate(ax):
        axis.yaxis.set_ticks_position('both')
        axis.xaxis.set_ticks_position('both')
        if i==nPanels - 1:
            axis.yaxis.set_ticks_position('both')
            lim = axis.get_ylim()
            axis.set_yticks([250,500,750,1000])
            # axis.set_ylim(400,1500)
            # axis.set_ylim(lim[0],1000)

        # turn off xtick labels for all but bottom panel
        if i < nPanels - 1:
            axis.set_xticklabels([])
        else:
            # axis.set_xlim([extent[0], extent[1]])
            axis.set_xlim(extent[1]*params["xCutFact"][0], extent[1]*params["xCutFact"][1])

            if params["xAxis"] == "distance":
                axis.set_xlabel("Along-Track Distance (km)")
            elif params["xAxis"] == "trace":
                axis.set_xlabel("Trace")

    if params["yAxis"] == "sample":
        ax0.set_ylabel("Sample", labelpad=22)
    elif params["yAxis"] == "time":
        ax0.set_ylabel(r'Two-Way Travel Time ($\mu$s)', labelpad=22)
            
    try:
        fig.tight_layout()
    except Exception:
        pass
    plt.subplots_adjust(hspace=0.13)#,left=0.1,right=0.9)
    # plt.tight_layout()
    # save figure
    fig.savefig(outPath, dpi=400, facecolor='white', transparent=True, bbox_inches="tight")
    # clear the figure
    plt.clf()
    plt.close("all")
    print(f'Radargram saved to:\t{outPath}')
    return


def main():

    # Set up CLI
    parser = argparse.ArgumentParser(
    description="""description: program for creating radargrams for a list of datafiles with corresponding pick files\n\nexample call: $python ragu_rgram_batch.py -f IRARES1B_20180819-215227.h5 -datpath /home/user/data/radar/ -pkpath /home/user/data/radar/pick/ -outpath /home/user/pres/ -elev -sim""",
    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", dest = "fileList", help="List of files for which to create radargrams, alternately a single data file name", nargs="+")
    parser.add_argument("-datpath", dest = "datPath", help="Path to radar datafiles")
    parser.add_argument("-pkpath", dest = "pickPath", help="Path to RAGU radar pick files", nargs="?")
    parser.add_argument("-outpath", dest = "outPath", help="Path to to output generated radargrams")
    parser.add_argument("-elev", help="Flag: Include elevation profile", default=False, action="store_true")
    parser.add_argument("-sim", help="Flag: Include clutter simulation if present", default=False, action="store_true")
    args = parser.parse_args()
    pkPath = args.pickPath

    # plotting parameters
    params = {}
    params["cmap"] = "Greys_r"                              # matplotlib.pyplot.imshow color map
    params["pnlHgt"] = 1                                    # panel height in inches for each panel in the generated radargram
    params["pnlWidth"] = 3                                # panel width in inches for each panel in the generated radargram
    params["yCutFact"] = 1/3                                # factor by which to trim the bottom portion of the radargram (.5 will preserve the upper half of the samples across the radargram)
    params["xCutFact"] = (.1,.55)                              # tuple, (left,right) factors by which to trim the radargram (0,1) will keep all traces
    params["yAxis"] = "time"                                # y axis label unit ("sample" or "time")
    params["xAxis"] = "distance"                            # x axis label unit ("trace" or "distance")
    params["lgd_pos"] = "upper left"                        # matplotlib legend position in subplot

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

        # parse appropriate year
        yr = f.split('_')[1][:4]

        rdata, horizons_dict, elev_dict = load(
            f,
            args.datPath + yr+ '/hdf5/',
            pkPath,
            args.elev
        )            
        outPath = args.outPath + f + '.png'

        radargram(
            rdata,
            horizons_dict,
            elev_dict,
            params,
            args.sim,
            outPath,
        )

# execute if run as a script
if __name__ == "__main__":
    main()