import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyproj as prj
import h5py
import glob
import argparse
import os
import ntpath
import struct
from datetime import datetime


def buildHd(file, ghog, z=np.array([-9999])):
    hd = {}

    step = float(ghog["restack"]["rx0"].attrs["stack_interval"].split(" ")[0])
    hd["NUMBER OF TRACES"] = ghog["restack"]["rx0"].shape[1]
    hd["NUMBER OF PTS/TRC"] = ghog["restack"]["rx0"].shape[0]
    hd["TIMEZERO AT POINT"] = "1"  # Figure out and fix this?
    hd["TOTAL TIME WINDOW"] = 1e9*hd["NUMBER OF PTS/TRC"]/ghog["raw"]["rx0"].attrs["fs"]
    hd["STARTING POSITION"] = 0
    hd["FINAL POSITION"] = ghog["restack"]["rx0"].shape[1]*step
    hd["STEP SIZE USED"] = step
    hd["POSITION UNITS"] = "m"
    hd["NOMINAL FREQUENCY"] = "1"
    hd["ANTENNA SEPARATION"] = np.median(np.round(ghog["restack"]["sep"]))
    hd["PULSER VOLTAGE (V)"] = "4000"
    hd["NUMBER OF STACKS"] = 1
    hd["SURVEY MODE"] = "Reflection"
    hd["STACKING TYPE"] = "F1, P" + str(hd["NUMBER OF STACKS"]) + ", DynaQ OFF"
    hd["TRIGGER MODE"] = "Manual"  # Fix to "Timed" or whatever the right word is
    hd["DATA TYPE"] = "F*4"
    hd["AMPLITUDE WINDOW (mV)"] = "FILLER"  # Put something here
    hd["TRACEHEADERDEF_26"] = "ORIENA"
    hd["GPR SERIAL#"] = "NA"
    hd["RX SERIAL#"] = "NA"
    hd["DVL SERIAL#"] = "NA"
    hd["TX SERIAL#"] = "NA"

    elevstr = "ELEVATION DATA ENTERED : MAX = %.2f MIN = %.2f" % (z.min(), z.max())

    fd = open(file, 'w')
    fd.write("420\nRuth - Gruondhog\n2020-Apr-23\n")
    for k, v in hd.items():
        fd.write("%s = %s\n" % (k, v))
        if(k == "STACKING TYPE"):
            fd.write(elevstr + "\n")
            fd.write("X Y Z POSITIONS ADDED - LatLong\n")

    fd.close()

    return 1


def buildDt1(file, ghog, atDist=[-9999], z=[-9999]):
    fd = open(file, 'wb')

    data = np.array(ghog["restack"]["rx0"]).astype(np.float32)

    data[np.isnan(data)] = 0

    #print("Skipping atDist check in buildDt1")
    if(atDist[0] == -9999):
        atDist = np.arange(data.shape[1]) * float(ghog["restack"]["rx0"].attrs["stack_interval"].split(" ")[0])

    #print("Skipping z check in buildDt1")
    if(z[0] == -9999):
        z = np.zeros(data.shape[1])

    # Scale data
    data = data * (50/(np.max(np.abs(data))))

    for i in range(data.shape[1]):
        # Write header
        hdr = b""
        hdr += struct.pack("f", i+1)  # Trace number
        hdr += struct.pack("f", atDist[i])  # Position along traverse
        hdr += struct.pack("f", data.shape[0])  # Num samples
        hdr += struct.pack("f", z[i])  # Elevation
        hdr += struct.pack("f", 0)  # Unassigned
        hdr += struct.pack("f", 4)  # Bytes per sample
        hdr += struct.pack("f", i+1)  # Aux trace number
        hdr += struct.pack("f", 1)  # Num stacks
        hdr += struct.pack("f", 1e9*data.shape[0]/ghog["raw"]["rx0"].attrs["fs"])  # Time window in ns
        for j in range(11):
            hdr += struct.pack("f", 0)  # Unassigned
        hdr += struct.pack("f", 0)  # Time 0 adjust
        hdr += struct.pack("f", 0)  # Data is zeros flag
        hdr += struct.pack("f", 0)  # Unassigned
        hdr += struct.pack("f", i)  # Time since midnight (fix)
        hdr += struct.pack("f", 0)  # Comment flag
        for j in range(7):
            hdr += struct.pack("f", 0)  # The comment

        fd.write(hdr)

        fd.write(data[:, i].astype(np.float32))

    fd.close()

    return 0


def buildGp2(file, cor):
    fd = open(file, 'w')

    # Header
    fd.write(";GPS@@@\n")
    fd.write(";Ver=1.1.0\n")
    fd.write(";DIP=2009-00152-00\n")
    fd.write(";DATE=" + cor["date"][0] + " " + cor["time"][0] + "\n")
    fd.write(";----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
    fd.write("traces,odo_tick,pos(m),time_elapsed(s),GPS\n")

    # Data
    for i, row in cor.iterrows():
        # Build GPGGA string
        dt = datetime.strptime(row.date+row.time, "%Y-%m-%d%H:%M:%S:%f")
        sow = ((dt.weekday() + 1) % 7) * 86400 + dt.hour*3600 + dt.minute*60 + dt.second + dt.microsecond*1e-6
        gpgga = "$GPGGA," + \
                str(sow) + ',' + \
                str(row.lat*100) + ',' + \
                row.latCard + ',' + \
                str(np.abs(row.lon)*100) + ',' + \
                row.lonCard + ',' + \
                "1,08,1.03," + \
                str(row.hgt) + ',' + \
                "M,0,M,00000*7B"


        fd.write(str(row.trace) + ',0,0.000000,0,')
        fd.write("\"" + gpgga + "\"\n")

    fd.close()

    return 0


def main():
    # Set up CLI
    parser = argparse.ArgumentParser(
        description="Program for converting a Groundhog HDF5 file to Sensors and Software dt1/hd format"
    )
    parser.add_argument("inpath", help="directory containing Groundhog HDF5 data files")
    parser.add_argument("outpath", help="directory for sensoft converted Groundhog data files")
    args = parser.parse_args()

    # loop through inpath
    files = glob.glob(args.inpath + '/*.h5')

    # Check that groundhog data files exists
    if len(files) == 0:
        print("No Groundhog HDF5 files foud at:\t", args.inpath)
        return 1
    if not os.path.isdir(args.outpath):
        print("Output path does not exist:\t", args.outpath)
        return 1

    for f in files:
        head, tail = ntpath.split(f)
        hdPath = args.outpath + tail.replace(".h5", ".hd")
        dt1Path = args.outpath + tail.replace(".h5", ".dt1")
        gp2Path = args.outpath + tail.replace(".h5", ".gp2")

        try:
            # Open groundhog file
            ghog = h5py.File(f)

            # Build Sensors & Software files
            print("Cutting Z out of buildHd")
            buildHd(hdPath, ghog) #, z)
            print("Cutting atDist and z out of Dt1")
            buildDt1(dt1Path, ghog) #, atDist, z)

            # print("Skipping buildGp2")
            #buildGp2(gp2Path, cor)
        except Exception as err:
            print("groundhog2sensoft error encountered for file:\t", tail, "\n", str(err))
            pass


main()
