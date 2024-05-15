"""
Restructure radar sounding data from Laurent Mingo's Blue Systems Air-IPR.
"""
### imports ###
import xml.etree.ElementTree as ET
import h5py
import numpy as np
import sys, glob, os, argparse, re
from datetime import datetime, timedelta


# Define a custom sorting function
def sort_by_numeric_value(s):
    return int(s.split('_')[1])  # Split the string by '_' and return the second part as an integer


def dm2dd(dms):
    return ((dms - dms % 100) / 100 + (dms % 100) / 60)


def _2dt(date_obj, ms):
    # Convert milliseconds to seconds
    sec = ms / 1000
    # Calculate the datetime object for the given UTC time
    utc_datetime = datetime.combine(date_obj, datetime.min.time()) + timedelta(seconds=sec)
    return utc_datetime


def _xmlGetVal(xml, name):
    """
    Borrowed from dlilien's impDAR: Look up a value in an XML fragment. Mod from Nat Wilson's irlib.
    """
    m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>'.format(
        name.replace(' ', r'\s')), xml, flags=re.IGNORECASE)
    if m is not None:
        tail = xml[m.span()[1]:]
        end = tail.find('</Val')
        return tail[:end]
    else:
        return None


def parseHeader(dig_xml):
    header = {}
    header["stack"] = int(_xmlGetVal(dig_xml, 'Stacking'))
    header["spt"] = int(_xmlGetVal(dig_xml, 'Record Length'))
    header["fs"] = float(_xmlGetVal(dig_xml, ' Sample Rate'))
    header["pre_trig"] = 0
    header["prf"] = np.nan

    return header


def buildH5(header, rx, gps, outfile):
    """
    buildhdf5 format borrowed from Groundhog radar system, courtesy of Michael Christoffersen
    """
    try:
        print("Saving ", outfile)

        fd = h5py.File(outfile, "w")

        raw = fd.create_group("raw")
        rx0 = raw.create_dataset("rx0", data=rx)

        if gps is not None:
            raw.create_dataset("gps0", data=gps)

        for k, v in header.items():
            rx0.attrs[k] = v

        fd.close()
    except Exception as e:
        print("Failure in buildH5")
        print(e)
        return -1

    return 0


def parse(fpath=''):

    # read in .h5 file
    with h5py.File(fpath, "r")  as fd:
        # get date string from file as datetime obj
        datestr = fpath.split('/')[-1].split(' ')[-1][:-3]
        date_obj = datetime.strptime(datestr, "%b.%d.%Y") 
        lns = sorted(fd.keys(), key=sort_by_numeric_value)
        for ln in lns:
            line = fd[ln]
            # get header info 
            header = parseHeader(line[list(line.keys())[0]]["datacapture_0/echogram_0"].attrs['Digitizer-MetaData_xml'])
            # initialize data array
            data = np.zeros((header["spt"], len(line)))

            # initialize latitude and longitude lists variables to store values
            lats = []
            lons = []
            hgts = []
            times = []

            # loop through traces
            for key in line.keys():

                i = int(key.split("_")[1])
                # fill in data array
                data[:,i] = line[key]["datacapture_0/echogram_0"][:]

                # parse gps xml data
                gps_xml = line[key]["datacapture_0/echogram_0"].attrs['GPS Cluster- MetaData_xml']
                lats.append(dm2dd(float(_xmlGetVal(gps_xml, 'Lat_N'))))
                lons.append(dm2dd(float(_xmlGetVal(gps_xml, 'Long_ W'))) * -1)
                hgts.append(float(_xmlGetVal(gps_xml, 'Alt_asl_m')))
                times.append(np.datetime64(_2dt(date_obj, float(_xmlGetVal(gps_xml, 'GPS_timestamp_UTC')))))

            fix =  {"lons": lons, "lats": lats, "hgts": hgts, "times": times}
            tnum = data.shape[1]

            # get fix to right datattype for hdf5
            if fix is not None:
                gps_t = np.dtype(
                    [("lon", "f8"), ("lat", "f8"), ("hgt", "f8"), ("utc", "S26")]
                )
                gps = [None] * tnum
                for i in range(tnum):
                    gps[i] = tuple(
                        [
                            fix["lons"][i],
                            fix["lats"][i],
                            fix["hgts"][i],
                            np.datetime_as_string(fix["times"][i]),
                        ]
                    )
                gps = np.array(gps, dtype=gps_t)

            # new file name
            d, fn = os.path.split(fpath)
            file = d + '/' + fn[:-3] + '_' + ln + '.h5'

            # rebuild h5 format, following same structure as groundhog for simplicity
            if buildH5(header, data, gps, file) == -1:
                print("%s - Failed to build HDF5." % file)
                continue

            print()


def main():

    # Set up CLI
    parser = argparse.ArgumentParser(
    description="Program for parsing and reformatting Blue Systems Air-IPR radar sounding data files"
    )
    parser.add_argument("-datapath", help="path to Blue Systems HDF5 data files")
    args = parser.parse_args()

    for fn in glob.glob(args.datapath + '/*.h5'):
        try:
            parse(fn)
        except Exception as err:
            print(f"Error on file {fn}: {err}")
            pass

# execute if run as a script
if __name__ == "__main__":
    main()