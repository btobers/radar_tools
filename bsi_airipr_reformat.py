"""
Restructure radar sounding data from Laurent Mingo's Blue Systems Air-IPR.
"""
### imports ###
import xml.etree.ElementTree as ET
import h5py
import numpy as np
import sys, glob, os, argparse, re
from datetime import datetime, timedelta

# Interpolate nans in array
def interp_xords(arr):
    nan_indices = np.isnan(arr)
    non_nan_indices = ~nan_indices
    non_nan_values = arr[non_nan_indices]
    indices = np.arange(len(arr))
    # linearly interpolate NaN values
    arr[nan_indices] = np.interp(indices[nan_indices], indices[non_nan_indices], non_nan_values)
    return arr

# Define a custom sorting function
def sort_by_numeric_value(s):
    return int(s.split('_')[1])  # Split the string by '_' and return the second part as an integer


# Convert degrees minutes seconds to decimal degrees
def dm2dd(dms):
    return int(dms/100) + (dms - 100*int(dms/100))/60

# Convert from GPS UTC timestamp to datetime object
def _2dt(date_obj, hhmmss):
    h = hhmmss[:2]
    m = hhmmss[2:4]
    s = hhmmss[4:]
    sec = int(h) * 3600 + int(m) * 60 + float(s)
    # Calculate the datetime object for the given UTC time
    utc_datetime = datetime.combine(date_obj, datetime.min.time()) + timedelta(seconds=sec)
    return utc_datetime

# Pull value from XML string given key
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

# Build header
def buildHeader(dig_xml):
    header = {}
    header["stack"] = int(_xmlGetVal(dig_xml, 'Stacking'))
    header["spt"] = int(_xmlGetVal(dig_xml, 'RecordLength'))
    header["fs"] = float(_xmlGetVal(dig_xml, 'SampleRate'))
    toffset = float(_xmlGetVal(dig_xml, 'RelativeInitialX')) # this is the time offset recorded by the digitizer - will need to account for this and shift data accordingly 
    header["pre_trig"] = int(toffset * header['fs'])
    header["prf"] = 512     # hard coded for AirIPR

    return header

# Build Groundhog-format hdf5 file
def buildH5(header, rx, gps, agl, outfile):
    """
    buildhdf5 format borrowed from Groundhog radar system, courtesy of Michael Christoffersen
    """
    try:
        print("Saving ", outfile)

        fd = h5py.File(outfile, "w")

        raw = fd.create_group("raw")
        rx0 = raw.create_dataset("rx0", data=rx)
        string_t = h5py.string_dtype(encoding='ascii')

        if gps is not None:
            raw.create_dataset("gps0", data=gps)

        if agl is not None:
            ext = fd.create_group("ext")
            agl0 = ext.create_dataset("agl0", data=agl, dtype=np.float32)
            agl0.attrs.create("unit", "meter", dtype=string_t)
            description_attr = "IPR platform height above ground level (AGL), as seen by onboard FMCW AGL sensor."
            agl0.attrs.create("description", description_attr, dtype=string_t)

        for k, v in header.items():
            rx0.attrs[k] = v

        string_t = h5py.string_dtype(encoding='ascii')
        fd.attrs.create("system", "Blue Systems Integration IceRadar", dtype=string_t)

        fd.close()
    except Exception as e:
        print("Failure in buildH5")
        print(e)
        return -1

    return 0

# Parse AirIPR hdf5 file
def parse(fpath='', outpath=''):

    # read in .h5 file
    with h5py.File(fpath, "r")  as fd:
        # get date string from file as datetime obj
        datestr = fpath.split('/')[-1].split('_')[-1].split('.')[0]
        date_obj = datetime.strptime(datestr, "%Y%m%d") 
        lns = sorted(fd.keys(), key=sort_by_numeric_value)
        for ln in lns:
            line = fd[ln]
            # get header info 
            header = buildHeader(line[list(line.keys())[0]]["datacapture_0/echogram_0"].attrs['DigitizerMetaData_xml'])
            # initialize data array
            data = np.zeros((header["spt"], len(line)))

            # loop through traces
            traces = sorted(line.keys(), key=sort_by_numeric_value)

            # initialize latitude and longitude lists variables to store values
            lats = []
            lons = []
            hgts = []
            times = []
            agls = []

            # loop through traces
            for t in traces:
                tmp = t.split("_")
                if len(tmp)==2:
                    i=int(tmp[-1])
                else:
                    data = data[:,:-1]
                    continue

                # fill in data array
                data[:,i] = line[t]["datacapture_0/echogram_0"][:]

                # parse gps xml data
                gps_xml = line[t]["datacapture_0/echogram_0"].attrs['GPSData_xml']
                try:
                    lats.append(dm2dd(float(_xmlGetVal(gps_xml, 'Lat'))))
                except:
                    lats.append(np.nan)
                try:
                    lons.append(dm2dd(float(_xmlGetVal(gps_xml, 'Long'))))
                except:
                    lons.append(np.nan)
                try:
                    hgts.append(float(_xmlGetVal(gps_xml, 'Alt_ASL_m')))
                except:
                    hgts.append(np.nan)
                try:
                    times.append(np.datetime64(_2dt(date_obj, str(_xmlGetVal(gps_xml, 'GPSTimestamp_UTC')))))
                except:
                    times.append(np.datetime64('NaT'))

                # parse agl xml data
                agl_xml = line[t]["datacapture_0/echogram_0"].attrs['AGLData_xml']
                try:
                    agls.append(float(_xmlGetVal(agl_xml, 'Data')))
                except:
                    agls.append(np.nan)

            #############################
            ### filter erroneous data ###
            #############################
            lons_filt = np.asarray(lons)
            lats_filt = np.asarray(lats)
            hgts_filt = np.asarray(hgts)

            idxs = np.logical_or((lons_filt >= -10) & (lons_filt <= 10),
                                 (lats_filt >= -10) & (lats_filt <= 10),
                                 )

            lats_filt[idxs] = np.nan
            lons_filt[idxs] = np.nan
            hgts_filt[idxs] = np.nan

            # interpolate over nans
            lats_filt = interp_xords(lats_filt)
            lons_filt = interp_xords(lons_filt)
            hgts_filt = interp_xords(hgts_filt)

            lons_filt = lons_filt.tolist()
            lats_filt = lats_filt.tolist()
            hgts_filt = hgts_filt.tolist()
            #############################
            #############################
            #############################

            fix =  {"lons": lons_filt, "lats": lats_filt, "hgts": hgts_filt, "times": times}
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
            agl = np.array(agls, dtype=np.float32)
            # new file name
            d, fn = os.path.split(fpath)
            if outpath and os.path.isdir(outpath):
                d = outpath
            # fn = fn.split('.')[]
            file = d + '/' + fn.split('.')[0] + '_' + ln + '.h5'

            # shift data to account for 'RelativeInitialX' attribute in the digitizer xml record - this should be a negative number
            if header['pre_trig'] < 0:
                data = np.roll(data, shift=header['pre_trig'], axis=0)
                header['pre_trig'] = 0

            # rebuild h5 format, following same structure as groundhog for simplicity
            if buildH5(header, data, gps, agl, file) == -1:
                print("%s - Failed to build HDF5." % file)
                continue

            print()


def main():

    # Set up CLI
    parser = argparse.ArgumentParser(
    description="Program for parsing and reformatting Blue Systems Air-IPR radar sounding data files"
    )
    parser.add_argument("-datapath", help="path to Blue Systems HDF5 data files")
    parser.add_argument("-outpath", help="path to save reformatted HDF5 data files")
    args = parser.parse_args()

    for fn in glob.glob(args.datapath + '/*.hdf5'):
        if '_line_' not in fn:
            try:
                parse(fn, args.outpath)
            except Exception as err:
                print(f"Error on file {fn}: {err}")
                pass

# execute if run as a script
if __name__ == "__main__":
    main()