"""
Restructure radar sounding data from Laurent Mingo's Blue Systems Air-IPR.
"""
### imports ###
import xml.etree.ElementTree as ET
import h5py
import numpy as np
import sys, glob, os, argparse


# Define a custom sorting function
def sort_by_numeric_value(s):
    return int(s.split('_')[1])  # Split the string by '_' and return the second part as an integer


def parseHeader():
    header = {}

    header["spt"] = 3200
    header["fs"] = 1
    header["pre_trig"] = 0
    # header["pre_trig"] = struct.unpack("q", data[12:20])[0]
    header["prf"] = 1000
    header["stack"] = 1
    # header["trig"] = struct.unpack("h", data[36:38])[0]
    # header["fs"] = struct.unpack("d", data[38:46])[0]

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
    
    # hardcode snum
    snum = 3200

    # read in .h5 file
    with h5py.File(fpath, "r")  as fd:
        lns = sorted(fd.keys(), key=sort_by_numeric_value)
        for ln in lns:
            line = fd[ln]
            # get header info 
            header = parseHeader()
            # initialize data array
            data = np.zeros((3200, len(line)))

            # initialize latitude and longitude lists variables to store values
            lats = []
            lons = []
            hgts = []
            times = []
            dates = []

            # loop through traces
            for key in line.keys():
                i = int(key.split("_")[1])
                # fill in data array
                data[:,i] = line[key]["datacapture_0/echogram_0"][:]

                # parse gps xml data
                gps_xml = ET.fromstring(line[key]["datacapture_0/echogram_0"].attrs['GPS Cluster- MetaData_xml'])
                # Iterate over all <String> elements and get lat/lon for each trace
                for string_elem in gps_xml.findall('./String'):
                    name_elem = string_elem.find('Name')
                    val_elem = string_elem.find('Val')
                    if name_elem is not None and val_elem is not None:
                        name = name_elem.text.strip()
                        val = val_elem.text.strip()
                        if name == 'Lat_N':
                            lats.append(float(val)/100)
                        elif name == 'Long_ W':  # Note the space in 'Long_ W'
                            lons.append(float(val)/100*-1)
                        elif name == 'Alt_asl_m': # Guessing this is in ref to WGS84 ellipsoid
                            hgts.append(float(val))
                        elif name == 'GPS_timestamp_UTC':
                            times.append(float(val))

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
                            np.nan,
                            # np.datetime_as_string(fix["times"][i]),
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
        # try:
        parse(fn)
        # except Exception as err:
        #     print(err)
        #     pass

# execute if run as a script
if __name__ == "__main__":
    main()