# script designed to plot elevation bed and gnd profile
# BST 06DEC19
# python

import numpy as np
import matplotlib.pyplot as plt
import sys

file = "/zippy/MARS/orig/supl/gis/shapefiles/UAF/radar/picks/2019/sep/20190925-210807_pk_elev_gnd_arcticDEM.csv"


dat = np.genfromtxt(file, delimiter=",", dtype = None, names = True)

c_ice = 3e8/np.sqrt(3.15)

# linearly interp over spike from bad DEM
dat["elev_gnd"][18305:18325] = dat["elev_gnd"][18326]


twtt_surf = 2*(dat["elev_air"]-dat["elev_gnd"])/3e8
thick = (dat["twtt_bed"]-twtt_surf)*c_ice/2


elev_bed = dat["elev_gnd"] - thick
lon =dat["lon"]
lat = dat["lat"]




dist = np.zeros(len(thick))
for _i in range(len(dist)):
    if _i>=1:
        dist[_i] = dist[_i-1] + np.sqrt((lon[_i] - lon[_i-1])**2 + (lat[_i] - lat[_i-1])**2)

dist = dist*100     # convert to km

plt.rcParams.update({'font.size': 22})
plt.plot(dist[39000:-500],dat["elev_gnd"][39000:-500])
plt.plot(dist[39000:-500],elev_bed[39000:-500])
plt.title('Bering Glacier Cross Flow Profile')
plt.xlabel('Distance  along track (km)')
plt.ylabel('Elevation above msl (m)')
plt.show()



# plt.plot(dist[5000:23000],dat["elev_gnd"][5000:23000])
# plt.plot(dist[5000:23000],elev_bed[5000:23000])
# plt.title('Bering Glacier Along Flow Profile')
# plt.xlabel('Distance  along track (km)')
# plt.ylabel('Elevation above msl (m)')
# plt.show()