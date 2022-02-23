# radar_tools
miscellaneous radar tools


This repository contains various useful radar scripts - mostly based on  [RAGU](https://github.com/btobers/RAGU) exports

### File Info
*OIB_track_query.py* Opens a NASA OIB radar track vector file and export a list of profiles which fall within a user specified lat-lon bounding area

*OIB_tracklist_cp.sh* Copies a tracklist produced by running *OIB_track_query.py* to a desired directory

*ragu_fence.ipynb* Creates fence diagrams and animations for radargrams using RAGU

*ragu_pick_crossover.ipynb* Performs a crossover analysis on RAGU radar picks

*ragu_picks_combine.py* Combines all RAGU pick files in a given directory to a single text and geopackage export

*ragu_rgram.ipynb* Creates compiled radargrams using RAGU
