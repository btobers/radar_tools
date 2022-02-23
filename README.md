# radar_tools
miscellaneous radar tools


This repository contains various useful radar scripts - mostly based on  [RAGU](https://github.com/btobers/RAGU) exports

### File Info
*OIB_track_query.py* opens a NASA OIB radar track vector file and export a list of profiles which fall within a user specified lat-lon bounding area

*OIB_tracklist_cp.sh* copies a tracklist produced by running *OIB_track_query.py* to a desired directory

*ragu_fence.ipynb* creates fence diagrams and animations for radargrams using RAGU

*ragu_pick_crossover.ipynb* performs a crossover analysis on RAGU radar picks

*ragu_picks_combine.py* combines all RAGU pick files in a given directory to a single text and geopackage export

*ragu_rgram.ipynb* creates compiled radargrams using RAGU
