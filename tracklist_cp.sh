#!/bin/bash
# copy OIB tracks from tracklist to directory
fname="/home/btober/Documents/data/radar/kennicott/kennicott_tracks.txt"
cpath="/home/btobers/Documents/data/radar/kennicott/hdf5"

rsync -av --no-relative --files-from=$fname /  $cpath