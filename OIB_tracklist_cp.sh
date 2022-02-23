#!/usr/bin/env sh
# copy OIB tracks from tracklist to directory of choice
# this bash script takes in a list of OIB tracks created using OIB_track_query.py

fname="/mnt/c/Users/btober/Desktop/ber-bag-mala/ber-bag-mala_2021_tracks.txt"

cpath="/mnt/c/Users/btober/Desktop/ber-bag-mala/hdf5"

rsync -av --no-relative --files-from=$fname /  $cpath