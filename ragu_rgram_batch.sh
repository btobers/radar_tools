#!/usr/bin/env sh
# wrapper for ragu_rgram_batch.py


for yyyy in {2013..2021..1}
do
    python ragu_rgram_batch.py -f F:/MARS/orig/xtra/OIB-AK/radar/$yyyy/pk_bst/flist.txt -datpath F:/MARS/targ/xped/OIB-AK/radar/$yyyy/hdf5/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/$yyyy/pk_bst/ -outpath F:/MARS/orig/xtra/OIB-AK/radar/$yyyy/pk_bst/fig/ -sim
done


python ragu_rgram_batch.py -f F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRUAFHF2/flist.txt -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRUAFHF2/ -outpath F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRUAFHF2/ -sim
python ragu_rgram_batch.py -f F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRARES2/flist.txt -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRARES2/ -outpath F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRARES2/ -sim