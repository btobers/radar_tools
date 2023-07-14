#!/usr/bin/env sh
# wrapper for ragu_rgram_batch.py


for yyyy in {2013..2021..1}
do
    python ragu_rgram_batch.py -f F:/MARS/orig/xtra/OIB-AK/radar/$yyyy/pk_bst/flist.txt -datpath F:/MARS/targ/xped/OIB-AK/radar/$yyyy/hdf5/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/$yyyy/pk_bst/ -outpath F:/MARS/orig/xtra/OIB-AK/radar/$yyyy/pk_bst/fig/ -sim
done


python ragu_rgram_batch.py -f F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRUAFHF2/flist.txt -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRUAFHF2/ -outpath F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRUAFHF2/ -sim
# python ragu_rgram_batch.py -f F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRARES2/flist.txt -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRARES2/ -outpath F:/MARS/orig/xtra/OIB-AK/radar/NSIDC/IRARES2/ -sim

python ragu_rgram_batch.py -f IRARES1B_20210513-021438.h5 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -outpath C:/Users/btober/OneDrive/Documents/pres/OIB-AK/pub/fig/ -sim





python ragu_rgram_batch.py -f IRARES2_20190925-234004 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -sim
python ragu_rgram_batch.py -f IRARES2_20210512-221454 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/ -elev

python ragu_rgram_batch.py -f IRARES2_20210513-183506 -datpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/repo/data/ -pkpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/repo/data/vector/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/bagley_power_reduction/ -elev


python ragu_rgram_batch.py -f IRARES2_20210513-183506 -datpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/repo/data/ -pkpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/repo/data/vector/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/bagley_power_reduction/ -elev



# commands for figures used in oib paper
python ragu_rgram_batch.py -f IRARES2_20190923-180738 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/ak_range/ -sim -elev
python ragu_rgram_batch.py -f IRARES2_20170530-225322 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/sargent_harding/ -sim -elev
python ragu_rgram_batch.py -f IRARES2_20170530-224155 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/sargent_harding/ -sim -elev
python ragu_rgram_batch.py -f IRARES2_20190924-195652 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/central/ -sim -elev
python ragu_rgram_batch.py -f IRARES2_20180529-185120 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/wrangells/ -sim -elev
python ragu_rgram_batch.py -f IRARES2_20190927-224307 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/jif/ -sim -elev
python ragu_rgram_batch.py -f IRARES2_20190928-221557 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/southeast_1/ -sim -elev
python ragu_rgram_batch.py -f IRUAFHF2_20140516-205224 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRUAFHF2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/southeast_2/ -sim -elev
python ragu_rgram_batch.py -f IRARES2_20200527-220522 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/southeast_2/ -sim -elev
python ragu_rgram_batch.py -f IRARES2_20210513-195651 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/logan/ -sim -elev
python ragu_rgram_batch.py -f IRARES2_20210510-182722 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/logan/ -sim -elev

IRARES2_20210513-021438

IRARES2_20210510-182722


python ragu_rgram_batch.py -f IRARES2_20170530-225322 -datpath F:/MARS/targ/xped/OIB-AK/radar/ -pkpath F:/MARS/orig/xtra/OIB-AK/radar/nsidc/IRARES2/ -outpath C:/Users/btober/OneDrive/Documents/MARS/note/pres/OIB-AK/pub/fig/picks/sargent_harding/ -sim -elev
