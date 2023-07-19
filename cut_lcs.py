from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sys.path.insert(0, '../')
from pyp2options import pip2options

# ##############################################################################

home = Path.home()
opt = pip2options(f"{home}/pip2/")
target = sys.argv[1]
runname = sys.argv[2]
use_detrended = int(sys.argv[3])
cut_start = sys.argv[4]
if len(sys.argv) == 6:
    cut_end = sys.argv[5]

ID = target.split('_')
field = f"{ID[0]}_{ID[1]}"
chip = ID[2]
num = ID[3].split('.')
lasttwo = num[0][-2:]

fc = f"{field}_{chip}"

if use_detrended == 1:
    tcol = 0
    magcol = 1
    errcol = 2
    lcext = 'lc.dtr'
elif use_detrended == 0:
    tcol=11
    magcol=14
    errcol=15
    lcext = 'lc'

plot=1

lcr = pd.read_csv(f"{opt['canddir']}/{field}/{runname}/{fc}/{target}/lightcurves/{target}.r.{lcext}",
                  sep='\s+', header=None)
lcz = pd.read_csv(f"{opt['canddir']}/{field}/{runname}/{fc}/{target}/lightcurves/{target}.z.{lcext}",
                  sep='\s+', header=None)

if (cut_start == 'list') or (cut_start == 'List') or (cut_start == 'LIST') or \
    (cut_start == 'l') or (cut_start == 'L'):
        
    cuts = np.loadtxt(f"{opt['canddir']}/{field}/{runname}/{fc}/{target}/{target}.cut_times",
                      delimiter=' ')
    lcr_cut=lcr.copy()
    lcz_cut = lcz.copy()
    for i in range(len(cuts)):
        lcr_cut = lcr_cut[(lcr_cut.iloc[:,tcol] < cuts[i,0]) |
                          (lcr_cut.iloc[:,tcol] > cuts[i,1])]
        lcz_cut = lcz_cut[(lcz_cut.iloc[:,tcol] < cuts[i,0]) |
                          (lcz_cut.iloc[:,tcol] > cuts[i,1])]
        print(lcr_cut.shape)

else:
    cuts = np.array([cut_start, cut_end], dtype='float')
    lcr_cut = lcr[lcr.iloc[:,tcol] < cuts[0] |
                  lcr.iloc[:,tcol] > cuts[1]]
    lcz_cut = lcz[lcz.iloc[:,tcol] < cuts[0] |
                  lcz.iloc[:,tcol] > cuts[1]]

lcr_cut.to_csv(f"{opt['canddir']}/{field}/{runname}/{fc}/{target}/lightcurves/{target}.r.{lcext}.cut",
               sep=' ', header=False, index=False)
lcz_cut.to_csv(f"{opt['canddir']}/{field}/{runname}/{fc}/{target}/lightcurves/{target}.z.{lcext}.cut",
               sep=' ', header=False, index=False)


if plot == 1:
    diff = lcr.iloc[:,magcol].median() - lcz.iloc[:,magcol].median()
    plt.errorbar(lcr_cut.iloc[:,tcol]-2458664, lcr_cut.iloc[:,magcol], 
                     yerr=lcr_cut.iloc[:,errcol], fmt='o', color='tab:blue', ms=3, 
                     rasterized=True, label='r-band')
    plt.errorbar(lcz_cut.iloc[:,tcol]-2458664, lcz_cut.iloc[:,magcol]+diff, 
                 yerr=lcz_cut.iloc[:,errcol], fmt='o', color='lightcoral', ms=3,
                 rasterized=True, label='z-band')
    plt.gca().invert_yaxis()
    plt.show()
