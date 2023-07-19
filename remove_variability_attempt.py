import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from other_funcs import load_data
from pathlib import Path
import sys
import os


# +
# filename_r='MISHAPS_F1_N19_01027538.r.lc.dtr'
# filename_z='MISHAPS_F1_N19_01027538.z.lc.dtr'
# #filename_r='MISHAPS_F1_N19_01045451.r.lc.dtr'
# #filename_z='MISHAPS_F1_N19_01045451.z.lc.dtr'
# #period = 0.4913
# #filename_r='MISHAPS_F1_N19_01030151.r.lc.dtr'
# #filename_z='MISHAPS_F1_N19_01030151.z.lc.dtr'
# filename_data='MISHAPS_F1_N19_r_2019refim.search_output'
# #period = 2.8723

# starid=filename_r.split('.')[0]
# period = 1.9888
# ecc = 0.0
# w = 0.0
# rldc_r = [0.6502, 0.0872] #limb darkening
# rldc_z = [0.4678, 0.1118]
# exclude = None #[2,36,48,738,770]
# exclude_transit= None #[740]
# include_transit= None #[33]
# exclude = [2, 36,48,738,770]
# exclude_transit=None #[739]
# include_transit=None #[2]

# +
#load search output data
#star_data, transit_data, transit_nights, all_nights = load_data(filename_data, starid, exclude=exclude, exclude_transit=exclude_transit, include_transit=include_transit)
# -

def remove_variability_all(filename_r, filename_z, transit_data, transit_nights, all_nights, starid):
    #starid=filename_r.split('.')[0]
    if os.path.exists('remove_variability_results')==False:
        os.mkdir('remove_variability_results')
    else:
        for f in os.listdir('remove_variability_results'):
            os.remove(os.path.join('remove_variability_results', f))
            
    target = filename_r.split('.')[0][5:]
    if not os.path.exists('cut_times'):
            os.mkdir('cut_times')
    with open(os.path.join('cut_times', 'cut_times_'+str(target)+'.txt'), 'w') as f:
        for night in transit_nights:
            specific_data=transit_data[transit_data['int(night)']==night]
            specific_data=specific_data.head(1)
            duration=specific_data['best_duration[n]']
            tcbest=specific_data['tcbest[n]']
            cut_time_f=float(0.8*duration+tcbest+2458664)
            cut_time_i=float(tcbest-0.8*duration+2458664)
            f.write(str(cut_time_i))
            f.write(' ')
            f.write(str(cut_time_f))
            f.write('\n')
    f.close()

    cut_lcs(filename_r, filename_z, all_nights)
    remove_variability(all_nights, transit_nights, filename_r, filename_z)


# -

# +
""
def cut_lcs(filename_r, filename_z, nights):
    target = filename_r.split('.')[0][5:]
    ID = target.split('_')
    field = f"{ID[0]}_{ID[1]}"
    chip = ID[2]
    num = ID[3].split('.')
    lasttwo = num[0][-2:]
    fc = f"{field}_{chip}"
    tcol = 0
    magcol = 1
    errcol = 2
    lcext = 'lc.dtr'

    lcr = pd.read_csv(filename_r, sep='\s+', header=None)
    lcz = pd.read_csv(filename_z, sep='\s+', header=None)

    cuts = np.loadtxt(os.path.join('cut_times', 'cut_times_'+str(target)+'.txt'),
                          delimiter=' ')
    #print(len(lcr))
    #print(len(lcz))
    tcol=0
    magcol=1
    errcol=2
    bjdzp = 2458664.0
    
    
    lcr_cut=lcr.copy()
    lcz_cut=lcz.copy()
    for i in range(len(cuts)):
        lcr_cut = lcr_cut[(lcr_cut.iloc[:,tcol] < cuts[i,0]) |
                              (lcr_cut.iloc[:,tcol] > cuts[i,1])]
        lcz_cut = lcz_cut[(lcz_cut.iloc[:,tcol] < cuts[i,0]) |
                              (lcz_cut.iloc[:,tcol] > cuts[i,1])]
        #print(lcr_cut.shape)
    #print(len(lcr_cut))
    #print(len(lcz_cut))
    lcr_cut.to_csv(filename_r+".cut",
                   sep=' ', header=False, index=False)
    lcz_cut.to_csv(filename_z+".cut",
                   sep=' ', header=False, index=False)
    
    
    lcr.iloc[:,tcol] -= bjdzp
    lcz.iloc[:,tcol] -= bjdzp
    
    full_medr = lcr.iloc[:,magcol].median()
    full_medz = lcz.iloc[:,magcol].median()
    full_diff  = full_medr- full_medz
    fig = plt.figure(figsize=(2*len(nights),26))
    for n, night in enumerate(nights):
        dr = lcr[(lcr.iloc[:,tcol]>night) & (lcr.iloc[:,tcol]<night+1.0)]
        dz = lcz[(lcz.iloc[:,tcol]>night) & (lcz.iloc[:,tcol]<night+1.0)]
        #dr_cut = lcr_cut[(lcr_cut.iloc[:,tcol]>night) & (lcr_cut.iloc[:,tcol]<night+100.0)]
        #dz_cut = lcz_cut[(lcz_cut.iloc[:,tcol]>night) & (lcz_cut.iloc[:,tcol]<night+1.0)]
        

        ax = plt.subplot(5, int(np.ceil(len(nights)/5.0)), int(n+1))

        ax.errorbar(dr.iloc[:,tcol], dr.iloc[:,magcol], yerr=dr.iloc[:,errcol],
                    fmt='o', color='blue', ms=1)
        ax.errorbar(dz.iloc[:,tcol], dz.iloc[:,magcol]+full_diff, yerr=dz.iloc[:,errcol],
                      fmt='o', color='red', ms=1)
        #ax.set_xlim(night, night+1.0)

#         ax.errorbar(dr.iloc[:,tcol], dr.iloc[:,magcol]-ntrend_r[:]+med_diff_r, yerr=dr.iloc[:,errcol],
#                     fmt='o', color='tab:blue', ms=3)
#         ax.errorbar(dz.iloc[:,tcol], dz.iloc[:,magcol]-ntrend_z[:]+full_diff+med_diff_z, yerr=dz.iloc[:,errcol],
#                     fmt='o', color='lightcoral', ms=3)
        # ax.plot(dr.iloc[:,tcol], ntrend_r+nmed_r)
        # ax.axhline(y=lcr.iloc[:,magcol].median())
        # ax.axhline(y=2*full_medr-nmed_r)
        ax.invert_yaxis()
    fig.suptitle('Before Variability Removal', fontsize=45.0)
    plt.savefig(os.path.join('remove_variability_results','var_before_removal_check.pdf'))
    plt.close()
    
    
    
    
    
    full_medr = lcr.iloc[:,magcol].median()
    full_medz = lcz.iloc[:,magcol].median()
    full_diff  = full_medr- full_medz
    
    
    lcr_cut.iloc[:,tcol] -= bjdzp
    lcz_cut.iloc[:,tcol] -= bjdzp
    fig = plt.figure(figsize=(2*len(nights),26))
    for n, night in enumerate(nights):
        #dr = lcr_cut[(lcr_cut.iloc[:,tcol]>night) & (lcr_cut.iloc[:,tcol]<night+1.0)]
        #dz = lcz_cut[(lcz_cut.iloc[:,tcol]>night) & (lcz_cut.iloc[:,tcol]<night+1.0)]
        dr_cut = lcr_cut[(lcr_cut.iloc[:,tcol]>night) & (lcr_cut.iloc[:,tcol]<night+1.0)]
        dz_cut = lcz_cut[(lcz_cut.iloc[:,tcol]>night) & (lcz_cut.iloc[:,tcol]<night+1.0)]
        
#         nmed_r = dr_cut.iloc[:,magcol].median()
#         nmed_z = dz_cut.iloc[:,magcol].median()

#         med_diff_r = full_medr - nmed_r
#         med_diff_z = full_medz - nmed_z

        #ropt, rcov = curve_fit(fit_func, dr_cut.iloc[:,tcol], dr_cut.iloc[:,magcol]-nmed_r)
        #ntrend_r = fit_func(dr.iloc[:,tcol], *ropt)

        #zopt, zcov = curve_fit(fit_func, dz_cut.iloc[:,tcol], dz_cut.iloc[:,magcol]-nmed_z)
        #ntrend_z = fit_func(dz.iloc[:,tcol], *zopt)    

        ax = plt.subplot(5, int(np.ceil(len(nights)/5.0)), int(n+1))

        ax.errorbar(dr_cut.iloc[:,tcol], dr_cut.iloc[:,magcol], yerr=dr_cut.iloc[:,errcol],
                    fmt='o', color='blue', ms=1)
        ax.errorbar(dz_cut.iloc[:,tcol], dz_cut.iloc[:,magcol]+full_diff, yerr=dz_cut.iloc[:,errcol],
                      fmt='o', color='red', ms=1)
        ax.set_xlim(night+0.45, night+1.0)

#         ax.errorbar(dr.iloc[:,tcol], dr.iloc[:,magcol]-ntrend_r[:]+med_diff_r, yerr=dr.iloc[:,errcol],
#                     fmt='o', color='tab:blue', ms=3)
#         ax.errorbar(dz.iloc[:,tcol], dz.iloc[:,magcol]-ntrend_z[:]+full_diff+med_diff_z, yerr=dz.iloc[:,errcol],
#                     fmt='o', color='lightcoral', ms=3)
        # ax.plot(dr.iloc[:,tcol], ntrend_r+nmed_r)
        # ax.axhline(y=lcr.iloc[:,magcol].median())
        # ax.axhline(y=2*full_medr-nmed_r)
        ax.invert_yaxis()
    fig.suptitle('Lightcurves With Transits Cut Out', fontsize=45.0)
    plt.savefig(os.path.join('remove_variability_results','var_cut_lc_check.pdf'))
    plt.close()

# +


#Ali's code from GitHub

from pathlib import Path
import sys

from scipy.stats import median_abs_deviation as MAD
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

"""
##############################################################################
##############################################################################
"""


def fit_func(x, m, b):
        return m*x + b
    
    
def remove_variability(nights, t_nights, filename_r, filename_z):
    save_new=1


    tcol=0
    magcol=1
    errcol=2

    bjdzp = 2458664.0

    # Full lightcurves
    lcr = pd.read_csv(filename_r, sep='\s+',
                      header=None)
    lcz = pd.read_csv(filename_z, sep='\s+',
                      header=None)
    # Lightcurves with transits removed
    lcr_cut = pd.read_csv(filename_r+'.cut', sep='\s+',
                      header=None)
    lcz_cut = pd.read_csv(filename_z+'.cut', sep='\s+',
                      header=None)
    #print(len(lcr))
    #print(len(lcz))
    lcr.iloc[:,tcol] -= bjdzp
    lcz.iloc[:,tcol] -= bjdzp
    lcr_cut.iloc[:,tcol] -= bjdzp
    lcz_cut.iloc[:,tcol] -= bjdzp

    full_medr = lcr.iloc[:,magcol].median()
    full_medz = lcz.iloc[:,magcol].median()
    full_diff  = full_medr- full_medz

    # Filter for transit datapoints
    transits_r = lcr[(lcr.iloc[:,tcol].isin(lcr_cut.iloc[:,tcol])==False)]
    transits_z = lcz[(lcz.iloc[:,tcol].isin(lcz_cut.iloc[:,tcol])==False)]

    # Filter for out of transit datapoints
    oot_r = lcr[(lcr.iloc[:,tcol].isin(lcr_cut.iloc[:,tcol])==True)]
    oot_z = lcz[(lcz.iloc[:,tcol].isin(lcz_cut.iloc[:,tcol])==True)]

    ootmed_r = oot_r.iloc[:,magcol].median()
    ootmed_z = oot_z.iloc[:,magcol].median()

    ###############################################################################
    ###############################################################################
    ###############################################################################
    # Use cut lightcurves to find the variability trend in each year? Week? Day?
    #plt.style.use('dark_background')
    

    fig = plt.figure(figsize=(2*len(nights),26))
    corrected_r = []
    corrected_z = []
    for n, night in enumerate(nights):
        dr = lcr[(lcr.iloc[:,tcol]>night) & (lcr.iloc[:,tcol]<night+1.0)]
        dz = lcz[(lcz.iloc[:,tcol]>night) & (lcz.iloc[:,tcol]<night+1.0)]
        dr_cut = lcr_cut[(lcr_cut.iloc[:,tcol]>night) & (lcr_cut.iloc[:,tcol]<night+1.0)]
        dz_cut = lcz_cut[(lcz_cut.iloc[:,tcol]>night) & (lcz_cut.iloc[:,tcol]<night+1.0)]

        nmed_r = dr_cut.iloc[:,magcol].median()
        nmed_z = dz_cut.iloc[:,magcol].median()

        med_diff_r = full_medr - nmed_r
        med_diff_z = full_medz - nmed_z

        ropt, rcov = curve_fit(fit_func, dr_cut.iloc[:,tcol], dr_cut.iloc[:,magcol]-nmed_r)
        ntrend_r = fit_func(dr.iloc[:,tcol], *ropt)

        zopt, zcov = curve_fit(fit_func, dz_cut.iloc[:,tcol], dz_cut.iloc[:,magcol]-nmed_z)
        ntrend_z = fit_func(dz.iloc[:,tcol], *zopt)    

        ax = plt.subplot(5, int(np.ceil(len(nights)/5.0)), int(n+1))

        ax.errorbar(dr.iloc[:,tcol], dr.iloc[:,magcol], yerr=dr.iloc[:,errcol],
                    fmt='o', color='wheat', ms=3)
        ax.errorbar(dz.iloc[:,tcol], dz.iloc[:,magcol]+full_diff, yerr=dz.iloc[:,errcol],
                      fmt='o', color='pink', ms=3)

        ax.errorbar(dr.iloc[:,tcol], dr.iloc[:,magcol]-ntrend_r[:]+med_diff_r, yerr=dr.iloc[:,errcol],
                    fmt='o', color='blue', ms=3)
        ax.errorbar(dz.iloc[:,tcol], dz.iloc[:,magcol]-ntrend_z[:]+full_diff+med_diff_z, yerr=dz.iloc[:,errcol],
                    fmt='o', color='red', ms=3)
        # ax.plot(dr.iloc[:,tcol], ntrend_r+nmed_r)
        # ax.axhline(y=lcr.iloc[:,magcol].median())
        # ax.axhline(y=2*full_medr-nmed_r)
        ax.invert_yaxis()

        corrected_r.append((dr.iloc[:,magcol]-ntrend_r[:]+med_diff_r).to_numpy())
        corrected_z.append((dz.iloc[:,magcol]-ntrend_z[:]+full_diff+med_diff_z).to_numpy())
    fig.suptitle('After Variability Removal', fontsize=45.0)
    plt.savefig(os.path.join('remove_variability_results','variability_removal_check.pdf'))
    plt.close()
    
    full_corrected_r = np.concatenate(corrected_r)
    full_corrected_z = np.concatenate(corrected_z)

    fcr = pd.DataFrame(list(zip(lcr.iloc[:,tcol]+bjdzp, full_corrected_r, lcr.iloc[:,errcol])))
    fcz = pd.DataFrame(list(zip(lcz.iloc[:,tcol]+bjdzp, full_corrected_z, lcz.iloc[:,errcol])))
    
    #print(len(fcr))
    #print(len(fcz))
    if save_new == 1:
        fcr.to_csv(filename_r+'.corrected2', sep=' ', 
                   header=False, index=False)
        fcz.to_csv(filename_z+'.corrected2', sep=' ' ,
                   header=False, index=False)




# +
#remove_variability_all(filename_r, filename_z, transit_data, transit_nights, all_nights, starid)
# +
#transit_nights
# -


