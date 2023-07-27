#!/usr/bin/env python
# coding: utf-8
# %%

# %%


import batman
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def make_total_plot(data_split, nights, best_periods, initial_period, avg_first_and_last, t0_mcmc, rldc_r, field, chip, star_id):
    tcol=0
    magcol=1
    errcol=2
    bjdzp = 2458664.0
    
    time_r_lists = list(data_split['time_r'])
    time_r = [item for sublist in time_r_lists for item in sublist]
    mag_r_lists = list(data_split['mag_r'])
    mag_r = [item for sublist in mag_r_lists for item in sublist]
    err_r_lists = list(data_split['err_r'])
    err_r = [item for sublist in err_r_lists for item in sublist]
    
    lcr = pd.DataFrame({'time_r':time_r, 'mag_r':mag_r, 'err_r':err_r})
    
    #lcr = pd.read_csv(filename_r, sep='\s+', header=None)
    #lcr=lc[['time_r', 'mag_r', 'err_r']]
    #lcr.iloc[:,tcol] -= bjdzp
    time_r = np.array(lcr.iloc[:,tcol])
    
    time_r = np.array(time_r)
    #print(lcr)
    
    periods=best_periods.head(5)
    periods_list = list(periods['Period'])
    t0s_list = list(periods['t0'])
    
    #periods_list.append(0.9752)
    idx=0
    for period in periods_list:
        params_r = batman.TransitParams()
        params_r.t0 = t0s_list[0]#results['50th'][0]                     #time of inferior conjunction
        params_r.per = period                 #orbital period
        params_r.rp = avg_first_and_last['rp_r']                   #planet radius (in units of stellar radii)
        params_r.a =  avg_first_and_last['a'] * (period/initial_period)                      #semi-major axis (in units of stellar radii)
        params_r.inc = avg_first_and_last['inc']                   #orbital inclination (in degrees) 
        params_r.ecc = 0.0                    #eccentricity
        params_r.w =  0.0                       #longitude of periastron (in degrees) 
        params_r.u =  rldc_r                #limb darkening coefficients [u1, u2] 
        params_r.limb_dark = "quadratic"
        m_r=batman.TransitModel(params_r, time_r)
        flux_r = m_r.light_curve(params_r)
        ans_r = avg_first_and_last['C_r']-2.5*np.log10(flux_r)
        lcr[period]=ans_r
        idx+=1
    
    fig = plt.figure(figsize=(2*len(nights),26))
    for n, night in enumerate(nights):
            dr = lcr[(lcr.iloc[:,tcol]>night) & (lcr.iloc[:,tcol]<night+1.0)]
            #dz = lcz[(lcz.iloc[:,tcol]>night) & (lcz.iloc[:,tcol]<night+1.0)]
            #dr_cut = lcr_cut[(lcr_cut.iloc[:,tcol]>night) & (lcr_cut.iloc[:,tcol]<night+100.0)]
            #dz_cut = lcz_cut[(lcz_cut.iloc[:,tcol]>night) & (lcz_cut.iloc[:,tcol]<night+1.0)]


            ax = plt.subplot(5, int(np.ceil(len(nights)/5.0)), int(n+1))

            ax.errorbar(dr.iloc[:,tcol], dr.iloc[:,magcol], yerr=dr.iloc[:,errcol],
                        fmt='o', color='black', ms=1, label='r-band data')
            #i=3
            #for period in periods_list:
            ax.plot(dr.iloc[:,tcol], dr.iloc[:,3], linewidth=4.0, label=str(round(periods_list[0],5)), zorder=5)
            ax.plot(dr.iloc[:,tcol], dr.iloc[:,4], linewidth=4.0, label=str(round(periods_list[1],5)), zorder=4)
            ax.plot(dr.iloc[:,tcol], dr.iloc[:,5], linewidth=4.0, label=str(round(periods_list[2],5)), zorder=3)
            ax.plot(dr.iloc[:,tcol], dr.iloc[:,6], linewidth=4.0, label=str(round(periods_list[3],5)), zorder=2)
            ax.plot(dr.iloc[:,tcol], dr.iloc[:,7], linewidth=4.0, label=str(round(periods_list[4],5)), zorder=1)

            #ax.plot(dr.iloc[:,tcol], dr.iloc[:,8], linewidth=4.0, label=str(round(periods_list[5],5)))
            ax.set_title('Night ' +str(night))
            #ax.set_ylabel('Magnitude')
            #ax.set_xlabel('Time (days)')
            ax.set_xlim(float(night)+0.45, float(night)+1.0)
            #ax.set_ylim(avg_first_and_last['C_r']-0.03, avg_first_and_last['rp_r']**2+0.03)
            ax.invert_yaxis()
            #ax.legend()
        
    labels=[str(round(periods_list[0],5)), str(round(periods_list[1],5)), str(round(periods_list[2],5)), str(round(periods_list[3],5)), str(round(periods_list[4],5))]
    fig.legend(labels=labels, loc='center right', fontsize=20)
    fig.suptitle('Best periods plotted with r-band data', fontsize=45.0)
    plt.savefig(os.path.join('period_fit_results','total_plot_'+str(star_id)+'.pdf'))
    plt.close()

