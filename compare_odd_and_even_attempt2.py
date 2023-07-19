#!/usr/bin/env python
# coding: utf-8
# %%

# %%


import batman
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import emcee
import corner
from PyAstronomy.pyasl import foldAt

from other_funcs import weighted_avg_and_std, difference, folded_lightcurve, odd_even_data
#from run_batman_inc_combined import run_batman
import mcmc_functions_inc_combined as mcmc_functions
from run_batman_folded import run_batman

def odd_even_data_new(data, t0_all, period, all_nights):
    t0_all=np.array(t0_all)
    t0_new = t0_all - t0_all[0]
    #print(t0_new)
    div=[round(t/period) for t in t0_new]
    #print(t0_all)
    #print(div)
    
    
    time_transits_r = list(data['time_r'])
    time_transits_z = list(data['time_z'])
    mag_transits_r = list(data['mag_r'])
    mag_transits_z = list(data['mag_z'])
    err_transits_r = list(data['err_r'])
    err_transits_z = list(data['err_z'])
    
#     times=[]
#     for time_list in time_transits_r:
#         times.append(int(time_list[0]))
    #initial_transit_time = t0_all[0]
#     n_max = int((t0_all[-1]-t0_all[0])/best_period)
#     if n_max%2!=0:
#         n_max+=1
#     even_ts=np.arange(t0_all[0], t0_all[0]+n_max*best_period, 2*period)
#     odd_ts=np.arange(t0_all[0]+period, t0_all[0]+(n_max+1)
    even_ts=[]
    t=t0_all[0]
    n=2
    while t<=time_transits_r[-1][-1]:
        even_ts.append(t)
        t=t0_all[0]+n*period
        n+=2
    t=t0_all[0]
    n=2
    while t>=0:
        t=t0_all[0]-n*period
        even_ts.append(t)
        n+=2
        
    odd_ts=[]
    t=t0_all[0]
    n=1
    while t<=time_transits_r[-1][-1]:
        t=t0_all[0]+n*period
        odd_ts.append(t)
        n+=2
    t=t0_all[0]
    n=1
    while t>=0:
        t=t0_all[0]-n*period
        odd_ts.append(t)
        n+=2
    
    even_ts=np.array(even_ts)
    odd_ts=np.array(odd_ts)
    
    time_r_total = np.array([item for sublist in time_transits_r for item in sublist])
    time_z_total = np.array([item for sublist in time_transits_z for item in sublist])
    mag_r_total = np.array([item for sublist in mag_transits_r for item in sublist])
    mag_z_total = np.array([item for sublist in mag_transits_z for item in sublist])
    err_r_total = np.array([item for sublist in err_transits_r for item in sublist])
    err_z_total = np.array([item for sublist in err_transits_z for item in sublist])
    
    data_all_r=pd.DataFrame({'time_r':time_r_total, 'mag_r':mag_r_total, 'err_r':err_r_total})
    data_all_z=pd.DataFrame({'time_z':time_z_total, 'mag_z':mag_z_total, 'err_z':err_z_total})
    data_even_r=pd.DataFrame({'time_r':[], 'mag_r':[], 'err_r':[]})
    data_odd_r=pd.DataFrame({'time_r':[], 'mag_r':[], 'err_r':[]})
    data_even_z=pd.DataFrame({'time_z':[], 'mag_z':[], 'err_z':[]})
    data_odd_z=pd.DataFrame({'time_z':[], 'mag_z':[], 'err_z':[]})
    idx_even_r=0
    idx_odd_r=0
    idx_even_z=0
    idx_odd_z=0
    for idx in range(len(data_all_r)):
        time=data_all_r.iloc[idx]['time_r']
        even_diffs=abs(even_ts-time)
        even_diff_min=even_diffs.min()
        odd_diffs=abs(odd_ts-time)
        odd_diff_min=odd_diffs.min()
        if even_diff_min < odd_diff_min:
            #data_even_r=data_even_r.append(data_all_r.iloc[idx])
            data_even_r.loc[idx_even_r] = data_all_r.iloc[idx]
            idx_even_r+=1
        else:
            #data_odd_r=data_odd_r.append(data_all_r.iloc[idx])
            data_odd_r.loc[idx_odd_r] = data_all_r.iloc[idx]
            idx_odd_r+=1
    for idx in range(len(data_all_z)):
        time=data_all_z.iloc[idx]['time_z']
        even_diffs=abs(even_ts-time)
        even_diff_min=even_diffs.min()
        odd_diffs=abs(odd_ts-time)
        odd_diff_min=odd_diffs.min()
        if even_diff_min < odd_diff_min:
            #data_even_z=data_even_z.append(data_all_z.iloc[idx])
            data_even_z.loc[idx_even_z] = data_all_z.iloc[idx]
            idx_even_z+=1
        else:
            #data_odd_z=data_odd_z.append(data_all_z.iloc[idx])
            data_odd_z.loc[idx_odd_z] = data_all_z.iloc[idx]
            idx_odd_z+=1
    return data_even_r, data_even_z, data_odd_r, data_odd_z
#     for n, night in enumerate(nights):
#         #dr = lcr_cut[(lcr_cut.iloc[:,tcol]>night) & (lcr_cut.iloc[:,tcol]<night+1.0)]
#         #dz = lcz_cut[(lcz_cut.iloc[:,tcol]>night) & (lcz_cut.iloc[:,tcol]<night+1.0)]
#         dr_cut = lcr_cut[(lcr_cut.iloc[:,tcol]>night) & (lcr_cut.iloc[:,tcol]<night+1.0)]
#         dz_cut = lcz_cut[(lcz_cut.iloc[:,tcol]>night) & (lcz_cut.iloc[:,tcol]<night+1.0)]
#     time_transits_even_r = []
#     mag_transits_even_r = []
#     err_transits_even_r = []
#     time_transits_odd_r = []
#     mag_transits_odd_r = []
#     err_transits_odd_r = []

#     time_transits_even_z = []
#     mag_transits_even_z = []
#     err_transits_even_z = []
#     time_transits_odd_z = []
#     mag_transits_odd_z = []
#     err_transits_odd_z = []
    
#     indices=np.arange(0, len(div))
#     for idx in range(len(time_transits_r)):
#         if div[idx]%2==0:
#             time_transits_even_r.extend(time_transits_r[idx])
#             mag_transits_even_r.extend(mag_transits_r[idx])
#             err_transits_even_r.extend(err_transits_r[idx])
#             time_transits_even_z.extend(time_transits_z[idx])
#             mag_transits_even_z.extend(mag_transits_z[idx])
#             err_transits_even_z.extend(err_transits_z[idx])
#         else:
#             time_transits_odd_r.extend(time_transits_r[idx])
#             mag_transits_odd_r.extend(mag_transits_r[idx])
#             err_transits_odd_r.extend(err_transits_r[idx])
#             time_transits_odd_z.extend(time_transits_z[idx])
#             mag_transits_odd_z.extend(mag_transits_z[idx])
#             err_transits_odd_z.extend(err_transits_z[idx])
            
#     data_even_z = pd.DataFrame({'time_z':time_transits_even_z, 'mag_z': mag_transits_even_z, 'err_z': err_transits_even_z})
#     data_even_r = pd.DataFrame({'time_r':time_transits_even_r, 'mag_r': mag_transits_even_r, 'err_r': err_transits_even_r})
#     data_odd_r = pd.DataFrame({'time_r':time_transits_odd_r, 'mag_r': mag_transits_odd_r, 'err_r': err_transits_odd_r})
#     data_odd_z = pd.DataFrame({'time_z':time_transits_odd_z, 'mag_z': mag_transits_odd_z, 'err_z': err_transits_odd_z})
#     return data_even_r, data_even_z, data_odd_r, data_odd_z

def compare_odd_and_even(data_new, t0s, p0_guess, best_period, ecc, w, rldc_r, rldc_z, all_nights, show_plot=False, ndim=7, nwalkers=20, nsteps=2000, nburn=500, nthin=125, exclude=None, even_guesses=None, odd_guesses=None, mcmc_sigmas=[10e-5, 10e-3, 10e-3, 10e-3, 10e-3, 10e-5, 10e-5]):

    no_even=False
    no_odd=False
    
    data_even_r, data_even_z, data_odd_r, data_odd_z = odd_even_data_new(data_new, t0s, best_period, all_nights)

    t0_1 = t0s[0]
    
    #data_even_r, data_even_z = folded_lightcurve(data_even_r, data_even_z, best_period, t0_1, show_plot)
    #data_odd_r, data_odd_z = folded_lightcurve(data_odd_r, data_odd_z, best_period, t0_1, show_plot)
    
    p0_guess[0]=t0_1

    np.random.seed(42) # set seed so results consistently converge
    labels = ["t0", "rp_r", "rp_z", "a", "inc", "C_r", "C_z"]

    if len(data_even_r)!=0:
        #param_even, cov_err_even, err_tot_adj_even, results_batman_even = run_batman(p0_guess, best_period, ecc, w, rldc_r, rldc_z, np.array(data_even_r['time_r']), np.array(data_even_z['time_z']), np.array(data_even_r['mag_r']), np.array(data_even_z['mag_z']), np.array(data_even_r['err_r']), np.array(data_even_z['err_z']), title='Transit Model for Even Data', folded=True, t0_1=t0_1)
        
        time_r_even = np.array(data_even_r['time_r']) - t0_1 + 0.25*best_period
        time_z_even = np.array(data_even_z['time_z']) - t0_1 + 0.25*best_period
        phases_r= foldAt(time_r_even, best_period, 0.0)
        phases_z= foldAt(time_z_even, best_period, 0.0)
        
        sortIndi = np.argsort(phases_r)
        phases_r = phases_r[sortIndi]
        mag_r = np.array(data_even_r['mag_r'])[sortIndi]
        err_r = np.array(data_even_r['err_r'])[sortIndi]
        time_r = np.array(data_even_r['time_r'])[sortIndi]

        sortIndi = np.argsort(phases_z)
        phases_z = phases_z[sortIndi]
        mag_z = np.array(data_even_z['mag_z'])[sortIndi]
        err_z = np.array(data_even_z['err_z'])[sortIndi]
        time_z = np.array(data_even_z['time_z'])[sortIndi]

        data_even_r['phase_r']=phases_r
        data_even_r['time_r']=time_r
        data_even_r['mag_r']=mag_r
        data_even_r['err_r']=err_r
        data_even_z['phase_z']=phases_z
        data_even_z['time_z']=time_z
        data_even_z['mag_z']=mag_z
        data_even_z['err_z']=err_z
    
        data_even_r_cut=data_even_r[(data_even_r['phase_r']>0.05) & (data_even_r['phase_r']<0.45)]
        data_even_z_cut=data_even_z[(data_even_z['phase_z']>0.05) & (data_even_z['phase_z']<0.45)]
        
        param_even, cov_err_even, err_tot_adj_even, results_batman_even = run_batman(p0_guess, best_period, ecc, w, rldc_r, rldc_z, data_even_r_cut, data_even_z_cut, data_even_r, data_even_z, title='Transit Model for Even Data', t0_1=t0_1)

        lower_bounds = [param_even[0]-0.1, 0.0, 0.0, 0.0, 70.0, 0., 0.]
        upper_bounds = [param_even[0]+0.1, 10.0, 10.0, 40.0, 130.0, 40., 40.]

        time_tot = np.append(np.array(data_even_r_cut['time_r']), np.array(data_even_z_cut['time_z'])+100000)
        mag_tot = np.append(np.array(data_even_r_cut['mag_r']), np.array(data_even_z_cut['mag_z']))
        err_tot = np.append(np.array(data_even_r_cut['err_r']), np.array(data_even_z_cut['err_z']))
        
        if even_guesses is None:
            even_guesses_new = param_even
        else:
            even_guesses_new = even_guesses

        variances = cov_err_even 
        sigmas = np.sqrt(variances)
        results_mcmc_even, results_mcmc_per_even, mcmc_chains_even = mcmc_functions.mcmc_all(ndim, nwalkers, nsteps, nburn, nthin, even_guesses_new, labels, time_tot, mag_tot, err_tot_adj_even, lower_bounds, upper_bounds, rldc_r, rldc_z, best_period, ecc, w, title='Transit Model for Even Data', rp_corner=True, sigmas=mcmc_sigmas)
    else:
        print('No even transits.')
        no_even=True
        
    if len(data_odd_r)!=0:
        #param_odd, cov_err_odd, err_tot_adj_odd, results_batman_odd = run_batman(p0_guess, best_period, ecc, w, rldc_r, rldc_z, np.array(data_odd_r['time_r']), np.array(data_odd_z['time_z']), np.array(data_odd_r['mag_r']), np.array(data_odd_z['mag_z']), np.array(data_odd_r['err_r']), np.array(data_odd_z['err_z']), title='Transit Model for Odd Data' , folded=True, t0_1=t0_1)
        
        time_r_odd = np.array(data_odd_r['time_r']) - t0_1 + 0.25*best_period
        time_z_odd = np.array(data_odd_z['time_z']) - t0_1 + 0.25*best_period
        phases_r= foldAt(time_r_odd, best_period, 0.0)
        phases_z= foldAt(time_z_odd, best_period, 0.0)
        
        sortIndi = np.argsort(phases_r)
        phases_r = phases_r[sortIndi]
        mag_r = np.array(data_odd_r['mag_r'])[sortIndi]
        err_r = np.array(data_odd_r['err_r'])[sortIndi]
        time_r = np.array(data_odd_r['time_r'])[sortIndi]

        sortIndi = np.argsort(phases_z)
        phases_z = phases_z[sortIndi]
        mag_z = np.array(data_odd_z['mag_z'])[sortIndi]
        err_z = np.array(data_odd_z['err_z'])[sortIndi]
        time_z = np.array(data_odd_z['time_z'])[sortIndi]

        data_odd_r['phase_r']=phases_r
        data_odd_r['time_r']=time_r
        data_odd_r['mag_r']=mag_r
        data_odd_r['err_r']=err_r
        data_odd_z['phase_z']=phases_z
        data_odd_z['time_z']=time_z
        data_odd_z['mag_z']=mag_z
        data_odd_z['err_z']=err_z
        
        data_odd_r_cut=data_odd_r[(data_odd_r['phase_r']>0.05) & (data_odd_r['phase_r']<0.45)]
        data_odd_z_cut=data_odd_z[(data_odd_z['phase_z']>0.05) & (data_odd_z['phase_z']<0.45)]
        
        param_odd, cov_err_odd, err_tot_adj_odd, results_batman_odd = run_batman(p0_guess, best_period, ecc, w, rldc_r, rldc_z, data_odd_r_cut, data_odd_z_cut, data_odd_r, data_odd_z, title='Transit Model for Odd Data' , t0_1=t0_1)

        lower_bounds = [param_odd[0]-0.1, 0.0, 0.0, 0.0, 70.0, 0., 0.]
        upper_bounds = [param_odd[0]+0.1, 10.0, 10.0, 40.0, 130.0, 40., 40.]
        
        time_tot = np.append(np.array(data_odd_r_cut['time_r']), np.array(data_odd_z_cut['time_z'])+100000)
        mag_tot = np.append(np.array(data_odd_r_cut['mag_r']), np.array(data_odd_z_cut['mag_z']))
        err_tot = np.append(np.array(data_odd_r_cut['err_r']), np.array(data_odd_z_cut['err_z']))
        
        if odd_guesses is None:
            odd_guesses_new = param_odd
        else:
            odd_guesses_new = odd_guesses
            
        variances = cov_err_odd 
        sigmas = np.sqrt(variances)
        results_mcmc_odd, results_mcmc_per_odd, mcmc_chains_odd = mcmc_functions.mcmc_all(ndim, nwalkers, nsteps, nburn, nthin, odd_guesses_new, labels, time_tot, mag_tot, err_tot_adj_odd, lower_bounds, upper_bounds, rldc_r, rldc_z, best_period, ecc, w, title='Transit model for Odd Data', rp_corner=True, sigmas=mcmc_sigmas)
    else:
        print('No odd transits.')
        no_odd=True
        
    if no_odd==True and no_even==False:
        return results_batman_even, None, results_mcmc_even, None, results_mcmc_per_even, None, even_guesses_new, None, mcmc_chains_even, None
    elif no_even==True and no_odd==False:
        return None, results_batman_odd, None, results_mcmc_odd, None, results_mcmc_per_odd, None, odd_guesses_new, None, mcmc_chains_odd
    elif no_even==True and no_odd==True:
        return None, None, None, None, None, None, None, None, None, None
    
    diff, diff_err = difference(results_batman_even['rp_r'], results_batman_even['rp_r_err'], results_batman_odd['rp_r'], results_batman_odd['rp_r_err'])
    #print('Difference between odd/even from batman for r-band: ' + str(round(diff, 5))+ ' Err: '+str(round(diff_err, 5)))
    diff, diff_err = difference(results_batman_even['rp_z'], results_batman_even['rp_z_err'], results_batman_odd['rp_z'], results_batman_odd['rp_z_err'])
    #print('Difference between odd/even from batman for z-band: ' + str(round(diff, 5))+ ' Err: '+str(round(diff_err, 5)))
    
    diff, diff_err = difference(results_mcmc_even['rp_r'], results_mcmc_even['rp_r_err'], results_mcmc_odd['rp_r'], results_mcmc_odd['rp_r_err'])
    #print('Difference between odd/even from mcmc for r-band: ' + str(round(diff, 5))+ ' Err: '+str(round(diff_err, 5)))
    diff, diff_err = difference(results_mcmc_even['rp_z'], results_mcmc_even['rp_z_err'], results_mcmc_odd['rp_z'], results_mcmc_odd['rp_z_err'])
    #print('Difference between odd/even from mcmc for z-band: ' + str(round(diff, 5))+ ' Err: '+str(round(diff_err, 5)))
    #plot as points with error bars, x axis as scenario, odd r band etc. y axis is values with error bars. add text to plot
    
    return results_batman_even, results_batman_odd, results_mcmc_even, results_mcmc_odd, results_mcmc_per_even, results_mcmc_per_odd, even_guesses_new, odd_guesses_new, mcmc_chains_even, mcmc_chains_odd

