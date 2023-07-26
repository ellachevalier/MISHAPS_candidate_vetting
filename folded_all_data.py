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

from other_funcs import weighted_avg_and_std, difference, folded_lightcurve
#from run_batman_inc_combined import run_batman
import mcmc_functions_inc_combined as mcmc_functions
from run_batman_folded import run_batman

def folded_all_data(data_new, t0s, p0_guess, best_period, ecc, w, rldc_r, rldc_z, show_plot=False, ndim=7, nwalkers=20, nsteps=2000, nburn=500, nthin=125, all_guesses=None, mcmc_sigmas=[10e-5, 10e-3, 10e-3, 10e-3, 10e-3, 10e-5, 10e-5]):
    
    time_transits_r = data_new['time_r']
    time_transits_z = data_new['time_z']
    mag_transits_r = data_new['mag_r']
    mag_transits_z = data_new['mag_z']
    err_transits_r = data_new['err_r']
    err_transits_z = data_new['err_z']
    
    time_transits_tot_r = []
    mag_transits_tot_r = []
    err_transits_tot_r = []
    time_transits_tot_z = []
    mag_transits_tot_z = []
    err_transits_tot_z = []

    t0_1 = t0s[0]
    for idx in range(len(time_transits_r)):
        time_transits_tot_r.extend(time_transits_r[idx])
        mag_transits_tot_r.extend(mag_transits_r[idx])
        err_transits_tot_r.extend(err_transits_r[idx])
        time_transits_tot_z.extend(time_transits_z[idx])
        mag_transits_tot_z.extend(mag_transits_z[idx])
        err_transits_tot_z.extend(err_transits_z[idx])

    time_transits_tot_r = np.array(time_transits_tot_r)
    mag_transits_tot_r = np.array(mag_transits_tot_r)
    err_r = np.array(err_transits_tot_r)
    time_transits_tot_z = np.array(time_transits_tot_z)
    mag_transits_tot_z = np.array(mag_transits_tot_z)
    err_z = np.array(err_transits_tot_z)

    data_r_all = pd.DataFrame({'time_r': time_transits_tot_r, 'mag_r': mag_transits_tot_r, 'err_r':err_r})
    data_z_all = pd.DataFrame({'time_z': time_transits_tot_z, 'mag_z': mag_transits_tot_z, 'err_z':err_z})
    #data_r_folded, data_z_folded = folded_lightcurve(data_r_all, data_z_all, best_period, t0_1, show_plot)
    
    #p0_guess[0]=0.25
    p0_guess[0]=t0_1
    #param_all, cov_err_all, err_tot_adj_all, results_batman_all = run_batman(p0_guess, best_period, ecc, w, rldc_r, rldc_z, np.array(data_r_all['time_r']), np.array(data_z_all['time_z']), np.array(data_r_all['mag_r']), np.array(data_z_all['mag_z']), np.array(data_r_all['err_r']), np.array(data_z_all['err_z']), title='Transit Model for All Data', folded=True, t0_1=t0_1)

    time_r = np.array(data_r_all['time_r']) - t0_1 + 0.25*best_period
    time_z = np.array(data_z_all['time_z']) - t0_1 + 0.25*best_period
    phases_r= foldAt(time_r, best_period, 0.0)
    phases_z= foldAt(time_z, best_period, 0.0)
    
    sortIndi = np.argsort(phases_r)
    phases_r = phases_r[sortIndi]
    mag_r = np.array(data_r_all['mag_r'])[sortIndi]
    err_r = np.array(data_r_all['err_r'])[sortIndi]
    time_r = np.array(data_r_all['time_r'])[sortIndi]
    
    sortIndi = np.argsort(phases_z)
    phases_z = phases_z[sortIndi]
    mag_z = np.array(data_z_all['mag_z'])[sortIndi]
    err_z = np.array(data_z_all['err_z'])[sortIndi]
    time_z = np.array(data_z_all['time_z'])[sortIndi]
    
    data_r_all['phase_r']=phases_r
    data_r_all['time_r']=time_r
    data_r_all['mag_r']=mag_r
    data_r_all['err_r']=err_r
    data_z_all['phase_z']=phases_z
    data_z_all['time_z']=time_z
    data_z_all['mag_z']=mag_z
    data_z_all['err_z']=err_z
    
    data_r_cut=data_r_all[(data_r_all['phase_r']>0.05) & (data_r_all['phase_r']<0.45)]
    data_z_cut=data_z_all[(data_z_all['phase_z']>0.05) & (data_z_all['phase_z']<0.45)]
        
    param_all, cov_err_all, err_tot_adj_all, results_batman_all = run_batman(p0_guess, best_period, ecc, w, rldc_r, rldc_z, data_r_cut, data_z_cut, data_r_all, data_z_all, title='Transit Model for All Data', t0_1=t0_1)

    time_tot = np.append(np.array(data_r_cut['time_r']), np.array(data_z_cut['time_z'])+100000)
    mag_tot = np.append(np.array(data_r_cut['mag_r']), np.array(data_z_cut['mag_z']))
    err_tot = np.append(np.array(data_r_cut['err_r']), np.array(data_z_cut['err_z']))
    
    labels = ["t0", "rp_r", "rp_z", "a", "inc", "C_r", "C_z"]

    lower_bounds = [param_all[0]-0.1, 0.0, 0.0, 0.0, 50.0, 0., 0.]
    upper_bounds = [param_all[0]+0.1, 10.0, 10.0, 40.0, 130.0, 40., 40.]
        
    if all_guesses is None:
        all_guesses_new = param_all
    else:
        all_guesses_new = all_guesses
        
    variances = cov_err_all 
    sigmas = np.sqrt(variances)
    results_mcmc_all, results_mcmc_per, mcmc_chains_all = mcmc_functions.mcmc_all(ndim, nwalkers, nsteps, nburn, nthin, all_guesses_new, labels, time_tot, mag_tot, err_tot_adj_all, lower_bounds, upper_bounds, rldc_r, rldc_z, best_period, ecc, w, title='Transit Model for All Data', rp_corner=True, sigmas=mcmc_sigmas)
    
    diff, diff_err = difference(results_batman_all['rp_r'], results_batman_all['rp_r_err'], results_batman_all['rp_z'], results_batman_all['rp_z_err'])
    #print('Difference between r- and z-bands from batman: ' + str(round(diff, 7))+ ' Err: '+str(round(diff_err, 7)))
    diff, diff_err = difference(results_mcmc_all['rp_r'], results_mcmc_all['rp_r_err'], results_mcmc_all['rp_z'], results_mcmc_all['rp_z_err'])
    #print('Difference between r- and z-bands from mcmc: ' + str(round(diff, 7))+ ' Err: '+str(round(diff_err, 7)))
    
    return results_batman_all, results_mcmc_all, results_mcmc_per, all_guesses_new, mcmc_chains_all

