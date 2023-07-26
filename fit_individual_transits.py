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
import os

from other_funcs import weighted_avg_and_std, difference
from run_batman_inc_combined import run_batman
import mcmc_functions_inc_combined as mcmc_functions

def fit_individual_transits(data_new, transit_nights, p0_guess, period, ecc, w, rldc_r, rldc_z, ndim=7, nwalkers=20, nsteps=2000, nburn=500, nthin=125):
    t0_all=[]
    rp_r_all=[]
    rp_z_all=[]

    t0_err_all=[]
    rp_r_err_all=[]
    rp_z_err_all=[]

    t0_mcmc=[]
    t0_err_mcmc=[]
    rp_r_mcmc=[]
    rp_r_err_mcmc=[]
    rp_z_mcmc=[]
    rp_z_err_mcmc=[]
    
    a_mcmc=[]
    inc_mcmc=[]
    C_r_mcmc=[]
    C_z_mcmc=[]

    time_transits_r = data_new['time_r']
    time_transits_z = data_new['time_z']
    mag_transits_r = data_new['mag_r']
    mag_transits_z = data_new['mag_z']
    err_transits_r = data_new['err_r']
    err_transits_z = data_new['err_z']
    
    t0_guesses = p0_guess[0]
    
    #indices=np.arange(0, len(time_transits_r))
    #if exclude != None:
        #indices=[idx for idx in indices if idx not in exclude]
    ans_r_list=[]
    ans_z_list=[]
    for idx in range(len(time_transits_r)):

        p0_guess[0]=t0_guesses[idx]
        param, cov_err, err_tot_adj, results_batman, ans_r, ans_z = run_batman(p0_guess, period, ecc, w, rldc_r, rldc_z, time_transits_r[idx], time_transits_z[idx], mag_transits_r[idx], mag_transits_z[idx], err_transits_r[idx], err_transits_z[idx], title='Individual Transit Model for Night '+ str(transit_nights[idx]), flag=True)
        
        ans_r_list.append(ans_r)
        ans_z_list.append(ans_z)
        
        time_tot = np.append(time_transits_r[idx], time_transits_z[idx]+100000)
        mag_tot = np.append(mag_transits_r[idx], mag_transits_z[idx])
        err_tot = np.append(err_transits_r[idx], err_transits_z[idx])

        guesses = param
        variances = cov_err 
        sigmas = np.sqrt(variances)

        np.random.seed(42) # set seed so results consistently converge
        labels = ["t0", "rp_r", "rp_z", "a", "inc", "C_r", "C_z"]

        lower_bounds = [param[0]-0.5, 0.0, 0.0, 0.0, 70.0, 0., 0.]
        upper_bounds = [param[0]+0.5, 10.0, 10.0, 40.0, 130.0, 40., 40.]

        results_mcmc, results_mcmc_per, mcmc_chains = mcmc_functions.mcmc_all(ndim, nwalkers, nsteps, nburn, nthin, guesses, labels, time_tot, mag_tot, err_tot_adj, lower_bounds, upper_bounds, rldc_r, rldc_z, period, ecc, w, title='Individual Transit Model for Night '+ str(transit_nights[idx]))

        t0_all.append(results_batman['t0'])
        t0_err_all.append(results_batman['t0_err'])
        rp_r_all.append(results_batman['rp_r'])
        rp_r_err_all.append(results_batman['rp_r_err'])
        rp_z_all.append(results_batman['rp_z'])
        rp_z_err_all.append(results_batman['rp_z_err'])

        t0_mcmc.append(results_mcmc['t0'])
        t0_err_mcmc.append(results_mcmc['t0_err'])
        rp_r_mcmc.append(results_mcmc['rp_r'])
        rp_r_err_mcmc.append(results_mcmc['rp_r_err'])
        rp_z_mcmc.append(results_mcmc['rp_z'])
        rp_z_err_mcmc.append(results_mcmc['rp_z_err'])
        
        a_mcmc.append(results_mcmc['a'])
        inc_mcmc.append(results_mcmc['inc'])
        C_r_mcmc.append(results_mcmc['C_r'])
        C_z_mcmc.append(results_mcmc['C_z'])
    
    rp_r_avg_batman = weighted_avg_and_std(rp_r_all, rp_r_err_all)
    rp_z_avg_batman = weighted_avg_and_std(rp_z_all, rp_z_err_all)
    diff, diff_err = difference(rp_r_avg_batman[0], rp_r_avg_batman[1], rp_z_avg_batman[0], rp_z_avg_batman[1])
    #print('Difference between r and z band from batman avg: '+str(round(diff,5))+' Error: '+str(round(diff_err,5)))

    rp_r_avg_mcmc = weighted_avg_and_std(rp_r_mcmc, rp_r_err_mcmc)
    rp_z_avg_mcmc = weighted_avg_and_std(rp_z_mcmc, rp_z_err_mcmc)
    diff, diff_err = difference(rp_r_avg_mcmc[0], rp_r_avg_mcmc[1], rp_z_avg_mcmc[0], rp_z_avg_mcmc[1])
    #print('Difference between r and z band from mcmc avg: '+str(round(diff,5))+' Error: '+str(round(diff_err,5)))
    
    avg_first_and_last={'rp_r':(rp_r_mcmc[0]+rp_r_mcmc[-1])/2, 'rp_z':(rp_z_mcmc[0]+rp_z_mcmc[-1])/2, 'a':(a_mcmc[0]+a_mcmc[-1])/2, 'inc':(inc_mcmc[0]+inc_mcmc[-1])/2, 'C_r':(C_r_mcmc[0]+C_r_mcmc[-1])/2, 'C_z':(C_z_mcmc[0]+C_z_mcmc[-1])/2}
    
    
    fig = plt.figure(figsize=(2*len(time_transits_r),26))
    for idx in range(len(time_transits_r)):
    #for n, night in enumerate(nights):
            #dr = lcr[(lcr.iloc[:,tcol]>night) & (lcr.iloc[:,tcol]<night+1.0)]
            #dz = lcz[(lcz.iloc[:,tcol]>night) & (lcz.iloc[:,tcol]<night+1.0)]
            #dr_cut = lcr_cut[(lcr_cut.iloc[:,tcol]>night) & (lcr_cut.iloc[:,tcol]<night+100.0)]
            #dz_cut = lcz_cut[(lcz_cut.iloc[:,tcol]>night) & (lcz_cut.iloc[:,tcol]<night+1.0)]

            ax = plt.subplot(5, int(np.ceil(len(time_transits_r)/5.0)), int(idx+1))

            ax.errorbar(time_transits_r[idx], mag_transits_r[idx], yerr=err_transits_r[idx],
                        fmt='o', color='lightskyblue', ms=1, label='r-band data')
            ax.plot(time_transits_r[idx], ans_r_list[idx], '-', color='blue')
            ax.errorbar(time_transits_z[idx], mag_transits_z[idx], yerr=err_transits_z[idx],
                        fmt='o', color='pink', ms=1, label='z-band data')
            ax.plot(time_transits_z[idx], ans_z_list[idx], '-', color='red')

            ax.set_title('Night ' +str(int(time_transits_r[idx][0])))
            ax.invert_yaxis()
            #ax.legend()
    labels=['r-band', 'z-band']
    fig.legend(labels=labels, loc='lower center', fontsize=20) 
    fig.suptitle('Individual transits with models', fontsize=45.0)
    #plt.savefig(os.path.join('poster_images','individual_transits_plot.png'), dpi=300)
    plt.savefig(os.path.join('figs','individual_transits_plot.pdf'))
    plt.close()
    
    return results_batman, results_mcmc, results_mcmc_per, t0_mcmc, t0_err_mcmc, avg_first_and_last

