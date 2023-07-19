#!/usr/bin/env python
# coding: utf-8
# %%
#separate r and z band flux parameters, fix time and vary baseline
#double check that primary and secondary cant be fit at the same time

# %%


#results_batman_all, results_mcmc_all, results_mcmc_per_all, all_guesses, mcmc_chains_all = folded_all_data(data_new, t0_mcmc, p0_guess, best_period, ecc, w, rldc_r, rldc_z, show_plot=True, nwalkers=nwalkers, nsteps=nsteps)
import batman
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from other_funcs import *
def fit_secondary_eclipse(data_new, results_mcmc_per_all, t0_mcmc, p0_guess, period, rldc, ecc=0.0, w=0.0):
    rldc_r=rldc[0]
    rldc_z=rldc[1]
    
    t0=results_mcmc_per_all['50th'][0]
    rp_r=results_mcmc_per_all['50th'][1]
    rp_z=results_mcmc_per_all['50th'][2]
    a=results_mcmc_per_all['50th'][3]
    inc=results_mcmc_per_all['50th'][4]
#     C_r=results_mcmc_per_all['50th'][5]
#     C_z=results_mcmc_per_all['50th'][6]

    t_sec=t0+0.5*period
    
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
    

    t0_1 = t0_mcmc[0]
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

    #def test(time_tot, t0, rp_r, rp_z, a, inc, C_r, C_z, fp, t_sec):
    def test(time_tot, fp_r, fp_z, C_r, C_z):
        time_z_try = time_tot[time_tot>100000]
        time_z_try = time_z_try - 100000
        time_r_try = time_tot[time_tot<100000]

        params_z = batman.TransitParams()
        params_z.t0 =  t0                  #time of inferior conjunction
        params_z.per = period                 #orbital period
        params_z.rp = rp_z                      #planet radius (in units of stellar radii)
        params_z.a =  a                       #semi-major axis (in units of stellar radii)
        params_z.inc = inc                    #orbital inclination (in degrees) 
        params_z.ecc = ecc                     #eccentricity
        params_z.w =  w                        #longitude of periastron (in degrees) 
        params_z.u =  rldc_z                #limb darkening coefficients [u1, u2] 
        params_z.fp = fp_z
        params_z.t_secondary = t_sec
        params_z.limb_dark = "quadratic"

        m_z=batman.TransitModel(params_z, time_z_try)
        flux_z = m_z.light_curve(params_z)
        mag_z_1 = C_z-2.5*np.log10(flux_z)

        params_r = batman.TransitParams()
        params_r.t0 =  t0                     #time of inferior conjunction
        params_r.per = period                 #orbital period
        params_r.rp = rp_r                    #planet radius (in units of stellar radii)
        params_r.a =  a                       #semi-major axis (in units of stellar radii)
        params_r.inc = inc                   #orbital inclination (in degrees) 
        params_r.ecc = ecc                    #eccentricity
        params_r.w =  w                       #longitude of periastron (in degrees) 
        params_r.u =  rldc_r
        params_r.fp = fp_r
        params_r.t_secondary = t_sec  
        params_r.limb_dark = "quadratic"

        m_r=batman.TransitModel(params_r, time_r_try)
        flux_r = m_r.light_curve(params_r)
        mag_r_2 = C_r-2.5*np.log10(flux_r)
        mag = np.append(mag_r_2, mag_z_1)

        return mag

        #if os.path.exists('figs')==False:
            #os.mkdir('figs')

    #p0=[p0_guess[0], p0_guess[1], p0_guess[2], p0_guess[3], p0_guess[4], p0_guess[5], p0_guess[6], 0.001, p0_guess[0]+0.5*period]
    p0=[p0_guess[0], p0_guess[1], p0_guess[2], p0_guess[3]]
    time_tot = np.append(time_transits_tot_r, time_transits_tot_z+100000)
    mag_tot = np.append(mag_transits_tot_r, mag_transits_tot_z)
    err_tot = np.append(err_transits_tot_r, err_transits_tot_z)

    #lower_bounds = np.array([-np.inf, 0.0, 0.0, 0.0, 70.0, 0., 0., -np.inf, -np.inf])
    #upper_bounds = np.array([np.inf, 10.0, 10.0, 40.0, 130.0, 40., 40., np.inf, np.inf])
    lower_bounds = np.array([-np.inf, -np.inf])
    upper_bounds = np.array([np.inf, np.inf])

    #param, param_cov = curve_fit(test, time_tot, mag_tot, p0=p0, sigma=err_tot, bounds=(lower_bounds, upper_bounds))
    param, param_cov = curve_fit(test, time_tot, mag_tot, p0=p0, sigma=err_tot)

    # t0 = param[0]
    # rp_r = param[1]
    # rp_z = param[2]
    # a = param[3]
    # inc = param[4]
    # C_r = param[5]
    # C_z = param[6]
    # fp = param[7]
    # t_sec = param[8]
    fp_r=param[0]
    fp_z=param[1]
    C_r=param[2]
    C_z=param[3]
    #t_sec=param[1]
    cov_err = np.sqrt(np.diag(param_cov))
    fp_r_err=cov_err[0]
    fp_z_err=cov_err[1]
    C_r_err=cov_err[2]
    C_z_err=cov_err[3]

    params1 = batman.TransitParams()
    params1.t0 = t0
    params1.per = period
    params1.rp = rp_r
    params1.a = a
    params1.inc = inc
    params1.ecc = ecc
    params1.w = w
    params1.u = rldc_r
    params1.fp = fp_r
    params1.t_secondary = t_sec 
    params1.limb_dark = "quadratic"  
    m1=batman.TransitModel(params1, time_transits_tot_r, transittype='secondary')
    flux1=m1.light_curve(params1)
    ans_r=C_r-2.5*np.log10(flux1)

    params2 = batman.TransitParams()
    params2.t0 = t0
    params2.per = period
    params2.rp = rp_z
    params2.a = a
    params2.inc = inc
    params2.ecc = ecc
    params2.w = w
    params2.u = rldc_z
    params2.fp = fp_z
    params2.t_secondary = t_sec 
    params2.limb_dark = "quadratic"  
    m2=batman.TransitModel(params2, time_transits_tot_z, transittype='secondary')
    flux2=m2.light_curve(params2)
    ans_z=C_z-2.5*np.log10(flux2)

    # if flag==True:
    #     const_r=np.ones(len(ans_r))*C_r
    #     const_z=np.ones(len(ans_z))*C_z
    #     if np.mean(ans_r)==C_r or np.mean(ans_z)==C_z:
    #         print('Fit failed for '+title)    
    
    data_r_all = pd.DataFrame({'time_r': time_transits_tot_r, 'mag_r': mag_transits_tot_r, 'err_r':err_transits_tot_r})
    data_z_all = pd.DataFrame({'time_z': time_transits_tot_z, 'mag_z': mag_transits_tot_z, 'err_z':err_transits_tot_z})
    initial_guess_r=np.array(test(time_tot, *p0)[:len(time_transits_tot_r)])
    initial_guess_z=np.array(test(time_tot, *p0)[len(time_transits_tot_r):])
    data_r_folded, data_z_folded, ans_r_folded, ans_z_folded, initial_guess_r_folded, initial_guess_z_folded = folded_lightcurve(data_r_all, data_z_all, period, t0_1, show_plot=False, ans_r=ans_r, ans_z=ans_z, initial_guess_r=initial_guess_r, initial_guess_z=initial_guess_z)
    plt.figure()
    plt.plot(data_r_folded['phases_r'], data_r_folded['mag_r'],'o', color='lightskyblue', label='data r', markersize=1.)
    plt.plot(data_z_folded['phases_z'], data_z_folded['mag_z'],'o', color='pink', label='data z', markersize=1.)
    #plt.plot(data_r_folded['phases_r'], initial_guess_r_folded, '--', color='blue', label='initial guess r')
    #plt.plot(data_z_folded['phases_z'], initial_guess_z_folded, '--', color='red', label='initial guess z')
    plt.gca().invert_yaxis()
    plt.plot(data_r_folded['phases_r'], ans_r_folded, color='blue', label='optimized fit r') #deepskyblue
    plt.plot(data_z_folded['phases_z'], ans_z_folded, color='red', label='optimized fit z') #blueviolet
    # plt.xlim(0,1)
    plt.legend()
    #plt.ylim(C_z+0.05, C_r-0.05)
    #plt.title(title+' for P='+str(round(period,4)))
    plt.xlabel('Phase')
    plt.ylabel('Magnitude')
    text_str='fp_r='+str(fp_r)+'$\pm$'+str(fp_r_err)+'\n fp_z='+str(fp_z)+'$\pm$'+str(fp_z_err)+'\n C_r='+str(C_r)+'$\pm$'+str(C_r_err)+'\n C_z='+str(C_z)+'$\pm$'+str(C_z_err)
    plt.text(0.4, C_z-0.03, text_str, fontsize=10, bbox=dict(facecolor='lightgray', alpha=0.5))
    plt.title('Secondary eclipse fit for P='+str(round(period,4)))
    plt.savefig(os.path.join('figs', 'secondary_eclipse_fit.pdf'))
    #plt.show()
    
    #Refit for all parameters
#     def test(time_tot, t0, rp_r, rp_z, a, inc, C_r, C_z, fp, t_sec):
#         time_z_try = time_tot[time_tot>100000]
#         time_z_try = time_z_try - 100000
#         time_r_try = time_tot[time_tot<100000]

#         params_z = batman.TransitParams()
#         params_z.t0 =  t0                  #time of inferior conjunction
#         params_z.per = period                 #orbital period
#         params_z.rp = rp_z                      #planet radius (in units of stellar radii)
#         params_z.a =  a                       #semi-major axis (in units of stellar radii)
#         params_z.inc = inc                    #orbital inclination (in degrees) 
#         params_z.ecc = ecc                     #eccentricity
#         params_z.w =  w                        #longitude of periastron (in degrees) 
#         params_z.u =  rldc_z                #limb darkening coefficients [u1, u2] 
#         params_z.fp = fp
#         params_z.t_secondary = t_sec
#         params_z.limb_dark = "quadratic"

#         m_z=batman.TransitModel(params_z, time_z_try)
#         flux_z = m_z.light_curve(params_z)
#         mag_z_1 = C_z-2.5*np.log10(flux_z)

#         params_r = batman.TransitParams()
#         params_r.t0 =  t0                     #time of inferior conjunction
#         params_r.per = period                 #orbital period
#         params_r.rp = rp_r                    #planet radius (in units of stellar radii)
#         params_r.a =  a                       #semi-major axis (in units of stellar radii)
#         params_r.inc = inc                   #orbital inclination (in degrees) 
#         params_r.ecc = ecc                    #eccentricity
#         params_r.w =  w                       #longitude of periastron (in degrees) 
#         params_r.u =  rldc_r
#         params_r.fp = fp
#         params_r.t_secondary = t_sec  
#         params_r.limb_dark = "quadratic"

#         m_r=batman.TransitModel(params_r, time_r_try)
#         flux_r = m_r.light_curve(params_r)
#         mag_r_2 = C_r-2.5*np.log10(flux_r)
#         mag = np.append(mag_r_2, mag_z_1)

#         return mag

#         #if os.path.exists('figs')==False:
#             #os.mkdir('figs')

#     #p0=[p0_guess[0], p0_guess[1], p0_guess[2], p0_guess[3], p0_guess[4], p0_guess[5], p0_guess[6], 0.001, p0_guess[0]+0.5*period]
#     p0=[t0, rp_r, rp_z, a, inc, C_r, C_z, p0_guess[0], p0_guess[1]]
#     time_tot = np.append(time_transits_tot_r, time_transits_tot_z+100000)
#     mag_tot = np.append(mag_transits_tot_r, mag_transits_tot_z)
#     err_tot = np.append(err_transits_tot_r, err_transits_tot_z)

#     lower_bounds = np.array([-np.inf, 0.0, 0.0, 0.0, 70.0, 0., 0., -np.inf, -np.inf])
#     upper_bounds = np.array([np.inf, 10.0, 10.0, 40.0, 130.0, 40., 40., np.inf, np.inf])
#     #lower_bounds = np.array([-np.inf, -np.inf])
#     #upper_bounds = np.array([np.inf, np.inf])

#     #param, param_cov = curve_fit(test, time_tot, mag_tot, p0=p0, sigma=err_tot, bounds=(lower_bounds, upper_bounds))
#     param_all, param_cov_all = curve_fit(test, time_tot, mag_tot, p0=p0, sigma=err_tot)

#     t0 = param_all[0]
#     rp_r = param_all[1]
#     rp_z = param_all[2]
#     a = param_all[3]
#     inc = param_all[4]
#     C_r = param_all[5]
#     C_z = param_all[6]
#     fp=param_all[7]
#     t_sec=param_all[8]
#     cov_err_all = np.sqrt(np.diag(param_cov_all))
#     fp_err=cov_err_all[7]
#     t_sec_err=cov_err_all[8]

#     params1 = batman.TransitParams()
#     params1.t0 = t0
#     params1.per = period
#     params1.rp = rp_r
#     params1.a = a
#     params1.inc = inc
#     params1.ecc = ecc
#     params1.w = w
#     params1.u = rldc_r
#     params1.fp = fp
#     params1.t_secondary = t_sec 
#     params1.limb_dark = "quadratic"  
#     m1=batman.TransitModel(params1, time_transits_tot_r, transittype='secondary')
#     flux1=m1.light_curve(params1)
#     ans_r=C_r-2.5*np.log10(flux1)

#     params2 = batman.TransitParams()
#     params2.t0 = t0
#     params2.per = period
#     params2.rp = rp_z
#     params2.a = a
#     params2.inc = inc
#     params2.ecc = ecc
#     params2.w = w
#     params2.u = rldc_z
#     params2.fp = fp
#     params2.t_secondary = t_sec 
#     params2.limb_dark = "quadratic"  
#     m2=batman.TransitModel(params2, time_transits_tot_z, transittype='secondary')
#     flux2=m2.light_curve(params2)
#     ans_z=C_z-2.5*np.log10(flux2)

#     # if flag==True:
#     #     const_r=np.ones(len(ans_r))*C_r
#     #     const_z=np.ones(len(ans_z))*C_z
#     #     if np.mean(ans_r)==C_r or np.mean(ans_z)==C_z:
#     #         print('Fit failed for '+title)    
    
#     data_r_all = pd.DataFrame({'time_r': time_transits_tot_r, 'mag_r': mag_transits_tot_r, 'err_r':err_transits_tot_r})
#     data_z_all = pd.DataFrame({'time_z': time_transits_tot_z, 'mag_z': mag_transits_tot_z, 'err_z':err_transits_tot_z})
#     initial_guess_r=np.array(test(time_tot, *p0)[:len(time_transits_tot_r)])
#     initial_guess_z=np.array(test(time_tot, *p0)[len(time_transits_tot_r):])
#     data_r_folded, data_z_folded, ans_r_folded, ans_z_folded, initial_guess_r_folded, initial_guess_z_folded = folded_lightcurve(data_r_all, data_z_all, period, t0_1, show_plot=False, ans_r=ans_r, ans_z=ans_z, initial_guess_r=initial_guess_r, initial_guess_z=initial_guess_z)
#     plt.figure()
#     plt.plot(data_r_folded['phases_r'], data_r_folded['mag_r'],'o', color='lightskyblue', label='data r', markersize=1.)
#     plt.plot(data_z_folded['phases_z'], data_z_folded['mag_z'],'o', color='pink', label='data z', markersize=1.)
#     #plt.plot(data_r_folded['phases_r'], initial_guess_r_folded, '--', color='blue', label='initial guess r')
#     #plt.plot(data_z_folded['phases_z'], initial_guess_z_folded, '--', color='red', label='initial guess z')
#     plt.gca().invert_yaxis()
#     plt.plot(data_r_folded['phases_r'], ans_r_folded, color='blue', label='optimized fit r') #deepskyblue
#     plt.plot(data_z_folded['phases_z'], ans_z_folded, color='red', label='optimized fit z') #blueviolet
#     # plt.xlim(0,1)
#     plt.legend()
#     #plt.ylim(C_z+0.05, C_r-0.05)
#     plt.title('Secondary eclipse fit for P='+str(round(period,4))+ ' using all parameters')
#     plt.xlabel('Phase')
#     plt.ylabel('Magnitude')
#     text_str='fp='+str(fp)+'$\pm$'+str(fp_err)+'\n t_sec='+str(t_sec)+'$\pm$'+str(t_sec_err)
#     plt.text(0.4, C_z-0.05, text_str, fontsize=10, bbox=dict(facecolor='lightgray', alpha=0.5))
#     plt.savefig(os.path.join('figs', 'secondary_eclipse_fit_all_parameters.pdf'))
    #plt.show()
        
    return param, param_cov #, param_all, param_cov_all

