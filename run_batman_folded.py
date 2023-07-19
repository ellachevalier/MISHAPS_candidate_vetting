#!/usr/bin/env python
# coding: utf-8
# %%

# %%

from PyAstronomy.pyasl import foldAt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import batman
import matplotlib.pyplot as plt
from sympy import symbols, solve, Eq
import os
from other_funcs import folded_lightcurve
import sys

def run_batman(p0_guess, period, ecc, w, rldc_r, rldc_z, data_cut_r, data_cut_z, data_r, data_z, title='Transit model', t0_1=None, flag=True):
    
    data_cut_r=data_cut_r.sort_values(by='time_r')
    data_cut_z=data_cut_z.sort_values(by='time_z')
    time_transits_r=np.array(data_cut_r['time_r'])
    mag_transits_r=np.array(data_cut_r['mag_r'])
    err_transits_r=np.array(data_cut_r['err_r'])
    time_transits_z=np.array(data_cut_z['time_z'])
    mag_transits_z=np.array(data_cut_z['mag_z'])
    err_transits_z=np.array(data_cut_z['err_z'])
    
    
    def test(time_tot, t0, rp_r, rp_z, a, inc, C_r, C_z):
        time_z_try = time_tot[time_tot>100000]
        time_z_try = time_z_try - 100000
        time_r_try = time_tot[time_tot<100000]

        params_z = batman.TransitParams()
        params_z.t0 =  t0                     #time of inferior conjunction
        params_z.per = period                 #orbital period
        params_z.rp = rp_z                      #planet radius (in units of stellar radii)
        params_z.a =  a                       #semi-major axis (in units of stellar radii)
        params_z.inc = inc                    #orbital inclination (in degrees) 
        params_z.ecc = ecc                     #eccentricity
        params_z.w =  w                        #longitude of periastron (in degrees) 
        params_z.u =  rldc_z                #limb darkening coefficients [u1, u2] 
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
        params_r.u =  rldc_r                #limb darkening coefficients [u1, u2] 
        params_r.limb_dark = "quadratic"

        m_r=batman.TransitModel(params_r, time_r_try)
        flux_r = m_r.light_curve(params_r)
        mag_r_2 = C_r-2.5*np.log10(flux_r)
        mag = np.append(mag_r_2, mag_z_1)

        return mag

    #if os.path.exists('figs')==False:
        #os.mkdir('figs')
#     if folded==True:
#         data_r=pd.DataFrame({'time_r':time_transits_r, 'mag_r':mag_transits_r, 'err_r':err_transits_r})
#         data_z=pd.DataFrame({'time_z':time_transits_z, 'mag_z':mag_transits_z, 'err_z':err_transits_z})
#         data_r['phase_r']= foldAt(time_transits_r, period, 0.0)
#         data_z['phase_z']= foldAt(time_transits_z, period, 0.0)
#         data_r_cut=data_r[(data_r['phase_r']>0.15) & (data_r['phase_r']<0.35)]
#         data_z_cut=data_z[(data_z['phase_z']>0.15) & (data_z['phase_z']<0.35)]
#         time_transits_r=np.array(data_r_cut['time_r'])
#         time_transits_z=np.array(data_z_cut['time_z'])
#         mag_transits_r=np.array(data_r_cut['mag_r'])
#         mag_transits_z=np.array(data_z_cut['mag_z'])
#         err_transits_r=np.array(data_r_cut['err_r'])
#         err_transits_z=np.array(data_z_cut['err_z'])
        
    
    p0=[p0_guess[0], p0_guess[1], p0_guess[2], p0_guess[3], p0_guess[4], p0_guess[5], p0_guess[6]]
    time_tot = np.append(time_transits_r, time_transits_z+100000)
    mag_tot = np.append(mag_transits_r, mag_transits_z)
    err_tot = np.append(err_transits_r, err_transits_z)

    lower_bounds = np.array([-np.inf, 0.0, 0.0, 0.0, 70.0, 0., 0.])
    upper_bounds = np.array([np.inf, 10.0, 10.0, 40.0, 130.0, 40., 40.])
    
    #param, param_cov = curve_fit(test, time_tot, mag_tot, p0=p0, sigma=err_tot, bounds=(lower_bounds, upper_bounds))
    param, param_cov = curve_fit(test, time_tot, mag_tot, p0=p0, sigma=err_tot)

    t0 = param[0]
    rp_r = param[1]
    rp_z = param[2]
    a = param[3]
    inc = param[4]
    C_r = param[5]
    C_z = param[6]

    params1 = batman.TransitParams()
    params1.t0 = t0
    params1.per = period
    params1.rp = rp_r
    params1.a = a
    params1.inc = inc
    params1.ecc = ecc
    params1.w = w
    params1.u = rldc_r               
    params1.limb_dark = "quadratic"  
    m1=batman.TransitModel(params1, time_transits_r)
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
    params2.limb_dark = "quadratic"  
    m2=batman.TransitModel(params2, time_transits_z)
    flux2=m2.light_curve(params2)
    ans_z=C_z-2.5*np.log10(flux2)
    
    if flag==True:
        const_r=np.ones(len(ans_r))*C_r
        const_z=np.ones(len(ans_z))*C_z
        if np.mean(ans_r)==C_r or np.mean(ans_z)==C_z:
            print('Fit failed for '+title)    
#Before changing to only fitting phases 0.15 to 0.35
#     if folded==True:
#         data_r_all = pd.DataFrame({'time_r': time_transits_r, 'mag_r': mag_transits_r, 'err_r':err_transits_r})
#         data_z_all = pd.DataFrame({'time_z': time_transits_z, 'mag_z': mag_transits_z, 'err_z':err_transits_z})
#         initial_guess_r=np.array(test(time_tot, *p0)[:len(time_transits_r)])
#         initial_guess_z=np.array(test(time_tot, *p0)[len(time_transits_r):])
#         data_r_folded, data_z_folded, ans_r_folded, ans_z_folded, initial_guess_r_folded, initial_guess_z_folded = folded_lightcurve(data_r_all, data_z_all, period, t0_1, show_plot=False, ans_r=ans_r, ans_z=ans_z, initial_guess_r=initial_guess_r, initial_guess_z=initial_guess_z)
#         plt.figure()
#         plt.plot(data_r_folded['phases_r'], data_r_folded['mag_r'],'o', color='blue', label='data r', markersize=1.)
#         plt.plot(data_z_folded['phases_z'], data_z_folded['mag_z'],'o', color='red', label='data z', markersize=1.)
#         plt.plot(data_r_folded['phases_r'], initial_guess_r_folded, '--', color='blue', label='initial guess r')
#         plt.plot(data_z_folded['phases_z'], initial_guess_z_folded, '--', color='red', label='initial guess z')
#         plt.gca().invert_yaxis()
#         plt.plot(data_r_folded['phases_r'], ans_r_folded, color='blue', label='optimized fit r') #deepskyblue
#         plt.plot(data_z_folded['phases_z'], ans_z_folded, color='red', label='optimized fit z') #blueviolet
#        # plt.xlim(0,1)
#         plt.legend()
#         #plt.ylim(C_z+0.05, C_r-0.05)
#         plt.title(title+' for P='+str(round(period,4)))
#         plt.xlabel('Phase')
#         plt.ylabel('Magnitude')
#         #plt.xlim(0,1)
#         plt.savefig(os.path.join('figs', title+'.png'))
#         #plt.show()
    
    
#     data_r_all = pd.DataFrame({'time_r': time_transits_r, 'mag_r': mag_transits_r, 'err_r':err_transits_r})
#     data_z_all = pd.DataFrame({'time_z': time_transits_z, 'mag_z': mag_transits_z, 'err_z':err_transits_z})
    initial_guess_r=np.array(test(time_tot, *p0)[:len(time_transits_r)])
    initial_guess_z=np.array(test(time_tot, *p0)[len(time_transits_r):])
    data_r_folded_cut, data_z_folded_cut, ans_r_folded, ans_z_folded, initial_guess_r_folded, initial_guess_z_folded = folded_lightcurve(data_cut_r, data_cut_z, period, t0_1, show_plot=False, ans_r=ans_r, ans_z=ans_z, initial_guess_r=initial_guess_r, initial_guess_z=initial_guess_z)
    data_r_folded, data_z_folded = folded_lightcurve(data_r, data_z, period, t0_1, show_plot=False) #ans_r=ans_r, ans_z=ans_z, initial_guess_r=initial_guess_r, initial_guess_z=initial_guess_z)
    plt.figure()
    plt.plot(data_r_folded['phases_r'], data_r_folded['mag_r'],'o', color='lightskyblue', label='data r', markersize=1.)
    plt.plot(data_z_folded['phases_z'], data_z_folded['mag_z'],'o', color='pink', label='data z', markersize=1.)
    plt.plot(data_r_folded_cut['phases_r'], initial_guess_r_folded, '--', color='blue', label='initial guess r')
    plt.plot(data_z_folded_cut['phases_z'], initial_guess_z_folded, '--', color='red', label='initial guess z')
    plt.gca().invert_yaxis()
    plt.plot(data_r_folded_cut['phases_r'], ans_r_folded, color='blue', label='optimized fit r') #deepskyblue
    plt.plot(data_z_folded_cut['phases_z'], ans_z_folded, color='red', label='optimized fit z') #blueviolet
       # plt.xlim(0,1)
    plt.legend()
        #plt.ylim(C_z+0.05, C_r-0.05)
    plt.title(title+' for P='+str(round(period,4)))
    plt.xlabel('Phase')
    plt.ylabel('Magnitude')
        #plt.xlim(0,1)
    plt.savefig(os.path.join('figs', title+'.png'))

    
#     else:   
#         plt.figure()
#         plt.plot(time_transits_r, mag_transits_r,'o', color='blue', label='data r', markersize=1.)
#         plt.plot(time_transits_z, mag_transits_z,'o', color='red', label='data z', markersize=1.)
#         #plt.errorbar(time_transits_r, mag_transits_r, err_transits_r, fmt='o', color='turquoise', label='data r', markersize=1.)
#         #plt.errorbar(time_transits_z, mag_transits_z, err_transits_z, fmt='o', color='mediumslateblue', label='data z', markersize=1.)
#         plt.plot(time_transits_r, test(time_tot, *p0)[:len(time_transits_r)], '--', color='lightblue', label='initial guess r')
#         plt.plot(time_transits_z, test(time_tot, *p0)[len(time_transits_r):], '--', color='lightcoral', label='initial guess z')
#         plt.gca().invert_yaxis()
#         plt.plot(time_transits_r, ans_r, color='blue', label='optimized fit r') #deepskyblue
#         plt.plot(time_transits_z, ans_z, color='red', label='optimized fit z') #blueviolet
#         plt.legend()
#         plt.title(title)
#         plt.xlabel('Time')
#         plt.ylabel('Magnitude')
#         plt.savefig(os.path.join('figs', title+'.png'))
#         #plt.show()
    
    
    
    cov_err = np.sqrt(np.diag(param_cov))
    t0_err = cov_err[0]
    rp_r_err = cov_err[1]
    rp_z_err = cov_err[2]
    a_err = cov_err[3]
    inc_err = cov_err[4]
    C_r_err = cov_err[5]
    C_z_err = cov_err[6]
    
    results = {'t0': t0, 't0_err':t0_err, 'rp_r':rp_r, 'rp_r_err': rp_r_err, 'rp_z':rp_z, 'rp_z_err': rp_z_err, 'a':a, 'a_err':a_err, 'inc':inc,'inc_err':inc_err, 'C_r':C_r,'C_r_err':C_r_err,'C_z':C_z,'C_z_err':C_z_err}
    #print('Curve fit parameters:\nt0 = '+str(round(t0,5))+' Err: '+str(round(t0_err,5))+'\nper = '+str(period)+'\nrp_r = '+str(round(rp_r,5))+' Err: '+str(round(rp_r_err,5))+'\nrp_z = '+str(round(rp_z,5))+' Err: '+str(round(rp_z_err,5))+'\na = '+str(round(a,5))+' Err: '+str(round(a_err,5))+'\ninc = '+str(round(inc,5))+' Err: '+str(round(inc_err))+'\nC_r = '+str(round(C_r,5))+' Err: '+str(round(C_r_err,5))+'\nC_z = '+str(round(C_z,5))+' Err: '+str(round(C_z_err,5)))
    #print('\nPlanet radius r-band: '+str(round(rp_r,5))+' Err: '+str(round(rp_r_err,5)))
    #print('Planet radius z-band: '+str(round(rp_z,5))+' Err: '+str(round(rp_z_err,5)))
    
    model_tot = np.append(ans_r, ans_z)

    N=len(mag_tot)
    x2 = np.sum(((mag_tot-model_tot)**2)/((err_tot)**2))
    C = symbols('C')
    eq = Eq((1/C**2)*x2*(1/N), 1)
    sol = solve(eq)
    const = float(sol[1])
    x2_new=(1/const**2)*x2
    err_tot_adj = err_tot*const
    
    return param, cov_err, err_tot_adj, results

