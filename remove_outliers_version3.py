#!/usr/bin/env python
# coding: utf-8
# %%

# %%
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import batman
import matplotlib.pyplot as plt

def weighted_avg_and_std(values, errors):
    """
    Return the weighted average and standard deviation.

    values, weights -- NumPy ndarrays with the same shape.
    """
    weights = 1/np.square(errors)
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

def remove_outliers(p0_guess, period, ecc, w, rldc_r, rldc_z, sigmas, data, transit_nights, show_plot=False):
    time_transits_r = data['time_r']
    time_transits_z = data['time_z']
    mag_transits_r = data['mag_r']
    mag_transits_z = data['mag_z']
    err_transits_r = data['err_r']
    err_transits_z = data['err_z']
    idx_transit=0
    for idx in range(len(time_transits_r)):
        if int(time_transits_r[idx][0]) in transit_nights:
            def test(time_tot, t0, rp_r, rp_z, a, inc, C_r, C_z):
                time_z_try = time_tot[time_tot>100000]
                time_z_try = time_z_try - 100000
                time_r_try = time_tot[time_tot<100000]

                params_z = batman.TransitParams()            
                params_z.t0 =  t0                     #time of inferior conjunction
                params_z.per = period                 #orbital period
                params_z.rp = rp_z                      #planet radius (in units of stellar radii)
                params_z.a =  a                       #semi-major axis (in units of stellar radii)
                params_z.inc = inc                   #orbital inclination (in degrees)                 
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
                params_r.ecc = ecc                     #eccentricity
                params_r.w =  w                       #longitude of periastron (in degrees) 
                params_r.u =  rldc_r                #limb darkening coefficients [u1, u2] 
                params_r.limb_dark = "quadratic"

                m_r=batman.TransitModel(params_r, time_r_try)
                flux_r = m_r.light_curve(params_r)
                mag_r_2 = C_r-2.5*np.log10(flux_r)
                mag = np.append(mag_r_2, mag_z_1)
                return mag


            p0=[p0_guess[0][idx_transit], p0_guess[1], p0_guess[2], p0_guess[3], p0_guess[4], p0_guess[5], p0_guess[6]]
            time_tot = np.append(time_transits_r[idx], time_transits_z[idx]+100000)
            mag_tot = np.append(mag_transits_r[idx], mag_transits_z[idx])
            err_tot = np.append(err_transits_r[idx], err_transits_z[idx])
            idx_transit+=1
            try:
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
                m1=batman.TransitModel(params1, time_transits_r[idx])
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
                m2=batman.TransitModel(params2, time_transits_z[idx])
                flux2=m2.light_curve(params2)
                ans_z=C_z-2.5*np.log10(flux2)

                if show_plot==True:
                    plt.figure()
                    plt.plot(time_transits_r[idx], mag_transits_r[idx],'o', color='blue', label='data r', markersize=1.)
                    plt.plot(time_transits_z[idx], mag_transits_z[idx],'o', color='red', label='data z', markersize=1.)
                    plt.plot(time_transits_r[idx], test(time_tot, *p0)[:len(time_transits_r[idx])], '--',label='initial guess r')
                    plt.plot(time_transits_z[idx], test(time_tot, *p0)[len(time_transits_r[idx]):], '--',label='initial guess z')
                    plt.gca().invert_yaxis()
                    plt.plot(time_transits_r[idx], ans_r, '-', color='blue', label='optimized fit r')
                    plt.plot(time_transits_z[idx], ans_z, '-', color='red', label='optimized fit z')
                    plt.legend()
                    plt.title('Transit model for r- and z- bands before outlier removal, night '+str(transit_nights[idx]))
                    plt.xlabel('Time')
                    plt.ylabel('Magnitude')
                    plt.show()

                res_r = mag_transits_r[idx] - ans_r
                res_z = mag_transits_z[idx] - ans_z
                res_r_mean = np.mean(res_r)
                res_z_mean = np.mean(res_z)

                const_r=np.ones(len(ans_r))*C_r
                const_z=np.ones(len(ans_z))*C_z
                flag=False
                ans_r=np.array(ans_r)
                ans_z=np.array(ans_z)
                #if np.mean(ans_r)==C_r or np.mean(ans_z)==C_z:
                if np.array_equal(ans_r, const_r) or np.array_equal(ans_z, const_z):
                    print('Fit failed.') 
                    flag=True

                sigma=sigmas[idx]
                data_r = pd.DataFrame({'time': time_transits_r[idx], 'mag': mag_transits_r[idx], 'err': err_transits_r[idx], 'res': res_r})
                data_z = pd.DataFrame({'time': time_transits_z[idx], 'mag': mag_transits_z[idx], 'err': err_transits_z[idx], 'res': res_z})
                before_r=len(data_r)
                before_z=len(data_z)
                data_r_final = data_r[(data_r['err']>=0.001)&(data_r['res']<(res_r_mean+sigma*3))&(data_r['res']>(res_r_mean-sigma*3))]
                data_z_final = data_z[(data_z['err']>=0.001)&(data_z['res']<(res_z_mean+sigma*3))&(data_z['res']>(res_z_mean-sigma*3))]
                after_r=len(data_r_final)
                after_z=len(data_z_final)
                print('Night: '+str(int(time_transits_r[idx][0])))
                print('r rejected:', before_r-after_r, ' z rejected:', before_z-after_z)
                frac_rejected_r = 1-(after_r/before_r)
                frac_rejected_z = 1-(after_z/before_z)
                #print('fraction rejected r:', frac_rejected_r, ' fraction rejected z:', frac_rejected_z)

                if frac_rejected_r <= 0.10 and frac_rejected_z <= 0.10 and flag==False:
                    remove_outliers = True
                else:
                    remove_outliers = False
                    print('Night: '+str(int(time_transits_r[idx][0]))+' -Outliers not removed.')
            except:
                remove_outliers = False
                
        else: #non transit nights
            weighted_mean_and_std_r = weighted_avg_and_std(mag_transits_r[idx], err_transits_r[idx])
            weighted_mean_r = weighted_mean_and_std_r[0] 
            weighted_std_r = weighted_mean_and_std_r[1]

            weighted_mean_and_std_z = weighted_avg_and_std(mag_transits_z[idx], err_transits_z[idx])
            weighted_mean_z = weighted_mean_and_std_z[0]
            weighted_std_z = weighted_mean_and_std_z[1]
            
            res_r = np.array(mag_transits_r[idx])-weighted_mean_r
            res_z = np.array(mag_transits_z[idx])-weighted_mean_z
            res_r_mean = np.mean(res_r)
            res_z_mean = np.mean(res_z)
            
            sigma=sigmas[idx]
            data_r = pd.DataFrame({'time': time_transits_r[idx], 'mag': mag_transits_r[idx], 'err': err_transits_r[idx], 'res': res_r})
            data_z = pd.DataFrame({'time': time_transits_z[idx], 'mag': mag_transits_z[idx], 'err': err_transits_z[idx], 'res': res_z})
            before_r=len(data_r)
            before_z=len(data_z)
            data_r_final = data_r[(data_r['err']>=0.001)&(data_r['res']<(res_r_mean+sigma*3))&(data_r['res']>(res_r_mean-sigma*3))]
            data_z_final = data_z[(data_z['err']>=0.001)&(data_z['res']<(res_z_mean+sigma*3))&(data_z['res']>(res_z_mean-sigma*3))]
            after_r=len(data_r_final)
            after_z=len(data_z_final)
            print('Night: '+str(int(time_transits_r[idx][0])))
            print('r rejected:', before_r-after_r, ' z rejected:', before_z-after_z)
            frac_rejected_r = 1-(after_r/before_r)
            frac_rejected_z = 1-(after_z/before_z)
            #print('fraction rejected r:', frac_rejected_r, ' fraction rejected z:', frac_rejected_z)

            if frac_rejected_r <= 0.10 and frac_rejected_z <= 0.10:
                remove_outliers = True
            else:
                remove_outliers = False
                print('Night: '+str(int(time_transits_r[idx][0]))+' -Outliers not removed.')
            
            
        if remove_outliers == True:
            time_transits_r[idx]=np.array(data_r_final['time'])
            time_transits_z[idx]=np.array(data_z_final['time'])
            mag_transits_r[idx]=np.array(data_r_final['mag'])
            mag_transits_z[idx]=np.array(data_z_final['mag'])
            err_transits_r[idx]=np.array(data_r_final['err'])
            err_transits_z[idx]=np.array(data_z_final['err'])
            print('Night: '+str(int(time_transits_r[idx][0]))+' -Outliers removed.')
        
        if show_plot==True and remove_outliers ==True:
            plt.figure()
            plt.plot(time_transits_r[idx], mag_transits_r[idx], 'o', color='blue', label='data r')
            plt.plot(time_transits_z[idx], mag_transits_z[idx], 'o', color='red', label='data z')
            plt.gca().invert_yaxis()
            plt.legend()
            plt.title('Data after outlier removal, night '+str(transit_nights[idx]))
            plt.xlabel('Time')
            plt.ylabel('Magnitude')
            plt.show()
    data['time_r'] = time_transits_r
    data['time_z'] = time_transits_z
    data['mag_r'] = mag_transits_r
    data['mag_z'] = mag_transits_z
    data['err_r'] = err_transits_r
    data['err_z'] = err_transits_z
    #data=data.sort_values(by='time_r')
    return data

