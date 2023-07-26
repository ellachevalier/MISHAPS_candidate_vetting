#!/usr/bin/env python
# coding: utf-8
# %%
import numpy as np
import pandas as pd
from other_funcs import get_split_data
import batman
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

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

#subtract n*period from calculated times to get phase*period. plot as func of time or n
def find_period(transit_data, data, avg_first_and_last, t, t_err, transit_nights, all_nights, period_guess, rldc, ecc=0.0, w=0.0, exclude=None, cutoff=100, final_cutoff=10):
    delta_t=t[-1]-t[0]
    delta_t_err=np.sqrt(t_err[-1]**2+t_err[0]**2)
    rldc_r=rldc[0]
    rldc_z=rldc[1]
    
    time_r_all = data['time_r'].tolist()
    time_z_all = data['time_z'].tolist()
    mag_r_all = data['mag_r'].tolist()
    mag_z_all = data['mag_z'].tolist()
    err_r_all = data['err_r'].tolist()
    err_z_all = data['err_z'].tolist()
    
    N_max = int(delta_t/0.25)
    N_min = int(delta_t/delta_t)
    N = range(N_min, N_max)
    x2_per_period={}
    periods=[]
    period_errs=[]
    x2_ts_list=[]
    x2_scaled=[]
    for n in N:
        period=delta_t/n #for range -err on period to +err and vary periods
        
        #period=1.9888
        #k=int(delta_t/period)
        #print(period)
        ts=[]
        for i in range(n+1):
            ts.append(t[0]+i*period) 

        ts_transit=[]
        for time in t: #switch to individual fits t0
            ts=np.array(ts)
            ts_diff=abs(ts-time)
            ts_new=pd.Series(ts)
            ts_diff=pd.Series(ts_diff)
            idx_min=ts_diff.idxmin()
            ts_transit.append(ts[idx_min])
        ts_transit=np.array(ts_transit)
        transit_nights=np.array(t)
        #print(ts_transit)
        #print(transit_nights)
        t_err=np.array(t_err)
        x2_ts=np.sum(((ts_transit-transit_nights)**2)/((t_err)**2)) #error from t0 from fits
        x2_ts_list.append(x2_ts)
        x2_scaled.append(x2_ts/n)
        periods.append(period)
        

    df=pd.DataFrame({'Period':periods, 'x2':x2_ts_list, 'x2_scaled':x2_scaled})    
    df['N']=delta_t/df['Period']
    df=df.sort_values(by='x2')
    
    if os.path.exists('period_fit_results')==False:
        os.mkdir('period_fit_results')
    
    plt.figure()
    plt.plot(df['N'], df['x2'], 'o')
    #plt.ylim(0,df['x2'].median())
    plt.ylabel('$\chi^2$')
    plt.xlabel('N')
    plt.title('$\chi^2$ vs. N before cutoff')
    plt.savefig(os.path.join('period_fit_results','n_vs_x2_before_cutoff.pdf'))
    #plt.show()
    #print(df)
    #plt.plot(df['n'], df['x2'])
    #df_new=df[df['x2_scaled']<len(transit_nights)]
    df_new=df.head(cutoff)
    
    plt.figure()
    plt.plot(df_new['N'], df_new['x2'], 'o')
    #plt.ylim(0,df['x2'].median()) 
    plt.ylabel('$\chi^2$')
    plt.xlabel('N')
    plt.title('$\chi^2$ vs. N after cutoff')
    plt.savefig(os.path.join('period_fit_results','n_vs_x2_after_cutoff.pdf'))
    #plt.show()
    
    periods_new=df_new['Period'].tolist()
    x2_per_periods=[]
    periods_final=[]
    periods_final_err=[]
    t0s_final=[]
    x2_transit_per_periods=[]
    x2_line_per_periods=[]
    
    #periods_new=[0.482011]
    show_plot=False
    for period in periods_new:
        #period=delta_t/n
        n=round(delta_t/period)
        period_err=delta_t_err/n
        #print(period)
        ts=[]
        ts_err=[]
        
        
        def test(time_tot, t0, period):
            time_z_try = time_tot[time_tot>100000]
            time_z_try = time_z_try - 100000
            time_r_try = time_tot[time_tot<100000]

            params_z = batman.TransitParams()
            params_z.t0 =  t0                     #time of inferior conjunction
            params_z.per = period                 #orbital period
            params_z.rp = rp_z                      #planet radius (in units of stellar radii)
            params_z.a =  a*period/period_guess                      #semi-major axis (in units of stellar radii)
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
            params_r.a =  a*period/period_guess                       #semi-major axis (in units of stellar radii)
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
        rp_r = avg_first_and_last['rp_r']  #0.1269746#0.13915168
        rp_z = avg_first_and_last['rp_z']  #0.1244935#0.14594736
        a = avg_first_and_last['a']     #7.284799#13.67223554
        inc = avg_first_and_last['inc']   #90.0022#87.0215763
        C_r = avg_first_and_last['C_r']   #16.94096#16.94083329
        C_z = avg_first_and_last['C_z']  
        
        p0=[t[0], period] #time_r_all
        time_r_total = np.array([item for sublist in time_r_all for item in sublist])
        time_z_total = np.array([item for sublist in time_z_all for item in sublist])
        mag_r_total = np.array([item for sublist in mag_r_all for item in sublist])
        mag_z_total = np.array([item for sublist in mag_z_all for item in sublist])
        err_r_total = np.array([item for sublist in err_r_all for item in sublist])
        err_z_total = np.array([item for sublist in err_z_all for item in sublist])

        time_tot = np.append(time_r_total, time_z_total+100000)
        mag_tot = np.append(mag_r_total, mag_z_total)
        err_tot = np.append(err_r_total, err_z_total)
        #print(t[0], period)
        #lower_bounds = np.array([t[0]-t_err[0],period-period_err])
        #upper_bounds = np.array([t[0]+t_err[0],period+period_err])

        #param, param_cov = curve_fit(test, time_tot, mag_tot, p0=p0, sigma=err_tot, bounds=(lower_bounds, upper_bounds))
        try:
            param, param_cov = curve_fit(test, time_tot, mag_tot, p0=p0, sigma=err_tot)
            cov_err = np.sqrt(np.diag(param_cov))

            t0_new=param[0]
            t0_new_err=cov_err[0]
            period_new=param[1]
            period_new_err=cov_err[1]
            
#             if period_new==period:
#                 print('fit returned same period')

            params1 = batman.TransitParams()
            params1.t0 = t0_new
            params1.per = period_new
            params1.rp = rp_r
            params1.a = a*period_new/period_guess
            params1.inc = inc
            params1.ecc = ecc
            params1.w = w
            params1.u = rldc_r               
            params1.limb_dark = "quadratic"  
            m1=batman.TransitModel(params1, time_r_total)
            flux1=m1.light_curve(params1)
            ans_r=C_r-2.5*np.log10(flux1)

            params2 = batman.TransitParams()
            params2.t0 = t0_new
            params2.per = period
            params2.rp = rp_z
            params2.a = a*period_new/period_guess
            params2.inc = inc
            params2.ecc = ecc
            params2.w = w
            params2.u = rldc_z                
            params2.limb_dark = "quadratic"  
            m2=batman.TransitModel(params2, time_z_total)
            flux2=m2.light_curve(params2)
            ans_z=C_z-2.5*np.log10(flux2)

            all_nights_df_r = pd.DataFrame({'time_r':time_r_total, 'mag_r':mag_r_total, 'err_r':err_r_total, 'ans_r':ans_r})
            all_nights_df_z = pd.DataFrame({'time_z':time_z_total, 'mag_z':mag_z_total, 'err_z':err_z_total, 'ans_z':ans_z})

            #start of old
            x2_transit_list = []
            x2_line_list = []
            for night in all_nights:
                night_df_r = all_nights_df_r[(all_nights_df_r.iloc[:,0]>night) & (all_nights_df_r.iloc[:,0]<night+1.0)]
                night_df_z = all_nights_df_z[(all_nights_df_z.iloc[:,0]>night) & (all_nights_df_z.iloc[:,0]<night+1.0)]

                time_r=np.array(night_df_r['time_r'])
                time_z=np.array(night_df_z['time_z'])
                mag_r=np.array(night_df_r['mag_r'])
                mag_z=np.array(night_df_z['mag_z'])
                err_r=np.array(night_df_r['err_r'])
                err_z=np.array(night_df_z['err_z'])
                ans_r=np.array(night_df_r['ans_r'])
                ans_z=np.array(night_df_z['ans_z'])

                res_mean_r=np.mean(mag_r-ans_r)
                res_mean_z=np.mean(mag_z-ans_z)
                ans_r=ans_r+res_mean_r
                ans_z=ans_z+res_mean_z

                model_tot = np.append(ans_r, ans_z)
                mag_tot = np.append(mag_r, mag_z)
                err_tot = np.append(err_r, err_z)

                x2_transit=np.sum(((mag_tot-model_tot)**2)/((err_tot)**2))
                x2_transit_list.append(x2_transit)

                weighted_mean_and_std_r = weighted_avg_and_std(mag_r, err_r)
                weighted_mean_r = weighted_mean_and_std_r[0] #shift by residuals (data-model) add mean of residuals to model
                weighted_std_r = weighted_mean_and_std_r[1]

                weighted_mean_and_std_z = weighted_avg_and_std(mag_z, err_z)
                weighted_mean_z = weighted_mean_and_std_z[0]
                weighted_std_z = weighted_mean_and_std_z[1]

                diff_r = mag_r-weighted_mean_r
                diff_z = mag_z-weighted_mean_z

                diff_tot = np.append(diff_r, diff_z)

                x2_line = np.sum(((diff_tot)**2)/((err_tot)**2))
                x2_line_list.append(x2_line)


                if show_plot==True:
                        plt.figure()
                        plt.plot(time_r, mag_r,'o', label='data r', markersize=1.)
                        plt.plot(time_z, mag_z,'o', label='data z', markersize=1.)
                        plt.plot(time_r, ans_r, '--', color='blue', label='model r')
                        plt.plot(time_z, ans_z, '--', color='red', label='model z')
                        plt.legend()
                        plt.title('Night: '+str(night)+' Period:'+str(period))
                        plt.xlabel('Time')
                        plt.ylabel('Magnitude')
                        plt.gca().invert_yaxis()
                        #plt.show()



            periods_final.append(period_new)
            periods_final_err.append(period_new_err)
            t0s_final.append(t0_new)
            #print(len(idx_x2.values()))
            #print(x2_diff_list)
            x2_transit_per_periods.append(sum(x2_transit_list))
            x2_line_per_periods.append(sum(x2_line_list))
        except:
            continue
        
        
    final_df=pd.DataFrame({'Period':periods_final, 'Period Err':periods_final_err, 'x2_transit':x2_transit_per_periods, 'x2_line':x2_line_per_periods, 't0':t0s_final})
    diffs = final_df['x2_transit']-final_df['x2_transit'].min()
    final_df['x2_diff']=diffs
    final_df=final_df.sort_values(by=['x2_transit'], ignore_index=True)
    final_df['N']=delta_t/final_df['Period']

#     min_x2_transit=final_df.iloc[0][2]
#     min_x2_transit
#     final_periods=final_df[final_df['x2_transit']<(min_x2_transit+final_cutoff)]
#     final_periods=final_periods['Period'].tolist()
#     final_periods
    
    #new graph stuff
#     plt.figure()
#     x2s=df['x2']-(df['x2'].min()-1)
#     x2_transits=final_df['x2_transit']-(final_df['x2_transit'].min()-1)
#     plt.plot(df['Period'], x2s, 'o', label='after transit time search')
#     line=final_df['x2_line'].iloc[0]
#     #lines=np.ones(len(df['Period']))*line
#     #print(lines)
#     #lines=lines-final_df['x2_transit'].min()
#     #plt.plot(np.log(df['Period']), np.log(lines))
#     plt.plot(final_df['Period'], x2_transits, 'o', label='after model search')
#     plt.xlim(0.22, 10)
#     plt.legend()
#     plt.xlabel('log(period)')
#     plt.ylabel('log($\chi^2$)')
#     plt.title('$\chi^2$ vs. Period')
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.savefig(os.path.join('period_fit_results','log_x2s_vs_log_periods.pdf'))
    #plt.show()
    
    df=df.sort_values(by='Period')
    df=df.reset_index(drop=True)
    min_periods=[]
    min_x2s=[]

    for idx in range(len(df)):
        if idx == 0:
            if df['x2'].iloc[idx]<df['x2'].iloc[idx+1]:
                min_periods.append(df['Period'].iloc[idx])
                min_x2s.append(df['x2'].iloc[idx])
        elif idx == len(df)-1:
            if df['x2'].iloc[idx]<df['x2'].iloc[idx-1]:
                min_periods.append(df['Period'].iloc[idx])
                min_x2s.append(df['x2'].iloc[idx])
        elif df['x2'].iloc[idx]<df['x2'].iloc[idx-1] and df['x2'].iloc[idx]<df['x2'].iloc[idx+1]:
            min_periods.append(df['Period'].iloc[idx])
            min_x2s.append(df['x2'].iloc[idx])
    
    df_peaks=pd.DataFrame({'Period':min_periods, 'x2':min_x2s})
    
    final_df=final_df.sort_values(by='Period')
    
    min_periods_final=[]
    min_x2s_final=[]
    final_df=final_df.reset_index(drop=True)
    
    for idx in range(len(final_df)):
        if idx == 0:   
            if final_df['x2_transit'].iloc[idx]<final_df['x2_transit'].iloc[idx+1]:
                min_periods_final.append(final_df['Period'].iloc[idx])
                min_x2s_final.append(final_df['x2_transit'].iloc[idx])
        elif idx == len(final_df)-1:
            if final_df['x2_transit'].iloc[idx]<final_df['x2_transit'].iloc[idx-1]:
                min_periods_final.append(final_df['Period'].iloc[idx])
                min_x2s_final.append(final_df['x2_transit'].iloc[idx]) 
        else:
            if final_df['x2_transit'].iloc[idx]<=final_df['x2_transit'].iloc[idx-1] and final_df['x2_transit'].iloc[idx]<=final_df['x2_transit'].iloc[idx+1]:
                min_periods_final.append(final_df['Period'].iloc[idx])
                min_x2s_final.append(final_df['x2_transit'].iloc[idx])
            #makes sure that minimums with only one point in the peak are still counted
            if abs(final_df['Period'].iloc[idx]-final_df['Period'].iloc[idx+1])>0.1 and abs(final_df['Period'].iloc[idx]-final_df['Period'].iloc[idx-1])>0.1:
                min_periods_final.append(final_df['Period'].iloc[idx])
                min_x2s_final.append(final_df['x2_transit'].iloc[idx])
    
    #print(min_periods_final)
    #print(min_x2s_final)
    final_df_peaks=pd.DataFrame({'Period':min_periods_final, 'x2':min_x2s_final})
    depth=transit_data['best_depth[n]'].mean()
    num=transit_data['best_npoints_intransit[n]'].mean()
    err=transit_data['std-all-data'].mean()
    chi2=((depth*num)**2)/((err)**2)
    final_df_peaks=final_df_peaks.sort_values(by='x2')
    #print(final_df_peaks)
    final_periods_old=list(final_df_peaks['Period'])
    #print(final_periods)
    final_df_peaks=final_df_peaks[final_df_peaks['x2']<(final_df_peaks['x2'].min()+4*chi2)]
    final_df_peaks=final_df_peaks.sort_values(by='x2')
    
    final_periods=list(final_df_peaks['Period'])
    #print('final periods')
    #print(final_periods)
    if len(final_periods)>5:
        final_periods=final_periods[:5]
    #if len(final_periods)==1:
        #final_periods.extend(final_periods_old[:2])
    
    plt.figure(figsize=(6.4,5.2)) #6.4 and 4.8 
    x2s=df['x2']-(df['x2'].min()-1)
    x2_transits=final_df['x2_transit']-(final_df['x2_transit'].min()-1)
    plt.plot(df['Period'], x2s, 'o', markersize=2.0, label='after transit time search')
    plt.plot(final_df['Period'], x2_transits, 'o', markersize=2.0, label='after model search')
    for period in list(final_df_peaks['Period']):
        plt.axvline(x=period, color='black', linewidth=0.5)
    plt.xlim(0.22, 10)
    plt.xlabel('log(period)')
    plt.ylabel('log($\chi^2$)')
    plt.title('$\chi^2$ vs. Period')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    #plt.savefig(os.path.join('poster_images','log_x2s_vs_log_periods_with_final_periods.png'), dpi=300)
    plt.savefig(os.path.join('period_fit_results','log_x2s_vs_log_periods_with_final_periods.pdf'))
    
    plt.figure()
    plt.plot(final_df['Period'],final_df['x2_diff'], 'o', markersize=2.0, label='$\Delta\chi^2$')
    line=final_df['x2_line'].iloc[0]
    line_change=line-final_df['x2_transit'].min()
    lines=np.ones(len(final_df['Period']))*line_change
    plt.plot(final_df['Period'],lines,label='$\chi^2$ for straight line')
    plt.title('$\Delta\chi^2$ vs. Period With Straight Line Comparison')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Period')
    plt.ylabel('$\Delta\chi^2$')
    plt.legend()
    plt.savefig(os.path.join('period_fit_results','deltax2_vs_period.pdf'))
    #plt.xlim(0.22, 10)
    
    df_new=df[df['Period']<2.5]
    x2s=df_new['x2']-df_new['x2'].min()
    x2_transits=final_df['x2_transit']-final_df['x2_transit'].min()
    #plt.plot(df_new['Period'], x2s, 'o')
    line=final_df['x2_line'].iloc[0]
    #lines=np.ones(len(df['Period']))*line
    #print(lines)
    #lines=lines-final_df['x2_transit'].min()
    #plt.plot(np.log(df['Period']), np.log(lines))
    #plt.plot(final_df['Period'], x2_transits, 'o')
    #for period in list(final_df_peaks['Period']):
        #plt.axvline(x=period, color='red', linewidth=0.5)
    #for period in list(df_peaks['Period']):
        #plt.axvline(x=period, color='black', linewidth=0.5)
    #plt.show()
    
    df=df.sort_values(by='x2')
    final_df=final_df.sort_values(by='x2_transit')
    return df, final_df, final_periods


# %%
# avg_first_and_last={'rp_r':0.1269746, 'rp_z':0.1244935, 'a':7.284799, 'inc':87.0, 'C_r':16.94096, 'C_z':16.94096}
# t_err=[0.0005142076083922166,
#  0.0009464611607341753,
#  0.0016085592302390967,
#  0.0007151504854749141,
#  0.0006273961957576624,
#  0.0006932796775990423,
#  0.000709058912092793,
#  0.004958962054308873,
#  0.0043480846363359695,
#  0.019423327685906416]
# t=[1.8097169851851527,
#  33.63632338300558,
#  35.629662950420844,
#  737.6819395435072,
#  739.6716860934674,
#  741.6601971528315,
#  743.6472319509671,
#  765.5215095952187,
#  767.5110501190792,
#  769.4822788167637]
# delta_t=t[-1]-t[0]
# delta_t_err=np.sqrt(t_err[-1]**2+t_err[0]**2)
# filename_r='MISHAPS_F1_N19_01027538.r.lc.dtr.corrected1'
# filename_z='MISHAPS_F1_N19_01027538.z.lc.dtr.corrected1'
# period_guess = 1.9888
# #filename_r='MISHAPS_F1_N19_01045451.r.lc.dtr'
# #filename_z='MISHAPS_F1_N19_01045451.z.lc.dtr'
# #period_guess=0.4913

# ecc = 0.0
# w = 0.0
# rldc_r = [0.6502, 0.0872] #limb darkening
# rldc_z = [0.4678, 0.1118]
# rldc=[rldc_r, rldc_z]
# columns=['starid', 'lcname', 'starflag', 'nightflag', 'int(night)', 'median-all-data', 'std-all-data', 'tcbest[n]',
# 'best_depth[n]', 'best_depth_sigma[n]', 'best_duration[n]', 'snrmax[n]', 'bestmedian_intransit[n]',
# 'best_std_intransit[n]', 'best_npoints_intransit[n]', 'best_median_outoftransit[n]', 'best_std_outoftransit[n]', 'best_npoints_outoftransit[n]', 'detections_so_far']
# data=pd.read_csv('MISHAPS_F1_N19_r_2019refim.search_output', sep='\s+',
#                   header=None)
# data.columns=columns

# starid=filename_r.split('.')[0]
# star_data = data[data['starid']==starid]
# transit_data = star_data[star_data['snrmax[n]']>5.0]
# transit_nights = transit_data['int(night)'].tolist()
# all_nights = star_data['int(night)'].tolist()
# tcbest=star_data['tcbest[n]'].tolist()
# #all_nights
# #data_split = get_split_data(filename_r, filename_z, all_nights)

# transit_nights = transit_data['int(night)'].tolist()
# #print(transit_nights)
# data_split = get_split_data(filename_r, filename_z, all_nights)
# #print(data_split.head())

# df, final_df, final_periods=find_period(transit_data, data_split, avg_first_and_last, t, t_err, transit_nights, all_nights, period_guess, rldc, ecc=0.0, w=0.0, exclude=None, cutoff=100, final_cutoff=10)

# %%
#Can copy and paste into command line to run function
#python final_script.py 'F1' 'N2' '01020506' 0.9752 --exclude 2 36 48 738 770 --include_transit 33 35 --Teff 4503 --logg 4.739 --add_period 0.9751834577 
#python final_script.py 'F1' 'N2' '01031115' 1.9413 --Teff 5540 --logg 4.499 --exclude 2 36 48 738 770
#python final_script.py 'F1' 'N19' '01045451' 0.4913 --exclude 2 36 48 738 770 --exclude_transit 740 --include_transit 33 --Teff 5607 --logg 4.481 --nsteps 20000 --even_guesses 0.59078904 0.46177587 0.47236939 2.20804344 55.45659689 17.99321475 17.923996
#python final_script.py 'F1' 'N19' '01027538' 1.9888 --exclude 2 36 48 738 770 --Teff 4398 --logg 4.629 --remove_variability True

# %%
#python final_script.py 'F1' 'N2' '01031115' 1.9413 --Teff 5540 --logg 4.499 --exclude 2 36 48 738 770
#python final_script.py 'F1' 'N19' '01045451' 0.4913 --exclude 2 36 48 738 770 --exclude_transit 740 --include_transit 33 --Teff 5607 --logg 4.481 --nsteps 20000 --even_guesses 0.59078904 0.46177587 0.47236939 2.20804344 55.45659689 17.99321475 17.923996

# %%
#python final_script.py 'F1' 'N2' '01020506' 0.9752 --exclude 2 36 48 738 770 --include_transit 33 35 --Teff 4503 --logg 4.739 --add_period 0.9751834577 
