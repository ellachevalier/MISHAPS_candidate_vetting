#!/usr/bin/env python
# coding: utf-8
# %%
#separate r and z band flux parameters, fix time and vary baseline
#double check that primary and secondary cant be fit at the same time

# %%


#results_batman_all, results_mcmc_all, results_mcmc_per_all, all_guesses, mcmc_chains_all = folded_all_data(data_new, t0_mcmc, p0_guess, best_period, ecc, w, rldc_r, rldc_z, show_plot=True, nwalkers=nwalkers, nsteps=nsteps)
import batman
import emcee
import corner
import os
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
    
    time_transits_r = data_new['time_r'].tolist()
    time_transits_z = data_new['time_z'].tolist()
    mag_transits_r = data_new['mag_r'].tolist()
    mag_transits_z = data_new['mag_z'].tolist()
    err_transits_r = data_new['err_r'].tolist()
    err_transits_z = data_new['err_z'].tolist()

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

    time_transits_tot_r = np.array(time_transits_tot_r, dtype='float64')
    mag_transits_tot_r = np.array(mag_transits_tot_r, dtype='float64')
    err_r = np.array(err_transits_tot_r, dtype='float64')
    time_transits_tot_z = np.array(time_transits_tot_z, dtype='float64')
    mag_transits_tot_z = np.array(mag_transits_tot_z, dtype='float64')
    err_z = np.array(err_transits_tot_z, dtype='float64')
    
    data_r_all = pd.DataFrame({'time_r': time_transits_tot_r, 'mag_r': mag_transits_tot_r, 'err_r':err_r})
    data_z_all = pd.DataFrame({'time_z': time_transits_tot_z, 'mag_z': mag_transits_tot_z, 'err_z':err_z})

    time_r = np.array(data_r_all['time_r']) - t0_1 + 0.25*period
    time_z = np.array(data_z_all['time_z']) - t0_1 + 0.25*period
    phases_r= foldAt(time_r, period, 0.0)
    phases_z= foldAt(time_z, period, 0.0)
    
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
    
    data_r_cut=data_r_all #[(data_r_all['phase_r']>0.55) & (data_r_all['phase_r']<0.95)]
    data_z_cut=data_z_all #[(data_z_all['phase_z']>0.55) & (data_z_all['phase_z']<0.95)]
    data_r_cut=data_r_cut.sort_values(by='time_r')
    data_z_cut=data_z_cut.sort_values(by='time_z')

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

        m_z=batman.TransitModel(params_z, time_z_try, transittype='secondary')
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

        m_r=batman.TransitModel(params_r, time_r_try, transittype='secondary')
        flux_r = m_r.light_curve(params_r)
        mag_r_2 = C_r-2.5*np.log10(flux_r)
        mag = np.append(mag_r_2, mag_z_1)

        return mag

        #if os.path.exists('figs')==False:
            #os.mkdir('figs')
    data_r_cut=data_r_cut.sort_values(by='time_r')
    data_z_cut=data_z_cut.sort_values(by='time_z')
    time_transits_r=np.array(data_r_cut['time_r'])
    mag_transits_r=np.array(data_r_cut['mag_r'])
    err_transits_r=np.array(data_r_cut['err_r'])
    time_transits_z=np.array(data_z_cut['time_z'])
    mag_transits_z=np.array(data_z_cut['mag_z'])
    err_transits_z=np.array(data_z_cut['err_z'])
    
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
    
    #data_r_all = pd.DataFrame({'time_r': time_transits_tot_r, 'mag_r': mag_transits_tot_r, 'err_r':err_transits_tot_r})
    #data_z_all = pd.DataFrame({'time_z': time_transits_tot_z, 'mag_z': mag_transits_tot_z, 'err_z':err_transits_tot_z})
    initial_guess_r=np.array(test(time_tot, *p0)[:len(time_transits_tot_r)])
    initial_guess_z=np.array(test(time_tot, *p0)[len(time_transits_tot_r):])
    data_r_folded_cut, data_z_folded_cut, ans_r_folded, ans_z_folded, initial_guess_r_folded, initial_guess_z_folded = folded_lightcurve(data_r_cut, data_z_cut, period, t0_1, show_plot=False, ans_r=ans_r, ans_z=ans_z, initial_guess_r=initial_guess_r, initial_guess_z=initial_guess_z)
    data_r_folded, data_z_folded = folded_lightcurve(data_r_all, data_z_all, period, t0_1, show_plot=False)
    plt.figure()
    plt.plot(data_r_folded['phases_r'], data_r_folded['mag_r'],'o', color='blue', alpha=0.25, label='data r', markersize=1.)
    plt.plot(data_z_folded['phases_z'], data_z_folded['mag_z'],'o', color='red', alpha=0.25, label='data z', markersize=1.)
    plt.plot(data_r_folded_cut['phases_r'], initial_guess_r_folded, '--', color='blue', label='initial guess r')
    plt.plot(data_z_folded_cut['phases_z'], initial_guess_z_folded, '--', color='red', label='initial guess z')
    plt.gca().invert_yaxis()
    plt.plot(data_r_folded_cut['phases_r'], ans_r_folded, color='blue', label='optimized fit r') #deepskyblue
    plt.plot(data_z_folded_cut['phases_z'], ans_z_folded, color='red', label='optimized fit z') #blueviolet
    # plt.xlim(0,1)
    #plt.legend()
    #plt.ylim(C_z+0.05, C_r-0.05)
    #plt.title(title+' for P='+str(round(period,4)))
    plt.xlabel('Phase')
    plt.ylabel('Magnitude')
    text_str='fp_r='+str(round(fp_r,5))+'$\pm$'+str(round(fp_r_err,5))+'\n fp_z='+str(round(fp_z,5))+'$\pm$'+str(round(fp_z_err,5))+'\n C_r='+str(round(C_r,5))+'$\pm$'+str(round(C_r_err,5))+'\n C_z='+str(round(C_z,5))+'$\pm$'+str(round(C_z_err,5))
    plt.text(0.4, C_z-0.03, text_str, fontsize=10, bbox=dict(facecolor='lightgray', alpha=0.5))
    plt.title('Secondary eclipse fit '+str(period))
    plt.savefig(os.path.join('figs', 'secondary_eclipse_fit.pdf'))
    #plt.show()
    
    def mcmc_model(theta, time, ld_coeff_r, ld_coeff_z, period, ecc, w, t0, rp_r, rp_z, a, inc, t_sec):
        fp_r, fp_z, C_r, C_z = theta

        time_z_try = time[time>100000]
        time_z_try = time_z_try - 100000
        time_r_try = time[time<100000]

        mcmc_params = batman.TransitParams()
        mcmc_params.t0 = t0
        mcmc_params.per = period
        mcmc_params.rp = rp_r
        mcmc_params.a = a
        mcmc_params.inc = inc
        mcmc_params.ecc = ecc
        mcmc_params.w = w
        mcmc_params.fp = fp_r
        mcmc_params.t_secondary = t_sec
        mcmc_params.limb_dark = 'quadratic'
        mcmc_params.u = ld_coeff_r
        mcmc_model = batman.TransitModel(mcmc_params, time_r_try, transittype='secondary')
        corrected1 = C_r - 2.5*np.log10(mcmc_model.light_curve(mcmc_params))

        mcmc_params = batman.TransitParams()
        mcmc_params.t0 = t0
        mcmc_params.per = period
        mcmc_params.rp = rp_z
        mcmc_params.a = a
        mcmc_params.inc = inc
        mcmc_params.ecc = ecc
        mcmc_params.w = w
        mcmc_params.fp = fp_z
        mcmc_params.t_secondary = t_sec
        mcmc_params.limb_dark = 'quadratic'
        mcmc_params.u = ld_coeff_z
        mcmc_model = batman.TransitModel(mcmc_params, time_z_try, transittype='secondary')
        corrected2 = C_z - 2.5*np.log10(mcmc_model.light_curve(mcmc_params))

        return np.append(corrected1, corrected2)

    def mcmc_likelihood(theta, time, mag, m_err, ld_coeff_r, ld_coeff_z, period, ecc, w, t0, rp_r, rp_z, a, inc, t_sec):
        fp_r, fp_z, C_r, C_z = theta
        model = mcmc_model(theta, time, ld_coeff_r, ld_coeff_z, period, ecc, w, t0, rp_r, rp_z, a, inc, t_sec)

        return -0.5*np.sum(((mag-model)**2)/(m_err**2))

    def mcmc_priors(theta, lower_bounds, upper_bounds):
        fp_r, fp_z, C_r, C_z = theta
        if lower_bounds[0] < t0 < upper_bounds[0] and \
            lower_bounds[1] < rp_r < upper_bounds[1] and \
            lower_bounds[2] < rp_z < upper_bounds[2] and \
            lower_bounds[3] < a < upper_bounds[3]:
                return 0.0
        return -np.inf

    def mcmc_prob(theta, time, mag, m_err, lower_bounds, upper_bounds, ld_coeff_r, ld_coeff_z, period, ecc, w, t0, rp_r, rp_z, a, inc, t_sec):
        lp = mcmc_priors(theta, lower_bounds, upper_bounds)
        if not np.isfinite(lp):
            return -np.inf
        probval = lp + mcmc_likelihood(theta, time, mag, m_err, ld_coeff_r, ld_coeff_z, period, ecc, w, t0, rp_r, rp_z, a, inc, t_sec)
        return probval
    
    def mcmc_all(ndim, nwalkers, nsteps, nburn, nthin, guesses, labels, time, mag, err, lower_bounds, upper_bounds, rldc_r, rldc_z, period, ecc, w, t0, rp_r, rp_z, a, inc, t_sec, title='mcmc model', rp_corner=False, sigmas=[10e-3, 10e-3, 10e-3, 10e-3]):

        pos = [guesses[0] + sigmas[0]*np.random.randn(nwalkers),
               guesses[1] + sigmas[1]*np.random.randn(nwalkers),
               guesses[2] + sigmas[2]*np.random.randn(nwalkers),
               guesses[3] + sigmas[3]*np.random.randn(nwalkers)]
        pos = np.vstack(pos).T
        sampler = emcee.EnsembleSampler(nwalkers, ndim, 
                                        mcmc_prob, 
                                        args=(time,
                                              mag,
                                              err,
                                              lower_bounds,
                                              upper_bounds,
                                              rldc_r, 
                                              rldc_z, 
                                              period, ecc, w, t0, rp_r, rp_z, a, inc, t_sec)
                                        )
        sampler.run_mcmc(pos, nsteps, progress=True)
        flat_samples = sampler.get_chain(discard=nburn, thin=nthin, flat=True)

        columns = ['fp_r', 'fp_z', 'C_r', 'C_z']
        df = pd.DataFrame(flat_samples, columns=columns)

        mcmc_fp_rs = flat_samples[:,0]
        mcmc_fp_zs = flat_samples[:,1]
        mcmc_C_rs = flat_samples[:,2]
        mcmc_C_zs = flat_samples[:,3]

        mcmc_chains = [list(mcmc_fp_rs), list(mcmc_fp_zs), list(mcmc_C_rs), list(mcmc_C_zs)]

        
        fp_r=np.nanmean(mcmc_fp_rs)
        fp_z=np.nanmean(mcmc_fp_zs)
        C_r=np.nanmean(mcmc_C_rs)
        C_z=np.nanmean(mcmc_C_zs)

        fp_r_err=np.nanstd(mcmc_fp_rs)
        fp_z_err=np.nanstd(mcmc_fp_zs)
        C_r_err=np.nanstd(mcmc_C_rs)
        C_z_err=np.nanstd(mcmc_C_zs)

        fp_r_median=np.median(mcmc_fp_rs)
        fp_z_median=np.median(mcmc_fp_zs)
        C_r_median=np.median(mcmc_C_rs)
        C_z_median=np.median(mcmc_C_zs)

        means=[fp_r, fp_z, C_r, C_z]
        stds=[fp_r_err, fp_z_err, C_r_err, C_z_err]
        medians=[fp_r_median, fp_z_median, C_r_median, C_z_median]

        #fig=plt.figure()
        fig=corner.corner(flat_samples, labels=labels, truths=guesses, num='corner',
                        rasterized=True, quartile=[0.16, 0.84])
        axes = np.array(fig.axes).reshape((ndim, ndim))
        fig.suptitle(title+' Corner Plot', fontsize=24)
        #fig.set_title(title+' corner plot')
    
        i=0

        from IPython.display import display, Math
        per50=[]
        per16=[]
        per84=[]
        for i in range(ndim):
            mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            per50.append(mcmc[1])
            per16.append(q[0])
            per84.append(q[1])
            txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
            txt = txt.format(mcmc[1], q[0], q[1], labels[i])
            #display(Math(txt))
            ax=axes[i, i]
            #ax.set_title(str(round(means[i],5))+' $\pm$ '+str(round(stds[i],5)), fontsize=15., fontdict={'family':'sans-serif'})
            ax.set_title(columns[i]+": "+r'$'+str(round(mcmc[1],4))+'_{'+str(round(q[0],4))+'}^{'+str(round(q[1],4))+'}$')
            i+=1
        #plt.show()
        fig.savefig(os.path.join('figs', 'secondary_corner_plot.pdf'))

        per_mcmc_results=pd.DataFrame({'50th':per50, '16th':per16, '84th':per84, 'Median':medians})
        #plot just rp corner plot


        fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True, rasterized=True)
        samples = sampler.get_chain()
        for i in range(ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("step number");
        #fig.title(title+' walkers')
        fig.suptitle(title+' Walkers', fontsize=14)
        fig.savefig(os.path.join('figs', 'secondary_walkers_plot.pdf'))

        print('mcmc parameters:')
        print('fp_r: ', round(np.nanmean(mcmc_fp_rs),5), 'Err:', round(np.nanstd(mcmc_fp_rs),5))
        print('fp_z: ', round(np.nanmean(mcmc_fp_zs),5), 'Err:', round(np.nanstd(mcmc_fp_zs),5))
        print('C_r:  ', round(np.nanmean(mcmc_C_rs),5), 'Err:', round(np.nanstd(mcmc_C_rs),5))
        print('C_z:  ', round(np.nanmean(mcmc_C_zs),5), 'Err:', round(np.nanstd(mcmc_C_zs),5))
    
    labels = ['fp_r', 'fp_z', 'C_r', 'C_z']
    
    lower_bounds = [-np.inf, -np.inf, -np.inf, -np.inf]
    upper_bounds = [np.inf, np.inf, np.inf, np.inf]
    
    ndim=4
    nwalkers=16
    nsteps=2000
    nburn=50
    nthin=125
    guesses=param

    mcmc_all(ndim, nwalkers, nsteps, nburn, nthin, guesses, labels, time_tot, mag_tot, err_tot, lower_bounds, upper_bounds, rldc_r, rldc_z, period, ecc, w, t0, rp_r, rp_z, a, inc, t_sec, title='mcmc model', rp_corner=False, sigmas=[10e-4, 10e-4, 10e-4, 10e-4])
    
    fig_list=os.listdir('figs')
    fig_list=[item for item in fig_list if ('walkers' in item or 'corner_plot' in item or 'fit' in item) and 'secondary' in item]
    
    pdfs = fig_list

    merger = PdfMerger()

    for pdf in pdfs:
        merger.append(os.path.join('figs', pdf))

    if os.path.exists('results')==False:
        os.mkdir('results')
    merger.write(os.path.join('results',"secondary_results"+"_period_"+str(round(period,4))+".pdf"))
    merger.close()
    
    for item in fig_list:
        os.remove(os.path.join('figs', item))
        
    return param, param_cov #, param_all, param_cov_all

