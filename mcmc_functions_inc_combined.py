#!/usr/bin/env python
# coding: utf-8
# %%
import numpy as np
import pandas as pd
import emcee
from scipy.optimize import curve_fit
import batman
import corner
import matplotlib.pyplot as plt
import os

def mcmc_model(theta, time, ld_coeff_r, ld_coeff_z, period, ecc, w):
    t0, rp_r, rp_z, a, inc, C_r, C_z = theta

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
    mcmc_params.limb_dark = 'quadratic'
    mcmc_params.u = ld_coeff_r
    mcmc_model = batman.TransitModel(mcmc_params, time_r_try)
    corrected1 = C_r - 2.5*np.log10(mcmc_model.light_curve(mcmc_params))

    mcmc_params = batman.TransitParams()
    mcmc_params.t0 = t0
    mcmc_params.per = period
    mcmc_params.rp = rp_z
    mcmc_params.a = a
    mcmc_params.inc = inc
    mcmc_params.ecc = ecc
    mcmc_params.w = w
    mcmc_params.limb_dark = 'quadratic'
    mcmc_params.u = ld_coeff_z
    mcmc_model = batman.TransitModel(mcmc_params, time_z_try)
    corrected2 = C_z - 2.5*np.log10(mcmc_model.light_curve(mcmc_params))

    return np.append(corrected1, corrected2)

def mcmc_likelihood(theta, time, mag, m_err, ld_coeff_r, ld_coeff_z, period, ecc, w):
    t0, rp_r, rp_z, a, inc, C_r, C_z = theta
    model = mcmc_model(theta, time, ld_coeff_r, ld_coeff_z, period, ecc, w)
    return -0.5*np.sum(((mag-model)**2)/(m_err**2))

def mcmc_priors(theta, lower_bounds, upper_bounds):
    t0, rp_r, rp_z, a, inc, C_r, C_z = theta
    if lower_bounds[0] < t0 < upper_bounds[0] and \
        lower_bounds[1] < rp_r < upper_bounds[1] and \
        lower_bounds[2] < rp_z < upper_bounds[2] and \
        lower_bounds[3] < a < upper_bounds[3] and \
        lower_bounds[4] < inc < upper_bounds[4] and \
        lower_bounds[5] < C_r < upper_bounds[5] and \
        lower_bounds[6] < C_z < upper_bounds[6]:
            return 0.0
    return -np.inf

def mcmc_prob(theta, time, mag, m_err, lower_bounds, upper_bounds, ld_coeff_r, ld_coeff_z, period, ecc, w):
    lp = mcmc_priors(theta, lower_bounds, upper_bounds)
    if not np.isfinite(lp):
        return -np.inf
    probval = lp + mcmc_likelihood(theta, time, mag, m_err, ld_coeff_r, ld_coeff_z, period, ecc, w)
    return probval
    
def mcmc_all(ndim, nwalkers, nsteps, nburn, nthin, guesses, labels, time, mag, err, lower_bounds, upper_bounds, rldc_r, rldc_z, period, ecc, w, title='mcmc model', rp_corner=False, sigmas=[10e-5, 10e-3, 10e-3, 10e-3, 10e-3, 10e-5, 10e-5]):

    pos = [guesses[0] + sigmas[0]*np.random.randn(nwalkers),
           guesses[1] + sigmas[1]*np.random.randn(nwalkers),
           guesses[2] + sigmas[2]*np.random.randn(nwalkers),
           guesses[3] + sigmas[3]*np.random.randn(nwalkers), 
           guesses[4] + sigmas[4]*np.random.randn(nwalkers),
           guesses[5] + sigmas[5]*np.random.randn(nwalkers),
           guesses[6] + sigmas[6]*np.random.randn(nwalkers)]
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
                                          period, ecc, w)
                                    )
    sampler.run_mcmc(pos, nsteps, progress=True)
    flat_samples = sampler.get_chain(discard=nburn, thin=nthin, flat=True)
    
    columns = ['t0', 'rp_r','rp_z', 'a', 'inc', 'C_r', 'C_z']
    df = pd.DataFrame(flat_samples, columns=columns)
    
    mcmc_t0s = flat_samples[:,0]
    mcmc_rp_rs = flat_samples[:,1]
    mcmc_rp_zs = flat_samples[:,2]
    mcmc_as = flat_samples[:,3]
    mcmc_incs = flat_samples[:,4]
    mcmc_C_rs = flat_samples[:,5]
    mcmc_C_zs = flat_samples[:,6]
    
    mcmc_chains = [list(mcmc_t0s), list(mcmc_rp_rs), list(mcmc_rp_zs), list(mcmc_as), list(mcmc_incs), list(mcmc_C_rs), list(mcmc_C_zs)]
    
    t0=np.nanmean(mcmc_t0s)
    rp_r=np.nanmean(mcmc_rp_rs)
    rp_z=np.nanmean(mcmc_rp_zs)
    a=np.nanmean(mcmc_as)
    inc=np.nanmean(mcmc_incs)
    C_r=np.nanmean(mcmc_C_rs)
    C_z=np.nanmean(mcmc_C_zs)
    
    t0_err=np.nanstd(mcmc_t0s)
    rp_r_err=np.nanstd(mcmc_rp_rs)
    rp_z_err=np.nanstd(mcmc_rp_zs)
    a_err=np.nanstd(mcmc_as)
    inc_err=np.nanstd(mcmc_incs)
    C_r_err=np.nanstd(mcmc_C_rs)
    C_z_err=np.nanstd(mcmc_C_zs)
    
    t0_median=np.median(mcmc_t0s)
    rp_r_median=np.median(mcmc_rp_rs)
    rp_z_median=np.median(mcmc_rp_zs)
    a_median=np.median(mcmc_as)
    inc_median=np.median(mcmc_incs)
    C_r_median=np.median(mcmc_C_rs)
    C_z_median=np.median(mcmc_C_zs)
    
    means=[t0, rp_r, rp_z, a, inc, C_r, C_z]
    stds=[t0_err, rp_r_err, rp_z_err, a_err, inc_err, C_r_err, C_z_err]
    medians=[t0_median, rp_r_median, rp_z_median, a_median, inc_median, C_r_median, C_z_median]
    
    #fig=plt.figure()
    fig=corner.corner(flat_samples, labels=labels, truths=guesses, num='corner',
                    rasterized=True, quartile=[0.16, 0.84])
    axes = np.array(fig.axes).reshape((ndim, ndim))
    fig.suptitle(title+' Corner Plot', fontsize=24)
    #fig.set_title(title+' corner plot')
    x=np.linspace(0.0, 10.0, 100)
    y=x
    ax=axes[2,1]
    ax.plot(x,y, 'r')
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
    fig.savefig(os.path.join('figs', 'corner_plot '+str(title)+'.pdf'))
    
    per_mcmc_results=pd.DataFrame({'50th':per50, '16th':per16, '84th':per84, 'Median':medians})
    #plot just rp corner plot
    
    if rp_corner==True:
        partial = np.atleast_2d(df.loc[:, ['rp_r', 'rp_z']])
        fig = corner.corner(partial, labels=['rp_r', 'rp_z'])
        axes = np.array(fig.axes).reshape((2, 2))
        ax=axes[0,0]
        ax.set_title('rp_r: '+r'$'+str(round(per50[1],4))+'_{'+str(round(per16[1],4))+'}^{'+str(round(per84[1],4))+'}$')
        ax=axes[1,1]
        ax.set_title('rp_z: '+r'$'+str(round(per50[2],4))+'_{'+str(round(per16[2],4))+'}^{'+str(round(per84[2],4))+'}$')
        x=np.linspace(0.0, 10.0, 100)
        y=x
        ax=axes[1,0]
        ax.plot(x,y, 'r')
        #rp_r_1=[np.nanmean(mcmc_rp_rs), np.nanmean(mcmc_rp_rs)]
        #rp_z_1=[np.nanmean(mcmc_rp_zs), np.nanmean(mcmc_rp_zs)]
        #print(rp_z)
        #corner.overplot_lines(fig, rp_r_1, color="C1")
        #corner.overplot_lines(fig, rp_z_1, color="C2")
        plt.savefig(os.path.join('figs', 'rp_comparision_plot '+str(title)+'.png'))

    
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
    fig.savefig(os.path.join('figs', 'walkers_plot '+str(title)+'.pdf'))
    
    print('mcmc parameters:')
    print('t0:   ', round(np.nanmean(mcmc_t0s),5), 'Err:', round(np.nanstd(mcmc_t0s),5))
    print('rp_r: ', round(np.nanmean(mcmc_rp_rs),5), 'Err:', round(np.nanstd(mcmc_rp_rs),5))
    print('rp_z: ', round(np.nanmean(mcmc_rp_zs),5), 'Err:', round(np.nanstd(mcmc_rp_zs),5))
    print('a:    ', round(np.nanmean(mcmc_as),5), 'Err:', round(np.nanstd(mcmc_as),5))
    print('inc:', round(np.nanmean(mcmc_incs),5), 'Err:', round(np.nanstd(mcmc_incs),5))
    print('C_r:  ', round(np.nanmean(mcmc_C_rs),5), 'Err:', round(np.nanstd(mcmc_C_rs),5))
    print('C_z:  ', round(np.nanmean(mcmc_C_zs),5), 'Err:', round(np.nanstd(mcmc_C_zs),5))
    
    results = {'t0': t0, 't0_err':t0_err, 'rp_r':rp_r, 'rp_r_err': rp_r_err, 'rp_z':rp_z, 'rp_z_err': rp_z_err, 'a':a, 'a_err':a_err, 'inc':inc,'inc_err':inc_err, 'C_r':C_r,'C_r_err':C_r_err,'C_z':C_z,'C_z_err':C_z_err}
    
    return results, per_mcmc_results, mcmc_chains



# %%
#python final_script_function.py 'F1' 'N13' '01041179' 1.3513 --exclude 2 36 48 738 770 --include_transit 0 --nsteps 20000 --initial_guesses 0.09580130224443757 0.09580130224443757 8.2390971854410395 89.0 17.439085927835052 17.554842295081972 --even_guesses 6.31562584e-01 9.48546699e-02 9.67997073e-02 5.33901897e+00 9.67720758e+01 1.74363890e+01 1.75525549e+01 --odd_guesses 6.31562584e-01 9.48546699e-02 9.67997073e-02 5.33901897e+00 9.67720758e+01 1.74363890e+01 1.75525549e+01
