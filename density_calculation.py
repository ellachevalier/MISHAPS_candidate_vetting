#!/usr/bin/env python
# coding: utf-8
# %%

# %%


# From Ali's code
import numpy as np

def impact_param_eqn(a_over_rstar, inclination):
    '''Solve for impact parameter b'''
    inclination_rad = inclination*(np.pi/180.0)
    impact_param = a_over_rstar*np.cos(inclination_rad)
    return impact_param

def duration_eqn(period, rp_over_rstar, a_over_rstar, impact_param, 
                 inclination):
    '''Solve for fit duration using Seager 2003 eqn 3'''
    inclination_rad = inclination*(np.pi/180.0)
    numerator = (1 + rp_over_rstar)**2 - impact_param**2
    denominator = 1 - (np.cos(inclination_rad))**2
    duration = (period/np.pi)*np.arcsin((1/a_over_rstar)*np.sqrt(numerator/denominator))
    return duration

def density_eqn(period, depth, impact_param, duration):
    '''Get stellar density using Seager 2003 eqn 9'''
    coeff = ((365.25**2)/(215**3))*(1/period**2)
    numerator = (1+np.sqrt(depth))**2 - (impact_param**2)*(1-(np.sin(duration*np.pi/period))**2)
    denominator = (np.sin(duration*np.pi/period))**2
    density = coeff*(numerator/denominator)**(3/2)
    return density

def basic_density_calc(period, aoverr, sig_aoverr, sig_period=0):
    '''Basic estimate from Sandford 2017'''
    r_sol = 0.00465047 #AU, radius of the sun
    numerator = 3*np.pi*(r_sol**3)*(aoverr**3)
    sig_num = 3*np.pi*(r_sol**3)*sig_aoverr*3*(aoverr**2)
    denominator = 4*np.pi*(period/365.25)**2
    sig_denom = 4*np.pi*2*period*sig_period*(1/365.25)**2
    density = numerator/denominator
    sig_density = numerator/denominator*np.sqrt((sig_num/numerator)**2 + (sig_denom/denominator)**2)
    return density, sig_density, numerator, sig_num, denominator, sig_denom

def density_calculation(period, mcmc_chains, band='r'):
    if band=='r':
        rp=np.array(mcmc_chains[1])
    else:
        rp=np.array(mcmc_chains[2])
    a=np.array(mcmc_chains[3])
    inc=np.array(mcmc_chains[4])
    fit_depth = rp**2
    fit_b = impact_param_eqn(a, inc)
    fit_dur = duration_eqn(period, rp, a, fit_b, inc)
    #errs=np.array(errs)
    # Density
    density = density_eqn(period, fit_depth, fit_b, fit_dur)
    
    density_val = np.mean(density)
    density_err = np.std(density)
    #density_val = basic_density_calc(period, a, errs[2])
    #basic_density = density_val[0]
    #basic_density_err = density_val[1]
    
    return density_val, density_err

# %%

#density, density_val, density_val_err = density_calculation(1.822, 0.147, 6.774, 83.937, 16.939, 1.98881, [0.0001, 0.004, 0.376, 0.695, 0.0001])
#print(density)
#print(density_val)
#print(density_val_err)
#print(density_try(6.774, 1.9881))
#print(density_try(2.419, 0.66294))
#0.7908*1.41


# %%


#density, density_val, density_val_err = density_calculation(1.822, 0.146, 2.419, 73.144, 16.939, 0.66294, [0.0001, 0.005, 0.146, 2.440, 0.0001])
#print(density)
#print(density_val, density_val_err)


# %%


#density, density_val, density_val_err = density_calculation(1.822, 0.138, 2.419, 73.144, 16.939, 0.66294, [0.0001, 0.005, 0.146, 2.440, 0.0001])
#print(density)
#print(density_val, density_val_err)


# %%





# %%




