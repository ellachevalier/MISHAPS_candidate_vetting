#!/usr/bin/env python
# coding: utf-8
# %%
import argparse
import os
import numpy as np
from other_funcs import *
from remove_outliers_version3 import remove_outliers
from fit_individual_transits import fit_individual_transits
from compare_odd_and_even_attempt2 import compare_odd_and_even
from folded_all_data import folded_all_data
from find_period_attempt7 import find_period
from remove_variability_attempt import remove_variability_all
from limb_darkening_calculation import ld_calc
from density_calculation import density_calculation
from total_plot import make_total_plot
from fit_secondary_eclipse_attempt2 import fit_secondary_eclipse


# %%
description='Compares the planet radii of a transiting body between the r- and z- bands and between the odd and even transits.'

parser=argparse.ArgumentParser(description=description)
parser.add_argument("field", metavar='field', type=str, help='field number (ex. F1)')
parser.add_argument("chip", metavar='chip', type=str, help='chip number (ex. N19)')
parser.add_argument("star_id",metavar='star_id', type=str, help='8 digit star id number')
parser.add_argument("period", metavar='period', type=float, help='enter initial guess for period')
parser.add_argument("--Teff", metavar='Teff', type=float, help='effective temperature', default=4872.743365833783)
parser.add_argument("--logg", metavar='logg', type=float, help='logg', default=4.562797429195741)
parser.add_argument("--ecc",metavar='ecc', type=float, help='eccentricity, default is 0.0', default=0.0)
parser.add_argument("--w", metavar='w',type=float, help='longitude of periastron, default is 0.0', default=0.0)
parser.add_argument("--exclude", nargs='+', metavar='exclude', type=int, help='list of nights to exclude', default=None)
parser.add_argument("--exclude_transit", nargs='+', metavar='exclude_transit',type=int, help='list of transit nights to exclude', default=None)
parser.add_argument("--include_transit", nargs='+',metavar='include_transit',type=int, help='list of transit nights to include', default=None)
parser.add_argument("--remove_variability", metavar='remove_variability',type=bool, help='True to remove variability, False otherwise', default=False)
parser.add_argument("--nwalkers_individual", metavar='nwalkers_individual',type=int, help='nwalkers for individual transit mcmc', default=20)
parser.add_argument("--nsteps_individual", metavar='nsteps_individual', type=int, help='nsteps for individual transit mcmc', default=2000)
parser.add_argument("--nwalkers", metavar='nwalkers',type=int, help='nwalkers for folded lc mcmc', default=20)
parser.add_argument("--nsteps", metavar='nsteps',type=int, help='nsteps for folded lc mcmc', default=2000)
parser.add_argument("--even_guesses", nargs='+', metavar='even_guesses', type=float, help='mcmc guesses for even data', default=None)
parser.add_argument("--odd_guesses", nargs='+', metavar='odd_guesses', type=float, help='mcmc guesses for odd data', default=None)
parser.add_argument("--all_guesses", nargs='+', metavar='all_guesses', type=float, help='mcmc guesses for all data', default=None)
parser.add_argument("--initial_guesses", nargs='+', metavar='initial_guesses', help='initial_guesses without t0 (rp_r rp_z a inc C_r C_z)', type=float, default=None)
parser.add_argument("--mcmc_sigmas", nargs='+', metavar='mcmc_sigmas', type=float, help='initial mcmc position errors', default=[10e-5, 10e-3, 10e-3, 10e-3, 10e-3, 10e-5, 10e-5])
parser.add_argument("--cutoff", metavar='cutoff', type=int, help='number of periods to keep and run with models', default=100)
parser.add_argument("--final_cutoff", metavar='final_cutoff', type=int, help='interval of x2 to include abov the lowest x2 period', default=10)
parser.add_argument("--add_period", metavar='add_period', nargs='+', type=float, help='list periods to add to analysis', default=None)
parser.add_argument("--fp_guess", nargs='+', metavar='fp_guess', type=float, help='initial guess for planet-to-star flux ratios [fp_guess_r, fp_guess_z]', default=[0.001, 0.001])
parser.add_argument("--period_fit_parameters", nargs='+', metavar='period_fit_parameters', type=float, help='transit parameters for model used in period fit, [rp_r, rp_z, a, inc, C_r, C_z]', default=None)
parser.add_argument("--C_guess_sec", nargs='+', metavar='C_guess_sec', type=float, help='initial guess for baseline C for secondary eclipse fit [C_guess_r, C_guess_z]', default=None)
parser.add_argument("--ignore_nontransit", nargs='+', metavar='ignore_nontransit', type=bool, help='True: folded lightcurves use only transit night data. False: folded lightcurves use all data.', default=False)
args=parser.parse_args()


field=args.field
chip=args.chip
star_id=args.star_id
period=args.period
Teff=args.Teff
logg=args.logg
ecc=args.ecc
w=args.w
exclude=args.exclude
exclude_transit=args.exclude_transit
include_transit=args.include_transit
remove_variability_test=args.remove_variability
nwalkers_individual=args.nwalkers_individual
nsteps_individual=args.nsteps_individual
nwalkers=args.nwalkers
nsteps=args.nsteps
even_guesses=args.even_guesses
odd_guesses=args.odd_guesses
all_guesses=args.all_guesses
initial_guess_test=args.initial_guesses
mcmc_sigmas=args.mcmc_sigmas
cutoff=args.cutoff
final_cutoff=args.final_cutoff
add_period=args.add_period
fp_guess=args.fp_guess
period_fit_parameters=args.period_fit_parameters
C_guess_sec=args.C_guess_sec
ignore_nontransit=args.ignore_nontransit

# %%
filename_r = os.path.join('data','MISHAPS_'+str(field)+"_"+str(chip)+"_"+str(star_id)+".r.lc.dtr")
filename_z = os.path.join('data','MISHAPS_'+str(field)+"_"+str(chip)+"_"+str(star_id)+".z.lc.dtr")
filename_data = os.path.join('search_outputs','MISHAPS_'+str(field)+"_"+str(chip)+"_r_2019refim.search_output")

qc, qe = ld_calc(Teff, logg)

rldc_r = qc[0].tolist()
rldc_z = qc[1].tolist()
    
#defines full starid including MISHAPS field and chip number
starid='MISHAPS_'+str(field)+'_'+str(chip)+'_'+str(star_id)


# %%
#creates or clears figs directory
if os.path.exists('figs')==False:
        os.mkdir('figs')
else:
    for f in os.listdir('figs'):
        os.remove(os.path.join('figs', f))


#remove variability if needed
if remove_variability_test==True:
    #load search output data
    star_data, transit_data, transit_nights, all_nights = load_data(filename_data, starid, exclude_transit=exclude_transit, include_transit=include_transit)
    remove_variability_all(filename_r, filename_z, transit_data, transit_nights, all_nights, starid)
    filename_r=filename_r+'.corrected2'
    filename_z=filename_z+'.corrected2'
    
star_data, transit_data, transit_nights, all_nights = load_data(filename_data, starid, exclude=exclude, exclude_transit=exclude_transit, include_transit=include_transit)
print(transit_nights)
#split file by transit
data_split = get_split_data(filename_r, filename_z)

data_split_transit = get_transit_data(data_split, transit_nights)
#get initial parameter guesses
p0_guess_new = initial_guesses(period, transit_data, data_split_transit)
if initial_guess_test is None:
    p0_guess = p0_guess_new
else:
    first=[p0_guess_new[0]]
    first.extend(initial_guess_test)
    p0_guess = first
print(p0_guess)

#remove outliers
sigmas=star_data['best_std_outoftransit[n]'].tolist()
data_new = remove_outliers(p0_guess, period, ecc, w, rldc_r, rldc_z, sigmas, data_split, transit_nights)
data_new_transit = get_transit_data(data_new, transit_nights)
#fit individual transits
results_batman, results_mcmc, results_mcmc_per, t0_mcmc, t0_err_mcmc, avg_first_and_last = fit_individual_transits(data_new_transit, transit_nights, p0_guess, period, ecc, w, rldc_r, rldc_z, nwalkers=nwalkers_individual, nsteps=nsteps_individual)

#find best period estimate
#data_period_fit = get_split_data(filename_r, filename_z, all_nights)
if period_fit_parameters is None:
    #period_fit_parameters_final=avg_first_and_last
    period_fit_parameters_final=results_mcmc
else:
    period_fit_parameters_final={}
    period_fit_parameters_final['rp_r']=period_fit_parameters[0]
    period_fit_parameters_final['rp_z']=period_fit_parameters[1]
    period_fit_parameters_final['a']=period_fit_parameters[2]
    period_fit_parameters_final['inc']=period_fit_parameters[3]
    period_fit_parameters_final['C_r']=period_fit_parameters[4]
    period_fit_parameters_final['C_z']=period_fit_parameters[5]
    
initial_best_periods, best_periods, final_periods=find_period(transit_data, data_new, period_fit_parameters_final, t0_mcmc, t0_err_mcmc, all_nights, period, [rldc_r, rldc_z], cutoff=cutoff, final_cutoff=final_cutoff)
make_chi2_csv(initial_best_periods, best_periods, field, chip, star_id)

if add_period is not None:
    final_periods.extend(add_period)

#folded lightcurve comparison between odd and even
p0_guess[0]=t0_mcmc[0]
if os.path.exists('chains')==False:
        os.mkdir('chains')
contents=[]
densities_r={}
densities_z={}
chains={}

if ignore_nontransit==True:
    final_data=data_new_transit
else:
    final_data=data_new


for best_period in final_periods:
    results_batman_even, results_batman_odd, results_mcmc_even, results_mcmc_odd, results_mcmc_per_even, results_mcmc_per_odd, even_guesses, odd_guesses, mcmc_chains_even, mcmc_chains_odd = compare_odd_and_even(final_data, t0_mcmc, p0_guess, best_period, ecc, w, rldc_r, rldc_z, all_nights, show_plot=True, nwalkers=nwalkers, nsteps=nsteps, even_guesses=even_guesses, odd_guesses=odd_guesses, mcmc_sigmas=mcmc_sigmas)
    
    results_batman_all, results_mcmc_all, results_mcmc_per_all, all_guesses, mcmc_chains_all = folded_all_data(final_data, t0_mcmc, p0_guess, best_period, ecc, w, rldc_r, rldc_z, show_plot=True, nwalkers=nwalkers, nsteps=nsteps, mcmc_sigmas=mcmc_sigmas, all_guesses=all_guesses)
    
    comparison_plot(results_mcmc_per_even, results_mcmc_per_odd, results_mcmc_per_all, results_batman_even, results_batman_odd, results_batman_all, period=best_period)
    
    make_pdf(best_period, starid)
    
    content=make_latex_tables(results_batman_all, results_batman_even, results_batman_odd, results_mcmc_per_all, results_mcmc_per_even, results_mcmc_per_odd, results_mcmc_all, results_mcmc_even, results_mcmc_odd, best_period)
    
    contents.extend(content)
    
    #errs_r=[results_mcmc_all['t0_err'], results_mcmc_all['rp_r_err'], results_mcmc_all['a_err'], results_mcmc_all['inc'], results_mcmc_all['C_r']]
    density_r, density_r_err = density_calculation(best_period, mcmc_chains_all, band='r')
    #errs_z=[results_mcmc_all['t0_err'], results_mcmc_all['rp_z_err'], results_mcmc_all['a_err'], results_mcmc_all['inc'], results_mcmc_all['C_z']]
    density_z, density_z_err = density_calculation(best_period, mcmc_chains_all, band='z')
    
    densities_r[best_period]=[density_r, density_r_err]
    densities_z[best_period]=[density_z, density_z_err]
    
    chains[best_period]=[mcmc_chains_all, mcmc_chains_even, mcmc_chains_odd]
    
    #fit secondary eclipse
    if C_guess_sec is None:
        p0_sec=[fp_guess[0], fp_guess[1], results_mcmc_per_all['50th'][5], results_mcmc_per_all['50th'][6]]
    else:
        p0_sec=[fp_guess[0], fp_guess[1], C_guess_sec[0], C_guess_sec[1]]
    #p0_sec=[0.001, 0.001, 16.94, 16.94]
    rldc=[rldc_r, rldc_z]
    param_sec, param_sec_cov = fit_secondary_eclipse(data_new, results_mcmc_per_all, t0_mcmc, p0_sec, best_period, rldc)

#make_total_plot(data_new, all_nights, best_periods, period, avg_first_and_last, t0_mcmc, rldc_r, field, chip, star_id)
make_total_plot(data_new, all_nights, best_periods, period, results_mcmc, t0_mcmc, rldc_r, field, chip, star_id)
make_json_file(chains, final_periods, field, chip, star_id)
make_final_pdf(contents, field, chip, star_id, p0_guess, even_guesses, odd_guesses, all_guesses, final_periods, densities_r, densities_z, best_periods, transit_nights, remove_variability=remove_variability_test)


# %%
# Can copy and paste into command line to run function
#python final_script.py 'F1' 'N2' '01020506' 0.9752 --exclude 2 36 48 738 770 --include_transit 33 35 --Teff 4503 --logg 4.739 --add_period 0.9751834577 --initial_guesses 0.12 0.12 9.0 89.0 18.444946 18.0233109
#python final_script.py 'F1' 'N2' '01031115' 1.9413 --Teff 5540 --logg 4.499 --exclude 2 36 48 738 770
#python final_script.py 'F1' 'N19' '01045451' 0.4913 --exclude 2 36 48 738 770 --exclude_transit 740 --include_transit 33 --Teff 5607 --logg 4.481 --nsteps 20000 --remove_variability True --initial_guesses 0.14581235179821908 0.14581235179821908 3.5063794852166175 89.0 17.99 17.99
#python final_script.py 'F1' 'N19' '01027538' 1.9888 --exclude 2 36 48 738 770 --Teff 4398 --logg 4.629 --remove_variability True --nsteps 20000
#python final_script.py 'F1' 'N6' '01050399' 1.0087 --exclude 2 36 48 738 744 770 --nsteps 20000 --remove_variability True --Teff 5895 --logg 4.621 --exclude_transit 768


#works but folded lightcurves don't look good so not good results
#python final_script.py 'F1' 'N13' '01041179' 1.3513 --exclude 2 36 48 738 770 --include_transit 0 --Teff 5508 --logg 4.411 --remove_variability True

#python final_script.py 'F1' 'N18' '01061878' 2.8073 --exclude 2 36 48 738 770 --include_transit 742 --Teff 4196 --logg 4.719 --remove_variability True --period_fit_parameters 0.1 0.1 8.8 90.7 18.25267 18.25121

# %%
