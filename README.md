# MISHAPS_candidate_vetting
 Code to compare transit depths, find valid periods, etc.


In the MISHAPS_candidate_vetting directory, create a folder called 'data' and a folder called 'search_outputs'. The data files and search output files should be stored in these folders. The detrended data files are used, and the name should be in the format MISHAPS_F1_N19_01027538.r.lc.dtr. In the command line, go to the MISHAPS_candidate_vetting directory. Type python final_script.py followed by the appropriate arguments. There are four required arguments: the field number (ex. 'F1'), the chip number (ex. 'N19'), the starid (ex. '01027538'), and the initial guess for the period (ex. 1.9888). The starid must match the starid in the name of the data file. 

There are many optional arguments that can be used to make the code more accurate.

# Optional parameters:

--Teff and --logg can be used to add the effective temperature in Kelvin and the logg value in g/cm^3.

--ecc and --w can be used to change the eccentricity and the longitude of periastron from their default values of 0.0

--exclude can be used to exclude nights with bad data. List the nights to exclude. (ex. --exclude 2 36 48 738 770)

--exclude_transit and --include transit can be used to exclude or include nights as transits. The program will automatically use what the search output determines as transit nights, but nights can be removed or added. Note: excluding a night from the transit nights does not exclude it from the data completely. (ex. --include_transit 0)

--remove_variability (True or False), default is False

--nwalkers_individual and --nsteps_individual can be used to set the nwalkers and nsteps for the mcmc when running the individual transit fits. Default values are nwalkers_individual=20 and nsteps_individual=2000.

--nwalkers and --nsteps can be used to set the nwalkers and nsteps for the mcmc when running the folded lightcurve fits. Default values are nwalkers=20 and nsteps=2000.

--even_guesses, --odd_guesses, and --all_guesses can be used to set the starting guesses for the mcmc for the even, odd, and all folded lightcurves, respectively. List as follows: --even_guesses t0 rp_r rp_z a inc C_r C_z. These arguments are useful if the mcmc for some of the folded lightcurves doesn't work. Look at the first page of the pdf to find what the even, odd, and all guesses were for each run.

--initial_guesses can be used to manually set the initial guesses for the batman fits (excluding t0). List as follows: --initial_guesses rp_r rp_z a inc C_r C_z. 

--mcmc_sigmas can be used to adjust how far away from the initial guesses the mcmc samples. Default values are 10e-5 10e-3 10e-3 10e-3 10e-3 10e-5 10e-5. If mcmc walkers are too narrow, these values can be decreased.

--cutoff and --final_cutoff can be used to adjust the thresholds for keeping periods in the period fit. These values don't need to be changed now because the period fit function was changed so that it does not need these. However, they are kept as arguments in case the period fit is reverted to old methods. 

--add_period can be used to add a period that was not found by the period search. The code will also create folded lightcurves for this period, along with the periods found from the period search.

--fp_guess can be used to adjust the initial guesses for the planet-to-star flux ratio for the secondary eclipse fit. List in this order: --fp_guess fp_r fp_z. Default is 0.001 0.001. 

--C_guess_sec can be used to adjust the initial guesses for the baseline magnitudes for the secondary eclipse fit. List in this order: --C_guess_sec C_r C_z

--period_fit_parameters can be used to adjust the parameters used to create the model in the period search. Currently, the model is created using the average of the parameters found from the individual transit fits. However, if the individual fits do not work well, it could be helpful to manually enter parameters that create a good model. List in the following order: --period_fit_parameters rp_r rp_z a inc C_r C_z

--ignore_nontransit (True or False) can be used to create the folded lightcurves with only the transit night data, ignoring the nontransit night data. Default is False. It could help with the curvefit to use only the transit nights, but if ignoring the nontransit night data is needed than the period is probably not a good fit.



<!-- #region -->
# Other notes:
Below are a few command line arguments that I used when I ran the code. You can copy and paste these if you want. These got me the best results for each candidate, but there is still much room for improvement. I have not thoroughly figured out the best starting initial guesses or even/odd/all guesses. (Make sure to adjust nsteps to your liking, keep at 2000 (default) so it runs faster but change to 20000 for better mcmc results)

Works okay:

python final_script.py 'F1' 'N2' '01020506' 0.9752 --exclude 2 36 48 738 770 --include_transit 33 35 --Teff 4503 --logg 4.739 --add_period 0.9751834577 --initial_guesses 0.12 0.12 9.0 89.0 18.444946 18.0233109

python final_script.py 'F1' 'N2' '01031115' 1.9413 --Teff 5540 --logg 4.499 --exclude 2 36 48 738 770

python final_script.py 'F1' 'N19' '01045451' 0.4913 --exclude 2 36 48 738 770 --exclude_transit 740 --include_transit 33 --Teff 5607 --logg 4.481 --nsteps 20000 --remove_variability True --initial_guesses 0.14581235179821908 0.14581235179821908 3.5063794852166175 89.0 17.99 17.99

python final_script.py 'F1' 'N19' '01027538' 1.9888 --exclude 2 36 48 738 770 --Teff 4398 --logg 4.629 --remove_variability True --nsteps 20000

python final_script.py 'F1' 'N6' '01050399' 1.0087 --exclude 2 36 48 738 744 770 --nsteps 20000 --remove_variability True --Teff 5895 --logg 4.621 --exclude_transit 768


Works but folded lightcurves don't look good, not good results:

python final_script.py 'F1' 'N13' '01041179' 1.3513 --exclude 2 36 48 738 770 --include_transit 0 --Teff 5508 --logg 4.411 --remove_variability True

python final_script.py 'F1' 'N18' '01061878' 2.8073 --exclude 2 36 48 738 770 --include_transit 742 --Teff 4196 --logg 4.719 --remove_variability True --period_fit_parameters 0.1 0.1 8.8 90.7 18.25267 18.25121

<!-- #endregion -->

```python

```
