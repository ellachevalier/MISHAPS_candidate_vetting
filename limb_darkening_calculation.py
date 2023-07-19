#!/usr/bin/env python
# coding: utf-8
# %%

# %%

import seaborn as sb
sb.set_context('notebook')
sb.set_style('ticks')
from ldtk import LDPSetCreator, BoxcarFilter, TabulatedFilter
from ldtk.filters import sdss_g, sdss_r, sdss_i, sdss_z

def ld_calc(Teff, logg, Teff_err=100, logg_err= 0.1):
#Teff_final=4398
#logg_final=4.629
    sc = LDPSetCreator(teff=(Teff, Teff_err), logg=(logg, logg_err), z=(0.0,0.05),
                       filters=[sdss_r, sdss_z], dataset='vis-lowres')

    ps = sc.create_profiles(nsamples=2000)
    ps.resample_linear_z(300)

    qc, qe = ps.coeffs_qd()
    
    return qc, qe


