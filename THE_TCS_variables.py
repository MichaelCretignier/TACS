import os

import numpy as np
import pandas as pd

cwd = os.getcwd()

HARPN_standards = ['HD4628','HD166620','HD122064','HD127334','HD144579','HD151288','HD10476']

NEID_standards = ['HD4628','HD9407','HD10700','HD127334','HD185144','HD211354',
                  'HD179957','HD116442','HD95735','HD166620','HD86728','HD68017',
                  'HD51419']

YARARA = pd.read_csv(cwd+'/TACS_Material/THE_YARARA.csv',index_col=0)

NEID_Gupta_25 = pd.read_csv(cwd+'/TACS_Material/THE_NEID.csv',index_col=0)
NEID_catalog = NEID_Gupta_25.loc[NEID_Gupta_25['Rank_NEID']>0]
NEID_Gupta_25 = NEID_Gupta_25.loc[NEID_Gupta_25['nobs_NEID']!=0]
YARARA_catalog = YARARA.loc[YARARA['nobs_DB']!=0]

cutoff_presurvey = {
    'teff_mean<':6000,
    'logg>':4.2,
    'VSINI<':5,
    'Fe/H>':-0.4,
    'ruwe_GAIA<':1.2,
    'multi_peak_GAIA<':1,
    'sky_contam_VIZIER<':0.1,
    'rv_trend_kms_DACE<':0.1,
    'eff_nights_1.75>':180,
    'season_length_1.75>':240,
    'HJ<':0.5,
    'BDW<':0.5,
    'RHK<':-4.7,
    'dist>':0,
    'gmag<':7.5,
    }

cutoff_tim = {
    'NIGHTS>':240,
    'Teff_phot<':6000,
    'Teff_spec<':6000,
    'Fe/H_spec>':-0.4,
    'parallax>':25,
    'log_ruwe<':np.round(np.log10(1.2),3),
    'vsini_spec<':3,
    'logRHK<':-4.8,
    'logg>':2.0,
    'gmag<':7,
    }

cutoff_jean = {
    'eff_nights_1.5>':160,
    'Teff_gspphot<':6000,
    'Teff_gspphot>':4500,
    'Fe/H>':-2.0,
    'parallax>':0,
    'log_ruwe<':np.round(np.log10(1.2),3),
    'vsini<':6,
    'logRHK<':-4.5,
    'logg_gspphot>':4,
    'gmag<':7,
    }

cutoff_sam = {
    'NIGHTS>':210,
    'Teff_phot<':5900,
    'Teff_phot>':4500,
    'Fe/H>':-2.0,
    'parallax>':0,
    'log_ruwe>':-2,
    'vsini<':50,
    'logRHK<':-4.0,
    'logg>':0.0,
    'gmag<':7.2,
}

cutoff_sam2 = {
    'NIGHTS>':240,
    'Teff_phot<':5900,
    'Teff_phot>':4500,
    'Fe/H>':-2.0,
    'parallax>':0,
    'log_ruwe>':0.041,
    'vsini<':10,
    'logRHK<':-4.0,
    'logg>':4.0,
    'gmag<':7.2,
}

cutoff_mick = {
    'eff_nights_1.5>':170,
    'teff_mean<':5900,
    'Fe/H>':-0.6,
    'parallax>':0,
    'log_ruwe<':np.round(np.log10(1.2),3),
    'vsini<':5,
    'logRHK<':-4.6,
    'logg>':4.2,
    'vmag<':7.3,
    'HZ_mp_min_osc+gr_texp15<':10,
}

cutoff_william1 = {
    'eff_nights_1.5>':180,
    'Teff_phot<':5980,
    'Teff_phot>':3865,
    'Teff_spec<':5980,
    'Teff_spec>':3865,
    'gmag<':6,
    }

cutoff_william2 = {
    'Teff_phot<':5980,
    'Teff_phot>':3865,
    'Teff_spec<':5980,
    'Teff_spec>':3865,
    'eff_nights_1.5>':180,
    'NIGHTS>':240,
    'Fe/H_spec>':-0.5,
    'parallax>':0,
    'log_ruwe<':np.round(np.log10(1.5),3),
    'logRHK<':-4.6,
    'logg>':0,
    'vsini_spec<':5,
    'gmag<':6.5,
    }

cutoff_stefano = {
 'dec_j2000>':-20,
 'teff_mean<':6000,
 'logg>':4.0,
 'vsini<':5,
 'Fe/H>':-0.4,
 'log_ruwe<':np.round(np.log10(1.2),3),
 'eff_nights_1.5>':160,
 'NIGHTS>':240,
 'HZ_period_inf>':250,
 'gmag<':7.5,
 'HZ_mp_min_osc+gr_texp15>':0.0,
 }



