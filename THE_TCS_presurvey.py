import matplotlib.pylab as plt
import numpy as np
import pandas as pd

import THE_TCS_classes as tcsc
import THE_TCS_variables as tcsv

#let's select -12 for the twilght (in best case, extratime for very bright targets)

tutorial = tcsc.tcs(sun_elevation=-12) 

#let's start with rough cutoff physically motivated

tutorial.func_cutoff(tagname='start',cutoff={
    'teff_mean<': 6000,
    'logg>': 4.2,
    'vsini<': 8,
    'Fe/H>': -0.4,
    'log_ruwe<':0.079,
    'HJ<': 0.5,
    'BDW<': 0.5,
    'RHK<': -4.7,
    'gmag<':7.5,
    'HZ_mp_min_osc+gr_texp15>':0
    }) # produce 250 stars in the pre-survey

# let's compute the average season length of this sample

min_obs_per_year = np.mean(tutorial.info_TA_stars_selected['start'].data['season_length_1.5']) 

# let's assume we want at least 1 observation 1 night over 2

print(min_obs_per_year*0.5) #this is equal to 120 measurement per year

# let's have a look on what is the exposure time to get 120 measurement for the final sample (~40 stars)

tutorial.plot_survey_stars(Nb_star=40) 

#to get 120 measurement per year (1 night over 2), texp_max = 20 minutes
#to get 240 measurement per year (every night), texp_max = 10 minutes

tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_phot', selection='start')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_osc', selection='start')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc', selection='start')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc+gr', selection='start')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.55, budget='_arve_phot+osc+gr', selection='start')

tutorial.compute_optimal_texp(snr=250, sig_rv=0.30, budget='_arve_phot+osc', texp_crit=20, selection='start')
tutorial.compute_optimal_texp(snr=250, sig_rv=0.55, budget='_arve_phot+osc+gr', texp_crit=20, selection='start')

tutorial.func_cutoff(cutoff=tcsc.mod_cutoff(tutorial.info_TA_cutoff['start'],
    {'snr_C22_texp15>':250, 
    'sig_rv_phot_texp20<':0.30}),
    par_space='ra_j2000&dec_j2000',par_crit='HWO==1')

tutorial.func_cutoff(tagname='final',
    cutoff=tcsc.mod_cutoff(tutorial.info_TA_cutoff['start'],
    {'snr_C22_texp15>':250, 
    'sig_rv_phot_texp20<':0.30,
    'season_length_1.75>':240,
    'HZ_mp_min_osc+gr_texp15<':16}),
    par_space='ra_j2000&dec_j2000',par_crit='HWO==1')

tutorial.func_cutoff(
    cutoff=tutorial.info_TA_cutoff['final'],
    par_space='teff_mean&snr_C22_texp15')

tutorial.plot_survey_stars(Nb_star=76) 

#### OPEN QUESTIONS:

# 1
# Why bright stars are mostly binaries ?!
# According to Andres, RUWE could be wrong for stars brighter than mv < 5

plt.figure(figsize=(18,8))
plt.subplot(1,2,1)
plt.scatter(tcsc.gr8_raw['snr_C22_texp15'],tcsc.gr8_raw['ruwe'],c=tcsc.gr8_raw['teff_mean'],cmap='jet',vmin=5000,vmax=6000) ; plt.colorbar()
plt.axhline(y=1.2,color='k',ls=':')
plt.yscale('log')
plt.ylabel('RUWE')
plt.xlabel('SNR_continuum')
plt.grid()
plt.subplot(1,2,2)
plt.scatter(tcsc.gr8_raw['vmag'],tcsc.gr8_raw['ruwe'],c=tcsc.gr8_raw['teff_mean'],cmap='jet',vmin=5000,vmax=6000) ; plt.colorbar()
plt.axhline(y=1.2,color='k',ls=':')
plt.yscale('log')
plt.ylabel('RUWE')
plt.xlabel('SNR_continuum')
plt.grid()

# 2 
# How long would it take to observe all the pre-survey to extract homogeneous RHK, vsini?

tutorial = tcsc.tcs(sun_elevation=-6) 
tutorial.func_cutoff(tagname='start',cutoff={
    'teff_mean<': 6000,
    'logg>': 4.2,
    'vsini<': 8,
    'Fe/H>': -0.4,
    'log_ruwe<':0.079,
    'HJ<': 0.5,
    'BDW<': 0.5,
    'RHK<': -4.7,
    'gmag<':7.5,
    'HZ_mp_min_osc+gr_texp15>':0
    }) # produce 250 stars in the pre-survey

tutorial.plot_survey_stars(Texp=10)

# by assuming a PRE-pre-survey of ~250 stars, it takes 365/25 = ~15 days
# in two weeks, we could have homogeneous Atmos + LOGRHK + vsini
# confirmation by a more mathematical compoutation

tutorial.compute_nb_nights_required(selection='start', texp=10,month=1)
tutorial.compute_optimal_texp(snr_crit=250, sig_rv_crit=0.00, budget='_phot', texp_crit=20, selection='start')
tutorial.compute_nb_nights_required(selection='start',texp='optimal',month=1)

# 3 
# Checking of the overlap between the selection of the TaCS members

question3 = tcsc.tcs()
question3.cutoff_ST()
stars = []
for members in list(question3.info_TA_stars_selected.keys())[1:]:
    loc = np.array(question3.info_TA_stars_selected[members].data.index)
    stars.append(loc)
stars = pd.DataFrame({'index':np.hstack(stars)})['index'].value_counts()
statistic = np.array(stars)
plt.close()

plt.figure()
plt.plot(statistic)
plt.xlabel('Nb stars',fontsize=15)
plt.ylabel('Nb of times selected',fontsize=15)
plt.ylim(0,9)
plt.xlim(0,300)
for j in range(1,9):
    loc = np.where(statistic==j)[0][-1]
    plt.scatter(loc,j,color='k')
    plt.text(loc,j,'%.0f'%(loc),ha='center',va='bottom')

