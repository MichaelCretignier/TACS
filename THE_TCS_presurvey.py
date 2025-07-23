import matplotlib.pylab as plt
import numpy as np

import THE_TCS_classes as tcsc

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

tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_phot', selection=tutorial.info_TA_stars_selected['start'].data)
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_osc', selection=tutorial.info_TA_stars_selected['start'].data)
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc', selection=tutorial.info_TA_stars_selected['start'].data)
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc+gr', selection=tutorial.info_TA_stars_selected['start'].data)
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.55, budget='_arve_phot+osc+gr', selection=tutorial.info_TA_stars_selected['start'].data)

tutorial.compute_optimal_texp(snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc', texp_crit=20, selection=tutorial.info_TA_stars_selected['start'].data)
tutorial.compute_optimal_texp(snr_crit=250, sig_rv_crit=0.55, budget='_arve_phot+osc+gr', texp_crit=20, selection=tutorial.info_TA_stars_selected['start'].data)

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

#toi4499
table = pd.DataFrame({
    'period':[3.490212,5.538772,145.7191],
    'k':[1.05,0.25,129.292],
    'ecc':[0,0,0.355],
    'peri':[0,0,295.548],
    't0':[2460338.925724,2460335.262371,2455510.079815],
    'mass':[1.0,1.0,1000]})

toi4499 = tcsc.tcs(instrument='HARPS3') 
toi4499.compute_night_length(sun_elevation=-12) #twilight -> [0-6] : civil ; [6-12] : nautical ; [12-18] : astronomical
toi4499.set_star(ra=18,dec=39) # change for a declination ra input (HD166620)
toi4499.plot_night_length()
toi4499.create_timeseries(airmass_max=1.75, nb_year=1, texp=30, weather=False)

toi4499.compute_exoplanet_rv_signal(keplerian_par=table, y0=2025)
toi4499.plot_keplerians(axhline=0)

plt.figure()
toi4499.info_XY_keplerian[1].plot()
selection = toi4499.info_XY_keplerian[1].night_subset(3,random=False)
toi4499.info_XY_keplerian[1].subset.plot()
selection = toi4499.info_XY_keplerian[1].night_subset(3,random=True)
toi4499.info_XY_keplerian[1].subset.plot()

#open question:

# 1
# Why bright stars are mostly binaries ?!
# According to Andres, RUWE could be wrong for stars brighter than mv ~ 5

plt.figure(figsize=(18,8))
plt.subplot(1,2,1)
plt.scatter(tcsc.gr8['snr_C22_texp15'],tcsc.gr8['ruwe'],c=tcsc.gr8['teff_mean'],cmap='jet',vmin=5000,vmax=6000) ; plt.colorbar()
plt.axhline(y=1.2,color='k',ls=':')
plt.yscale('log')
plt.ylabel('RUWE')
plt.xlabel('SNR_continuum')
plt.grid()
plt.subplot(1,2,2)
plt.scatter(tcsc.gr8['vmag'],tcsc.gr8['ruwe'],c=tcsc.gr8['teff_mean'],cmap='jet',vmin=5000,vmax=6000) ; plt.colorbar()
plt.axhline(y=1.2,color='k',ls=':')
plt.yscale('log')
plt.ylabel('RUWE')
plt.xlabel('SNR_continuum')
plt.grid()

