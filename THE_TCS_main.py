import matplotlib.pylab as plt

import THE_TCS_classes as tcsc

#### VISUALIZATION OF KNOWN EXOPLANETS #####

summary = tcsc.plot_exoplanets2(cutoff={})
summary = tcsc.plot_exoplanets2(cutoff={'MIST Teff<':6000,'ruwe<':1.2,'MIST logg>':4.0})

# START
tutorial = tcsc.tcs() 
#twilight -> [0-6] : civil ; [6-12] : nautical ; [12-18] : astronomical
tutorial.compute_night_length(sun_elevation=-12)
tutorial.plot_survey_stars(Texp=15)
tutorial.plot_survey_stars(Texp=9)
tutorial.plot_survey_stars(Nb_star=40)

### ------- EXAMPLES ------- ###

##### COMPUTE SEASON AND NIGHT LENGTH #####
star = tcsc.tcs(sun_elevation=-12, starname='HD19019') 

#analysis
plt.figure(figsize=(12,12))

plt.subplot(2,2,1) ; star.compute_nights(airmass_max=11, weather=False, plot=True)
plt.subplot(2,2,2) ; star.compute_nights(airmass_max=1.5, weather=False, plot=True)

star.set_star(ra=18,dec=38) # change for a declination ra input (HD166620)

plt.subplot(2,2,3) ; star.compute_nights(airmass_max=1.5, weather=False, plot=True)
plt.subplot(2,2,4) ; star.compute_nights(airmass_max=1.5, weather=True, plot=True)

star.plot_night_length()

##### COMPUTE 10 YEARS TIME-SERIES #####

star2 = tcsc.tcs(sun_elevation=-12, starname='HD217014')
star2.plot_exoplanets_db(y_var='k')

star2.create_timeseries(airmass_max=1.75, nb_year=1, texp=15, weather=False)
star2.compute_exoplanet_rv_signal(y0=2025) #Nov comissioning
star2.plot_keplerians()

star2.set_star(starname='HD75732')
star2.create_timeseries(airmass_max=1.75, nb_year=1, texp=15, weather=False)
star2.compute_exoplanet_rv_signal(y0=2026) #start in 2026
star2.plot_keplerians()


##### SG CALENDAR #####
star3 = tcsc.tcs()
star3.compute_SG_calendar(sun_elevation=-6, airmass_max=1.75, alpha_step=1, dec_step=5)
#change the default cutoff if necessary
star3.compute_SG_calendar(sun_elevation=-6, airmass_max=1.75, alpha_step=1, dec_step=5,
                          cutoff={'ruwe<':1.2,
                                  'MIST Teff<':6000,
                                  'MIST logg>':4.2, 
                                  'Vmag<':5.5,
                                  'vsini<':8})

star3.commpute_SG_month(month=1,
                        cutoff=star3.info_TA_cutoff)

star3.commpute_SG_month(month=1,
                        cutoff={'ruwe<':1.2,
                                  'MIST Teff<':6000,
                                  'MIST logg>':4.2, 
                                  'Vmag<':6.5,
                                  'vsini<':8})


#TESS LIGHTCURVES
star4 = tcsc.tcs(sun_elevation=-12, starname='HD99492')
star4.show_lightcurve(rm_gap=True)

#CUTOFF
star5 = tcsc.tcs(sun_elevation=-12, starname='HD99492')
star5.func_cutoff(cutoff=None)
star5.func_cutoff(cutoff=None,par_space='ra_j2000&dec_j2000',par_crit='HWO==1')
star5.func_cutoff(cutoff=None,par_space='Teff&dist',par_box=['4500->5300','0->30'])


#

#computation of the maximum exposure time

tutorial = tcsc.tcs(sun_elevation=-12) 
tutorial.func_cutoff(cutoff={
    'Teff<': 6000,
    'logg>': 4.2,
    'vsini<': 8,
    'Fe/H>': -0.4,
    'log_ruwe<':0.079,
    'HJ<': 0.5,
    'BDW<': 0.5,
    'RHK<': -4.6})
min_obs_per_year = np.mean(tutorial.info_TA_stars_selected['season_length_1.5'])*0.5 # at least 1 observation 1 night over 2
tutorial.plot_survey_stars(Nb_star=40) #final sample
#to get 120 measurement per year (1 night over 2), texp_max = 20 minutes
#to get 240 measurement per year (every night), texp_max = 10 minutes

tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_phot')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_osc')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_osc+gr')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve')

tutorial.func_cutoff(cutoff={
    'Teff<': 6000,
    'logg>': 4.2,
    'vsini<': 8,
    'Fe/H>': -0.4,
    'log_ruwe<': 0.079,
    'HJ<': 0.5,
    'BDW<': 0.5,
    'RHK<': -4.7,
    'snr_texp20>':250, 
    'sig_rv_phot_texp20<':0.30, 
    },par_space='ra_j2000&dec_j2000',par_crit='HWO==1')

tutorial.func_cutoff(cutoff={
    'Teff<': 6000,
    'logg>': 4.2,
    'vsini<': 8,
    'Fe/H>': -0.4,
    'log_ruwe<': 0.079,
    'HJ<': 0.5,
    'BDW<': 0.5,
    'RHK<': -4.7,
    'snr_texp20>':250, 
    'sig_rv_phot_texp20<':0.30, 
    'season_length_1.75>':240,
    'vmag<':7.5,
    },par_space='ra_j2000&dec_j2000',par_crit='HWO==1')

tutorial.plot_survey_stars(Nb_star=113) 

tutorial.plot_survey_snr_texp(tutorial.info_TA_stars_selected, snr_crit=0, sig_rv_crit=1.0, texp=10)

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


plt.scatter(toi4499.info_XY_keplerian[1].x,toi4499.info_XY_keplerian[1].y,c=toi4499.info_XY_keplerian[1].x.astype('int')%20,cmap='tab20')
