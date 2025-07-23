import matplotlib.pylab as plt

import THE_TCS_classes as tcsc

#### VISUALIZATION OF KNOWN EXOPLANETS #####

summary = tcsc.plot_exoplanets2(cutoff={})
summary = tcsc.plot_exoplanets2(cutoff={'teff_mean<':6000,'ruwe<':1.2,'logg>':4.2})

# START
tutorial = tcsc.tcs(sun_elevation=-12) 
# Twilight -> [0-6] : civil ; [6-12] : nautical ; [12-18] : astronomical
# As a recall, Total_time = Nstar * Texp * Nb_obs
tutorial.plot_survey_stars(Texp=20)
tutorial.plot_survey_stars(Texp=10) # NB, the number are not twice because overhead = 1min
tutorial.plot_survey_stars(Nb_star=40)
tutorial.plot_survey_stars(Nb_obs_per_year=75) # 150 measurements for a 2-year pre-survey

### ------- EXAMPLES ------- ###

##### COMPUTE SEASON AND NIGHT LENGTH #####
star = tcsc.tcs(sun_elevation=-12, starname='HD217014') #using starname

#analysis visibility
plt.figure(figsize=(12,12))

plt.subplot(2,2,1) ; star.compute_nights(airmass_max=11, weather=False, plot=True)
plt.subplot(2,2,2) ; star.compute_nights(airmass_max=1.5, weather=False, plot=True)

star.set_star(ra=18,dec=38) # change for a DEC vs RA input (HD166620)

plt.subplot(2,2,3) ; star.compute_nights(airmass_max=1.5, weather=False, plot=True)
plt.subplot(2,2,4) ; star.compute_nights(airmass_max=1.5, weather=True, plot=True)

#night duration
star.plot_night_length()

##### COMPUTE 10 YEARS TIME-SERIES #####

star2 = tcsc.tcs(sun_elevation=-12, starname='HD217014')
star2.plot_exoplanets_db(y_var='k')

star2.create_timeseries(airmass_max=1.75, nb_year=1, texp=15, weather=False)
star2.compute_exoplanet_rv_signal(y0=2025) #Nov. comissioning
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
                                  'teff_mean<':6000,
                                  'logg>':4.2, 
                                  'vmag<':5.5,
                                  'vsini<':8})

star3.compute_SG_month(month=1)

star3.compute_SG_month(month=1,
                        cutoff={'ruwe<':1.2,
                                  'MIST Teff<':6000,
                                  'MIST logg>':4.2, 
                                  'Vmag<':6.5,
                                  'vsini<':8})


#OPTIMAL EXPOSURE TIME

star3.compute_optimal_texp(snr_crit=250, sig_rv_crit=0.30, budget='_phot', texp_crit=15)

#TESS LIGHTCURVES
star4 = tcsc.tcs(sun_elevation=-12, starname='HD99492')
star4.show_lightcurve(rm_gap=True)

#CUTOFF
star5 = tcsc.tcs(sun_elevation=-12)
star5.func_cutoff(cutoff=None)
star5.func_cutoff(cutoff=None,par_space='ra_j2000&dec_j2000',par_crit='HWO==1')
star5.func_cutoff(cutoff=None,par_space='teff_mean&dist',par_box=['4500->5300','0->30'])

star6 = tcsc.tcs(sun_elevation=-6)
star6.cutoff_ST()

