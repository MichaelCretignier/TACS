import matplotlib.pylab as plt

import THE_TCS_classes as tcsc
import THE_TCS_variables as tcsv

#### PRESURVEY CUTOFFF ####

presurvey = tcsc.tcs(sun_elevation=-12, instrument='HARPS3') #HARPS3 is the default
presurvey.func_cutoff(cutoff=tcsv.cutoff_presurvey, tagname='presurvey') #this line is already running by default in tcsc.tsc()

#testing if a star is in a list (otherwise why not)
presurvey.which_cutoff('HD166620', tagname='presurvey') #can also take cutoff={} input
presurvey.which_cutoff('HD16160', tagname='presurvey')
presurvey.which_cutoff('51Peg', tagname='presurvey')
 
#following the K sample
presurvey.func_cutoff(cutoff=tcsv.cutoff_presurvey, show_sample='K', tagname='dustbin')

#### EXAMPLE OF CUTOFF VISUALISATION ####

example1 = tcsc.tcs(sun_elevation=-12)
#following a binary flag
example1.func_cutoff(par_space='ra_j2000&dec_j2000', par_crit='HWO==1', cutoff=tcsv.cutoff_presurvey, tagname='dustbin')

#following a parameter space box
example1.func_cutoff(par_space='teff_mean&dist', par_box=['4500->5300','0->30'], cutoff=tcsv.cutoff_presurvey, tagname='dustbin')

#### VISUALIZATION OF KNOWN EXOPLANETS #####

summary = tcsc.plot_exoplanets2(cutoff={})
summary = tcsc.plot_exoplanets2(cutoff={'teff_mean<':6000,'ruwe<':1.2,'logg>':4.2})

### ------- EXAMPLES ------- ###

# Twilight -> [0-6] : civil ; [6-12] : nautical ; [12-18] : astronomical
survey = tcsc.tcs(sun_elevation=-12) #HARPS3 is the default
survey.func_cutoff(tagname='bright!',cutoff={'gmag<':6,'teff_mean<':6000})

##### COMPUTE SEASON AND NIGHT LENGTH #####
star = tcsc.tcs(sun_elevation=-12, starname='HD127334') #using starname

plt.figure(figsize=(12,12))

plt.subplot(1,2,1) ; star.compute_nights(airmass_max=11, weather=False, plot=True)
plt.subplot(1,2,2) ; star.compute_nights(airmass_max=1.5, weather=False, plot=True)

plt.figure(figsize=(18,5))
for n,ins in enumerate(['HARPS3','HARPS','NEID','ESPRESSO','KPF']):
    star = tcsc.tcs(sun_elevation=-12, instrument=ins)
    star.set_star(ra=18,dec=20) # change for a DEC vs RA input
    plt.subplot(1,5,n+1) ; star.compute_nights(airmass_max=1.5, weather=False, plot=True) ; plt.title(ins)
plt.subplots_adjust(left=0.05,right=0.96)

#night duration
star.plot_night_length()

##### COMPUTE 10 YEARS TIME-SERIES #####

star2 = tcsc.tcs(sun_elevation=-12, starname='HD217014')
star2.plot_exoplanets_db(y_var='mass')

star2.create_timeseries(airmass_max=1.75, nb_year=1, texp=15, weather=False)
star2.compute_exoplanet_rv_signal(y0=2025) #Nov. comissioning
star2.plot_keplerians()

star2.set_star(starname='HD75732')
star2.create_timeseries(airmass_max=1.75, nb_year=2, texp=15, weather=True)
star2.compute_exoplanet_rv_signal(y0=2026) #start in 2026
star2.plot_keplerians()

##### SG CALENDAR #####
star3 = tcsc.tcs()
star3.compute_SG_calendar(sun_elevation=-6, airmass_max=1.75, alpha_step=1, dec_step=5)

star3.compute_SG_month(month=1,plot=True)

#OPTIMAL EXPOSURE TIME

# START
tutorial = tcsc.tcs(sun_elevation=-12) 
# As a recall, Total_time = Nstar * (Texp + overhead) * Nb_obs
tutorial.plot_survey_stars(Nb_star=100)

tutorial.plot_survey_stars(Texp=20)
tutorial.plot_survey_stars(Texp=10) # NB, the number are not twice because overhead = 1min
tutorial.plot_survey_stars(Nb_star=40)
tutorial.plot_survey_stars(Nb_obs_per_year=60) # 120 measurements for a 2-year pre-survey

tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_phot', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_osc', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc+gr', selection='presurvey')

tutorial.compute_optimal_texp(snr=150, sig_rv=0.30, budget='_arve_phot+osc', texp_crit=15, selection='presurvey')
tutorial.compute_optimal_texp(snr=200, sig_rv=0.30, budget='_arve_phot+osc', texp_crit=15, selection='presurvey')

tutorial.plot_survey_stars(Texp=15,selection='presurvey',color='green') 
tutorial.plot_survey_stars(Texp=None,selection='presurvey',ranking='HZ_mp_min_osc+gr_texp15',color='C1') 
tutorial.plot_survey_stars(Texp=None,selection='presurvey',ranking='texp_optimal',color='C1') 

tutorial.create_table_scheduler(
    selection='presurvey',
    year=2026,
    texp=1000,
    n_obs=50,
    month_obs_baseline=3
    )

tutorial.create_table_scheduler(
    selection=tutorial.info_TA_stars_selected['presurvey'].data.sort_values(by='HZ_mp_min_osc+gr_texp15')[0:40],
    year=2026,
    texp=700,
    freq_obs=1,
    ranking=None,
    month_obs_baseline=12,
    tagname='40stars'
    )


#TESS LIGHTCURVES
star4 = tcsc.tcs(sun_elevation=-12, starname='HD99492')
star4.show_lightcurve(rm_gap=True)

#Investigate a pre-determined star list (cross-matched with GR8)

neid = tcsc.tcs(sun_elevation=-12)
neid.create_star_selection(tcsv.NEID_catalog['HD'],tagname='NEID')
neid.create_star_selection(tcsv.NEID_standards,tagname='NEID_standards')

neid.info_TA_stars_selected['NEID_standards'].plot(y='dec_j2000',x='ra_j2000')
neid.info_TA_stars_selected['presurvey'].plot(y='dec_j2000',x='ra_j2000',c='k',GUI=False)

neid.compute_SG_calendar(
    sun_elevation = -6, 
    airmass_max = 1.75, 
    alpha_step = 1, 
    dec_step = 5,
    selection='NEID_standards')

neid.compute_SG_month(month=1,plot=True,selection='NEID_standards')

##
test = tcsc.tcs(sun_elevation=-12, instrument='HARPS3') #HARPS3 is the default

for selection in ['GR8','presurvey']:
    rhk1 = test.info_TA_stars_selected[selection].data['logRHK']
    rhk2 = test.info_TA_stars_selected[selection].data['logRHK_BoroSaika+18']
    rhk3 = test.info_TA_stars_selected[selection].data['logRHK_DACE']
    rhk4 = test.info_TA_stars_selected[selection].data['logRHK_YARARA']
    rhk = np.array([rhk1,rhk2,rhk3,rhk4]).T

    plt.figure()
    for n in [1,4]:
        r = np.nanmean(rhk[:,0:n],axis=1)
        myf.hist(r[r==r],bins=np.arange(-6,-4,0.1),label='%.0f%%(%.0f)'%(np.sum(r==r)*100/len(r),np.sum(r==r)),color='C%.0f'%(n-1))
    plt.legend()