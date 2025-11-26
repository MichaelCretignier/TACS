import matplotlib.pylab as plt

import THE_TCS_classes as tcsc
import THE_TCS_variables as tcsv


#USEFUL FUNCTIONS
starnames = tcsc.get_info_starname('HD16160')
tcsc.get_info_prot('HD16160')
tcsc.get_info_binary('HD16160')
tcsc.plot_binary('HD16160',seeing=0.75,source='COMPOSITE')

#### PRESURVEY CUTOFFF ####

# Twilight -> [0-6] : civil ; [6-12] : nautical ; [12-18] : astronomical

presurvey = tcsc.tcs(version='5.0') #last catalog = 5.0 (default value) old = [1.0,2.0]
presurvey.print_sp_stat()

#these lines are already run by default in tcsc.tsc()
presurvey.func_cutoff(cutoff=tcsv.cutoff_RVopti, tagname='RVopti') 
presurvey.func_cutoff(cutoff=tcsv.cutoff_megan, tagname='solartwins', protection=False) 
samples = presurvey.union('RVopti','solartwins',union_name='presurvey',Xmarker={'under_review>':0.5},ordering='HD')

#plot sample

plt.figure(figsize=(12,9))
plt.subplot(2,2,1) ; presurvey.info_TA_stars_selected['GR8'].plot_space_mission(newfig=False) ; plt.title('GR8 (%.0f)'%(len(presurvey.info_TA_stars_selected['GR8'].data)))
plt.subplot(2,2,2) ; presurvey.info_TA_stars_selected['RVopti'].plot_space_mission(newfig=False) ; plt.title('RV opti (%.0f)'%(len(presurvey.info_TA_stars_selected['RVopti'].data)))
plt.subplot(2,2,3) ; presurvey.info_TA_stars_selected['solartwins'].plot_space_mission(newfig=False) ; plt.title('Solar-twins (%.0f)'%(len(presurvey.info_TA_stars_selected['solartwins'].data)))
plt.subplot(2,2,4) ; presurvey.info_TA_stars_selected['presurvey'].plot_space_mission(newfig=False) ; plt.title('Union (%.0f)'%(len(presurvey.info_TA_stars_selected['presurvey'].data)))
plt.subplots_adjust(left=0.08,right=0.96,top=0.95,bottom=0.09)

#testing if a star is in a list (otherwise why not)
presurvey.which_cutoff('51Peg', tagname='solartwins')
presurvey.which_cutoff('51Peg', tagname='RVopti')
presurvey.which_cutoff('HD219134', tagname='RVopti') 
presurvey.which_cutoff('HD22049', tagname='RVopti')

presurvey.which_cutoff(['HD166620','HD16160','51Peg','61CygB'], tagname='RVopti')
presurvey.which_cutoff(presurvey.info_TA_stars_selected['minimal'].data.sort_values(by='vmag')['HD'][0:20], tagname='RVopti',plot=True)

presurvey.which_cutoff(tcsv.catalog_NEID['HD'], tagname='RVopti', plot=True)
presurvey.which_cutoff(tcsv.catalog_2ES['GAIA'], tagname='RVopti',plot=True)

#following the K sample
presurvey.func_cutoff(cutoff=tcsv.cutoff_RVopti, show_sample='K', tagname='dustbin')

#### EXAMPLE OF CUTOFF VISUALISATION ####

example1 = tcsc.tcs()
#following a binary flag
example1.func_cutoff(par_space='ra_j2000&dec_j2000', par_crit='HWO==1', cutoff=tcsv.cutoff_RVopti, tagname='dustbin')
example1.func_cutoff(par_space='ra_j2000&dec_j2000', par_crit='PLATO==1', cutoff=tcsv.cutoff_RVopti, tagname='dustbin')

#following a parameter space box
example1.func_cutoff(par_space='teff&distance', par_box=['4500->5300','0->30'], cutoff=tcsv.cutoff_RVopti, tagname='dustbin')

#### VISUALIZATION OF KNOWN EXOPLANETS #####

summary = tcsc.tcs() 
summary1 = tcsc.plot_exoplanets2(summary.info_TA_stars_selected['GR8'].data)
summary2 = tcsc.plot_exoplanets2(summary.info_TA_stars_selected['presurvey'].data)

### ------- EXAMPLES ------- ###

# Twilight -> [0-6] : civil ; [6-12] : nautical ; [12-18] : astronomical
survey = tcsc.tcs(sun_elevation=-12) #HARPS3 is the default
survey.func_cutoff(tagname='bright!',cutoff={'gmag<':5.5,'teff<':6000})
print(survey.info_TA_stars_selected['bright!'].data)

##### COMPUTE SEASON AND NIGHT LENGTH #####

star = tcsc.tcs(sun_elevation=-18, instrument='HARPS3') #HARPS3 is the default
star.set_star(starname='HD146233')

plt.figure(figsize=(12,8))

plt.subplot(1,2,1) ; star.compute_nights(airmass_max=11, weather=False, plot=True)
plt.subplot(1,2,2) ; star.compute_nights(airmass_max=1.5, weather=False, plot=True)

#night duration
star.plot_night_length()

#other instruments
plt.figure(figsize=(18,5))
for n,ins in enumerate(['HARPS3','HARPS','NEID','ESPRESSO','KPF','SOPHIE']):
    star = tcsc.tcs(sun_elevation=-12, instrument=ins)
    star.set_star(ra=8,dec=5) # change for a DEC vs RA input
    plt.subplot(1,6,n+1) ; star.compute_nights(airmass_max=1.5, weather=False, plot=True) ; plt.title(ins)
plt.subplots_adjust(left=0.05,right=0.96)

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
star3 = tcsc.tcs(sun_elevation=-6, instrument='HARPS3')
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

selection = 'presurvey'

tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_phot', selection=selection)
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_osc', selection=selection)
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc', selection=selection)
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc+gr', selection=selection)

tutorial.compute_optimal_texp(snr=200, sig_rv=0.30, budget='_arve_phot+osc', texp_crit=15, selection=selection)
tutorial.compute_optimal_texp(snr=250, sig_rv=0.30, budget='_arve_phot+osc', texp_crit=15, selection=selection)

tutorial.plot_survey_stars(Texp=15,selection=selection,color='green') 
tutorial.plot_survey_stars(Texp=None,selection=selection,ranking='HZ_mp_min_osc+gr_texp15',color='C1') 
tutorial.plot_survey_stars(Texp=None,selection=selection,ranking='texp_optimal',color='C1') 

tutorial.create_table_scheduler(
    selection=selection,
    year=2026,
    texp=900,
    n_obs=65,
    ranking=None,
    month_obs_baseline=3,
    )

tutorial.create_table_scheduler(
    selection=tutorial.info_TA_stars_selected[selection].data.sort_values(by='HZ_mp_min_osc+gr_texp15')[0:40],
    year=2026,
    texp=660,
    freq_obs=1,
    ranking=None,
    month_obs_baseline=12,
    tagname='_40stars'
    )

tutorial.info_TA_stars_selected[selection].plot(y='dec_j2000',x='ra_j2000')

#TESS LIGHTCURVES
star4 = tcsc.tcs(sun_elevation=-12, starname='HD4628')
plt.figure()
star4.show_lightcurve(rm_gap=True)

#Investigate a pre-determined star list (cross-matched with GR8)

neid = tcsc.tcs(sun_elevation=-12)
neid.create_star_selection(tcsv.catalog_NEID['HD'],tagname='NEID')
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
