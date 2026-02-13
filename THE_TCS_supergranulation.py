import matplotlib.pylab as plt

import numpy as np
import tacs

# standard stars based on NEID
# HD86728, HD68017, HD51419, HD127334
standards = tacs.inner_gr8(tacs.NEID_standards)[1:] #remove HD4628 that is not quiet
star1 = tacs.tcs(sun_elevation=-6) 
plt.figure(figsize=(16,8))
for n,s in enumerate(standards): 
    star1.set_star(starname=s,verbose=False)
    s1 = plt.subplot(int(np.ceil(len(standards)/3)),3,n+1)
    plt.title('%s'%(star1.info_SC_starname['HD']))
    star1.plot_night_length(figure=s1,legend=False,color=None) #peak in April
    plt.ylim(-1,10)
plt.subplots_adjust(hspace=0.45,top=0.95,bottom=0.10)

#standards1 = tacs.inner_gr8(['HD4628','HD69830','HD146233','HD186408']) #remove HD4628 that is not quiet
standards1 = tacs.THE_standards 
star1 = tacs.tcs(sun_elevation=-6) 
plt.figure(figsize=(16,8))
s1 = plt.subplot(1,1,1)
for n,s in enumerate(standards1): 
    star1.set_star(starname=s,verbose=False)
    star1.plot_night_length(figure=s1,legend=False,airmass_max=[1.5],sun_elevation=[-12, -18],color='C%.0f'%(n),showname=True) #peak in April
    plt.ylim(-1,10)
plt.subplots_adjust(hspace=0.45,top=0.95,bottom=0.10)


#create a timesampling for star1
star1.create_timeseries(airmass_max=1.75, nb_year=1, texp=10, weather=False)
dustbin = star1.info_XY_timestamps.night_subset(obs_per_night=3,random=False,replace=False)

plt.figure()
star1.info_XY_timestamps.plot()
star1.info_XY_timestamps.subset.plot()

#compute the sky night length over the year
star3 = tacs.tcs()
star3.compute_SG_calendar(
    sun_elevation = -6, 
    airmass_max = 1.75, 
    alpha_step = 10, 
    dec_step = 1,
    selection = 'presurvey')

star3.compute_SG_month(month=5, plot=True, selection='presurvey')
star3.info_TA_stars_selected['minimal'].plot('vmag','night_length_May',print_names=False,GUI=True,alpha=0.2)
star3.info_TA_stars_selected['presurvey'].plot('vmag','night_length_May',print_names=True,GUI=False)

star3.compute_SG_month(month=6, plot=False, selection='SG')
star3.info_TA_stars_selected['minimal'].plot('vmag','night_length_Jun',print_names=False,GUI=True,alpha=0.2)
star3.info_TA_stars_selected['SG'].plot('vmag','night_length_Jun',print_names=True,GUI=False)

# you can also start with your own hardcoded list of stars
star4 = tacs.tcs(sun_elevation=-6, starname='HD55575') 

starnames = ['HD55575','HD89269','HD56124','HD90839','HD95128']
star4.create_star_selection(starnames,tagname='my_selection')

star4.compute_SG_calendar(
    sun_elevation = -6, 
    airmass_max = 1.75, 
    alpha_step = 10, 
    dec_step = 1,
    selection='my_selection')

star4.compute_SG_month(month=6, plot=True, selection='my_selection') #june
star4.info_TA_stars_selected['my_selection'].plot('vmag','night_length_Jun',print_names=True)

# the code also work for other instruments and not only HARPS3
for n,ins in enumerate(['HARPS3','NEID','KPF','EXPRES','SOPHIE']):
    star5 = tacs.tcs(sun_elevation=-6, starname='HD55575',instrument=ins,verbose=False) 
    #plt.subplot(2,2,1+n) ; star5.compute_nights(airmass_max=11, weather=False, plot=True)
    star5.create_timeseries(airmass_max=1.75, nb_year=1, month=1, texp=10, weather=False)
    star5.info_XY_timestamps.plot(label=ins)
plt.xlim(0,30)
plt.legend()

