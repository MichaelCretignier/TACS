import matplotlib.pylab as plt

import THE_TCS_classes as tcsc
import THE_TCS_variables as tcsv

# standard stars based on NEID
# HD86728, HD68017, HD51419, HD127334
standards = tcsc.inner_gr8(tcsv.NEID_standards)[1:] #remove HD4628 that is not quiet
star1 = tcsc.tcs(sun_elevation=-6) 
plt.figure(figsize=(16,8))
for n,s in enumerate(standards): 
    star1.set_star(starname=s,verbose=False)
    s1 = plt.subplot(len(standards)/3,3,n+1)
    plt.title('%s'%(star1.info_SC_starname['HD']))
    star1.plot_night_length(figure=s1,legend=False) #peak in April
    plt.ylim(-1,10)
plt.subplots_adjust(hspace=0.45,top=0.95,bottom=0.10)


#create a timesampling for star1
star1.create_timeseries(airmass_max=1.75, nb_year=1, texp=10, weather=False)
dustbin = star1.info_XY_timestamps.night_subset(obs_per_night=3,random=False,replace=False)

plt.figure()
star1.info_XY_timestamps.plot()
star1.info_XY_timestamps.subset.plot()

#SG calendar
star3 = tcsc.tcs()
presurvey = star3.info_TA_cutoff['presurvey']

#we can add or modify cutoff selection if needed
cutoff = tcsc.mod_cutoff(presurvey,{'gmag<':7.5,'logRHK_known<':-4.8,'vsini_known<':5}) # 'known' means we want an existing value in the DB

#compute the sky night length over the year
star3.compute_SG_calendar(
    sun_elevation = -6, 
    airmass_max = 1.75, 
    alpha_step = 0.5, 
    dec_step = 1,
    cutoff = cutoff)

star3.compute_SG_month(month=1, plot=False, selection='SG') #january
star3.info_TA_stars_selected['SG'].plot('vmag','night_length_Jan',print_names=True)
star3.info_TA_stars_selected['minimal'].plot('vmag','night_length_Jan',print_names=False,GUI=False,alpha=0.2)

star3.compute_SG_month(month=2, plot=False, selection='SG') #february
star3.info_TA_stars_selected['SG'].plot('vmag','night_length_Feb',print_names=True)
star3.info_TA_stars_selected['minimal'].plot('vmag','night_length_Feb',print_names=False,GUI=False,alpha=0.2)

# you can also start with your own hardcoded list of stars
star4 = tcsc.tcs(sun_elevation=-6, starname='HD55575') 

starnames = ['HD55575','HD89269','HD56124','HD90839','HD95128']
star4.create_star_selection(starnames,tagname='my_selection')

star4.compute_SG_calendar(
    sun_elevation = -6, 
    airmass_max = 1.75, 
    alpha_step = 0.5, 
    dec_step = 1,
    selection='my_selection')

star4.compute_SG_month(month=1, plot=False, selection='my_selection') #january
star4.info_TA_stars_selected['my_selection'].plot('vmag','night_length_Jan',print_names=True)

# the code also work for other instruments and not only HARPS3
for n,ins in enumerate(['HARPS3','NEID','KPF','EXPRES']):
    star5 = tcsc.tcs(sun_elevation=-6, starname='HD55575',instrument=ins) 
    #plt.subplot(2,2,1+n) ; star5.compute_nights(airmass_max=11, weather=False, plot=True)
    star5.create_timeseries(airmass_max=1.75, nb_year=1, month=1, texp=10, weather=False)
    star5.info_XY_timestamps.plot(label=ins)
plt.xlim(0,30)
plt.legend()

