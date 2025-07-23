import matplotlib.pylab as plt

import THE_TCS_classes as tcsc

#standard stars based on HARPN

#first and best standard star
star1 = tcsc.tcs(sun_elevation=-6, starname='HD127334') 
s1 = plt.subplot(2,1,1)
star1.plot_night_length(figure=s1) #peak in April

#second standard star
star2 = tcsc.tcs(sun_elevation=-6, starname='HD144579') 
s2 = plt.subplot(2,1,2)
star2.plot_night_length(figure=s2) #peak in May

#create a timesampling for star1
star1.create_timeseries(airmass_max=1.75, nb_year=1, texp=10, weather=False)
dustbin = star1.info_XY_timestamps.night_subset(obs_per_night=2,random=True,replace=False)

plt.figure()
star1.info_XY_timestamps.plot(ytext=3)
star1.info_XY_timestamps.subset.plot()

#SG calendar
star3 = tcsc.tcs()
presurvey = star3.info_TA_cutoff['presurvey']
cutoff = tcsc.mod_cutoff(presurvey,{'gmag<':7.5,'RHK_known<':4.8,'vsini_known<':8}) # 'known' means we want an existing value in the DB

#compute the sky night length over the year
star3.compute_SG_calendar(
    sun_elevation = -6, 
    airmass_max = 1.75, 
    alpha_step = 0.5, 
    dec_step = 1,
    cutoff = cutoff)

star3.compute_SG_month(month=1, plot=False) #january
star3.compute_SG_month(month=2, plot=False) #february

plt.figure(figsize=(10,10)) ; star3.info_TA_stars_selected['SG'].plot('vmag','night_length_Jan',print_names=True)
plt.figure(figsize=(10,10)) ; star3.info_TA_stars_selected['SG'].plot('vmag','night_length_Feb',print_names=True)


