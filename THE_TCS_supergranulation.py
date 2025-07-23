import matplotlib.pylab as plt

import THE_TCS_classes as tcsc

#let's select -12 for the twilght (in best case, extratime for very bright targets)

# standard stars based on HARPN

#first and best standard star
star1 = tcsc.tcs(sun_elevation=-6, starname='HD127334') 
s1 = plt.subplot(2,1,1)
star1.plot_night_length(figure=s1)
star1.create_timeseries(airmass_max=1.75, nb_year=1, texp=15, weather=False)

#second standard star
star2 = tcsc.tcs(sun_elevation=-6, starname='HD144579') 
s2 = plt.subplot(2,1,2)
star2.plot_night_length(figure=s2)
star1.create_timeseries(airmass_max=1.75, nb_year=1, texp=15, weather=False)

#SG calender
star3 = tcsc.tcs()
star3.compute_SG_calendar(sun_elevation=-6, airmass_max=1.75, alpha_step=0.5, dec_step=1)
presurvey = star3.info_TA_cutoff['presurvey']
star3.compute_SG_month(month=1,cutoff=tcsc.mod_cutoff(presurvey,{'gmag<':6.5,'RHK_known<':4.8,'vsini_known<':8}))


