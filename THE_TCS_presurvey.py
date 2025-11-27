import matplotlib.pylab as plt
import numpy as np
import pandas as pd

import THE_TCS_classes as tcsc
import THE_TCS_variables as tcsv

#let's select -12 for the twilght (in best case, extratime for very bright targets)

tutorial = tcsc.tcs(sun_elevation=-12) 

# let's compute the average season length of the presurvey

min_obs_per_year = np.mean(tutorial.info_TA_stars_selected['presurvey'].data['season_length_1.5']) 

# let's assume we want at least 1 observation 1 night over 2

print(min_obs_per_year*0.5) #this is equal to 130 measurement per year

# let's have a look on what is the exposure time to get 130 measurement for the final sample (~40 stars)

tutorial.plot_survey_stars(Nb_star=40) 

#to get 130 measurements per year (1 night over 2), texp_max = 18 minutes
#to get 260 measurements per year (every night), texp_max = 9 minutes

tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_phot', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_osc', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc+gr', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=20, snr_crit=250, sig_rv_crit=0.55, budget='_arve_phot+osc+gr', selection='presurvey')

tutorial.compute_optimal_texp(snr=250, sig_rv=0.30, budget='_arve_phot+osc', texp_crit=20, selection='presurvey')
print(int(60*tutorial.info_SC_nb_hours_per_yr_eff/(1+np.mean(tutorial.info_TA_stars_selected['presurvey'].data['texp_optimal']))/len(tutorial.info_TA_stars_selected['presurvey'].data)))

tutorial.create_table_scheduler(
    selection='presurvey',
    year=2026,
    texp='optimal',
    n_obs=90,
    ranking=None,
    month_obs_baseline=3,
    )

tutorial.create_table_scheduler(
    selection='presurvey',
    year=2026,
    texp=687,
    n_obs=90,
    ranking=None,
    month_obs_baseline=3,
    )

tutorial.compute_optimal_texp(snr=200, sig_rv=0.30, budget='_arve_phot+osc', texp_crit=20, selection='presurvey')
print(int(60*tutorial.info_SC_nb_hours_per_yr_eff/(np.mean(tutorial.info_TA_stars_selected['presurvey'].data['texp_optimal'])+1)/len(tutorial.info_TA_stars_selected['presurvey'].data)))

tutorial.create_table_scheduler(
    selection='presurvey',
    year=2026,
    texp='optimal',
    n_obs=132,
    ranking=None,
    month_obs_baseline=5,
    )


tutorial.func_cutoff(cutoff=tcsc.mod_cutoff(tutorial.info_TA_cutoff['presurvey'],
    {'snr_C22_texp15>':250, 
    'sig_rv_phot_texp20<':0.30}),
    par_space='ra_j2000&dec_j2000',par_crit='HWO==1')

tutorial.func_cutoff(tagname='final',
    cutoff=tcsc.mod_cutoff(tutorial.info_TA_cutoff['presurvey'],
    {'snr_C22_texp15>':250, 
    'sig_rv_phot_texp20<':0.30,
    'season_length_1.75>':240,
    'HZ_mp_min_osc+gr_texp15<':16}),
    par_space='ra_j2000&dec_j2000',par_crit='HWO==1')

tutorial.func_cutoff(
    cutoff=tutorial.info_TA_cutoff['final'],
    par_space='teff&snr_C22_texp15')

tutorial.plot_survey_stars(Nb_star=76) 

#### OPEN QUESTIONS:

# 1
# Why bright stars are mostly binaries ?!
# According to Andres, RUWE could be wrong for stars brighter than mv < 5

plt.figure(figsize=(18,8))
plt.subplot(1,2,1)
plt.scatter(tcsc.gr8_raw['snr_C22_texp15'],tcsc.gr8_raw['ruwe'],c=tcsc.gr8_raw['teff'],cmap='jet',vmin=5000,vmax=6000) ; plt.colorbar()
plt.axhline(y=1.2,color='k',ls=':')
plt.yscale('log')
plt.ylabel('RUWE')
plt.xlabel('SNR_continuum')
plt.grid()
plt.subplot(1,2,2)
plt.scatter(tcsc.gr8_raw['vmag'],tcsc.gr8_raw['ruwe'],c=tcsc.gr8_raw['teff'],cmap='jet',vmin=5000,vmax=6000) ; plt.colorbar()
plt.axhline(y=1.2,color='k',ls=':')
plt.yscale('log')
plt.ylabel('RUWE')
plt.xlabel('SNR_continuum')
plt.grid()

# 2 
# How long would it take to observe all the pre-survey to extract homogeneous RHK, vsini?

tutorial = tcsc.tcs(sun_elevation=-6) 

tutorial.plot_survey_stars(Texp=10)

# by assuming a PRE-pre-survey of ~150 stars, it takes 365/25 = ~15 days
# in two weeks, we could have homogeneous Atmos + LOGRHK + vsini
# confirmation by a more mathematical compoutation

tutorial.compute_nb_nights_required(selection='presurvey', texp=10,month=1)
tutorial.compute_optimal_texp(snr=250, sig_rv=0.00, budget='_phot', texp_crit=20, selection='start')
tutorial.compute_nb_nights_required(selection='presurvey',texp='optimal',month=1)

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

#4 

#5 some peculiar population to check:

import os

from PIL import ImageGrab

tutorial = tcsc.tcs() 

#TESS SAMPLE
#2 over 7 rejected
cutoff = tcsc.mod_cutoff(tutorial.info_TA_cutoff['minimal'],{'TESS>':0.5})
dust = tutorial.func_cutoff(tagname='TESS',cutoff=cutoff,protection=False)
plt.close('cumulative')
tess = np.array(tutorial.info_TA_stars_selected['TESS'].data['HD'])
tutorial.which_cutoff(tess, tagname='presurvey')
for t in np.sort(tess):
    os.system('cls' if os.name == 'nt' else 'clear')
    print('======'*12)
    tutorial.which_cutoff(t, tagname='presurvey',display=['nobs_DB','prot','pmag','SG_NGT_len','Rank_THE'])
    if len(tutorial.info_TA_stars_missing)>0:
        bbox = (970*2, 100*2, 1630*2, 930*2)  # adjust coordinates
        screenshot = ImageGrab.grab(bbox)
        screenshot.save('/Users/cretignier/Documents/THE/TCS/STARS_TO_CHECK/%s.png'%(t))

#BRIGHT SAMPLE
#20 over 28 rejected
cutoff = tcsc.mod_cutoff(tutorial.info_TA_cutoff['minimal'],{'gmag<':5.5})
dust = tutorial.func_cutoff(tagname='bright', cutoff=cutoff, protection=False)
plt.close('cumulative')
bright = np.array(tutorial.info_TA_stars_selected['bright'].data['HD'])
tutorial.which_cutoff(bright, tagname='presurvey')
for b in np.sort(bright):
    os.system('cls' if os.name == 'nt' else 'clear')
    print('======'*12)
    tutorial.which_cutoff(b, tagname='presurvey',display=['nobs_DB','prot','pmag','SG_NGT_len','Rank_THE'])
    tutorial.print_sp_stat(tcsc.gr8['2.0'].loc[tcsc.gr8['2.0']['HD']==b,'SPclass'].values[0])
    if len(tutorial.info_TA_stars_missing)>0:
        bbox = (970*2, 100*2, 1630*2, 930*2)  # adjust coordinates
        screenshot = ImageGrab.grab(bbox)
        screenshot.save('/Users/cretignier/Documents/THE/TCS/STARS_TO_CHECK/%s.png'%(b))
 
#HWO SAMPLE
#33 over 46 rejected
cutoff = tcsc.mod_cutoff(tutorial.info_TA_cutoff['minimal'],{'HWO>':0.5})
dust = tutorial.func_cutoff(tagname='HWO',cutoff=cutoff,protection=False)
plt.close('cumulative')
hwo = np.array(tutorial.info_TA_stars_selected['HWO'].data['HD'])
tutorial.which_cutoff(hwo, tagname='presurvey')
for h in np.sort(hwo):
    os.system('cls' if os.name == 'nt' else 'clear')
    print('======'*12)
    tutorial.which_cutoff(h, tagname='presurvey',display=['nobs_DB','prot','pmag','SG_NGT_len','Rank_THE'])
    tutorial.print_sp_stat(tcsc.gr8['2.0'].loc[tcsc.gr8['2.0']['HD']==h,'SPclass'].values[0])
    if len(tutorial.info_TA_stars_missing)>0:
        bbox = (970*2, 100*2, 1630*2, 930*2)  # adjust coordinates
        screenshot = ImageGrab.grab(bbox)
        screenshot.save('/Users/cretignier/Documents/THE/TCS/STARS_TO_CHECK/%s.png'%(h))

#HIGH DB measurement
#13 over 26 rejected
cutoff = tcsc.mod_cutoff(tutorial.info_TA_cutoff['minimal'],{'nobs_DB>':200})
dust = tutorial.func_cutoff(tagname='HDB',cutoff=cutoff,protection=False)
plt.close('cumulative')
hdb = np.array(tutorial.info_TA_stars_selected['HDB'].data['HD'])
tutorial.which_cutoff(hdb, tagname='presurvey')
for h in np.sort(hdb):
    os.system('cls' if os.name == 'nt' else 'clear')
    print('======'*12)
    tutorial.which_cutoff(h, tagname='presurvey',display=['nobs_DB','prot','pmag','SG_NGT_len','Rank_THE'])
    tutorial.print_sp_stat(tcsc.gr8['2.0'].loc[tcsc.gr8['2.0']['HD']==h,'SPclass'].values[0])
    if len(tutorial.info_TA_stars_missing)>0:
        bbox = (970*2, 100*2, 1630*2, 930*2)  # adjust coordinates
        screenshot = ImageGrab.grab(bbox)
        screenshot.save('/Users/cretignier/Documents/THE/TCS/STARS_TO_CHECK/%s.png'%(h))

#LOW DB measurement
#14 over 20 rejected
cutoff = tcsc.mod_cutoff(tutorial.info_TA_cutoff['minimal'],{'nobs_DB<':1,'logRHK_known<':4.8})
dust = tutorial.func_cutoff(tagname='LDB',cutoff=cutoff,protection=False)
plt.close('cumulative')
ldb = tutorial.info_TA_stars_selected['LDB'].data.sort_values(by=['HZ_mp_min_osc+gr_texp15'])
ldb = np.array(ldb['HD'])[0:20]
ldb = ldb[ldb!='-']
tutorial.which_cutoff(ldb, tagname='presurvey')
for h in np.sort(ldb):
    os.system('cls' if os.name == 'nt' else 'clear')
    print('======'*12)
    tutorial.which_cutoff(h, tagname='presurvey',display=['nobs_DB','prot','pmag','SG_NGT_len','Rank_THE'])
    tutorial.print_sp_stat(tcsc.gr8['2.0'].loc[tcsc.gr8['2.0']['HD']==h,'SPclass'].values[0])
    if len(tutorial.info_TA_stars_missing)>0:
        bbox = (970*2, 100*2, 1630*2, 930*2)  # adjust coordinates
        screenshot = ImageGrab.grab(bbox)
        screenshot.save('/Users/cretignier/Documents/THE/TCS/STARS_TO_CHECK/%s.png'%(h))

#SG January-February stars

cutoff = tcsc.mod_cutoff(tutorial.info_TA_cutoff['minimal'],{'logRHK_known<':-4.8,'vsini_known<':5}) # 'known' means we want an existing value in the DB
dust = tutorial.func_cutoff(tagname='SG',cutoff=cutoff,protection=False)

sg = tutorial.info_TA_stars_selected['SG'].data
sg = sg.sort_values(by=['SG_NGT_len'],ascending=False)[['nobs_DB','HD','SPclass','SG_NGT_len','vmag']][0:30]
sg = np.array(sg['HD'])
tutorial.which_cutoff(sg, tagname='presurvey')
for h in np.sort(sg):
    os.system('cls' if os.name == 'nt' else 'clear')
    print('======'*12)
    tutorial.which_cutoff(h, tagname='presurvey',display=['nobs_DB','prot','pmag','SG_NGT_len','Rank_THE'])
    tutorial.print_sp_stat(tcsc.gr8['2.0'].loc[tcsc.gr8['2.0']['HD']==h,'SPclass'].values[0])
    if len(tutorial.info_TA_stars_missing)>0:
        bbox = (970*2, 100*2, 1630*2, 930*2)  # adjust coordinates
        screenshot = ImageGrab.grab(bbox)
        screenshot.save('/Users/cretignier/Documents/THE/TCS/STARS_TO_CHECK/%s.png'%(h))

