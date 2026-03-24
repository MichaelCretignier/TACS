import matplotlib.pylab as plt
import numpy as np
import pandas as pd

import tacs

#let's select -12 for the twilght (in best case, extratime for very bright targets)

tutorial = tacs.tcs(sun_elevation=-12) 
star = 'HD146233' # try also with HD111398
tacs.plot_rv_texp(star,budget='osc',use_vsini=False)

# let's compute the average season length of the presurvey

min_obs_per_year = np.mean(tutorial.info_TA_stars_selected['presurvey'].data['season_length_1.5']) 

# let's assume we want at least 1 observation 1 night over 2

print(min_obs_per_year*0.5) #this is equal to 130 measurement per year

# let's have a look on what is the exposure time to get 130 measurement for the final sample (~40 stars)

tutorial.plot_survey_stars(Nb_star=40) 

#to get 130 measurements per year (1 night over 2), texp_max = 18 minutes
#to get 260 measurements per year (every night), texp_max = 9 minutes

tutorial.plot_survey_snr_texp(texp=15, snr_crit=250, sig_rv_crit=0.30, budget='_phot', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=15, snr_crit=250, sig_rv_crit=0.30, budget='_arve_osc', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=15, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=15, snr_crit=250, sig_rv_crit=0.30, budget='_arve_phot+osc+gr', selection='presurvey')
tutorial.plot_survey_snr_texp(texp=15, snr_crit=250, sig_rv_crit=0.55, budget='_arve_phot+osc+gr', selection='presurvey')

tutorial.compute_optimal_texp(snr=1, sig_rv=0.30, budget='_extempo_phot+osc', texp_crit=50, selection='presurvey', use_vsini=True)

tutorial.compute_optimal_texp(snr=250, sig_rv=0.30, budget='_extempo_phot+osc', texp_crit=50, selection='presurvey', use_vsini=True)

tutorial.compute_optimal_texp(snr=200, sig_rv=0.30, budget='_extempo_phot+osc', texp_crit=50, selection='presurvey', use_vsini=True)
tutorial.compare_obs_strategy('presurvey',budget='_arve_phot+osc',color='C0', figname='texp')

tutorial.compute_optimal_texp(snr=200, sig_rv=0.30, budget='_arve_phot+osc', texp_crit=50, selection='presurvey')
tutorial.compare_obs_strategy('presurvey',budget='_arve_phot+osc',color='C1', figname='texp')

tutorial = tacs.tcs(sun_elevation=-12) 
tutorial.compute_optimal_texp(snr=200, sig_rv=0.30, budget='_arve_phot+osc', texp_crit=50, selection='presurvey')

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
nmax = int(60*tutorial.info_SC_nb_hours_per_yr/(np.mean(tutorial.info_TA_stars_selected['presurvey'].data['texp_optimal'])+1)/len(tutorial.info_TA_stars_selected['presurvey'].data))
print(nmax)

tutorial.create_table_scheduler(
    selection='presurvey',
    year=2026,
    texp='optimal',
    n_obs=nmax,
    ranking=None,
    month_obs_baseline=5,
    )


tutorial.func_cutoff(cutoff=tacs.mod_cutoff(tutorial.info_TA_cutoff['RVopti'],
    {'snr_C22_texp15>':250, 
    'sig_rv_phot_texp20<':0.30}),
    par_space='ra_j2000&dec_j2000',par_crit='HWO!=0')

tutorial.func_cutoff(tagname='final',
    cutoff=tacs.mod_cutoff(tutorial.info_TA_cutoff['RVopti'],
    {'snr_C22_texp15>':250, 
    'sig_rv_phot_texp20<':0.30,
    'season_length_1.75>':240,
    'HZ_mp_min_osc+gr_texp15<':16}),
    par_space='ra_j2000&dec_j2000',par_crit='HWO!=0')

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
plt.scatter(tacs.gr8_raw['5.2']['snr_C22_texp15'],tacs.gr8_raw['5.2']['ruwe_GAIA'],c=tacs.gr8_raw['5.2']['teff'],cmap='jet',vmin=5000,vmax=6000) ; plt.colorbar()
plt.axhline(y=1.2,color='k',ls=':')
plt.yscale('log')
plt.ylabel('RUWE')
plt.xlabel('SNR_continuum')
plt.grid()
plt.subplot(1,2,2)
plt.scatter(tacs.gr8_raw['5.2']['vmag'],tacs.gr8_raw['5.2']['ruwe_GAIA'],c=tacs.gr8_raw['5.2']['teff'],cmap='jet',vmin=5000,vmax=6000) ; plt.colorbar()
plt.axhline(y=1.2,color='k',ls=':')
plt.yscale('log')
plt.ylabel('RUWE')
plt.xlabel('SNR_continuum')
plt.grid()

# 2 
# How long would it take to observe all the pre-survey to extract homogeneous RHK, vsini?

tutorial = tacs.tcs(sun_elevation=-6) 

tutorial.plot_survey_stars(Texp=10)

# by assuming a PRE-pre-survey of ~150 stars, it takes 365/25 = ~15 days
# in two weeks, we could have homogeneous Atmos + LOGRHK + vsini
# confirmation by a more mathematical compoutation

tutorial.compute_nb_nights_required(selection='presurvey', texp=10,month=1)
tutorial.compute_optimal_texp(snr=250, sig_rv=0.00, budget='_phot', texp_crit=20, selection='presurvey',use_vsini=True)
tutorial.compute_nb_nights_required(selection='presurvey',texp='optimal',month=1)


#5 some peculiar population to check:

import os

from PIL import ImageGrab

tutorial = tacs.tcs() 

#TESS SAMPLE
#2 over 7 rejected
cutoff = tacs.mod_cutoff(tutorial.info_TA_cutoff['wide'],{'TESS>':0.5})
dust = tutorial.func_cutoff(tagname='TESS',cutoff=cutoff,protection=False)
plt.close('cumulative')
tess = np.array(tutorial.info_TA_stars_selected['TESS'].data['HD'])

#BRIGHT SAMPLE
#20 over 28 rejected
cutoff = tacs.mod_cutoff(tutorial.info_TA_cutoff['wide'],{'gmag<':5.5})
dust = tutorial.func_cutoff(tagname='bright', cutoff=cutoff, protection=False)
plt.close('cumulative')
bright = np.array(tutorial.info_TA_stars_selected['bright'].data['HD'])

#HWO SAMPLE
#33 over 46 rejected
cutoff = tacs.mod_cutoff(tutorial.info_TA_cutoff['wide'],{'HWO>':0.5})
dust = tutorial.func_cutoff(tagname='HWO',cutoff=cutoff,protection=False)
plt.close('cumulative')
hwo = np.array(tutorial.info_TA_stars_selected['HWO'].data['HD'])

#HIGH DB measurement
#13 over 26 rejected
cutoff = tacs.mod_cutoff(tutorial.info_TA_cutoff['wide'],{'nobs_DB>':200})
dust = tutorial.func_cutoff(tagname='HDB',cutoff=cutoff,protection=False)
plt.close('cumulative')
hdb = np.array(tutorial.info_TA_stars_selected['HDB'].data['HD'])

#LOW DB measurement
#14 over 20 rejected
cutoff = tacs.mod_cutoff(tutorial.info_TA_cutoff['wide'],{'nobs_DB<':1,'logRHK_known<':4.8})
dust = tutorial.func_cutoff(tagname='LDB',cutoff=cutoff,protection=False)
plt.close('cumulative')
ldb = tutorial.info_TA_stars_selected['LDB'].data.sort_values(by=['HZ_mp_min_osc+gr_texp15'])
ldb = np.array(ldb['HD'])[0:20]
ldb = ldb[ldb!='-']

#SG January-February stars

cutoff = tacs.mod_cutoff(tutorial.info_TA_cutoff['wide'],{'logRHK_known<':-4.8,'vsini_known<':5}) # 'known' means we want an existing value in the DB
dust = tutorial.func_cutoff(tagname='SG',cutoff=cutoff,protection=False)

sg = tutorial.info_TA_stars_selected['SG'].data
sg = sg.sort_values(by=['SG_NGT_len'],ascending=False)[['nobs_DB','HD','SPclass','SG_NGT_len','vmag']][0:30]
sg = np.array(sg['HD'])
tacs.which_cutoff(sg, tutorial.info_TA_cutoff['RVopti'])
for h in np.sort(sg):
    os.system('cls' if os.name == 'nt' else 'clear')
    print('======'*12)
    tacs.which_cutoff(h, tutorial.info_TA_cutoff['RVopti'],display=['nobs_DB','prot','pmag','SG_NGT_len','Rank_THE'])
    tutorial.print_sp_stat(tacs.gr8['2.0'].loc[tacs.gr8['2.0']['HD']==h,'SPclass'].values[0])
    if len(tutorial.info_TA_stars_missing)>0:
        bbox = (970*2, 100*2, 1630*2, 930*2)  # adjust coordinates
        screenshot = ImageGrab.grab(bbox)
        screenshot.save('/Users/cretignier/Documents/THE/TCS/STARS_TO_CHECK/%s.png'%(h))

catalog_version = tacs.last_catalog #last catalog = 5.2 
presurvey = tacs.tcs(version=catalog_version)
table = presurvey.info_TA_stars_selected['presurvey'].data
cutoff = presurvey.info_TA_cutoff['RVopti']
gr8 = presurvey.info_TA_stars_selected['GR8'].data
os.system('rm -f /Users/cretignier/Documents/THE/figures/All_summary/*.png')
os.system('rm -f /Users/cretignier/Documents/THE/figures/All_summary_PRIVATE/*.png')
for index in gr8.index:
    hd = gr8.loc[index,'HD']
    tacs.plot_summary(index, show_private=False,selection=table,cutoff=cutoff)
    plt.savefig('/Users/cretignier/Documents/THE/figures/All_summary/THE%s_%s.png'%(str(index).zfill(4),hd))
    plt.close('all')
    tacs.plot_summary(index, show_private=True,selection=table,cutoff=cutoff)
    plt.savefig('/Users/cretignier/Documents/THE/figures/All_summary_PRIVATE/THE%s_%s.png'%(str(index).zfill(4),hd))
    plt.close('all')

catalog_version = tacs.last_catalog #last catalog = 5.2 
presurvey = tacs.tcs(version=catalog_version)
table = presurvey.info_TA_stars_selected['presurvey'].data
cutoff = presurvey.info_TA_cutoff['RVopti']
gr8 = presurvey.info_TA_stars_selected['GR8'].data
os.system('rm -f /Users/cretignier/Documents/THE/figures/Presurvey_summary/*.png')
os.system('rm -f /Users/cretignier/Documents/THE/figures/Presurvey_summary_PRIVATE/*.png')
for index in table.index:
    hd = gr8.loc[index,'HD']
    tacs.plot_summary(hd, show_private=False,selection=table,cutoff=cutoff)
    plt.savefig('/Users/cretignier/Documents/THE/figures/Presurvey_summary/%s_THE%s.png'%(hd,str(index).zfill(4)))
    plt.close('all')
    tacs.plot_summary(hd, show_private=True,selection=table,cutoff=cutoff)
    plt.savefig('/Users/cretignier/Documents/THE/figures/Presurvey_summary_PRIVATE/%s_THE%s.png'%(hd,str(index).zfill(4)))
    plt.close('all')

table = presurvey.info_TA_stars_selected['wide'].data
table = table.loc[~np.in1d(table['GAIA'],presurvey.info_TA_stars_selected['presurvey'].data['GAIA'])]
table = table.sort_values(by='vmag')[0:30]
cutoff = presurvey.info_TA_cutoff['RVopti']
gr8 = presurvey.info_TA_stars_selected['GR8'].data
os.system('rm -f /Users/cretignier/Documents/THE/figures/Bright_summary_PRIVATE/*.png')
for index in table.index:
    hd = gr8.loc[index,'HD']
    tacs.plot_summary(index, show_private=True,cutoff=cutoff)
    plt.savefig('/Users/cretignier/Documents/THE/figures/Bright_summary_PRIVATE/%s_THE%s.png'%(hd,str(index).zfill(4)))
    plt.close('all')
    cutoff = presurvey.info_TA_cutoff['RVopti'].copy()


######

#
survey = tacs.tcs(sun_elevation=-12) 

survey.compute_optimal_texp(
    snr = 250, 
    sig_rv = 0.30, 
    budget = '_arve_phot+osc', 
    texp_crit = 50, 
    texp_extra = 2, # +2min
    texp_min = 8,   #  8min
    use_vsini = True,
    selection = 'presurvey')

best = survey.compute_ranking(selection='presurvey', budget='arve_phot+osc+gr', use_vsini=True, texp='optimal')
survey.compare_obs_strategy('presurvey',budget='_arve_phot+osc',color='C1', figname='texp')

config = {
    '100':{'baseline':4,'texp_min':8,'nobs':90},
    '80':{'baseline':5,'texp_min':8,'nobs':130},
    '60':{'baseline':8,'texp_min':8,'nobs':200},
    '40':{'baseline':12,'texp_min':12,'nobs':250},
    }

for N in [100,80,60,40]:
    c = config[str(N)]
    texp_min = c['texp_min']
    baseline = c['baseline']
    n_obs = c['nobs']
    survey.compute_optimal_texp(
        snr=250, 
        sig_rv=0.30, 
        budget='_arve_phot+osc', 
        texp_crit = 25, 
        texp_extra = 2, # +2min
        texp_min = texp_min,   #  8min
        use_vsini=True,
        selection='presurvey')
    survey.compare_obs_strategy('presurvey',budget='_arve_phot+osc',color='C0', figname='texp')

    survey_tab = survey.info_TA_stars_selected['presurvey'].data.copy()
    survey_tab = survey_tab.loc[survey_tab['texp_optimal']<100]
    survey_tab = survey_tab.sort_values(by='gmag')[0:N]
    survey.create_table_scheduler(
        selection=survey_tab,
        year = 2026,
        texp = 'optimal',
        t_slew = 60,
        n_obs = n_obs,
        ranking = None,
        month_obs_baseline = baseline,
        standards = True
        )
    plt.savefig('/Users/cretignier/Documents/Analysis/N%.0f.pdf'%(N))


standards = tacs.inner_gr8(['HD4628','HD146233','HD186408']) 
star1 = tacs.tcs(sun_elevation=-6) 
plt.figure(figsize=(16,8))
s1 = plt.subplot(1,1,1)
for n,s in enumerate(standards): 
    star1.set_star(starname=s,verbose=False)
    star1.plot_night_length(figure=s1,legend=False,airmass_max=[1.5],sun_elevation=[-12]) #peak in April
    plt.ylim(-1,10)
plt.subplots_adjust(hspace=0.45,top=0.95,bottom=0.10)
