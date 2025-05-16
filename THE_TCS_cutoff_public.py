import os

import matplotlib.pylab as plt
import numpy as np
import pandas as pd

try:
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    module_warning = False
except ModuleNotFoundError:
    module_warning = True

##########################
####### IMPORT TAB #######
##########################

cwd = os.getcwd()

#table = pd.read_csv(cwd+'done_the_isochrones.csv',index_col=0) #no more used
table = pd.read_csv(cwd+'/Master_table.csv',index_col=0)

##########################
####### TAB EXTENT #######
##########################

table['Teff'] = table['Teff_spec']
table.loc[table['Teff']!=table['Teff'],'Teff'] = table.loc[table['Teff']!=table['Teff'],'Teff_gspphot']

table['logg'] = table['logg_spec']
table.loc[table['logg']!=table['logg'],'logg'] = table.loc[table['logg']!=table['logg'],'logg_gspphot']

table['dist'] = 1000/table['parallax']
table['log_ruwe'] = np.log10(table['ruwe'])
table.loc[table['logRHK']!=table['logRHK'],'logRHK'] = -6.0
#airmass = pd.read_pickle('/Users/cretignier/Documents/THE/Airmass_gr8.p')
table['index'] = np.arange(len(table))

#HWO flag
table['HWO'] = 0
table.loc[(table['Teff']>5300)&(table['Teff']<6000)&(table['dist']<20),'HWO'] = 1
table.loc[(table['Teff']>4500)&(table['Teff']<5300)&(table['dist']<12),'HWO'] = 1
table.loc[(table['Teff']<4500)&(table['Teff']<5300)&(table['dist']<5),'HWO'] = 1

##########################
####### FUNCTIONS #######
##########################

def find_cutoff(table,parameter,N):
    if N>len(table):
        N = len(table)
    
    if parameter[-1]=='<':
        return np.sort(table[parameter[:-1]])[N]
    else:
        return np.sort(table[parameter[:-1]])[-N]


def func_cutoff(table, cutoff, tagname='', plot=True, par_space='', par_box=['',''], par_crit=''):
    'par_space format : P1 & P2'
    'par_box format : P1_min -> P1_max & P2_min -> P2_max'

    table2 = table.copy()
    count=0
    nb_rows = (len(cutoff)-1)//4+1
    old_value = np.nan
    old_value2 = len(table)
    for kw in cutoff.keys():
        count+=1
        value = cutoff[kw]
        if kw[-1]=='<':
            mask = table2[kw[0:-1]]<value
        else:
            mask = table2[kw[0:-1]]>value
        
        if plot:
            plt.figure('cumulative'+tagname,figsize=(16,4*nb_rows))
            plt.subplot(nb_rows,4,count)
            plt.title(kw+str(value))
            plt.hist(table2[kw[0:-1]],cumulative=True,bins=100)
            plt.axvline(x=value,label='%.0f / %.0f'%(sum(mask),len(mask)),color='k')
            if kw[-1]=='<':
                plt.axvspan(xmin=value,xmax=np.nanmax(table2[kw[0:-1]]),alpha=0.2,color='k')
            else:
                plt.axvspan(xmax=value,xmin=np.nanmin(table2[kw[0:-1]]),alpha=0.2,color='k')
            plt.legend()
            plt.grid()
            plt.xlim(np.nanmin(table2[kw[0:-1]]),np.nanmax(table2[kw[0:-1]]))
            table2 = table2[mask].reset_index(drop=True)
            if par_space!='':
                p1 = par_space.split('&')[0].replace(' ','')
                p2 = par_space.split('&')[1].replace(' ','')
                plt.figure('para'+tagname,figsize=(18,4*nb_rows))
                if count==1:
                    plt.subplot(nb_rows,4,count)
                    ax1 = plt.gca()
                else:
                    plt.subplot(nb_rows,4,count,sharex=ax1,sharey=ax1)
                plt.scatter(table[p1],table[p2],color='k',alpha=0.1,marker='.')
                plt.scatter(table2[p1],table2[p2],color='r',ec='k',marker='.',label='%.0f (-%.0f)'%(len(table2),old_value2-len(table2)))
                old_value2 = len(table2)
                if (par_box[0]!='')|(par_box[1]!=''):
                    p1x = np.array(par_box[0].split('->')).astype('float')
                    p2y = np.array(par_box[1].split('->')).astype('float')
                    mask_box = (table2[p1]>p1x[0])&(table2[p1]<p1x[1])&(table2[p2]>p2y[0])&(table2[p2]<p2y[1])
                    plt.scatter(table2.loc[mask_box,p1],table2.loc[mask_box,p2],color='g',ec='k',marker='o',label='%.0f (-%.0f)'%(sum(mask_box),old_value-sum(mask_box)))
                    old_value = sum(mask_box)
                if par_crit!='':
                    p1c = par_crit.split('==')[0]
                    p1c_val = float(par_crit.split('==')[1])
                    mask_box = (table2[p1c].astype('float')==p1c_val)
                    plt.scatter(table2.loc[mask_box,p1],table2.loc[mask_box,p2],color='g',ec='k',marker='o',label='%.0f (-%.0f)'%(sum(mask_box),old_value-sum(mask_box)))
                    old_value = sum(mask_box)   

                plt.xlabel(p1)
                plt.ylabel(p2)
                plt.legend(loc=1)
                plt.title(kw+str(value))

    if plot:
        plt.figure('cumulative'+tagname,figsize=(18,4*nb_rows))
        plt.subplots_adjust(hspace=0.3,wspace=0.3,top=0.95,bottom=0.08,left=0.08,right=0.95)
        if par_space!='':
            plt.figure('para'+tagname,figsize=(18,4*nb_rows))
            plt.subplots_adjust(hspace=0.3,wspace=0.3,top=0.95,bottom=0.08,left=0.08,right=0.95)
        plt.show()
    table2 = table2.sort_values(by='HZ_mp_min_osc+gr_texp15')
    
    
    #printable table
    print_table = table2[0:30][['ra_j2000','dec_j2000','primary_name','gmag','eff_nights','dist','Teff','HZ_mp_min_osc+gr_texp15']]
    print_table['gmag'] = np.round(print_table['gmag'],2)
    print_table['dist'] = np.round(print_table['dist'],1)
    print_table['eff_nights'] = (np.round(print_table['eff_nights'],0)).astype('int')
    print_table['Teff'] = (np.round(print_table['Teff'],0)).astype('int')
    print_table['HZ_mp_min_osc+gr_texp15'] = np.round(print_table['HZ_mp_min_osc+gr_texp15'],2)
    print('\n [INFO] Here are the top 30-ranked stars of your THE list:\n')          
    print(print_table)
    
    return table2

def plot_sky(table, color='k', alpha=1.0, s=5, fig=None, ax=None, plato=False, kepler=False):

    if module_warning:
        print(' [ERROR] Astropy Python library is not installed')
    else:
        ra = np.array(table['ra_j2000'])
        dec = np.array(table['dec_j2000'])
        dist = np.array(table['dist'])

        stars = SkyCoord(ra,dec,dist, unit=(u.deg, u.deg, u.pc)) 
        sph = stars.spherical

        first_fig = False
        if fig is None:
            first_fig = True
            fig, ax = plt.subplots(1,1,figsize=(12, 6), 
                                subplot_kw=dict(projection="aitoff"))

        ax.scatter(-sph.lon.wrap_at(180*u.deg).radian,
                        sph.lat.radian,s=s,alpha=alpha,
                        color=color)
        
        if plato: #TO BE CHECKED
            plato = SkyCoord([60,100,100,60,60],[0,0,45,45,0],[10,10,10,10,10], unit=(u.deg, u.deg, u.pc)).spherical
            plt.plot(-plato.lon.wrap_at(180*u.deg).radian,plato.lat.radian,color='k',ls='-.')
            plato = SkyCoord([70,80,90,80,70],[23,15,23,32,23],[10,10,10,10,10], unit=(u.deg, u.deg, u.pc)).spherical
            plt.plot(-plato.lon.wrap_at(180*u.deg).radian,plato.lat.radian,color='k',ls='-')

        if kepler: #TO BE CHECKED
            kepler = SkyCoord([360-285,360-299,360-299,360-285,360-285],[38,38,50,50,38],[10,10,10,10,10], unit=(u.deg, u.deg, u.pc)).spherical
            plt.plot(-kepler.lon.wrap_at(180*u.deg).radian,kepler.lat.radian,color='b',ls='-')

        #pcolormesh

        if first_fig:
            plt.grid()
        plt.show()

        return fig, ax


##########################
####### MAIN CODE #######  <---------------------------------------------------------
##########################

"""
ra_j2000: --- Right ascension coordinate
dec_j2000: --- Declination coordinate
gmag: --- Gaia G magnitude
dist: --- Distance in parsec
Teff: --- Effective temperature of the star
logg: --- Surface gravity in cgs
Fe/H: --- Stellar mettalicity 
vsini: --- projected equatorial velocity in kms
log_ruwe: --- Log10 of the GAIA RUWE value (binaries cutoff)
airmass_min: --- minimum airmass along the year
tyr_rise: --- Decimal year start of the season
tyr_set: --- Decimal year end of the season
NIGHTS: --- geometrical season length at La Palma above airmass of 1.75 (Tim version)
season_length: --- geometrical season length at La Palma above airmass of 1.75
eff_nights: --- Effective number of nights (Product of weather forecast with geometrical season length)
eff_airmass_mean: --- Mean optimal airmass value (Product of weather forecast with visibility curve)
eff_seeing: --- Mean expected seeing value (Function of mean_airmass and seasonal seeing)
nb_subexp: --- Number of subexp required to not saturate the detector on 15 minutes exposures
snr_550: 
snr_550_texp550: --- Expected SNR at 550 nm for 15min exposure (Clark ETC + Atmos parameters)
sig_rv_texp15: --- Expected RV precision at 550 nm for 15min exposure (Clark ETC + Atmos parameters)
texp_snr_250: --- Required exposure time to reach a SNR=250 at 550 nm
sig_rv_osc_texp15: --- Expected RV residual due to p-mode with 15min exposure
sig_rv_osc+gr_texp15: --- Expected RV residual due to p-mode + granulation with 15min exposure
HZ_period_inf: --- Inner habitable zone period
HZ_period_inf: --- outer habitable zone period
HZ_amp_inf: --- Inner habitable zone RV K-semi amplitude (1 earth mass)
HZ_amp_sup: --- Outer habitable zone RV K-semi amplitude (1 earth mass)
HZ_mp_min_texp15: --- the minimum mass planet that can be detected in the inner HZ if observed every night for 8 years (photo noise)
HZ_mp_min_texp15: --- the minimum mass planet that can be detected in the inner HZ if observed every night for 8 years (photo noise)
HZ_mp_min_osc_texp15: --- the minimum mass planet that can be detected in the inner HZ if observed every night for 8 years (photo noise + p-mode)
HZ_mp_min_osc+gr_texp15: --- the minimum mass planet that can be detected in the inner HZ if observed every night (photo noise + p-mode + granulation)
"""

# example1 of classical cutoff parameters
# use a < or > sign after the parameter name
# the cutoff are applied in the list order of the cutoff0 dictionnary

cutoff0 = {
    'Teff<':6000,
    'logg>':4.0,
    'vsini<':8,
    'Fe/H>':-0.4,
    'airmass_min>':0.99,
    'eff_nights>':160,
    'NIGHTS>':240,
    'gmag<':7.0,
    'log_ruwe<':np.round(np.log10(1.2),3),
    'dist>':0,
    'sig_rv_texp15>':0,
    'sig_rv_osc_texp15>':0,
    'sig_rv_osc+gr_texp15>':0,
    'HZ_mp_min_texp15>':0,
    'HZ_mp_min_osc+gr_texp15<':20,
    }

table_filtered = func_cutoff(table,cutoff0,tagname='')

# example2 following a peculiar parameter space box statistic
# let's define a cutoff list

cutoff1 = {
    'Teff<':6000,
    'logg>':4.0,
    'vsini<':8,
    'Fe/H>':-0.4,
    'log_ruwe<':np.round(np.log10(1.2),3),
    'eff_nights>':160,
    'NIGHTS>':240,
    'dist>':0,
    }

# we want to follow the stastistic in a 'Teff' vs 'dist' box
# let's focus on the K-dwarf sample of nearby star
# par_space = parameter1&parameter2 (name convention)
# par_box = ['lim1->lim2','lim1->lim2'] (name convention)

table_filtered = func_cutoff(table,cutoff1,par_space='Teff&dist',par_box=['4500->5300','0->30'])

# we can also follow a binary flag column e.g the target list of HWO

table_filtered = func_cutoff(table,cutoff1,par_space='Teff&dist',par_crit='HWO==1')

# we can also visualize them in the sky

table_filtered = func_cutoff(table,cutoff1,par_space='ra_j2000&dec_j2000',par_crit='HWO==1')


# example of the Tim 33 stars sample
# Note there are now 47 stars since those without RHK values are no more rejected

cutoff_tim = {
    'gmag<':7,
    'NIGHTS>':240,
    'Teff_phot<':6000,
    'Teff_spec<':6000,
    'Fe/H_spec>':-0.4,
    'parallax>':25,
    'log_ruwe<':np.round(np.log10(1.2),3),
    'vsini_spec<':3,
    'logRHK<':-4.8,
    'logg>':2.0,
    'eff_nights>':0,
    'HZ_mp_min_osc+gr_texp15>':0,
    }

table_filtered = func_cutoff(table,cutoff_tim,tagname='_Tim',par_space='ra_j2000&dec_j2000',par_crit='HWO==1')

