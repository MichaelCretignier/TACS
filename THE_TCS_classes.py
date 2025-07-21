import os

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import THE_TCS_functions as tcsf

#IMPORT MAIN TABLES

cwd = os.getcwd()

gr8 = pd.read_csv(cwd+'/TACS_Material/Master_table.csv',index_col=0)
db_starname = pd.read_csv(cwd+'/TACS_Material/simbad_name.csv',index_col=0)
table_time = pd.read_csv(cwd+'/TACS_Material/time_conversion.csv',index_col=0)

seeing = interp1d(table_time['deci'].astype('float')-2026,table_time['seeing'].astype('float'), kind='cubic', bounds_error=False, fill_value='extrapolate')(np.linspace(0,1,365))
downtime = interp1d(table_time['deci'].astype('float')-2026,table_time['downtime'].astype('float'), kind='cubic', bounds_error=False, fill_value='extrapolate')(np.linspace(0,1,365))

def crossmatch_names(tab1,kw):
    names1 = np.array(tab1[kw])
    index = [np.nan]*len(names1)
    for column in db_starname.keys():
        for m,n in enumerate(names1):
            loc = np.where(np.array(db_starname[column]==n))[0]
            if len(loc):
                index[m] = loc[0]
    index = np.array(index)
    for column in list(db_starname.keys())[:-1]:
        tab1[column] = np.array(db_starname.loc[index,column])
    return tab1

lightcurves = pd.read_pickle(cwd+'/TACS_Material/Lightcurves.p')
db_exoplanets = pd.read_csv(cwd+'/TACS_Material/exoplanets_db.csv',index_col=0)
db_exoplanets.loc[db_exoplanets['k']!=db_exoplanets['k'],'k'] = 0.1
db_exoplanets = crossmatch_names(db_exoplanets,'host')
db_exoplanets = db_exoplanets.sort_values(by='PRIMARY').reset_index(drop=True)

gr8['Teff'] = gr8['Teff_spec']
gr8.loc[gr8['Teff']!=gr8['Teff'],'Teff'] = gr8.loc[gr8['Teff']!=gr8['Teff'],'Teff_gspphot']

gr8['logg'] = gr8['MIST logg']
gr8.loc[gr8['logg']!=gr8['logg'],'logg'] = gr8.loc[gr8['logg']!=gr8['logg'],'logg_spec']

gr8.loc[gr8['Fe/H']!=gr8['Fe/H'],'Fe/H'] = 2

gr8['dist'] = 1000/gr8['parallax']
gr8['log_ruwe'] = np.log10(gr8['ruwe'])
gr8.loc[gr8['logRHK']!=gr8['logRHK'],'logRHK'] = -6.0
gr8.loc[gr8['RHK']!=gr8['RHK'],'RHK'] = -6.0
gr8.loc[gr8['vsini']!=gr8['vsini'],'vsini'] = 0.0
gr8.loc[gr8['vsini']>50,'vsini'] = 50

gr8['HWO'] = 0
gr8.loc[(gr8['Teff']>5300)&(gr8['Teff']<6000)&(gr8['dist']<20),'HWO'] = 1
gr8.loc[(gr8['Teff']>4500)&(gr8['Teff']<5300)&(gr8['dist']<12),'HWO'] = 1
gr8.loc[(gr8['Teff']<4500)&(gr8['Teff']<5300)&(gr8['dist']<5),'HWO'] = 1


#FUNCTIONS

def auto_format(values):
    maxi = np.nanmax(values)
    if maxi!=0:
        maxi = np.round(np.log10(maxi),0)
    digit = int(4-maxi)
    return np.round(values,digit)


def dropout_year():
    down = interp1d(table_time['deci'].astype('float')-2026,table_time['downtime'].astype('float'), kind='cubic', bounds_error=False, fill_value='extrapolate')(np.linspace(0,1,365))
    save = []
    for d in down:
        a = np.random.choice([0]*int(d*100)+[1]*(10000-int(d*100)),1)
        save.append(a)
    save = np.array(save).astype('bool')
    return save

def query_table(ra,dec,table):
    dist = abs(table['dec_j2000']-dec) + abs(table['ra_j2000']-ra)
    return table.loc[np.argmin(dist)]

def resolve_starname(name,verbose=True):
    button = 0
    for columns in list(db_starname.keys()):
       loc =  np.where(np.array(db_starname[columns])==name)[0]
       if len(loc):
            loc = loc[0]
            button = 1
            break
       else:
           loc =  np.where(np.array(db_starname[columns])==name.replace(' ',''))[0]
           if len(loc):
                loc = loc[0]
                button = 1
                break
            
    if button:
        output = db_starname.iloc[loc]
        output['INDEX'] = loc
        if verbose:
            print('\n [INFO] Starnames found:')
            print(output,'\n')
        return output
    else:
        print(' [INFO] Starname has not been found')
        return None


def get_starname(entry):
    gaia_ID = np.where(db_starname['GAIA']=='Gaia DR3 '+str(entry['gaiaedr3_source_id']))[0]
    selected = db_starname.loc[gaia_ID]
    for order in ['HD','HIP','GJ','CSTL','PRIMARY']:
        if selected[order].values[0]!='-':
            break
    return (selected[order].values[0], gaia_ID)

def star_info(entry, format='v1'):
    name, ID = get_starname(entry)
    if format=='v1':
        info = ' ID : %.0f \n Star : %s   Mv = %.2f \n Ra = %.2f    Dec = %.2f \n Teff = %.0f   Logg = %.2f \n FeH = %.2f    RHK = %.2f   Vsini = %.1f'%(ID,name,entry['Vmag'], entry['ra_j2000'], entry['dec_j2000'], entry['Teff'], entry['MIST logg'], entry['Fe/H'], entry['RHK'], entry['vsini'])
    else:
        info = ' ID : %.0f   Star : %s   Mv = %.2f \n Ra = %.2f    Dec = %.2f \n Teff = %.0f   Logg = %.2f    FeH = %.2f    RHK = %.2f   Vsini = %.1f'%(ID,name,entry['Vmag'], entry['ra_j2000'], entry['dec_j2000'], entry['Teff'], entry['MIST logg'], entry['Fe/H'], entry['RHK'], entry['vsini'])
    return info

def plot_TESS_CVZ():
    theta = np.linspace(0,2*np.pi,100)
    plt.plot(np.cos(theta)*12/360*24+18,np.sin(theta)*12+66,lw=1,ls='-.',color='k')
    plt.text(18,66,'TESS',ha='center',va='center')

def plot_KEPLER_CVZ():
    theta = np.linspace(0,2*np.pi,100)
    plt.plot(np.cos(theta)*8/360*24+19.5,np.sin(theta)*8+44.5,lw=1,ls=':',color='k')
    plt.text(19.5,44.5,'KEPLER',ha='center',va='center')

def plot_exoplanets(y_var='k'):
    fig = plt.figure(figsize=(8,8))

    plt.scatter(db_exoplanets['period'],db_exoplanets[y_var],color='k',marker='.')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Period [days]',fontsize=14)
    if y_var=='k':
        plt.ylabel('K semi-amplitude [m/s]',fontsize=14)
        plt.ylim(0.1,1000)
    else:
        plt.ylabel(r'Planetary mass [$M_{\oplus}$]',fontsize=14)       
        plt.ylim(1,10000)
    plt.axvspan(xmin=60,xmax=400,ymin=0,ymax=0.25,color='g',alpha=0.2) 
    plt.grid()
    plt.scatter(db_exoplanets['period'],db_exoplanets[y_var],c=db_exoplanets['radius'],cmap='brg',vmin=1,vmax=10)
    ax = plt.colorbar(pad=0)
    ax.ax.set_ylabel(r'Planetary radius [$R_{\oplus}$]',fontsize=14)
    plt.subplots_adjust(left=0.10,right=0.95)
    return fig 

def plot_exoplanets2(cutoff={'MIST Teff<':6000},mcrit_sup=4000,mcrit_inf=100):
    table_gr8 = gr8.copy() 
    for kw in cutoff.keys():
        if kw[-1]=='<':
            table_gr8 = table_gr8.loc[table_gr8[kw[:-1]]<cutoff[kw]]
        else:
            table_gr8 = table_gr8.loc[table_gr8[kw[:-1]]>cutoff[kw]]

    selected = np.in1d(np.array([i.split(' ')[-1] for i in np.array(db_exoplanets['host'])]),np.array(table_gr8['gaiaedr3_source_id']))
    db_exoplanets['pre_survey'] = selected.astype('int')

    fig = plt.figure(figsize=(18,12))
    count=0
    markersize = np.array(db_exoplanets['mass']>0)*10
    markersize += np.array(db_exoplanets['mass']>10)*20
    markersize += np.array(db_exoplanets['mass']>30)*50
    markersize += np.array(db_exoplanets['mass']>100)*150
    db_exoplanets['marker'] = np.sqrt(db_exoplanets['mass'])*10
    pmin = db_exoplanets['period']*(1-db_exoplanets['ecc'])/np.sqrt(1+db_exoplanets['ecc'])
    pmax = db_exoplanets['period']*(1+db_exoplanets['ecc'])/np.sqrt(1-db_exoplanets['ecc'])

    db_exoplanets['p_eccmin'] = pmin
    db_exoplanets['p_eccmax'] = pmax

    for j in range(3):
        plt.subplot(1,3,j+1) ; plt.xscale('log')
        plt.axvspan(xmin=60,xmax=400,color='g',alpha=0.2)
        plt.tick_params(labelleft=False)
        plt.xlim(0.5,80000)
        plt.xlabel('Period [days]',fontsize=13)
    
    db = np.sort(np.unique(db_exoplanets['PRIMARY']))
    db = np.array([np.where(np.array(db_starname['PRIMARY'])==d)[0][0] for d in db])

    summary = []
    for system in db_starname.loc[db,'GAIA']:
        plt.subplot(1,3,(abs(count)//30)+1)
        count-=1
        mask = np.array(db_exoplanets['host']==system)
        syst = db_exoplanets.loc[mask]
        plt.plot([1,50000],[count,count],lw=1,color='k',ls='-',alpha=0.2)
        plt.scatter(syst['period'],syst['period']*0+count,s=syst['marker'],color='k',zorder=10)
        plt.scatter(syst['period'],syst['period']*0+count,s=syst['marker'],c=syst['radius'],cmap='brg',vmin=1,vmax=8,zorder=10,edgecolor='k')
        condition1 = np.sum((syst['p_eccmin']<400)&(syst['mass']>mcrit_inf)).astype('bool')
        condition2 = np.sum((syst['mass']>mcrit_sup)).astype('bool')
        condition = condition1|condition2
        summary.append([system,int(condition1),int(condition2)])       
        color_condition = ['k','r'][int(condition)]
        indicator = ['x','•'][np.array(syst['pre_survey'])[0]]
        plt.text(100000,count,indicator+' '+db_starname.loc[db_starname['GAIA']==system,'PRIMARY'].values[0],va='center',ha='left',color=color_condition,alpha=[0.25,1][np.array(syst['pre_survey'])[0]])
        for p1,p2 in zip(syst['p_eccmin'],syst['p_eccmax']):
            plt.plot([p1,p2],[count,count],color='k',lw=3)
        for mass,period,p1 in np.array(syst[['mass','period','p_eccmin']]):
            if mass>mcrit_sup:
                plt.text(period,count,'%.0f'%(np.round(mass/95,0)),color='r',va='center',ha='center',zorder=1000)
            elif mass>mcrit_inf:
                plt.text(period,count,'%.0f'%(np.round(mass/95,0)),color=['white','r'][int(p1<400)],va='center',ha='center',zorder=1000)
    summary = np.array(summary)
    summary = pd.DataFrame(summary,columns=['GAIA','HJ','BDW'])
    plt.subplots_adjust(left=0.03,right=0.93,wspace=0.30)
    return summary

class tableXY(object):
    def __init__(self, x=None, y=None, xlabel='', ylabel='', ls='o'):
        if x is None:
            x = np.arange(len(y))

        if y is None:
            y = np.zeros(len(x))

        self.x = x
        self.y = y
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.ls = ls

    def interpolate(self,new_grid):
        self.y = interp1d(self.x,self.y, kind='cubic', bounds_error=False, fill_value='extrapolate')(new_grid)
        self.x = new_grid

    def monthly_average(self):
        borders = [0,31,59,90,120,151,181,212,243,273,304,334,365]
        statistic = []
        for b1,b2 in zip(borders[0:-1],borders[1:]):
            mask = (self.x>=b1)&(self.x<=b2)
            statistic.append(np.mean(self.y[mask]))
        statistic = np.array(statistic)
        return auto_format(statistic)

    def plot(self, alpha=0.5, label=None, ytext=60, figure=None, ls=None, subset=None):
        if figure is not None:
            plt.figure(figure)
        if ls is None:
            ls = self.ls
        if subset is None:
            subset = np.arange(0,len(self.x)).astype('int')
        
        if ls=='o':
            plt.scatter(self.x[subset],self.y[subset],alpha=alpha,label=label)
        else:
            plt.plot(self.x,self.y,alpha=alpha,label=label,ls=ls,color='k')
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        ax = plt.gca()
        if self.xlabel=='Nights [days]':
            for j,t in zip([31,59,90,120,151,181,212,243,273,304,334,365],['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
                plt.axvline(x=j,color='k',ls='-',alpha=0.5)
                plt.axvline(x=j,color='white',ls=':')
                plt.text(j-15,ytext,t,color='k',ha='center')
            plt.xlim(0,365)

    def create_subset(self,subset):
        self.subset = tableXY(x=self.x[subset], y=self.y[subset], xlabel=self.xlabel, ylabel=self.ylabel)

    def night_subset(self,obs_per_night,random=False,replace=False):
        
        if obs_per_night is not None:
            nights = self.x.astype('int')
            selection = []
            for n in np.unique(nights):
                loc = np.where(nights==n)[0]
                m = np.min([len(loc),obs_per_night])
                if random:
                    selection.append(np.random.choice(loc,m,replace=False))
                else:
                    index = np.unique(np.round(np.linspace(0+np.min(loc),np.max(loc),obs_per_night),0)).astype('int')
                    selection.append(index)
            selection = np.sort(np.hstack(selection))
        else:
            selection = np.arange(len(self.x))

        if replace:
            self.x = self.x[selection]
            self.y = self.y[selection]
            return None
        else:
            self.create_subset(selection)
            return selection


class image(object):
    def __init__(self,data,xlabel='',ylabel='',zlabel=''):
        self.data = data
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel

    def colorbar(self):
        plt.colorbar()
        ax = plt.gca()
        #ax.set_ylabel(self.zlabel)

    def plot(self,alpha=1.0,colorbar=False,vmin=None,vmax=None,ytext=60):
        plt.imshow(self.data,alpha=alpha,aspect='auto',vmin=vmin,vmax=vmax)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        if colorbar:
            self.colorbar()
        if self.ylabel=='UT Time []':
            plt.axhline(y=500,color='k',ls=':',lw=1)
            for j in np.arange(0,1441,120):
                plt.axhline(y=j,lw=1,color='white',alpha=0.5)
            plt.yticks(np.arange(0,1441,120),['%.0f'%(i) for i in np.arange(-12,13,2)])

        if self.xlabel=='Nights [days]':
            for j,t in zip([31,59,90,120,151,181,212,243,273,304,334,365],['J','F','M','A','M','J','J','A','S','O','N','D']):
                plt.axvline(x=j,color='k',ls='-',alpha=0.5)
                plt.axvline(x=j,color='white',ls=':')
                plt.text(j-15,ytext,t,color='w')
            plt.xlim(0,365)

class tcs(object):
    
    def __init__(self, sun_elevation=None, starname=None, instrument='HARPS3'):    
        self.info_XY_telescope_open = tableXY(y=np.ones(365))
        self.info_XY_downtime = tableXY(x=np.arange(365),y=table_time['downtime'])
        self.simu_SG_calendar = None
        self.simu_counter_survey = 0
        self.info_SC_starname = None
        self.info_SC_instrument = instrument
        self.info_TA_cutoff = {
            'Teff<':6000,
            'logg>':4.2,
            'vsini<':8,
            'Fe/H>':-0.4,
            'eff_nights_1.75>':180,
            'season_length_1.75>':240,
            'log_ruwe<':np.round(np.log10(1.2),3),
            'HJ<':0.5,
            'BDW<':0.5,
            'RHK<':-4.8,
            'dist>':0,
            'gmag<':7.0,
            }

        if sun_elevation is not None:
            self.compute_night_length(sun_elevation=sun_elevation)
        if starname is not None:
            self.set_star(starname=starname)
            self.random_weather()

    def compute_night_length(self, sun_elevation=-12, verbose=True):
        almanac_table = pd.read_pickle(cwd+'/TACS_Material/almanac.p')[self.info_SC_instrument]
        almanac = almanac_table['sun']
        self.info_SC_night_def = sun_elevation
        self.info_IM_night = image(
            ((almanac>sun_elevation).T).astype('float'),
            xlabel='Nights [days]',
            ylabel='UT time []',
            zlabel='Night time')
        
        if (self.info_SC_instrument=='HARPN')|(self.info_SC_instrument=='HARPS3'):
            table_time = pd.read_csv(cwd+'/TACS_Material/time_conversion.csv',index_col=0)
            down = interp1d(table_time['deci'].astype('float')-2026,table_time['downtime'].astype('float'), kind='cubic', bounds_error=False, fill_value='extrapolate')(np.linspace(0,1,365))
            hours = np.sum(1-self.info_IM_night.data,axis=0)/60
            nb_hours = np.sum(hours)*0.60 #60% of GTO
            gto = {'HARPS3':0.60,'HARPN':0.50}[self.info_SC_instrument]
            nb_hours_weather = np.sum(hours*(100-down)/100)*gto #60% of GTO
            self.info_SC_nb_hours_per_yr = np.round(nb_hours,1)
            self.info_SC_nb_hours_per_yr_eff = np.round(nb_hours_weather,1)
            if verbose:
                print(' [INFO] Nb of hours per year : %.0f (%.0f)'%(self.info_SC_nb_hours_per_yr_eff,self.info_SC_nb_hours_per_yr))

        if self.info_SC_starname is not None:
            self.compute_night_star()

    def random_weather(self):
        output = np.ravel(dropout_year().astype('int'))
        print(' [INFO] Number of bad/good nights = %.0f/%.0f'%(len(output)-sum(output),sum(output)))
        self.info_XY_telescope_open = tableXY(
            y=output,
            xlabel='Nights [days]',
            ylabel='Telescope open')

    def set_star(self,ra=0,dec=0,starname=None,id=None,verbose=True):
        self.info_TA_starnames = None
        if id is not None:
            starname = gr8.iloc[id]['primary_name']

        if starname is not None:
            starname = resolve_starname(starname,verbose)
            self.info_TA_starnames = starname
            if starname is not None:
                ra = np.array(gr8.loc[gr8['primary_name']==starname['PRIMARY'],'ra_j2000'])[0]/360*24
                dec = np.array(gr8.loc[gr8['primary_name']==starname['PRIMARY'],'dec_j2000'])[0]
            else:
                print('[ERROR] STARNAME NOT FOUND')
        self.info_SC_ra = ra
        self.info_SC_dec = dec
        self.info_SC_starname = starname
        self.compute_night_star()

    def compute_night_star(self):

        alpha = self.info_SC_ra
        dec = self.info_SC_dec
        
        hours,airmass = tcsf.star_observability(alpha, dec, instrument=self.info_SC_instrument)

        save2 = np.array([np.roll(airmass,-rolling) for rolling in (np.arange(0,365)/365*len(hours)).astype('int')])
        self.info_IM_airmass = image(
            save2.T,
            xlabel='Nights [days]',
            ylabel='UT Time []',
            zlabel='Airmass')

        self.info_IM_airmass_night = image(
            save2.T*(1-self.info_IM_night.data),
            xlabel='Nights [days]',
            ylabel='UT Time []',
            zlabel='Airmass')

    def compute_statistic_star(self, min_airmass=1.75, texp=15):
        airmass = self.info_IM_airmass_night.data.copy()
        airmass[airmass>min_airmass] = np.nan
        airmass[airmass<1.0] = np.nan
        airmass[:,np.sum(airmass==airmass,axis=0)<texp] = np.nan
        min_airmass = np.nanmin(airmass,axis=0)
        weight = 1-downtime.copy()/100
        weight[min_airmass!=min_airmass] = np.nan
        season_length = int(np.round(np.sum(min_airmass==min_airmass),0))
        eff_nights = int(np.round(np.nansum(weight,0)))
        eff_airmass = np.nansum(min_airmass*weight)/np.nansum(weight)
        eff_seeing = np.nansum(seeing*(min_airmass**0.6)*weight)/np.nansum(weight)
        gap = np.where(min_airmass!=min_airmass)[0]
        if gap[0]!=0:
            tyr_set = gap[0]+365
            tyr_rise = gap[-1]
        else:
            observable = np.where(min_airmass==min_airmass)[0]
            if len(observable):
                tyr_rise = observable[0]
                tyr_set = observable[-1]
            else:
                tyr_rise = 0
                tyr_set = 0

        print(' Season length = %.0f \n Number of eff. nights = %.0f \n Mean eff. airmass = %.3f \n Mean eff. seeing = %.3f \n tyr_rise = %.0f \n tyr_set = %.0f'%(season_length,eff_nights,eff_airmass,eff_seeing,tyr_rise,tyr_set))
        return [season_length, eff_nights, eff_airmass, eff_seeing, tyr_rise, tyr_set]

    def show_lightcurve(self,rm_gap=True,plot=True):
        lc = lightcurves['QLP'][self.info_TA_starnames['GAIA']]
        
        if lc is not None:
            t = lc[0].copy()
            f = lc[1].copy()
        else:
            t = []
            f = []

        if rm_gap:
            loc = np.where(np.diff(t)>30)[0]
            for l in loc:
                t[l+1:] = t[l+1:] - (t[l+1] - t[l]) + 30
        lc = tableXY(x=t,y=f*100,xlabel='Time [days]',ylabel='Flux [%]')
        self.info_XY_lc_qlp = lc
        if plot:
            lc.plot()

    def compute_nights(self,airmass_max=1.5,weather=False, plot=False):
        
        map_night = 1-self.info_IM_night.data
        map_star = self.info_IM_airmass.data<airmass_max
        map_weather = self.info_XY_telescope_open.y*np.ones(1441)[:,np.newaxis]

        total_map = map_night*map_star
        self.info_XY_season = tableXY(
            y = np.sum(total_map,axis=0)>10,
            xlabel='Nights [days]',
            ylabel='Visible'
            ) #15-min exposure
        self.info_XY_night_duration = []
        for t in total_map.T:
            loc = np.where(t==1)[0]
            if len(loc):
                duration = (np.max(loc)-np.min(loc))*24/1441
            else:
                duration = 0
            self.info_XY_night_duration.append(duration)
        self.info_XY_night_duration = tableXY(
            y = np.array(self.info_XY_night_duration),
            xlabel='Nights [days]',
            ylabel='Night duration [hours]'
            )

        if weather:
            total_map = total_map*map_weather

        self.info_IM_observable = image(
            total_map,
            xlabel='Nights [days]',
            ylabel='UT Time []',
            zlabel='Observable')
        
        if plot:
            self.info_IM_observable.plot()
        
    def create_timeseries(self, airmass_max=1.5, nb_year=10, texp=15, weather=True, nb_subexp=None):
        dt = 24*60/1440 #dt in minutes
        N = int(np.round(texp/dt,0))
        j0 = 0.0#61041.0
        jdb = []
        for year in range(nb_year):
            if weather:
                self.random_weather()
            self.compute_nights(airmass_max=airmass_max, 
                                weather=weather, 
                                plot=False)
            epochs = []
            for d,n in enumerate(self.info_IM_observable.data.T):
                loc = np.where(n==1)[0]                
                if len(loc)>N:
                    loc = loc[int(N/2):-int(N/2)]
                    mini = loc[0]
                    maxi = loc[-1]
                    t0 = np.mean(loc)+np.arange(-100,100,1)*N
                    t0 = t0[(t0>mini)&(t0<maxi)]
                    epochs.append(j0+d+0.5+(t0/1441-0.5)+365*year)
            jdb.append(epochs)
        jdb = np.hstack(jdb)
        jdb = np.hstack(jdb)
        self.info_XY_timestamps = tableXY(x=jdb,y=np.random.randn(len(jdb)),xlabel='Nights [days]')
        
        if nb_subexp is not None:
            for n in np.unique(jdb//1):
                l2 = np.where((jdb//1)==n)[0]

    def compute_SG_calendar(
            self,
            sun_elevation=-12,
            airmass_max=1.8,
            alpha_step=1, 
            dec_step=1,
            cutoff=None):
        
        if cutoff is None:
            cutoff = self.info_TA_cutoff

        self.compute_night_length(sun_elevation=sun_elevation) 
        
        button = 1
        if self.simu_SG_calendar is not None:
            if (self.simu_SG_calendar['param1']==sun_elevation)&(self.simu_SG_calendar['param2']==airmass_max)&(self.simu_SG_calendar['param3']==alpha_step)&(self.simu_SG_calendar['param4']==dec_step):
                button = 0
                print(' [INFO] Old simulation found.')
        
        if button==1:
            output = []
            params = []
            RA, DEC = np.meshgrid(np.arange(0,30,alpha_step),np.arange(-30,90,dec_step))
            loading = np.round(len(np.ravel(RA))*np.arange(0,101,10)/100,0).astype('int')
            counter=0
            for i,j in zip(np.ravel(RA),np.ravel(DEC)):
                if counter in loading:
                    print(' [INFO] Progress... [%.0f%%]'%(np.where(loading==counter)[0][0]*10))
                self.set_star(ra=i,dec=j) 
                self.compute_nights(airmass_max=airmass_max, weather=False, plot=False)
                params.append([i,j])
                output.append(self.info_XY_night_duration.monthly_average())
                counter+=1
            print(' [INFO] Finished!')
            print(' [INFO] Stacking tables...')
            params = np.array(params)
            output = np.array(output)

            self.simu_SG_calendar = {
                'param1':sun_elevation,
                'param2':airmass_max,
                'param3':alpha_step,
                'param4':dec_step,
                'outputs':(params,output,RA,DEC)}
        else:
            params,output,RA,DEC = self.simu_SG_calendar['outputs']
        
        print(' [INFO] Producing Calendar plot... Wait...')
        fig = plt.figure(figsize=(18,12))
        fig.suptitle('Sun elevation = %.0f | Airmass max = %.2f'%(sun_elevation,airmass_max))
        plt.subplots_adjust(left=0.05,right=0.98,top=0.94,bottom=0.05,hspace=0.30,wspace=0.30)

        downtime = np.array([32., 31., 31., 27., 15.,  4.,  3.,  8., 19., 30., 37., 38., 32.])

        table_gr8 = gr8.copy() 
        for kw in cutoff.keys():
            if kw[-1]=='<':
                table_gr8 = table_gr8.loc[table_gr8[kw[:-1]]<cutoff[kw]]
            else:
                table_gr8 = table_gr8.loc[table_gr8[kw[:-1]]>cutoff[kw]]

        for j in range(12):
            plt.subplot(3,4,j+1)
            plt.title(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'][j]+' (Bad weather = %.0f%%)'%(downtime[j]))
            cp = plt.contour(RA,DEC,np.reshape(output[:,j],np.shape(RA)),levels=[6,7,8,9,10])
            plt.clabel(cp, inline=True, fontsize=8,fmt="%.0f")
            plt.grid()
            plt.xlim(0,24)
            plt.ylim(-30,90)
            plt.xlabel('RA [hours]')
            plt.ylabel('Dec [deg]')
            plt.scatter(gr8['ra_j2000']/360*24,gr8['dec_j2000'],s=(7.5-gr8['Vmag'])*30,alpha=0.15,c=gr8['Teff'],cmap='jet_r',vmin=5000,vmax=6000)
            plt.scatter(gr8['ra_j2000']/360*24+24,gr8['dec_j2000'],s=(7.5-gr8['Vmag'])*30,alpha=0.15,c=gr8['Teff'],cmap='jet_r',vmin=5000,vmax=6000)
            plt.scatter(table_gr8['ra_j2000']/360*24,table_gr8['dec_j2000'],s=(7.5-table_gr8['Vmag'])*30,c=table_gr8['Teff'],cmap='jet_r',vmin=5000,vmax=6000,ec='k')
            plt.scatter(table_gr8['ra_j2000']/360*24+24,table_gr8['dec_j2000'],s=(7.5-table_gr8['Vmag'])*30,c=table_gr8['Teff'],cmap='jet_r',vmin=5000,vmax=6000,ec='k')
            plot_TESS_CVZ()
            plot_KEPLER_CVZ()

    def commpute_SG_month(self,month=1,cutoff=None):
        if cutoff is None:
            cutoff = self.info_TA_cutoff
        
        table_gr8 = gr8.copy() 
        for kw in cutoff.keys():
            if kw[-1]=='<':
                table_gr8 = table_gr8.loc[table_gr8[kw[:-1]]<cutoff[kw]]
            else:
                table_gr8 = table_gr8.loc[table_gr8[kw[:-1]]>cutoff[kw]]
        table_gr8 = table_gr8.reset_index(drop=True)
        
        params,output,RA,DEC = self.simu_SG_calendar['outputs']
        sun_elevation = self.simu_SG_calendar['param1']
        airmass_max = self.simu_SG_calendar['param2']

        month_tag = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'][month-1]

        fig = plt.figure(figsize=(18,12))
        fig.suptitle('Sun elevation = %.0f \nAirmass max = %.2f \nMonth = %s'%(sun_elevation,airmass_max,month_tag))
        #plt.title()
        cp = plt.contour(RA,DEC,np.reshape(output[:,month-1],np.shape(RA)),levels=[6,7,8,9,10])
        plt.clabel(cp, inline=True, fontsize=8,fmt="%.0f")
        plt.grid()
        plt.xlim(0,24)
        plt.ylim(-30,90)
        plt.xlabel('RA [hours]')
        plt.ylabel('Dec [deg]')
        plot_TESS_CVZ()
        plot_KEPLER_CVZ()
        plt.scatter(gr8['ra_j2000']/360*24,gr8['dec_j2000'],s=(7.5-gr8['Vmag'])*30,alpha=0.15,c=gr8['Teff'],cmap='jet_r',vmin=5000,vmax=6000)
        plt.scatter(gr8['ra_j2000']/360*24+24,gr8['dec_j2000'],s=(7.5-gr8['Vmag'])*30,alpha=0.15,c=gr8['Teff'],cmap='jet_r',vmin=5000,vmax=6000)
        plt.scatter(table_gr8['ra_j2000']/360*24,table_gr8['dec_j2000'],s=(7.5-table_gr8['Vmag'])*30,c=table_gr8['Teff'],cmap='jet_r',vmin=5000,vmax=6000,ec='k')
        plt.scatter(table_gr8['ra_j2000']/360*24+24,table_gr8['dec_j2000'],s=(7.5-table_gr8['Vmag'])*30,c=table_gr8['Teff'],cmap='jet_r',vmin=5000,vmax=6000,ec='k')
        plt.colorbar(pad=0)  
        plt.subplots_adjust(left=0.05,right=1.1)          

        info_text = plt.text(0,107,'Info Star',fontsize=13,ha='left',va='top')
        l, = plt.plot([-5],[130],marker='x',color='k',markersize=10)

        class Index(object):
            def __init__(self):
                self.info_text = ''
                self.marker = None
            def update(self,newx,newy):
                new_star = query_table(newx/24*360,newy,table_gr8)
                text_fmt = star_info(new_star)
                self.info_text.set_text(text_fmt)
                self.marker.set_data([new_star['ra_j2000']/360*24,new_star['dec_j2000']])
                
                plt.draw()
                fig.canvas.draw_idle()
                
        t = Index()
        t.info_text = info_text
        t.marker = l
        
        def onclick(event):
            if event.dblclick:
                t.update(event.xdata,event.ydata)

        plt.gcf().canvas.mpl_connect('button_press_event', onclick)


    def plot_exoplanets_db(self,y_var='k'):
        fig = plot_exoplanets(y_var=y_var)

        gaia_name = self.info_TA_starnames['GAIA']
        found = np.where(np.array(db_exoplanets['host'])==gaia_name)[0]
        if len(found):
            init_text = star_info(gr8.iloc[self.info_TA_starnames['INDEX']],format='v2')
            l, = plt.plot(
                np.array(db_exoplanets.loc[found,'period']),
                np.array(db_exoplanets.loc[found,y_var]),
                marker='x',color='k',markersize=10)
        else:
            init_text = 'Info Star'
            l, = plt.plot([0.3],[1000],marker='x',color='k',markersize=10)
        
        info_text = plt.text(0.7,{'k':3000,'mass':30000}[y_var],init_text,fontsize=13,ha='left',va='top')

        class Index(object):
            def __init__(self):
                self.info_text = ''
                self.marker = None
            def update(self,newx,newy):
                loc = np.argmin(abs(np.log10(db_exoplanets[y_var])-np.log10(newy))+abs(np.log10(db_exoplanets['period'])-np.log10(newx)))
                info = resolve_starname(db_exoplanets.loc[loc,'host'],verbose=False)
                new_star = gr8.iloc[info['INDEX']]
                text_fmt = star_info(new_star,format='v2')
                self.info_text.set_text(text_fmt)
                exo = db_exoplanets.loc[db_exoplanets['host']==info['GAIA']]
                self.marker.set_data([exo['period'],exo[y_var]])
                plt.draw()
                fig.canvas.draw_idle()
                
        t = Index()
        t.info_text = info_text
        t.marker = l
        
        def onclick(event):
            if event.dblclick:
                t.update(event.xdata,event.ydata)

        plt.gcf().canvas.mpl_connect('button_press_event', onclick)

    def compute_exoplanet_rv_signal(self, jdb=None, keplerian_par={}, y0=2026):
        if len(keplerian_par)==0:
            gaia_name = self.info_TA_starnames['GAIA']
            mask = np.array(db_exoplanets['host']==gaia_name)
        else:
            mask = np.array(db_exoplanets['host']=='')
        if jdb is None:
            jdb = self.info_XY_timestamps.x

        j0 = 0
        if np.mean(jdb)<40000:
            j0 = 2461041.5+365*(y0-2026) #2026-01-01
            xlabel = 'Nights [days]'
        else:
            xlabel=''
        jdb = jdb + j0

        if (sum(mask)!=0)|(len(keplerian_par)!=0):
            self.info_XY_keplerian = [[]]
            self.info_XY_keplerian_model = [[]]
            if len(keplerian_par)==0:
                syst = db_exoplanets.loc[mask]
                self.info_TA_exoplanets_known = syst[['name','period','k','mass','radius','ecc','peri','t0']]
            else:
                syst = keplerian_par.copy()

            jdb_model = np.arange(np.min(jdb),np.max(jdb),np.min(syst['period'])/20)
            for P,K,e,omega,t0 in np.array(syst[['period','k','ecc','peri','t0']]):
                signal = tcsf.Keplerian_rv(jdb, P, K, e, omega, t0)
                self.info_XY_keplerian.append(tableXY(x=jdb-j0, y=signal, xlabel=xlabel,ls='o', ylabel='RV [m/s]'))
                signal2 = tcsf.Keplerian_rv(jdb_model, P, K, e, omega, t0)
                self.info_XY_keplerian_model.append(tableXY(x=jdb_model-j0, y=signal2, xlabel=xlabel,ls='-', ylabel='RV [m/s]'))
            rv_tot = np.sum([k.y for k in self.info_XY_keplerian[1:]],axis=0)
            self.info_XY_keplerian[0] = tableXY(x=jdb-j0, y=rv_tot, xlabel=xlabel,ls='o', ylabel='RV tot [m/s]')
            rv_tot2 = np.sum([k.y for k in self.info_XY_keplerian_model[1:]],axis=0)
            self.info_XY_keplerian_model[0] = tableXY(x=jdb_model-j0, y=rv_tot2, xlabel=xlabel,ls='-', ylabel='RV tot [m/s]')
        
        else:
            self.info_TA_exoplanets_known = []
            print(' [INFO] No exoplanets found for this target.')

    def plot_keplerians(self,axhline=None, obs_per_night=None, random=False):
        nb = len(self.info_XY_keplerian)
        plt.figure(figsize=(18,10))
        for j in np.arange(nb):
            if not j:
                plt.subplot(nb,1,j+1) ; ax = plt.gca()
            else:
                plt.subplot(nb,1,j+1,sharex=ax)
            ymax = np.max(self.info_XY_keplerian_model[j].y)
            ymin = np.min(self.info_XY_keplerian_model[j].y)
            
            self.info_XY_keplerian_model[j].plot(ytext=ymax*1.2)
            if j==0:
                selection = self.info_XY_keplerian[j].night_subset(obs_per_night,random=random)
            self.info_XY_keplerian[j].create_subset(selection)
            self.info_XY_keplerian[j].subset.plot(ytext=ymax*1.2)
            plt.ylim(1.5*ymin,1.5*ymax)
            if axhline is not None:
                plt.axhline(y=axhline,alpha=0.4,color='k',lw=1)


    def plot_night_length(self):
        backup = np.array([self.info_SC_night_def]).copy()
        
        self.compute_night_length(sun_elevation=-12, verbose=False) 
        self.compute_nights(airmass_max=1.5, weather=False, plot=False)
        self.info_XY_night_duration.plot(figure='NightLength',label='Z=1.5 | S=-12',ytext=-0.5) 
        
        self.compute_nights(airmass_max=1.8, weather=False, plot=False)
        self.info_XY_night_duration.plot(figure='NightLength',label='Z=1.8 | S=-12') 
        
        self.compute_night_length(sun_elevation=-6, verbose=False) #change the sunset/rise parameter
        self.compute_nights(airmass_max=1.8, weather=False, plot=False)
        self.info_XY_night_duration.plot(figure='NightLength',label='Z=1.8 | S=-6')
        plt.legend()
        self.compute_night_length(sun_elevation=backup[0], verbose=False) 
        plt.ylim(-1)

    def plot_survey_snr_texp(self, selection=None, texp=None, snr_crit=250, sig_rv_crit=0.30, budget='_phot'):
        
        if selection is None:
            selection = gr8.copy()
        
        snr_texp15 = 0.5*(selection['snr_420_texp15']+selection['snr_550_texp15']) #Cretignier et al. +22
        #snr_texp15 = gr8['snr_550_texp15']
        texp_snr250 = selection['texp_snr_250']
        if budget=='_phot':
            sig_rv_texp15 = selection['sig_rv'+budget+'_texp15']
            sig_rv = sig_rv_texp15/(np.sqrt(texp/15))
        else:
            texp = tcsf.find_nearest(np.array([1,5,10,12,15,20]),texp)[1][0]
            print(' [INFO] Closest Texp tabulated = %.0f min'%(texp))
            sig_rv = selection['sig_rv'+budget+'_texp%.0f'%(texp)]

        snr = snr_texp15*(np.sqrt(texp/15))
        mask = (snr>snr_crit)&(sig_rv<sig_rv_crit)
        plt.figure('SNR_SIGRV_TEXP%.0fMIN'%(texp),figsize=(8,8))
        plt.axes([0.12,0.12,0.7,0.7])
        plt.scatter(snr,sig_rv,color='k',label='%.0f'%(len(mask)))
        plt.scatter(snr[mask],sig_rv[mask],color='g',marker='.',label='%.0f'%(sum(mask)))
        plt.axvline(x=snr_crit,color='k',ls=':')
        plt.axhline(y=sig_rv_crit,color='k',ls=':')
        plt.axvspan(xmin=0,xmax=snr_crit,color='k',alpha=0.2)
        plt.axhspan(ymin=sig_rv_crit,ymax=1.0,color='k',alpha=0.2)
        plt.xlabel('SNR continuum',fontsize=14)
        plt.ylabel(r'$\sigma_{{\gamma}}$ $RV$ [m/s]',fontsize=14)
        plt.legend()
        plt.xlim(0,1000)
        plt.ylim(0,0.75)
        ax = plt.gca()
        plt.axes([0.82,0.12,0.10,0.7],sharey=ax) ; plt.tick_params(labelleft=False,right=True,labelright=True)
        plt.hist(sig_rv,bins=np.arange(0,0.75,0.01),color='k',orientation='horizontal',alpha=0.4)
        plt.hist(sig_rv[mask],bins=np.arange(0,0.75,0.01),color='g',orientation='horizontal',alpha=0.4)
        plt.axes([0.12,0.82,0.7,0.10],sharex=ax) ; plt.tick_params(labelbottom=False,top=True,labeltop=True)
        plt.hist(snr,bins=np.arange(0,1000,25),color='k',alpha=0.4) 
        plt.hist(snr[mask],bins=np.arange(0,1000,25),color='g',alpha=0.4) 

    def plot_survey_stars(self,weather=True, Nb_star=None, Texp=None, overhead=1):
        if weather:
            nb_hours = self.info_SC_nb_hours_per_yr_eff
        else:
            nb_hours = self.info_SC_nb_hours_per_yr

        texp = np.array([5,10,12,15,20])
        nb_exp = nb_hours*60/texp + overhead 

        ti = np.arange(5,20.1,0.1)
        ni = np.arange(20,200,1)
        texp,nstars = np.meshgrid(ti,ni)

        ls = ['-','--','-.',':'][self.simu_counter_survey%4]
        self.simu_counter_survey +=1

        plt.figure('Survey Strategy',figsize=(10,10))
        plt.axes([0.08,0.08,0.86,0.6])
        c = plt.contour(texp,nstars,nb_hours*60/((texp+overhead)*nstars),levels=[25,50,75,100,150,200,250,300])
        plt.clabel(c,fmt='%.0f')
        plt.xlabel('Texp [min]',fontsize=14)
        plt.ylabel('Nb stars []',fontsize=14)
        if self.simu_counter_survey==1:
            plt.grid()
        if Texp is not None:
            plt.axvline(x=Texp,color='k',ls=ls)
            plt.axes([0.08,0.73,0.86,0.25])
            plt.plot(ni,nb_hours*60/((Texp+overhead)*ni),color='k',label='Texp = %.0f min'%(Texp),ls=ls)
            plt.scatter(np.arange(20,201,20),nb_hours*60/((Texp+overhead)*np.arange(20,201,20)),color='k')
            for i in np.arange(20,201,20):
                plt.text(i,nb_hours*60/((Texp+overhead)*i)+10,'%.0f'%(nb_hours*60/((Texp+overhead)*i)),color='k',ha='center')
            if self.simu_counter_survey==1:
                plt.grid()
            plt.ylabel('Nb yearly obs. per *')
            plt.xlabel('Nb star []')
            plt.xlim(20,200) ; plt.ylim(0,365)
            plt.legend()
        if Nb_star is not None:
            plt.axhline(y=Nb_star,color='k',ls=ls)
            plt.axes([0.08,0.73,0.86,0.25])
            plt.plot(ti,nb_hours*60/((ti+overhead)*Nb_star),ls=ls,color='k')
            plt.scatter(np.arange(5,21,2.5),nb_hours*60/((np.arange(5,21,2.5)+overhead)*Nb_star),color='k')
            for i in np.arange(5,21,2.5):
                plt.text(i,nb_hours*60/((i+overhead)*Nb_star)+10,'%.0f'%(nb_hours*60/((i+overhead)*Nb_star)),color='k',ha='center')
            if self.simu_counter_survey==1:
                plt.grid()
            plt.ylabel('Nb yearly obs. per *')
            plt.xlabel('Texp [min]')
            plt.xlim(4,20)

    def func_cutoff(self, cutoff=None, par_space='', par_box=['',''], par_crit=''):
        """example : table_filtered = func_cutoff(table,cutoff1,par_space='Teff&dist',par_box=['4500->5300','0->30'])"""
        GR8 = gr8.copy()
        if cutoff is None:
            cutoff = self.info_TA_cutoff
        for c in cutoff.keys():
            if c[0:8]=='snr_texp':
                texp=float(c[:-1].split('exp')[1])
                snr_texp15 = 0.5*(gr8['snr_420_texp15']+gr8['snr_550_texp15']) #Cretignier et al. +22
                snr = snr_texp15*(np.sqrt(texp/15))
                GR8[c[:-1]] = snr
            if c[0:11]=='sig_rv_texp':
                texp=float(c[:-1].split('exp')[1])
                sig_rv_texp15 = gr8['sig_rv_texp15']
                sig_rv = sig_rv_texp15/(np.sqrt(texp/15))
                GR8[c[:-1]] = sig_rv
            if c[0:16]=='sig_rv_phot_texp':
                texp=float(c[:-1].split('exp')[1])
                sig_rv_texp15 = gr8['sig_rv_phot_texp15']
                sig_rv = sig_rv_texp15/(np.sqrt(texp/15))
                GR8[c[:-1]] = sig_rv

        table_filtered = tcsf.func_cutoff(GR8,cutoff,tagname='',par_space=par_space, par_box=par_box, par_crit=par_crit)
        self.info_TA_my_cutoff = cutoff
        self.info_TA_stars_selected = table_filtered
