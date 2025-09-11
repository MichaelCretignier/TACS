import os

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from colorama import Fore
from scipy.interpolate import interp1d

import THE_TCS_functions as tcsf
import THE_TCS_variables as tcsv

#IMPORT MAIN TABLES

print("""\n[INFO USER] READ ME CAREFULLY 
[INFO USER] The RUWE is currently disabled for stars brighter than mv<5 
[INFO USER] An issue or an upgrade? Contact me at:  michael.cretignier@physics.ox.ac.uk
      """)

cwd = os.getcwd()

Teff_var = 'teff'

db_starname = pd.read_csv(cwd+'/TACS_Material/THE_SIMBAD.csv',index_col=0)
table_time = pd.read_csv(cwd+'/TACS_Material/time_conversion.csv',index_col=0)

seeing = interp1d(table_time['deci'].astype('float')-2026,table_time['seeing'].astype('float'), kind='cubic', bounds_error=False, fill_value='extrapolate')(np.linspace(0,1,365))
downtime = interp1d(table_time['deci'].astype('float')-2026,table_time['downtime'].astype('float'), kind='cubic', bounds_error=False, fill_value='extrapolate')(np.linspace(0,1,365))
month_year = np.ones(365)
for i in [31,59,90,120,151,181,212,243,273,304,334,365]:
    month_year[np.arange(1,366)>i] = month_year[np.arange(1,366)>i]+1

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
db_exoplanets = crossmatch_names(db_exoplanets,'GAIA')
db_exoplanets = db_exoplanets.sort_values(by='PRIMARY').reset_index(drop=True)

db_tess_candidates = pd.read_csv(cwd+'/TACS_Material/TESS_candidates.csv',index_col=0)
db_tess_candidates['name'] = 'TOI'+db_tess_candidates['TOI'].astype('str')
db_tess_candidates['method'] = 'Transit'
db_exoplanets = db_exoplanets.merge(db_tess_candidates,how='outer')

def produce_gr8(version='2.0',verbose=True):
    #GR8 TABLE FORMATION

    gr8_raw = pd.read_csv(cwd+'/TACS_Material/THE_Master_table_v'+version+'.csv',index_col=0)
    gr8 = gr8_raw.copy()

    if verbose:
        if len(gr8)==1112:
            print('[INFO USER] You are using the [PRIVATE] version of the code, do NOT share it outside THE members\n')
        else:
            print('[INFO USER] You are using the [PUBLIC] version of the code. If you are THE member, read the README.\n')

    gr8['vsini_known'] = np.array(gr8['vsini'])
    gr8.loc[gr8['vsini_known']!=gr8['vsini_known'],'vsini_known'] = 30
    gr8.loc[gr8['vsini']!=gr8['vsini'],'vsini'] = 0.0

    gr8['logRHK_known'] = np.array(gr8['logRHK'])
    gr8.loc[gr8['logRHK_known']!=gr8['logRHK_known'],'logRHK_known'] = -4.0
    gr8.loc[gr8['logRHK']!=gr8['logRHK'],'logRHK'] = -6.0

    #based on Andreas comment, would be better to check RVs databases drift
    gr8.loc[gr8['logg']<3.5,'logg'] = 3.5
    gr8.loc[gr8['vmag']<5,'ruwe_GAIA'] = 0.1
    gr8.loc[gr8['ruwe_GAIA']>4,'ruwe_GAIA'] = 4
    gr8.loc[gr8['rv_trend_kms_DACE']>1,'rv_trend_kms_DACE'] = 1
    gr8.loc[gr8['sky_contam_VIZIER']>1,'sky_contam_VIZIER'] = 1
    gr8.loc[gr8['multi_peak_GAIA']>10,'multi_peak_GAIA'] = 10

    gr8.loc[gr8['HWO']==1,'ruwe_GAIA'] = 0.0
    gr8.loc[gr8['HWO']==1,'rv_trend_kms_DACE'] = 0.0

    gr8['SPclass'] = '-'
    gr8.loc[(gr8[Teff_var]>6000),'SPclass'] = 'F'
    gr8.loc[(gr8[Teff_var]>5600)&(gr8[Teff_var]<=6000),'SPclass'] = 'S'
    gr8.loc[(gr8[Teff_var]>5200)&(gr8[Teff_var]<=5600),'SPclass'] = 'G'
    gr8.loc[(gr8[Teff_var]<=5200),'SPclass'] = 'K'
    return gr8, gr8_raw

v1 = produce_gr8('1.0')
v2 = produce_gr8('2.0',verbose=False)

gr8 = {'1.0':v1[0],'2.0':v2[0]}
gr8_raw = {'1.0':v1[1],'2.0':v2[1]}

#FUNCTIONS

def mod_cutoff(cutoff,new_cutoff):
    cutoff2 = cutoff.copy()
    for kw in new_cutoff:
        cutoff2[kw] = new_cutoff[kw]
    return cutoff2

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
    table = table.reset_index(drop=True)
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
            print('\n[INFO] Starnames found:\n')
            print(output,'\n')
            gaia_name = output['GAIA']
            planets = db_exoplanets.loc[db_exoplanets['GAIA']==gaia_name]
            if len(planets):
                print('\n[INFO] Planets found:\n')
                display = planets[['name','method','year','period','mass','radius']]
                display.index = display['name']
                display['period'] = np.round(display['period'],1)
                display['mass'] = np.round(display['mass'],0)
                display['radius'] = np.round(display['radius'],1)
                del display['name']
                print(display,'\n')
        return output
    else:
        print('[WARNING] Starname %s has not been found'%(name))
        return None

def starname_resolver(stars):
    flag = []
    output = []
    for s in np.array(stars):
        loc = resolve_starname(s,verbose=False)
        if loc is not None:
            flag.append(1)
            output.append(loc['INDEX'])
        else:
            flag.append(np.nan)
            output.append(0)
    output = np.array(output)
    flag = np.array(flag)

    names = db_starname.loc[output]
    names[flag!=flag] = np.nan
    return names

def create_gr8like(table,name_col,gr8_column='HD'):
    table = table.dropna(subset=[name_col]).reset_index(drop=True)
    match = starname_resolver(table[name_col])
    table[gr8_column] = np.array(match[gr8_column])
    table = table.dropna(subset=[gr8_column])
    table = table.drop_duplicates(subset=[gr8_column])
    stat_after = len(table)
    output = pd.merge(left=db_starname,right=table,on=gr8_column,how='left')
    print('%.0f stars matched from the 1112 GR8 initial list'%(stat_after))
    return output

def inner_gr8(list_names):
    inner = []
    for s in list_names:
        a = resolve_starname(s,verbose=False)
        if a is not None:
            inner.append(s)
    return inner

def get_starname(entry):
    gaia_ID = np.where(db_starname['GAIA']==entry['GAIA'])[0]
    selected = db_starname.loc[gaia_ID]
    for order in ['HD','HIP','GJ','CSTL','PRIMARY']:
        if selected[order].values[0]!='-':
            break
    return (selected[order].values[0], gaia_ID)

def star_info(entry, format='v1'):
    name, ID = get_starname(entry)
    if format=='v1':
        info = ' ID : %.0f \n Star : %s   Mv = %.2f   Ra = %.2f    Dec = %.2f \n Teff = %.0f   Logg = %.2f   FeH = %.2f    RHK = %.2f   Vsini = %.1f \n RUWE = %.2f   HJ = %.0f   BDW = %.0f   GZ = %.0f   NEP = %.0f   SE = %.0f'%(ID,name,entry['vmag'], entry['ra_j2000'], entry['dec_j2000'], entry[Teff_var], entry['logg'], entry['feh'], entry['logRHK'], entry['vsini'], entry['ruwe_GAIA'], entry['HJ'], entry['BDW'], entry['GZ'], entry['NEP'], entry['SE'])
    else:
        info = ' ID : %.0f   Star : %s   Mv = %.2f   Ra = %.2f    Dec = %.2f \n Teff = %.0f   Logg = %.2f    FeH = %.2f    RHK = %.2f   Vsini = %.1f \n RUWE = %.2f   HJ = %.0f   BDW = %.0f   GZ = %.0f   NEP = %.0f   SE = %.0f'%(ID,name,entry['vmag'], entry['ra_j2000'], entry['dec_j2000'], entry[Teff_var], entry['logg'], entry['feh'], entry['logRHK'], entry['vsini'], entry['ruwe_GAIA'], entry['HJ'], entry['BDW'], entry['GZ'], entry['NEP'], entry['SE'])
    return info

def plot_TESS_CVZ(ra_unit='hours'):
    theta = np.linspace(0,2*np.pi,100)
    if ra_unit=='hours':
        plt.plot(np.cos(theta)*12/360*24+18,np.sin(theta)*12+66,lw=1,ls='-.',color='k')
    else:
        plt.plot(np.cos(theta)*12+18*360/24,np.sin(theta)*12+66,lw=1,ls='-.',color='k')
    plt.text(18,66,'TESS',ha='center',va='center')

def plot_KEPLER_CVZ(ra_unit='hours'):
    theta = np.linspace(0,2*np.pi,100)
    if ra_unit=='hours':
        plt.plot(np.cos(theta)*8/360*24+19.5,np.sin(theta)*8+44.5,lw=1,ls=':',color='k')
    else:
        plt.plot(np.cos(theta)*8+19.5*360/24,np.sin(theta)*8+44.5,lw=1,ls=':',color='k')
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
        for j in [10,30,4000]:
            plt.axhline(y=j,color='C0',alpha=0.5)
        plt.ylabel(r'Planetary mass [$M_{\oplus}$]',fontsize=14)       
        plt.ylim(1,10000)
    plt.axvspan(xmin=60,xmax=400,ymin=0,ymax=0.25,color='g',alpha=0.2) 
    plt.grid()
    plt.scatter(db_exoplanets['period'],db_exoplanets[y_var],c=db_exoplanets['radius'],cmap='brg',vmin=1,vmax=10)
    ax = plt.colorbar(pad=0)
    ax.ax.set_ylabel(r'Planetary radius [$R_{\oplus}$]',fontsize=14)
    plt.subplots_adjust(left=0.10,right=0.95)
    return fig 

def plot_exoplanets2(selection,cutoff={'teff<':6000},mcrit_sup=4000,mcrit_inf=50):
    
    table_gr8 = selection.copy()
    
    for kw in cutoff.keys():
        if kw[-1]=='<':
            table_gr8 = table_gr8.loc[table_gr8[kw[:-1]]<cutoff[kw]]
        else:
            table_gr8 = table_gr8.loc[table_gr8[kw[:-1]]>cutoff[kw]]

    selected = np.in1d(np.array(db_exoplanets['GAIA']),np.array(table_gr8['GAIA']))
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

    nraw = 33
    ncol = len(np.unique(db_exoplanets['GAIA']))//nraw+1
    for j in range(ncol):
        plt.subplot(1,ncol,j+1) ; plt.xscale('log')
        plt.axvspan(xmin=60,xmax=400,color='g',alpha=0.2)
        plt.tick_params(labelleft=False)
        plt.xlim(0.5,99000)
        plt.xlabel('Period [days]',fontsize=13)
    
    db = np.sort(np.unique(db_exoplanets['PRIMARY']))
    db = np.array([np.where(np.array(db_starname['PRIMARY'])==d)[0][0] for d in db])

    print(db_exoplanets)

    summary = []
    for system in db_starname.loc[db,'GAIA']:
        plt.subplot(1,ncol,(abs(count)//nraw)+1)
        plt.axvline(x=10*365.25,color='k',lw=0.5,ls='--')
        count-=1
        mask = np.array(db_exoplanets['GAIA']==system)
        syst = db_exoplanets.loc[mask]
        plt.plot([1,50000],[count,count],lw=1,color='k',ls='-',alpha=0.2)
        plt.scatter(syst['period'],syst['period']*0+count,s=syst['marker'],color='k',zorder=10)
        plt.scatter(syst['period'],syst['period']*0+count,s=syst['marker'],c=syst['radius'],cmap='brg',vmin=1,vmax=8,zorder=10,edgecolor='k')
        condition1 = np.sum((syst['p_eccmin']<400)&(syst['mass']>mcrit_inf)).astype('bool')
        condition2 = np.sum((syst['mass']>mcrit_sup)).astype('bool')
        condition_GZ = np.sum((syst['mass']>30)&(syst['mass']<=mcrit_sup)).astype('int')
        condition_NE = np.sum((syst['mass']>10)&(syst['mass']<=30)).astype('int')
        condition_SE = np.sum((syst['mass']<=10)).astype('int')
        condition_transit = np.sum(syst['radius']==syst['radius'])
        condition_imaging = np.sum(syst['method']=='Imaging')
        if (condition_transit!=0)|(condition_imaging!=0):
            plt.arrow(0.5,count,0.2,0,color='k',head_width=0.1)
        if condition_imaging!=0:
            plt.arrow(0.6,count-0.25,0.0,0.5,color='k',head_width=0.1)

        summary.append([system,int(condition1),int(condition2),int(condition_GZ),int(condition_NE),int(condition_SE),int(condition_transit)])       

        condition_rejected = condition1|condition2
        color_condition = ['k','r'][int(condition_rejected)]
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
    summary = pd.DataFrame(summary,columns=['GAIA','HJ','BDW','GZ','NEP','SE','TRNS'])
    plt.subplots_adjust(left=0.03,right=0.93,wspace=0.30,top=0.96,bottom=0.09)
    return summary

def plot_tess_candidates():
    pass

class tableXY(object):
    def __init__(self, x=None, y=None, yerr=None, xlabel='', ylabel='', ls='o'):
        if x is None:
            x = np.arange(len(y))

        if y is None:
            y = np.zeros(len(x))
        
        if yerr is None:
            yerr = np.zeros(len(y))

        self.x = x
        self.y = y
        self.yerr = yerr
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.ls = ls

    def interpolate(self,new_grid):
        self.y = interp1d(self.x,self.y, kind='cubic', bounds_error=False, fill_value='extrapolate')(new_grid)
        self.x = new_grid

    def monthly_average(self):
        borders = tcsv.month_border
        statistic = []
        for b1,b2 in zip(borders[0:-1],borders[1:]):
            mask = (self.x>=b1)&(self.x<=b2)
            statistic.append(np.mean(self.y[mask]))
        statistic = np.array(statistic)
        return auto_format(statistic)

    def plot(self, alpha=0.5, label=None, ytext=0, figure=None, ls=None, subset=None, color=None):
        if figure is not None:
            if type(figure)==str:
                plt.figure(figure)
        if ls is None:
            ls = self.ls
        if subset is None:
            subset = np.arange(0,len(self.x)).astype('int')
        
        if ls=='o':
            if np.sum(abs(self.yerr))!=0:
                plt.errorbar(self.x[subset],self.y[subset],yerr=self.yerr, capsize=0, marker='o', ls='', alpha=alpha, label=label, color=color)
            else:
                plt.scatter(self.x[subset],self.y[subset],alpha=alpha,label=label,color=color)
        else:
            plt.plot(self.x,self.y,alpha=alpha,label=label,ls=ls,color='k')
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        ax = plt.gca()
        if self.xlabel=='Nights [days]':
            for j,t in zip(tcsv.month_border,['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
                plt.axvline(x=j,color='k',ls='-',alpha=0.5)
                plt.axvline(x=j,color='white',ls=':')
                plt.text(j+15,ytext,t,color='k',ha='center')
            plt.xlim(0,365)

    def create_subset(self,subset):
        self.subset = tableXY(x=self.x[subset], y=self.y[subset], yerr=self.yerr[subset], xlabel=self.xlabel, ylabel=self.ylabel)

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
            self.yerr = self.yerr[selection]
            return None
        else:
            self.create_subset(selection)
            return selection
        
    def split_chunck(self, nb_chunck=1):
        nb = len(self.x)
        liste = np.array_split(np.arange(nb),nb_chunck)
        self.chunck = [[]]*nb_chunck
        for n,subset in enumerate(liste):
            self.chunck[n] = tableXY(x=self.x[subset], y=self.y[subset], yerr=self.yerr[subset], xlabel=self.xlabel, ylabel=self.ylabel)

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

class table_star(object):
    def __init__(self,data):
        self.data = data

    def print(self,columns=[]):
        if len(columns)==0:
            print(self.data)
        else:
            print(self.data[['ra_j2000','dec_j2000','primary_name','vmag']+columns])

    def print_columns(self):
        print(list(self.data.keys()))

    def plot(self, x, y, c=None, s=None, print_names=False, GUI=True, alpha=1.0):
        
        if GUI:
            fig = plt.figure(figsize=(10,10))
            plt.axes([0.1,0.1,0.85,0.75])
        
        dataframe = self.data
        xval = np.array(dataframe[x])
        yval = np.array(dataframe[y])
        if c is None:
            for spflag,color,marker in zip(['K','G','S','F'],['r','C1','gold','cyan'],['o','s','*','x']):
                mask3 = np.array(dataframe['SPclass']==spflag)
                plt.scatter(xval[mask3],yval[mask3],color=color,marker=marker,s=30,zorder=10,alpha=alpha)   
        else:
            if len(c)<3:
                plt.scatter(xval,yval,color=c,marker='o',s=10,zorder=1, alpha=alpha)   

        mask3 = np.array(dataframe['protected']==1)
        if sum(mask3):
            plt.scatter(xval[mask3],yval[mask3],facecolors="none",edgecolors='k',marker='o',s=100,zorder=10,alpha=alpha,)   

        if print_names:
            index = np.array(list(dataframe.index))
            n = np.array(db_starname.loc[index,'HD'])
            for xi,yi,ti in zip(xval,yval,n):
                plt.text(xi,yi,ti,ha='center',va='bottom',fontsize=7)
        plt.xlabel(x,fontsize=14)
        plt.ylabel(y,fontsize=14)

        if GUI:
            ax = plt.gca()
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            info_text = plt.text(xlim[0]+0.03*(xlim[1]-xlim[0]),ylim[1]+0.15*(ylim[1]-ylim[0]),'Double click somewhere',fontsize=13,ha='left',va='top')
            l, = plt.plot([xlim[0]],[ylim[0]],marker='x',color='k',markersize=10,zorder=100)

            class Index(object):
                def __init__(self):
                    self.info_text = ''
                    self.marker = None
                def update(self,newx,newy):
                    dist = abs((xval-newx)/np.nanstd(xval))+abs((yval-newy)/np.nanstd(yval))
                    loc = np.array(dataframe.index)[np.argmin(dist)]
                    new_star = dataframe.loc[loc]
                    text_fmt = star_info(new_star)
                    self.info_text.set_text(text_fmt)
                    self.marker.set_data([new_star[x],new_star[y]])
                    
                    plt.draw()
                    fig.canvas.draw_idle()
                    
            t = Index()
            t.info_text = info_text
            t.marker = l
            
            def onclick(event):
                if event.dblclick:
                    t.update(event.xdata,event.ydata)

            plt.gcf().canvas.mpl_connect('button_press_event', onclick)

class tcs(object):
    
    def __init__(self, sun_elevation=None, starname=None, instrument='HARPS3', verbose=True, method='fast', version='2.0'):    
        self.info_XY_telescope_open = []
        self.info_XY_downtime = tableXY(x=np.arange(365),y=downtime)
        self.simu_SG_calendar = None
        self.simu_counter_survey = 0
        self.simu_tag_survey = ''
        self.info_SC_starname = None
        self.info_SC_instrument = instrument

        GR8 = gr8[version]
        protected = starname_resolver(tcsv.stars_protected)
        protected = protected.dropna(subset=['PRIMARY'])
        GR8['protected'] = 0
        GR8.loc[np.array(protected.index),'protected'] = 1
        self.info_TA_stars_selected = {'GR8':table_star(GR8)}
        
        self.info_TA_cutoff = {}

        self.info_TA_cutoff['presurvey'] = tcsv.cutoff_presurvey
        self.info_TA_cutoff['minimal'] = tcsv.cutoff_minimal

        self.func_cutoff(tagname='minimal',cutoff=tcsv.cutoff_minimal, verbose=False)
        plt.close('cumulative')

        self.func_cutoff(tagname='presurvey',cutoff=tcsv.cutoff_presurvey, verbose=False)
        plt.close('cumulative')

        if type(verbose)!=list:
            verbose = [verbose]*3

        self.random_weather(verbose=verbose[0])

        if sun_elevation is not None:
            self.compute_night_length(sun_elevation=sun_elevation, verbose=verbose[1])

        if starname is not None:
            self.set_star(starname=starname, verbose=verbose[2], method=method)

        sig_ins = 0.0
        if self.info_SC_instrument!='HARPS3':
            sig_ins = 0.50
        self.info_SC_instrument_noise = sig_ins
        self.info_XY_noise_instrument = self.create_ins_noise(sig_ins)

    def create_ins_noise(self, sig_ins=0.00, nb_years=10):
        ins_noise = tableXY(x=np.arange(365*nb_years),y=np.random.randn(365*nb_years)*sig_ins) #generate a 10-years instrument stability
        return ins_noise
    
    def create_star_selection(self,starnames,tagname='my_selection'):
        gr8 = self.info_TA_stars_selected['GR8'].data.copy()
        if type(starnames) is not list:
                starnames = list(starnames)
        gaia_names = []
        for s in starnames:
            output = resolve_starname(s,verbose=False)
            if output is not None:
                gaia_names.append(output['GAIA'])
        gaia_names = np.array(gaia_names)
        mask_star = np.in1d(np.array(gr8['GAIA']),gaia_names)
        self.info_TA_stars_selected[tagname] = table_star(gr8.loc[mask_star])
    

    def compute_night_length(self, sun_elevation=-12, verbose=True, method='approximation'):
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
                print('[INFO] Nb of hours per year : %.0f (%.0f)'%(self.info_SC_nb_hours_per_yr_eff,self.info_SC_nb_hours_per_yr))

        if self.info_SC_starname is not None:
            self.compute_night_star(method=method)

    def random_weather(self,verbose=True):
        output = np.ravel(dropout_year().astype('int'))
        if verbose:
            print('[INFO] Number of bad/good nights = %.0f/%.0f'%(len(output)-sum(output),sum(output)))
        self.info_XY_telescope_open.append(tableXY(y=output,xlabel='Nights [days]',ylabel='Telescope open'))

    def set_star(self,ra=0,dec=0,starname=None,id=None,verbose=True,method='fast'):
        "ra in hours, dec in degree"
        self.info_TA_starnames = None
        gr8 = self.info_TA_stars_selected['GR8'].data.copy()
        if id is not None:
            starname = gr8.iloc[id]['PRIMARY']

        if starname is not None:
            starname = resolve_starname(starname,verbose)
            self.info_TA_starnames = starname
            if starname is not None:
                ra = np.array(gr8.loc[gr8['PRIMARY']==starname['PRIMARY'],'ra_j2000'])[0]/360*24
                dec = np.array(gr8.loc[gr8['PRIMARY']==starname['PRIMARY'],'dec_j2000'])[0]
            else:
                print('[ERROR] STARNAME NOT FOUND')
        self.info_SC_ra = ra
        self.info_SC_dec = dec
        self.info_SC_starname = starname
        self.compute_night_star(method=method)

    def compute_night_star(self, method='fast'):

        alpha = self.info_SC_ra
        dec = self.info_SC_dec
        
        if method=='fast':
            hours,airmass = tcsf.star_observability(alpha, dec, instrument=self.info_SC_instrument, month=1, day=1)
            save2 = np.array([np.roll(airmass,-rolling) for rolling in (np.arange(0,365)/365*len(hours)).astype('int')])
        else:            
            save2 = []
            for j in range(1,13):
                for i in np.arange(1,tcsv.month_len[j-1]+1):
                    hours,airmass = tcsf.star_observability(alpha, dec, instrument=self.info_SC_instrument, month=j, day=i)            
                    save2.append(airmass)
            save2 = np.array(save2)
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

    def compute_nights(self, airmass_max=1.5, weather=False, plot=False, year=0):
        
        map_night = 1-self.info_IM_night.data
        map_star = self.info_IM_airmass.data<airmass_max
        if weather:
            open_dome = self.info_XY_telescope_open[year].y
        else:
            open_dome = np.ones(365)
        map_weather = open_dome*np.ones(1441)[:,np.newaxis]

        total_map = map_night*map_star
        self.info_XY_season = tableXY(
            y = np.sum(total_map,axis=0)>10,
            xlabel='Nights [days]',
            ylabel='Visible'
            ) #15-min exposure
        self.info_XY_night_duration = tableXY(
            y = np.sum(total_map,axis=0)/60,
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
        
    def create_timeseries(self, airmass_max=1.5, nb_year=10, month=None, texp=15, weather=True):

        if weather:
            nb_yr_bn = len(self.info_XY_telescope_open)
            extra = nb_year-nb_yr_bn
            if extra<0:
                extra = 0
            for i in np.arange(extra):
                self.random_weather()

        dt = 24*60/1440 #dt in minutes
        N = int(np.round(texp/dt,0))
        j0 = 0.0#61041.0
        jdb = []
        for year in range(nb_year):
            self.compute_nights(airmass_max=airmass_max, 
                                weather=weather, 
                                plot=False,
                                year=year)
            epochs = []
            l0 = np.median(np.where(self.info_IM_observable.data==1)[0])
            dt = int(l0-0.5*1440)
            maps = np.roll(self.info_IM_observable.data.T,-dt,axis=1)
            for d,n in enumerate(maps):
                loc = np.where(n==1)[0]  
                if len(loc)>N:                            
                    loc = loc[int(N/2):-int(N/2)]
                    mini = loc[0]
                    maxi = loc[-1]
                    t0 = np.mean(loc)+np.arange(-100,100,1)*N
                    t0 = t0[(t0>mini)&(t0<maxi)]
                    t0 += dt
                    epochs.append(j0+d+0.5+(t0/1441-0.5)+365*year)
            jdb.append(np.concatenate(epochs))
            #jdb.append(epochs)
        #jdb = np.hstack(jdb)
        jdb = np.hstack(jdb)
        rv = np.random.randn(len(jdb))
        
        if month is not None:
            days = np.where(month_year==month)[0]
            mask = np.in1d(jdb.astype('int'),days)
            jdb = jdb[mask]
            rv = rv[mask]

        self.info_XY_timestamps = tableXY(x=jdb, y=rv, xlabel='Nights [days]')
        

    def create_survey(self, selection, airmass_max=1.5, nb_year=10, texp=15, weather=True, nb_subexp=None):
        if selection not in self.info_TA_stars_selected.keys():
            print(' [ERROR] This selection %s is not part of the selections:',list(self.info_TA_stars_selected.keys()))
        else:
            table = self.info_TA_stars_selected[selection]
            timeseries = []
            for ID in np.array(table.index):
                self.set_star(id=ID)
                self.create_timeseries(airmass_max=airmass_max, nb_year=nb_year, texp=texp, weather=weather, nb_subexp=nb_subexp)
                timeseries.append(self.info_XY_timestamps)


    def compute_SG_calendar(
            self,
            sun_elevation=-12,
            airmass_max=1.8,
            alpha_step=1, 
            dec_step=1,
            cutoff=None,
            selection='presurvey'):
        
        gr8 = self.info_TA_stars_selected['GR8'].data.copy()

        if selection is not None:
            table_gr8 = self.info_TA_stars_selected[selection].data
            self.info_TA_cutoff['SG'] = cutoff
            self.info_TA_stars_selected['SG'] = table_star(table_gr8.copy())
        else:
            if cutoff is None:
                cutoff = self.info_TA_cutoff['presurvey']
            
            table_gr8 = self.info_TA_stars_selected['GR8'].data.copy()
            for kw in cutoff.keys():
                if kw[-1]=='<':
                    table_gr8 = table_gr8.loc[table_gr8[kw[:-1]]<cutoff[kw]]
                else:
                    table_gr8 = table_gr8.loc[table_gr8[kw[:-1]]>cutoff[kw]]

            self.info_TA_cutoff['SG'] = cutoff
            self.info_TA_stars_selected['SG'] = table_star(table_gr8.copy())

        self.compute_night_length(sun_elevation=sun_elevation) 
        
        button = 1
        if self.simu_SG_calendar is not None:
            if (self.simu_SG_calendar['param1']==sun_elevation)&(self.simu_SG_calendar['param2']==airmass_max)&(self.simu_SG_calendar['param3']==alpha_step)&(self.simu_SG_calendar['param4']==dec_step):
                button = 0
                print('[INFO] Old simulation found.')
        
        if button==1:
            output = []
            params = []
            RA, DEC = np.meshgrid(np.arange(0,30,alpha_step),np.arange(-30,90,dec_step))
            loading = np.round(len(np.ravel(RA))*np.arange(0,101,10)/100,0).astype('int')
            counter=0
            for i,j in zip(np.ravel(RA),np.ravel(DEC)):
                if counter in loading:
                    print('[INFO] Progress... [%.0f%%]'%(np.where(loading==counter)[0][0]*10))
                self.set_star(ra=i,dec=j) 
                self.compute_nights(airmass_max=airmass_max, weather=False, plot=False)
                params.append([i,j])
                output.append(self.info_XY_night_duration.monthly_average())
                counter+=1
            print('[INFO] Finished!')
            print('[INFO] Stacking tables...')
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
        
        print('[INFO] Producing Calendar plot... Wait...')
        fig = plt.figure(figsize=(16,10))
        fig.suptitle('Sun elevation = %.0f | Airmass max = %.2f'%(sun_elevation,airmass_max))
        plt.subplots_adjust(left=0.05,right=0.98,top=0.94,bottom=0.05,hspace=0.30,wspace=0.30)

        downtime = np.array([32., 31., 31., 27., 15.,  4.,  3.,  8., 19., 30., 37., 38., 32.])

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
            plt.scatter(gr8['ra_j2000']/360*24,gr8['dec_j2000'],s=(7.5-gr8['vmag'])*30,alpha=0.15,c=gr8[Teff_var],cmap='jet_r',vmin=5000,vmax=6000)
            plt.scatter(gr8['ra_j2000']/360*24+24,gr8['dec_j2000'],s=(7.5-gr8['vmag'])*30,alpha=0.15,c=gr8[Teff_var],cmap='jet_r',vmin=5000,vmax=6000)
            plt.scatter(table_gr8['ra_j2000']/360*24,table_gr8['dec_j2000'],s=(7.5-table_gr8['vmag'])*30,c=table_gr8[Teff_var],cmap='jet_r',vmin=5000,vmax=6000,ec='k')
            plt.scatter(table_gr8['ra_j2000']/360*24+24,table_gr8['dec_j2000'],s=(7.5-table_gr8['vmag'])*30,c=table_gr8[Teff_var],cmap='jet_r',vmin=5000,vmax=6000,ec='k')
            plot_TESS_CVZ()
            plot_KEPLER_CVZ()

    def compute_SG_month(self, month=1, plot=False, selection='SG'):
        
        params,output,RA,DEC = self.simu_SG_calendar['outputs']
        sun_elevation = self.simu_SG_calendar['param1']
        airmass_max = self.simu_SG_calendar['param2']

        ra = np.ravel(RA)
        dec = np.ravel(DEC)

        month_tag = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'][month-1]
        gr8 = self.info_TA_stars_selected['GR8'].data

        #add the info column
        for selections in ['presurvey','minimal',selection]:
            table = self.info_TA_stars_selected[selections].data
            dist = abs(np.array(table['ra_j2000'])/360*24-ra[:,np.newaxis])+abs(np.array(table['dec_j2000'])-dec[:,np.newaxis])
            loc = np.argmin(dist,axis=0)
            table['night_length_%s'%(month_tag)] = output[loc][:,month-1]
            self.info_TA_stars_selected[selections] = table_star(table.copy())

        if plot:
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
            plt.scatter(gr8['ra_j2000']/360*24,gr8['dec_j2000'],s=(7.5-gr8['vmag'])*30,alpha=0.15,c=gr8[Teff_var],cmap='jet_r',vmin=5000,vmax=6000)
            plt.scatter(gr8['ra_j2000']/360*24+24,gr8['dec_j2000'],s=(7.5-gr8['vmag'])*30,alpha=0.15,c=gr8[Teff_var],cmap='jet_r',vmin=5000,vmax=6000)
            plt.scatter(table['ra_j2000']/360*24,table['dec_j2000'],s=(7.5-table['vmag'])*30,c=table[Teff_var],cmap='jet_r',vmin=5000,vmax=6000,ec='k')
            plt.scatter(table['ra_j2000']/360*24+24,table['dec_j2000'],s=(7.5-table['vmag'])*30,c=table[Teff_var],cmap='jet_r',vmin=5000,vmax=6000,ec='k')
            plt.colorbar(pad=0)  
            plt.subplots_adjust(left=0.05,right=1.1)          

            info_text = plt.text(0,107,'Double click somewhere',fontsize=13,ha='left',va='top')
            l, = plt.plot([-5],[130],marker='x',color='k',markersize=10)

            class Index(object):
                def __init__(self):
                    self.info_text = ''
                    self.marker = None
                def update(self,newx,newy):
                    new_star = query_table(newx/24*360,newy,table)
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
        gr8 = self.info_TA_stars_selected['GR8'].data.copy()
        gaia_name = self.info_TA_starnames['GAIA']
        found = np.where(np.array(db_exoplanets['GAIA'])==gaia_name)[0]
        if len(found):
            init_text = star_info(gr8.iloc[self.info_TA_starnames['INDEX']],format='v2')
            l, = plt.plot(
                np.array(db_exoplanets.loc[found,'period']),
                np.array(db_exoplanets.loc[found,y_var]),
                marker='x',color='k',markersize=10)
        else:
            init_text = 'Double click somewhere'
            l, = plt.plot([0.3],[1000],marker='x',color='k',markersize=10)
        
        info_text = plt.text(0.5,{'k':3000,'mass':30000}[y_var],init_text,fontsize=13,ha='left',va='top')

        class Index(object):
            def __init__(self):
                self.info_text = ''
                self.marker = None
            def update(self,newx,newy):
                loc = np.argmin(abs(np.log10(db_exoplanets[y_var])-np.log10(newy))+abs(np.log10(db_exoplanets['period'])-np.log10(newx)))
                info = resolve_starname(db_exoplanets.loc[loc,'GAIA'],verbose=False)
                new_star = gr8.iloc[info['INDEX']]
                text_fmt = star_info(new_star,format='v2')
                self.info_text.set_text(text_fmt)
                exo = db_exoplanets.loc[db_exoplanets['GAIA']==info['GAIA']]
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

    def compute_exoplanet_rv_signal(self, jdb=None, keplerian_par={}, y0=2026, photon_noise=0.0):
        if len(keplerian_par)==0:
            gaia_name = self.info_TA_starnames['GAIA']
            mask = np.array(db_exoplanets['GAIA']==gaia_name)
        else:
            mask = np.array(db_exoplanets['GAIA']=='')
        if jdb is None:
            jdb = self.info_XY_timestamps.x

        self.info_XY_keplerian_noise = []

        phot_noise = np.random.randn(len(jdb))*photon_noise
        ins_noise = self.info_XY_noise_instrument.y[(jdb-int(np.min(jdb))).astype('int')]
        red_noise = np.zeros(len(jdb))

        self.info_XY_keplerian_noise.append(tableXY(x=jdb,y=phot_noise))
        self.info_XY_keplerian_noise.append(tableXY(x=jdb,y=ins_noise))
        self.info_XY_keplerian_noise.append(tableXY(x=jdb,y=red_noise))
        noise_tot = phot_noise + ins_noise + red_noise
        noise_err = np.sqrt(photon_noise**2+self.info_SC_instrument_noise**2)

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
            self.info_XY_keplerian[0] = tableXY(x=jdb-j0, y=rv_tot, xlabel=xlabel,ls='o', ylabel='RV Kep tot [m/s]')
            rv_tot2 = np.sum([k.y for k in self.info_XY_keplerian_model[1:]],axis=0)
            self.info_XY_keplerian_model[0] = tableXY(x=jdb_model-j0, y=rv_tot2, xlabel=xlabel,ls='-', ylabel='RV tot [m/s]')
        
            self.info_XY_keplerian.append(tableXY(x=jdb-j0, y=rv_tot+noise_tot, yerr=np.ones(len(jdb))*noise_err, xlabel=xlabel,ls='o', ylabel='RV tot [m/s]'))
            self.info_XY_keplerian_model.append(tableXY(x=jdb_model-j0, y=rv_tot2, xlabel=xlabel,ls='-', ylabel='RV tot [m/s]'))

        else:
            self.info_TA_exoplanets_known = []
            print('[INFO] No exoplanets found for this target.')

    def plot_keplerians(self,axhline=None, obs_per_night=None, random=False):
        nb = len(self.info_XY_keplerian)
        fig = plt.figure(figsize=(18,10))
        if self.info_TA_starnames is not None:
            fig.suptitle('%s | %s | %s | %s'%(self.info_TA_starnames['PRIMARY'], self.info_TA_starnames['HD'], self.info_TA_starnames['HIP'], self.info_TA_starnames['GAIA']))
        
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


    def plot_night_length(self,figure='NightLength',legend=True):
        backup = np.array([self.info_SC_night_def]).copy()
        
        self.compute_night_length(sun_elevation=-12, verbose=False) 
        self.compute_nights(airmass_max=1.5, weather=False, plot=False)
        self.info_XY_night_duration.plot(figure=figure,label='Z=1.5 | S=-12',ytext=-0.5) 
        
        self.compute_nights(airmass_max=1.8, weather=False, plot=False)
        self.info_XY_night_duration.plot(figure=figure,label='Z=1.8 | S=-12',ytext=-0.5) 
        
        self.compute_night_length(sun_elevation=-6, verbose=False) #change the sunset/rise parameter
        self.compute_nights(airmass_max=1.8, weather=False, plot=False)
        self.info_XY_night_duration.plot(figure=figure,label='Z=1.8 | S=-6',ytext=-0.5)
        if legend:
            plt.legend()
        self.compute_night_length(sun_elevation=backup[0], verbose=False) 
        plt.ylim(-1)

    def compute_optimal_texp(self, selection=None, snr=250, sig_rv=0.30, texp_crit=20, budget='_phot'):
        """ budget = '_arve_osc+gr' """
        
        if snr<1:
            snr=1
        
        if sig_rv<0.001:
            sig_rv=100
        
        if selection is None:
            selection = self.info_TA_stars_selected['GR8'].data.copy()
        else:
            selection = self.info_TA_stars_selected[selection].data

        snr_texp15 = np.array(0.5*(selection['snr_420_texp15']+selection['snr_550_texp15'])) #Cretignier et al. +22
        texp_snr_crit = 15*(snr/snr_texp15)**2

        if budget!='_phot':
            texp_tabulated = np.array([1,5,8,10,12,15,20,25,30])
            kws = ['sig_rv'+budget+'_texp%.0f'%(j) for j in texp_tabulated]
            sig_rvi = np.array(selection[kws])<=sig_rv
            texp_sig_rv_crit = []
            for s in sig_rvi:
                loc = np.where(s)[0]
                if len(loc):
                    texp_sig_rv_crit.append(texp_tabulated[loc[0]])
                else:
                    texp_sig_rv_crit.append(45)
            texp_sig_rv_crit = np.array(texp_sig_rv_crit)
        else:
            sig_rv_texp15 = np.array(selection['sig_rv_phot_texp15'])
            texp_sig_rv_crit = 15*(sig_rv_texp15/sig_rv)**2

        optimal_time = np.max([texp_snr_crit,texp_sig_rv_crit],axis=0)
        statistic = np.argmax([texp_snr_crit,texp_sig_rv_crit],axis=0)

        selection['texp_optimal'] = optimal_time

        plt.figure(figsize=(18,6))
        plt.axes([0.03,0.33,0.2,0.8])
        plt.pie([np.sum(statistic==0),np.sum(statistic==1)],labels=['SNR \nlimited','Sig RV \nlimited'], autopct='%.0f%%')

        plt.axes([0.32,0.1,0.65,0.8])
        plt.scatter(np.arange(len(selection)),texp_snr_crit,label=r'SNR > %.0f'%(snr))
        plt.scatter(np.arange(len(selection)),texp_sig_rv_crit,label=r'$\sigma_{RV}$ (%s) < %.2f'%(budget[1:],sig_rv))
        plt.scatter(np.arange(len(selection)),optimal_time,color='k',label='Optimal',marker='.')
        plt.legend()
        plt.ylabel('Texp [min]') ; plt.ylim(0,50)
        plt.xlabel('Star ID')
        plt.axhline(y=texp_crit,color='r',ls=':')
        plt.title('Nb stars valid = %.0f / %.0f'%(np.sum(optimal_time<=texp_crit),len(selection)))

        plt.axes([0.05,0.1,0.2,0.35])
        plt.hist(texp_snr_crit,bins=np.arange(0,46,1),color='C0',alpha=0.4,label=r'SNR > %.0f'%(snr))
        plt.hist(texp_sig_rv_crit,bins=np.arange(0,46,1),color='C1',alpha=0.4,label=r'$\sigma_{RV}$ (%s) < %.2f'%(budget[1:],sig_rv))
        plt.hist(optimal_time,bins=np.arange(0,46,1),color='k',alpha=0.4,label='Optimal')
        plt.xlabel('Texp [min]') ; plt.xlim(0,50)

    def compute_nb_nights_required(self, selection='', texp=15, month=1):
        """Only thing missing is to check the number of star with season gap"""
        tab = self.info_TA_stars_selected[selection].data



        nights = (1-self.info_IM_night.data)
        observing = 1 - downtime/100
        hours = np.sum(nights*observing,axis=0)/60*0.60 #60% of GTO

        day0 = np.where(month_year==month)[0][0]

        cumu_hours = np.cumsum(np.roll(hours,-day0))

        overhead = 1
        if texp!='optimal':
            total_time = (texp+overhead)*len(tab)/60
        else:
            total_time = (np.sum(tab['texp_optimal'])+len(tab)*overhead)/60
        
        nb_days = np.where(cumu_hours>total_time)[0][0]
        print('[INFO] Total time needed = %.1f hours'%(total_time))
        print('[INFO] Number of days needed : %.0f'%(nb_days))
        

    def plot_survey_snr_texp(self, selection=None, texp=None, snr_crit=250, sig_rv_crit=0.30, budget='_phot'):
        gr8 = self.info_TA_stars_selected['GR8'].data.copy()
        if selection is None:
            selection = self.info_TA_stars_selected['GR8'].data.copy()
        else:
            selection = self.info_TA_stars_selected[selection].data

        #snr_texp15 = gr8['snr_550_texp15']
        texp_snr250 = selection['texp_snr_250']
        if budget=='_phot':
            sig_rv_texp15 = selection['sig_rv'+budget+'_texp15']
            sig_rv = sig_rv_texp15/(np.sqrt(texp/15))
            sig_rv_gr8 = gr8['sig_rv'+budget+'_texp15']/(np.sqrt(texp/15))
        else:
            texp = tcsf.find_nearest(np.array([1,5,8,10,12,15,20,25,30]),texp)[1][0]
            print('[INFO] Closest Texp tabulated = %.0f min'%(texp))
            sig_rv = selection['sig_rv'+budget+'_texp%.0f'%(texp)]
            sig_rv_gr8 = gr8['sig_rv'+budget+'_texp%.0f'%(texp)]
        snr = (np.sqrt(texp/15))*selection['snr_C22_texp15'] #Cretignier et al. +22

        mask = (snr>snr_crit)&(sig_rv<sig_rv_crit)
        plt.figure('SNR_SIGRV_TEXP%.0fMIN_%s'%(texp,budget.split('_')[-1]),figsize=(8,8))
        plt.axes([0.12,0.12,0.7,0.7])
        for t1,t2,color,marker in zip([4000,5200,5600,6000],[5200,5600,6000,7000],['r','C1','gold','cyan'],['o','s','*','x']):
            mask2 = np.array((selection[Teff_var]>t1)&(selection[Teff_var]<=t2))
            mask3 = np.array((gr8[Teff_var]>t1)&(gr8[Teff_var]<=t2))
            plt.scatter((np.sqrt(texp/15))*gr8.loc[mask3,'snr_C22_texp15'],sig_rv_gr8.loc[mask3],color='gray',marker=marker,s=30,alpha=0.5)                
            plt.scatter(snr[mask2],sig_rv[mask2],color='k',marker=marker,s=30)    

        plt.scatter(snr[mask],sig_rv[mask],color='g',marker='o',s=7,label='All (%.0f / %.0f)'%(sum(mask),len(mask)))
        for t1,t2,color,marker in zip([4000,5200,5600,6000],[5200,5600,6000,7000],['r','C1','gold','cyan'],['o','s','*','x']):
            mask2 = np.array((selection[Teff_var]>t1)&(selection[Teff_var]<=t2))
            plt.scatter(snr[mask&mask2],sig_rv[mask&mask2],color=color,marker=marker,s=7,label='T = [%.0fK - %.0fK] | (%.0f / %.0f)'%(t1,t2,sum(mask&mask2),sum(mask2)))

        plt.axvline(x=snr_crit,color='k',ls=':')
        plt.axhline(y=sig_rv_crit,color='k',ls=':')
        plt.axvspan(xmin=0,xmax=snr_crit,color='k',alpha=0.2)
        plt.axhspan(ymin=sig_rv_crit,ymax=1.0,color='k',alpha=0.2)
        plt.xlabel('SNR continuum',fontsize=14)
        plt.ylabel(r'$\sigma_{{\gamma}}$ $RV$ [m/s]',fontsize=14)
        plt.legend(loc=1,markerscale=2.0)
        plt.xlim(0,1299)
        plt.ylim(0,1.00)
        ax = plt.gca()
        plt.axes([0.82,0.12,0.10,0.7],sharey=ax) ; plt.tick_params(labelleft=False,right=True,labelright=True)
        plt.hist(sig_rv,bins=np.arange(0,1.00,0.01),color='k',orientation='horizontal',alpha=0.4)
        plt.hist(sig_rv[mask],bins=np.arange(0,1.00,0.01),color='g',orientation='horizontal',alpha=0.4)
        plt.axes([0.12,0.82,0.7,0.10],sharex=ax) ; plt.tick_params(labelbottom=False,top=True,labeltop=True)
        plt.hist(snr,bins=np.arange(0,1000,10),color='k',alpha=0.4) 
        plt.hist(snr[mask],bins=np.arange(0,1000,10),color='g',alpha=0.4) 

    
    def plot_survey_stars(self,weather=True, Nb_star=None, Texp=None, Nb_obs_per_year=None, overhead=1, selection=None, ranking='HZ_mp_min_osc+gr_texp15', color='k'):
        if weather:
            nb_hours = self.info_SC_nb_hours_per_yr_eff
        else:
            nb_hours = self.info_SC_nb_hours_per_yr

        texp = np.array([5,8,10,12,15,20,25,30])
        nb_exp = nb_hours*60/texp + overhead 

        ti = np.arange(4,30.1,0.1)
        ni = np.arange(20,200,1)
        texp,nstars = np.meshgrid(ti,ni)

        code = ''
        if (Texp is not None)&(selection is None):
            code='_T'
        if Nb_star is not None:
            code='_N'

        if self.simu_tag_survey!=code:
            self.simu_counter_survey = 0
        self.simu_tag_survey = code

        ls = ['-','--','-.',':'][self.simu_counter_survey%4]

        self.simu_counter_survey = self.simu_counter_survey+1

        if selection is not None:
            table = self.info_TA_stars_selected[selection].data
            table = table.sort_values(by=ranking)
            total_time = self.info_SC_nb_hours_per_yr_eff*60 # total time per year in min
            if Texp is not None:
                texp = np.ones(len(table))*(Texp+overhead)
                label= 'Texp=%.0fmin'%(Texp)
            else:
                texp = np.array(table['texp_optimal']+overhead)
                label = ranking

            texp_cumu = np.cumsum(texp)

            nb_obs_per_year = total_time/texp_cumu
            plt.figure('Survey_Strategy')
            plt.title('Total number of stars in the selection = %.0f'%(len(texp_cumu)))
            plt.plot(np.arange(20,len(texp_cumu)),nb_obs_per_year[20:],color=color,ls=ls,label=label)
            for j in range(20,len(texp),10):
                plt.scatter(j,nb_obs_per_year[j],color='k',marker='o')
                plt.text(j,nb_obs_per_year[j],'%.0f'%(nb_obs_per_year[j]),ha='left',va='bottom')
                plt.xlabel('Number of stars')
                plt.ylabel('Number of obs. per year')
            plt.legend()
        else:
            plt.figure('Survey_Strategy'+code,figsize=(10,10))
            plt.axes([0.08,0.08,0.86,0.6])
            c = plt.contour(texp,nstars,nb_hours*60/((texp+overhead)*nstars),levels=[25,50,75,100,150,200,250,300])
            plt.clabel(c,fmt='%.0f')
            plt.xlabel('Texp [min]',fontsize=14)
            plt.ylabel('Nb stars []',fontsize=14)
            if self.simu_counter_survey==1:
                plt.grid()
            if Texp is not None:
                plt.axvline(x=Texp,color='k',ls=ls,lw=3)
                plt.axes([0.08,0.73,0.86,0.25])
                plt.plot(ni,nb_hours*60/((Texp+overhead)*ni),color=color,label='Texp = %.0f min'%(Texp),ls=ls)
                plt.scatter(np.arange(20,401,20),nb_hours*60/((Texp+overhead)*np.arange(20,401,20)),color=color)
                for i in np.arange(20,401,20):
                    plt.text(i,nb_hours*60/((Texp+overhead)*i)+10,'%.0f'%(nb_hours*60/((Texp+overhead)*i)),color='k',ha='center')
                if self.simu_counter_survey==1:
                    plt.grid()
                plt.ylabel('Nb yearly obs. per *')
                plt.xlabel('Nb star []')
                plt.xlim(20,200) ; plt.ylim(0,365)
                plt.legend()
                self.info_XY_survey_stat = tableXY(x=ni,y=nb_hours*60/((Texp+overhead)*ni))
            if Nb_star is not None:
                plt.axhline(y=Nb_star,color='k',ls=ls,lw=3)
                plt.axes([0.08,0.73,0.86,0.25])
                plt.plot(ti,nb_hours*60/((ti+overhead)*Nb_star),ls=ls,color=color)
                plt.scatter(np.arange(5,31,2.5),nb_hours*60/((np.arange(5,31,2.5)+overhead)*Nb_star),color=color)
                for i in np.arange(5,31,2.5):
                    plt.text(i,nb_hours*60/((i+overhead)*Nb_star)+10,'%.0f'%(nb_hours*60/((i+overhead)*Nb_star)),color='k',ha='center')
                if self.simu_counter_survey==1:
                    plt.grid()
                plt.ylabel('Nb yearly obs. per *')
                plt.xlabel('Texp [min]')
                plt.xlim(4,30)
                self.info_XY_survey_stat = tableXY(x=ti,y=nb_hours*60/((ti+overhead)*Nb_star))
            if Nb_obs_per_year is not None:
                c2 = plt.contour(texp,nstars,nb_hours*60/((texp+overhead)*nstars),levels=[0,Nb_obs_per_year],cmap='Greys',linestyles=ls,linewidths=3) 
                plt.clabel(c2,fmt='%.0f')
                plt.axes([0.08,0.73,0.86,0.25])
                plt.plot(ti,nb_hours*60/((ti+overhead)*Nb_obs_per_year),ls=ls,color=color)
                plt.scatter(np.arange(5,31,2.5),nb_hours*60/((np.arange(5,31,2.5)+overhead)*Nb_obs_per_year),color=color)
                for i in np.arange(5,31,2.5):
                    plt.text(i,nb_hours*60/((i+overhead)*Nb_obs_per_year)+10,'%.0f'%(nb_hours*60/((i+overhead)*Nb_obs_per_year)),color='k',ha='center')
                if self.simu_counter_survey==1:
                    plt.grid()
                plt.ylabel('Nb stars')
                plt.xlabel('Texp [min]')
                plt.xlim(4,30)
                self.info_XY_survey_stat = tableXY(x=ti,y=nb_hours*60/((ti+overhead)*Nb_obs_per_year))

    def which_cutoff(self,starname,cutoff=None,tagname=None,plot=False):
        gr8 = self.info_TA_stars_selected['GR8'].data.copy()
        if tagname is not None:
            try:
                cutoff = self.info_TA_cutoff[tagname]
            except:
                print('[ERROR] this tagname is not found, current list is: ',list(self.info_TA_cutoff.keys()))
        
        if type(starname)==str:
            starname = [starname]
        starname = np.array(starname)

        outputs = []
        for s in starname:
            index = resolve_starname(s,verbose=(len(starname)==1))
            output = []
            if index is not None:
                star = gr8.loc[index['INDEX']]
                for kws in cutoff.keys():
                    kw = kws[:-1]
                    condition = kws[-1]
                    value = cutoff[kws]
                    if condition=='>':
                        test = int(star[kw]>value)
                    else:
                        test = int(star[kw]<value)
                    highlight=0
                    if (kw=='protected')&(star[kw]==1):
                        highlight=1
                    if (kw=='gmag'):
                        highlight = np.sum(np.array([7.25, 5.75, 5.25, 4.5])>star[kw])-1
                    if (kw=='HWO')&(star[kw]==1.0):
                        highlight=1
                    if (kw=='nobs_DB'):
                        highlight = np.sum(np.array([-1, 130, 250, 400])<star[kw])-1
                    if (kw=='season_length_1.75'):
                        highlight = np.sum(np.array([225,275,290,305])<star[kw])-1
                    if (kw=='season_length_1.5'):
                        highlight = np.sum(np.array([210,260,275,285])<star[kw])-1
                    output.append([s,['--->',''][test],kw,condition,value,star[kw],['FALSE','TRUE'][test],['','❂','❂❂','❂❂❂','X'][highlight],['<---',''][test]])
                output = pd.DataFrame(output,columns=['starname','!','feature','condition','threshold','value','test','❂','!!'])
                output['value'] = np.round(output['value'],2)
                if len(starname)==1:
                    protection = int(np.array(output.loc[output['feature']=='protected','value']==1)[0])
                    if (sum(output['test']=='FALSE')!=0)&(protection==0):
                        print(Fore.RED+'[INFO] -- NO -- %s was rejected.\n'%(s)+Fore.RESET)
                    elif (sum(output['test']=='FALSE')!=0)&(protection==1):
                        print(Fore.YELLOW+'[INFO] -- NO -- %s was rejected (but is now protected!).\n'%(s)+Fore.RESET)                
                    else:
                        print(Fore.GREEN+'[INFO] -- YES -- %s is still selected.\n'%(s)+Fore.RESET)
                    print(output[output.columns[1:]])
            else:
                output.append([s,'--->','starname','!=','UNFOUND',s,'FALSE','<---'])
                output = pd.DataFrame(output,columns=['starname','!','feature','condition','threshold','value','test','!!'])
            if sum(output['test']=='FALSE'):
                outputs.append(output.loc[output['test']=='FALSE'])
        
        pd.set_option("display.max_rows", None)
        if len(starname)!=1:
            outputs = pd.concat(outputs)
            print(outputs)
            rejected = len(np.unique(outputs['starname']))
            total = len(starname)
            missing = len(np.unique(outputs.loc[outputs['feature']=='starname','starname']))
            print('\n[INFO] %.0f (%.0f%%) stars rejected from the %.0f stars in the list (%.0f stars not in GR8, %.0f rejected).\n'%(rejected,100*rejected/total,total,missing,rejected-missing))
            if plot:
                plt.figure()
                outputs['feature'].value_counts().plot.pie(autopct="%1.0f%%")
        pd.reset_option("display.max_rows")

        self.info_TA_stars_missing = outputs

    def func_cutoff(self, tagname='handmade', tagname_fig='', cutoff=None, par_space='', par_box=['',''], par_crit='', verbose=True, show_sample=None, protection=True):
        """example : table_filtered = func_cutoff(table,cutoff1,par_space='Teff&dist',par_box=['4500->5300','0->30'])"""
        GR8 = self.info_TA_stars_selected['GR8'].data.copy()
        if show_sample is not None:
            GR8 = GR8.loc[GR8['SPclass']==show_sample]

        if protection is False:
            GR8['protected'] = 0

        if cutoff is None:
            cutoff = self.info_TA_cutoff['presurvey']
            tagname = 'presurvey'
        
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

        if par_space!='':
            tagname_fig = par_space
        if par_box[0]!='':
            tagname_fig=par_box[0]
        
        table_filtered = tcsf.func_cutoff(GR8,cutoff,tagname=tagname_fig,par_space=par_space, par_box=par_box, par_crit=par_crit, verbose=verbose)
        if tagname!='dustbin':
            self.info_TA_cutoff[tagname] = cutoff
            self.info_TA_stars_selected[tagname] = table_star(table_filtered.copy())


    def cutoff_ST(self):
        for tagname,cutoff in zip(['Tim','Jean','Sam1','Sam2','Miku','William1','William2','Stefano'],[tcsv.cutoff_tim,tcsv.cutoff_jean,tcsv.cutoff_sam,tcsv.cutoff_sam2,tcsv.cutoff_mick,tcsv.cutoff_william1,tcsv.cutoff_william2,tcsv.cutoff_stefano]):
            self.func_cutoff(tagname=tagname,cutoff=cutoff)



    def create_table_scheduler(self, selection, year=2026, month_obs_baseline=12, texp=900, n_obs='auto', freq_obs=None, ranking='HZ_mp_min_osc+gr_texp15', tagname='', plot_ranking_priority=False, plot_real_ID=False, need_help=False):

        if type(selection)==str:
            table_scheduler = self.info_TA_stars_selected[selection].data.copy()
        else:
            table_scheduler = selection.copy()
        table_scheduler['GR8_ID'] = table_scheduler.index

        if tagname=='':
            tagname = '_'+selection

        table_scheduler = table_scheduler.sort_values(by='ra_j2000').reset_index(drop=True)

        table_scheduler['priority'] = 6
        if ranking is not None:
            ranking = table_scheduler[ranking]
            table_scheduler.loc[ranking<np.nanpercentile(ranking,33),'priority'] = 3
            table_scheduler.loc[ranking>np.nanpercentile(ranking,66),'priority'] = 9
        else:
            table_scheduler['priority'] = 9

        table_scheduler['expN'] = 1

        table_scheduler['groupEnableTime'] = np.nan
        table_scheduler['groupDisableTime'] = np.nan

        table_scheduler['schedulingMode'] = 'MONITORING'
        table_scheduler['t0'] = 'NULL'
        table_scheduler['period'] = 1
        table_scheduler['delta'] = 0.3
        table_scheduler['moonMaxFI'] = 0
        table_scheduler['minMoonDist1'] = 5
        table_scheduler['ifExceedsMinFI1'] = 15
        table_scheduler['minMoonDist2'] = 0
        table_scheduler['ifExceedsMinFI2'] = 0
        table_scheduler['rvMoonRange'] = 10
        table_scheduler['rvMoonMinFI'] = 0.5
        table_scheduler['rvMoonMinAlt'] = -5
        table_scheduler['twilight'] = -12
        table_scheduler['maxSeeing'] = 2
        table_scheduler['maxExtinction'] = 1.1
        table_scheduler['maxSkyBrightness'] = 15
        table_scheduler['minAltitude'] = 35
        table_scheduler['acqType'] = 'OBJECT'
        table_scheduler['GDR3_ID_number'] = np.array([i.split(' ')[-1] for i in table_scheduler['GAIA']])

        if plot_ranking_priority:
            table_scheduler = table_scheduler.sort_values(by=['priority','ra_j2000']).reset_index(drop=True)            
        else:
            table_scheduler = table_scheduler.sort_values(by=['ra_j2000']).reset_index(drop=True)

        tyr_rise = np.array(table_scheduler['tyr_rise_1.75'])-2025+year
        tyr_set = np.array(table_scheduler['tyr_set_1.75'])-2025+year

        tyr_mid = 0.5*(tyr_set+tyr_rise)
        month_obs_baseline = int(month_obs_baseline)
        if month_obs_baseline!=12:
            tyr_rise2 = tyr_mid-0.5*month_obs_baseline/12
            tyr_set2 = tyr_mid+0.5*month_obs_baseline/12
            tyr_rise[tyr_rise2>tyr_rise] = tyr_rise2[tyr_rise2>tyr_rise]
            tyr_set[tyr_set2<tyr_set] = tyr_set2[tyr_set2<tyr_set]
        tyr_set[tyr_rise>year+1] = tyr_set[tyr_rise>year+1]-1
        tyr_rise[tyr_rise>year+1] = tyr_rise[tyr_rise>year+1]-1

        tyr_rise = tcsf.conv_time(list(tyr_rise))
        tyr_set = tcsf.conv_time(list(tyr_set))
        season_length = tyr_set[0]-tyr_rise[0]

        if need_help:
            self.plot_survey_stars(Nb_star=len(table_scheduler))
            pouet

        if type(texp)==str:
            texp = 900

        loc = tcsf.find_nearest(self.info_XY_survey_stat.x*60,texp)[0][0]
        nobs_max = self.info_XY_survey_stat.y[loc]

        warning = 0
        if type(n_obs)==str:
            n_obs = nobs_max

        if freq_obs is not None:
            table_scheduler['obsN'] = (season_length*freq_obs).astype('int')
        else:
            table_scheduler['obsN'] = n_obs

        table_scheduler['expTime'] = texp

        total_time = np.sum(table_scheduler['obsN']*(table_scheduler['expTime']+60))/3600
        total_max = self.info_SC_nb_hours_per_yr
        total_max_eff = self.info_SC_nb_hours_per_yr_eff

        fraction = 100*total_time/total_max
        fraction_eff = 100*total_time/total_max_eff
        if fraction>125:
            print('\n---->[WARNING] You overfilled the scheduler! Carefully consider this option.\n')
            warning = 1
        print('[INFO] You filled the GTO time (60%% of the telescope time) at %.0f%%'%(fraction))
        print('[INFO] You filled the effective GTO time (60%% of the telescope time + weather) at %.0f%%'%(fraction_eff))
        
        if plot_real_ID:
            table_scheduler['ID_table'] = table_scheduler['GR8_ID']           
        else:
            table_scheduler['ID_table'] = np.arange(len(table_scheduler))+1

        table_scheduler2 = table_scheduler.copy()
        for n in range(len(table_scheduler)):
            t1 = tyr_rise[1][n]
            t2 = tyr_set[1][n]
            if t2>(year+1):
                frag1 = year+1-t1
                frag2 = t2-year-1
                f1 = frag1/(frag2+frag1)
                f2 = frag2/(frag2+frag1)

                table_scheduler.loc[n,'groupEnableTime'] = '%.0f-01-01T00:00:00.000'%(year)
                table_scheduler.loc[n,'groupDisableTime'] = str(year)+tyr_set[2][n][4:]
                table_scheduler.loc[n,'obsN'] = int(f2*table_scheduler.loc[n,'obsN'])

                table_scheduler2.loc[n,'groupEnableTime'] = tyr_rise[2][n]
                table_scheduler2.loc[n,'groupDisableTime'] = '%.0f-12-31T23:59:00.000'%(year)
                table_scheduler2.loc[n,'obsN'] = table_scheduler2.loc[n,'obsN'] - table_scheduler.loc[n,'obsN']
            else:
                table_scheduler.loc[n,'groupEnableTime'] = tyr_rise[2][n]
                table_scheduler.loc[n,'groupDisableTime'] = tyr_set[2][n] 

        table_scheduler_final = pd.concat([table_scheduler,table_scheduler2],axis=0)
        table_scheduler_final = table_scheduler_final.loc[table_scheduler_final['obsN']!=0]

        table_scheduler_final = table_scheduler_final.dropna(subset=['groupEnableTime']).reset_index(drop=True) 

        variables = ['priority','GDR3_ID_number','expTime','expN','obsN','groupEnableTime','groupDisableTime',
                     'acqType','schedulingMode','t0','period','delta',
                     'moonMaxFI','minMoonDist1','ifExceedsMinFI1','minMoonDist2','ifExceedsMinFI2','rvMoonRange','rvMoonMinFI','rvMoonMinAlt',
                     'twilight','maxSeeing','maxExtinction','maxSkyBrightness','minAltitude']

        gto_time = np.sum(1-self.info_IM_night.data,axis=0)

        t0 = tcsf.conv_time([str(year)+'-01-01T00:00:00.000'])[0]
        plt.figure(figsize=(10,10))
        plt.axes([0.08,0.1,0.85,0.65])
        jdb1 = tcsf.conv_time(list(table_scheduler_final['groupEnableTime']))[0]
        jdb2 = tcsf.conv_time(list(table_scheduler_final['groupDisableTime']))[0]
        obs = []
        for n in np.arange(len(table_scheduler_final)):
            N = np.array(table_scheduler_final.loc[n,'obsN'])
            ID = np.array(table_scheduler_final.loc[n,'ID_table'])
            rank = np.array(table_scheduler_final.loc[n,'priority'])
            texp = np.array(table_scheduler_final.loc[n,'expTime'])
            days = np.arange(jdb1[n],jdb2[n]+1,1)
            Draw = np.min([N,len(days)])
            obs.append([texp*np.ones(Draw), np.random.choice(days,Draw,replace=False)])
            plt.scatter(obs[-1][1],ID*np.ones(Draw),s=rank,alpha=rank/10,color='k')
        obs = np.hstack(obs)
        obs[1] = ((obs[1]-t0)/366*12).astype('int')+1
        stat = tcsv.months_specie.copy()
        for j in range(1,13):
            stat[stat==j] = np.sum(obs[0][obs[1]==j])/60/tcsv.month_len[j-1]

        plt.scatter([np.nan],[np.nan],color='k',s=9, alpha=0.9,label='Priotity=9')
        plt.scatter([np.nan],[np.nan],color='k',s=6,alpha=0.6,label='Priotity=6')
        plt.scatter([np.nan],[np.nan],color='k',s=3,alpha=0.3,label='Priotity=3')
        for j in t0+np.array(tcsv.month_border):
            plt.axvline(x=j,lw=1,alpha=0.4,color='k')
        plt.legend()
        plt.ylabel('Star ID')
        plt.xlabel('Jdb time - 2,400,000 [days]')
        ax = plt.gca()

        plt.axes([0.08,0.77,0.85,0.2],sharex=ax)
        plt.tick_params(top=True,labeltop=True,labelbottom=False)
        plt.plot(t0+self.info_XY_downtime.x,0.60*gto_time,label='GTO')
        plt.plot(t0+self.info_XY_downtime.x,0.60*gto_time*(1-self.info_XY_downtime.y/100),label='GTO + weather')
        plt.plot(t0+np.arange(0,365,1),stat,marker='.',label='Obs time')
        plt.legend()
        plt.ylabel('Time per night [min]')
        for j in t0+np.array(tcsv.month_border):
            plt.axvline(x=j,lw=1,alpha=0.4,color='k')

        table_scheduler_final = table_scheduler_final[variables]

        now = tcsf.now()[0:19].replace(':','-')

        plt.savefig(cwd+'/TACS_OUTPUT/TAB_SCHEDULER/scheduler_%s_%.0f_B%.0f%s.png'%(now,year,int(month_obs_baseline),tagname))
        table_scheduler_final.to_csv(cwd+'/TACS_OUTPUT/TAB_SCHEDULER/scheduler_%s_%.0f_B%.0f%s.csv'%(now,year,int(month_obs_baseline),tagname))

