from datetime import datetime, timezone

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from astropy import units as u


def find_nearest(array,value,dist_abs=True,closest='abs'):
    if type(array)!=np.ndarray:
        array = np.array(array)
    if type(value)!=np.ndarray:
        value = np.array([value])
    
    array[np.isnan(array)] = 1e16

    dist = array-value[:,np.newaxis]
    if closest=='low':
        dist[dist>0] = -np.inf
    elif closest=='high':
        dist[dist<0] = np.inf
    idx = np.argmin(np.abs(dist),axis=1)
    
    distance = abs(array[idx]-value) 
    if dist_abs==False:
        distance = array[idx]-value
    return idx, array[idx], distance

def HZ(Ms,mp):
    """Made by Baptiste Klein"""
    return 

def Keplerian_rv(jdb, P, K, e, omega, t0):
    # Keplerian parameters
    P = P * u.day             # Orbital period
    e = e                    # Eccentricity
    omega = omega * u.deg         # Argument of periastron
    K = K * u.m/u.s           # RV semi-amplitude
    t0 = t0 * u.day             # Time of conjunction passage
    gamma = 0 * u.m/u.s        # Systemic velocity

    # Observation times
    t = jdb * u.day

    # Compute true anomaly at conjunction
    fc = np.pi/2 - omega.to(u.rad).value  # Inferior conjunction

    # Eccentric anomaly at conjunction
    tan_half_E = np.tan(fc / 2) * np.sqrt((1 - e) / (1 + e))
    Ec = 2 * np.arctan(tan_half_E)

    # Normalize to [0, 2pi)
    Ec = np.mod(Ec, 2 * np.pi)

    # Mean anomaly at conjunction
    Mc = Ec - e * np.sin(Ec)

    # Compute T0 (periastron time) from Tc
    T0 = t0 - (P * Mc / (2 * np.pi))

    # Mean anomaly
    M = 2 * np.pi * ((t - T0) / P).decompose()  # unitless

    # Solve Kepler’s equation (E - e*sin(E) = M)
    def solve_kepler(M_val, e, tol=1e-10):
        E = M_val.copy()
        for _ in range(100):
            delta = (E - e * np.sin(E) - M_val) / (1 - e * np.cos(E))
            E -= delta
            if np.all(np.abs(delta) < tol):
                break
        return E

    E = solve_kepler(M.value, e)  # Pass in unitless mean anomaly

    # True anomaly (unitless)
    f = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2),
                    np.sqrt(1 - e) * np.cos(E / 2))

    # RV equation
    RV = gamma + K * (np.cos(f + omega.to(u.rad).value) + e * np.cos(omega.to(u.rad).value))
    return np.array(RV.to(u.m/u.s))


def Fisher_std(jdb, K=1, P=100, phi=0.0, sigma=0.5):
    jdb = jdb - np.min(jdb)
    # Dérivées
    omega_t = 2 * np.pi * jdb / P
    s = np.sin(omega_t + phi)
    c = np.cos(omega_t + phi)
    
    dA = s
    dphi = K * c
    dP = -K * (2 * np.pi * jdb / P**2) * c

    # Matrice de Fisher
    I = np.zeros((3, 3))
    I[0, 0] = np.sum(dA * dA)
    I[0, 1] = np.sum(dA * dphi)
    I[0, 2] = np.sum(dA * dP)
    I[1, 1] = np.sum(dphi * dphi)
    I[1, 2] = np.sum(dphi * dP)
    I[2, 2] = np.sum(dP * dP)

    # Symétriser
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]

    # Appliquer facteur 1/σ²
    I /= sigma**2

    # Erreurs (écart-types)
    cov = np.linalg.inv(I)
    errors = np.sqrt(np.diag(cov))

    print("\n1-sigma uncertainties (std dev):")
    print(f"Amplitude: {errors[0]:.4f}")
    print(f"Phase:     {errors[1]:.4f} rad")
    print(f"Période:   {errors[2]:.4f}")


# -------- FONCTIONS -------- #

def julian_date(dt):
    """Convertit une datetime UTC en date julienne"""
    timestamp = dt.timestamp()
    return timestamp / 86400.0 + 2440587.5

def greenwich_sidereal_time(jd):
    """Temps sidéral à Greenwich en degrés"""
    T = (jd - 2451545.0) / 36525.0
    gst = 280.46061837 + 360.98564736629 * (jd - 2451545) \
          + 0.000387933 * T**2 - T**3 / 38710000.0
    return gst % 360

def local_sidereal_time(gst_deg, longitude_deg):
    """Temps sidéral local en degrés"""
    return (gst_deg + longitude_deg) % 360

def hour_angle(lst_deg, ra_hours):
    """Angle horaire en degrés"""
    ra_deg = ra_hours * 15
    return (lst_deg - ra_deg + 360) % 360

def altitude(lat_deg, dec_deg, ha_deg):
    """Altitude de l'étoile (en degrés)"""
    lat = np.radians(lat_deg)
    dec = np.radians(dec_deg)
    ha = np.radians(ha_deg)
    sin_alt = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(ha)
    return np.degrees(np.arcsin(sin_alt))

def airmass_kasten_young(alt_deg):
    """Airmass selon la formule de Kasten & Young (1989)"""
    z = 90 - alt_deg
    airmass = 1 / (np.cos(np.radians(z)) + 0.50572 * (96.07995 - z)**-1.6364)
    airmass[alt_deg<0] = 10 
    airmass[airmass>10] = 10 
    return airmass

def star_observability(alpha_h, delta_deg, tstamp_min=1, Plot=False, instrument='HARPS3'):
    """tstamp_min being the sampling in minute"""
    
    lat_deg = {
        'HARPN':28.7540,
        'HARPS3':28.7540,
        'HARPS':-29.260972,
        'CORALIE':-29.260972,
        'ESPRESSO':-24.627622,
        'SOPHIE':43.930833,
        'EXPRES':34.74444,
        'NEID':31.9584,
        'KPF':19.8261,
        }[instrument]     
    
    lon_deg = {
        'HARPN':-17.8890,
        'HARPS3':-17.8890,
        'HARPS':-70.731694,
        'CORALIE':-70.731694,
        'ESPRESSO':-70.405075,
        'SOPHIE':5.713333,
        'EXPRES':-68.578056,
        'NEID':-111.5987,
        'KPF':-155.4700,
        }[instrument]  

    # Date et heure UTC de l'observation
    utc_datetime = datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc)

    jd = julian_date(utc_datetime)

    hours = np.linspace(-0.5,0.5,int(24*60/tstamp_min)+1)

    gst = greenwich_sidereal_time(jd+hours)
    lst = local_sidereal_time(gst, lon_deg)
    ha = hour_angle(lst, alpha_h)
    alt = altitude(lat_deg, delta_deg, ha)
    airmass = airmass_kasten_young(alt)

    if Plot:
        plt.plot(hours,airmass,marker='.')

    return hours, airmass


def func_cutoff(table, cutoff, tagname='', plot=True, par_space='', par_box=['',''], par_crit='', verbose=True):
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
            mask = table2[kw[0:-1]]<=value
        else:
            mask = table2[kw[0:-1]]>=value
        
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
        plt.subplots_adjust(hspace=0.45,wspace=0.3,top=0.95,bottom=0.08,left=0.08,right=0.95)
        if par_space!='':
            plt.figure('para'+tagname,figsize=(18,4*nb_rows))
            plt.subplots_adjust(hspace=0.45,wspace=0.3,top=0.95,bottom=0.08,left=0.08,right=0.95)
        plt.show()
    table2 = table2.sort_values(by='HZ_mp_min_osc+gr_texp15')
    
    
    #printable table
    print_table = table2[0:30][['ra_j2000','dec_j2000','primary_name','vmag','eff_nights_1.5','dist','teff_mean','HZ_mp_min_osc+gr_texp15']]
    print_table['vmag'] = np.round(print_table['vmag'],2)
    print_table['dist'] = np.round(print_table['dist'],1)
    print_table['eff_nights_1.5'] = (np.round(print_table['eff_nights_1.5'],0)).astype('int')
    print_table['teff_mean'] = (np.round(print_table['teff_mean'],0)).astype('int')
    print_table['HZ_mp_min_osc+gr_texp15'] = np.round(print_table['HZ_mp_min_osc+gr_texp15'],2)
        
    if verbose:
        print('\n [INFO] %.0f stars in the final sample'%(len(table2)))
        print('\n [INFO] Here are the top 30-ranked stars of your THE list:\n')  
        print(print_table)
    
    return table2

