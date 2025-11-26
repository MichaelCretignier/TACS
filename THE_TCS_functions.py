import matplotlib

try:
    matplotlib.use('Qt5Agg',force=True)
except:
    matplotlib.use('Agg',force=True)

import warnings
from datetime import datetime, timedelta, timezone

import astropy.time as Time
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from astropy import units as u
from IPython import get_ipython
from matplotlib import MatplotlibDeprecationWarning
from matplotlib.collections import LineCollection

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)
warnings.filterwarnings('ignore', message='ERFA function.*yielded.*')

ipython = get_ipython()
if ipython is not None:
    try:
        ipython.magic('%matplotlib qt')
    except:
        try:
            ipython.run_line_magic('matplotlib', 'qt')
        except:
            pass

def sort_hd(names):
    number = [n[2:] for n in names]
    for n in np.arange(len(number)):
        try:
            int(number[n][-1])
        except:
            number[n] = number[n][:-1]
    number = np.array(number).astype('int')
    sort = np.argsort(number)
    return names[sort],sort

def _flatten_dict(d, parent_key=""):
    items = {}
    for k, v in d.items():
        new_key = f"{parent_key}/{k}" if parent_key else k
        if isinstance(v, dict):
            items.update(_flatten_dict(v, new_key))
        else:
            items[new_key] = v
    return items

def _unflatten_dict(d):
    result = {}
    for key, value in d.items():
        parts = key.split("/")
        cur = result
        for p in parts[:-1]:
            cur = cur.setdefault(p, {})
        cur[parts[-1]] = value
    return result

def save_nested_dict_npz(d, filename):
    flat = _flatten_dict(d)
    np.savez(filename, **flat)

def load_nested_dict_npz(filename):
    data = np.load(filename, allow_pickle=True)
    return _unflatten_dict({k: data[k] for k in data.files})

def format_number(nb,digit=3):
    nb_log = int(np.round(np.log10(nb)-0.5,0))
    nb_digit = digit-nb_log
    if nb_digit<0:
        nb_digit=0
    
    output = ['%.0f'%(nb),'%.1f'%(nb),'%.2f'%(nb),'%.3f'%(nb),'%.4f'%(nb),'%.5f'%(nb),'%.6f'%(nb)][nb_digit]
    
    return output

def format_mass(mass_planet):
    symbol = '⊕'
    mass_printed = mass_planet
    if mass_planet>3300*6:
        symbol = '⊙'
        mass_printed = mass_printed/333030
    elif mass_planet>300:
        symbol = '♃'
        mass_printed = mass_printed/317.8   
    elif mass_planet>95:
        symbol = '♄'
        mass_printed = mass_printed/95.16                    
    elif mass_planet>15:
        symbol = '♆'
        mass_printed = mass_printed/17.15
    return '%.1f M%s'%(mass_printed,symbol)

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

def now():
    iso_time = datetime.now(timezone.utc).isoformat()
    return iso_time

def decimal_year_to_iso(decimal_years):
    iso_times = []
    for y in decimal_years:
        year = int(np.floor(y))
        remainder = y - year
        
        # Check if leap year
        days_in_year = 366 if ((year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)) else 365
        
        # Fractional year → days + seconds
        seconds_in_year = days_in_year * 24 * 3600
        delta_seconds = remainder * seconds_in_year
        
        dt = datetime(year, 1, 1) + timedelta(seconds=delta_seconds)
        iso_times.append(dt.isoformat())
    iso_times = pd.DataFrame(np.array([list(decimal_years),iso_times]).T,columns=['deciyear','iso'])
    iso_times['iso'] = iso_times['iso'].str[0:10]+'T00:00:00.000000'
    return iso_times

def julian_date(dt):
    """Convertit une datetime UTC en date julienne"""
    timestamp = dt.timestamp()
    return timestamp / 86400.0 + 2440587.5

def conv_time(time):    
    time = np.array(time)
    if (type(time[0])==np.float64)|(type(time[0])==np.int64):
        fmt='mjd'
        if time[0]<2030:
            fmt='decimalyear'
        elif np.mean(time)<20000:
            time+=50000
        if fmt=='mjd':
            t0 = time
            t1 = np.array([Time.Time(i, format=fmt).decimalyear for i in time])
            t2 = np.array([Time.Time(i, format=fmt).isot for i in time])
        else:
            t0 = np.array([Time.Time(i, format=fmt).mjd for i in time])
            t1 = time
            t2 = np.array([Time.Time(i, format=fmt).isot for i in time])            
    elif type(time[0])==np.str_:
        fmt='isot'
        t0 = np.array([Time.Time(i, format=fmt).jd-2400000 for i in time]) 
        t1 = np.array([Time.Time(i, format=fmt).decimalyear for i in time])
        t2 = time  
    return t0,t1,t2

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

def star_observability(alpha_h, delta_deg, tstamp_min=1, Plot=False, instrument='HARPS3', day=1, month=1):
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
        '2ES':-29.25786,
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
        '2ES':-70.73666,
        }[instrument]  

    # Date et heure UTC de l'observation
    utc_datetime = datetime(2025, month, day, 0, 0, 0, tzinfo=timezone.utc)

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
    nb_rows = (len(cutoff)-1)//7+1
    if par_crit!='':
        p1c = par_crit.split('==')[0]
        p1c_val = float(par_crit.split('==')[1])
        mask_box = (table2[p1c].astype('float')==p1c_val)
        old_value = sum(mask_box)
    else:
        old_value = np.nan
    old_value2 = len(table)
    for kw in cutoff.keys():
        count+=1
        value = cutoff[kw]
        if kw[-1]=='<':
            mask = (table2[kw[0:-1]]<value)|(table2['under_review']==1)
        else:
            mask = (table2[kw[0:-1]]>value)|(table2['under_review']==1)
            
        if plot:
            plt.figure('cumulative'+tagname,figsize=(16,3*nb_rows))
            plt.subplot(nb_rows,7,count)
            if (kw[0:-1]!='under_review')&(kw[0:-1]!='HWO')&(kw[0:-1]!='PLATO'):
                plt.title(kw+str(value))
            else:
                plt.title(kw[0:-1]+'(%.0f)'%(sum(table2[kw[0:-1]])))
            plt.hist(table2[kw[0:-1]],cumulative=True,bins=100)
            plt.axvline(x=value,label='%.0f / %.0f'%(sum(mask),len(mask)),color='k')
            if len(table2):
                xmax = np.nanmax(table2[kw[0:-1]])
                xmin = np.nanmin(table2[kw[0:-1]])
            else:
                xmax=value
                xmin=value
            if kw[-1]=='<':
                plt.axvspan(xmin=value,xmax=xmax,alpha=0.2,color='k')
            else:
                plt.axvspan(xmax=value,xmin=xmin,alpha=0.2,color='k')
            plt.legend(loc=4)
            plt.grid()
            plt.xlim(xmin,xmax)
            table2 = table2[mask]#.reset_index(drop=True)
            if par_space!='':
                p1 = par_space.split('&')[0].replace(' ','')
                p2 = par_space.split('&')[1].replace(' ','')
                plt.figure('para'+tagname,figsize=(18,3.5*nb_rows))
                if count==1:
                    plt.subplot(nb_rows,6,count)
                    ax1 = plt.gca()
                else:
                    plt.subplot(nb_rows,6,count,sharex=ax1,sharey=ax1)
                plt.scatter(table[p1],table[p2],color='k',alpha=0.1,marker='.')
                plt.scatter(table2[p1],table2[p2],color='r',ec='k',marker='.',label='%.0f (-%.0f)'%(len(table2),old_value2-len(table2)))
                old_value2 = len(table2)
                if (par_box[0]!='')|(par_box[1]!=''):
                    p1x = np.array(par_box[0].split('->')).astype('float')
                    p2y = np.array(par_box[1].split('->')).astype('float')
                    mask_box = (table2[p1]>p1x[0])&(table2[p1]<p1x[1])&(table2[p2]>p2y[0])&(table2[p2]<p2y[1])
                    mask_box = np.array(mask_box)
                    plt.scatter(np.array(table2[p1])[mask_box],np.array(table2[p2])[mask_box],color='g',ec='k',marker='o',label='%.0f (-%.0f)'%(sum(mask_box),old_value-sum(mask_box)))
                    old_value = sum(mask_box)
                if par_crit!='':
                    p1c = par_crit.split('==')[0]
                    p1c_val = float(par_crit.split('==')[1])
                    mask_box = (table2[p1c].astype('float')==p1c_val)
                    mask_box = np.array(mask_box)
                    plt.scatter(np.array(table2[p1])[mask_box],np.array(table2[p2])[mask_box],color='g',ec='k',marker='o',label='%.0f (-%.0f)'%(sum(mask_box),old_value-sum(mask_box)))
                    old_value = sum(mask_box)   

                plt.xlabel(p1)
                plt.ylabel(p2)
                plt.legend(loc=1)
                plt.title(kw+str(value))

    if plot:
        plt.figure('cumulative'+tagname,figsize=(18,4*nb_rows))
        plt.subplots_adjust(hspace=0.45,wspace=0.3,top=0.93,bottom=0.08,left=0.08,right=0.95)
        if par_space!='':
            plt.figure('para'+tagname,figsize=(18,4*nb_rows))
            plt.subplots_adjust(hspace=0.45,wspace=0.3,top=0.93,bottom=0.08,left=0.08,right=0.95)
        plt.show()
    ranking = 'HZ_mp_min_osc+gr_texp15'
    table2 = table2.sort_values(by=ranking)
    
    
    #printable table
    print_table = table2[0:30][['ra_j2000','dec_j2000','PRIMARY','vmag','eff_nights_1.5','distance','teff',ranking]]
    print_table['vmag'] = np.round(print_table['vmag'],2)
    print_table['distance'] = np.round(print_table['distance'],1)
    print_table['eff_nights_1.5'] = (np.round(print_table['eff_nights_1.5'],0)).astype('int')
    print_table['teff'] = (np.round(print_table['teff'],0)).astype('int')
    print_table[ranking] = np.round(print_table[ranking],2)
        
    if verbose:
        print('\n[INFO] %.0f stars in the final sample'%(len(table2)))
        print('\n[INFO] Here are the top 30-ranked stars of your THE list:\n')  
        print(print_table)
    
    return table2


# analytic transform (ICRS/J2000 approx) RA,Dec (deg) -> ecliptic latitude beta (deg)
def ra_dec_to_ecl_lat(ra_deg, dec_deg):
    eps = 23.439291111  # obliquity of the ecliptic (deg, J2000)
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    eps_r = np.deg2rad(eps)
    sinb = np.sin(dec) * np.cos(eps_r) - np.cos(dec) * np.sin(eps_r) * np.sin(ra)
    # numeric safety
    sinb = np.clip(sinb, -1.0, 1.0)
    return np.rad2deg(np.arcsin(sinb))

# mapping |beta| -> approximate number of sectors (linear from 54->1 to 90->13)
def beta_to_sectors(abs_beta_deg):
    s = 1.0 + (abs_beta_deg - 54.0) / (90.0 - 54.0) * 12.0
    s = np.clip(s, 1.0, 13.0)
    return np.rint(s).astype(int)

def plot_TESS_sectors(ra_unit='hours'):

    # grid
    n_ra = 360   # 1-deg steps
    n_dec = 1801  # 1-deg steps
    ra_vals = np.linspace(0, 360, n_ra, endpoint=False)
    dec_vals = np.linspace(-90, 90, n_dec)
    RA, DEC = np.meshgrid(ra_vals, dec_vals)

    # compute beta and sectors
    ra_flat = RA.ravel()
    dec_flat = DEC.ravel()
    beta_flat = ra_dec_to_ecl_lat(ra_flat, dec_flat)
    abs_beta = np.abs(beta_flat)
    sectors_flat = beta_to_sectors(abs_beta)
    sectors = sectors_flat.reshape(RA.shape)

    # plot
    img = plt.imshow(sectors, origin='lower', extent=(0,[360,24][int(ra_unit=='hours')],-90,90), aspect='auto',cmap='Reds')
    # contour lines

    cs = plt.contour(RA/[1,15][int(ra_unit=='hours')], DEC, sectors, levels=np.arange(1,14), linewidths=0.4, colors='k', alpha=0.8)
    plt.clabel(cs, inline=True, fontsize=8, fmt='%d')
    plt.show()

def _circle_patch(ra_c_deg, dec_c_deg, r_deg, ra_unit='hours', **plot_kw):
    t = np.linspace(0,2*np.pi,200)
    ra_off = (r_deg * np.cos(t)) / np.cos(np.radians(dec_c_deg))
    dec = dec_c_deg + r_deg * np.sin(t)
    ra_deg = ra_c_deg + ra_off
    if ra_unit == 'hours':
        ra_plot = ra_deg / 15.0
        ra_label = ra_c_deg / 15.0
    else:
        ra_plot = ra_deg
        ra_label = ra_c_deg
    plt.plot(ra_plot, dec, **plot_kw)
    return ra_label, dec_c_deg

def plot_TESS_CVZ(ra_unit='hours'):
    ra_label, dc = _circle_patch(18.0*15, 66.5, 12.0, ra_unit, lw=1, ls='-.', color='k')
    plt.text(ra_label, dc, 'TESS\nCVZ', color='k', ha='center', va='center',fontweight='bold')

def plot_KEPLER(ra_unit='hours'):
    ra_label, dc = _circle_patch(19.3777777778*15, 44.5, 8.0, ra_unit, lw=1, ls=':', color='k')
    plt.text(ra_label, dc, 'KEPLER', color='k', ha='center', va='center',fontweight='bold')

def plot_PLATO_North_LOP(ra_unit='hours'):
    ra_label, dc = _circle_patch(277.1, 52.9, 15.0, ra_unit, lw=1.2, ls='--', color='k')
    plt.text(ra_label, dc, 'PLATO\nNorth', color='k', ha='center', va='center',fontweight='bold')

def holman_stability(m1,m2,e,P):
    """Return Acrit in AU and Pcrit in years"""
    M = m1+m2
    mu = m2/(m1+m2)
    acrit = (0.464 - 0.380*mu - 0.631*e + 0.586*mu*e + 0.150*e**2 - 0.198*mu*e**2)*(M*P**2)**(1/3)
    pcrit = np.sqrt(acrit**3/m1)
    return acrit,pcrit

def kepler_E(M, e, tol=1e-12, maxiter=100):
    """Solve Kepler's equation M = E - e*sin(E), vectorized Newton-Raphson."""
    M = np.asarray(M, dtype=float)
    # bring M to [-pi, pi] for numerical stability
    M_wrapped = (M + np.pi) % (2*np.pi) - np.pi
    # initial guess
    E = M_wrapped + (e * np.sign(np.sin(M_wrapped))) * 0.85
    # fallback for small e
    E[np.isnan(E)] = M_wrapped[np.isnan(E)]
    for _ in range(maxiter):
        f = E - e * np.sin(E) - M_wrapped
        fp = 1 - e * np.cos(E)
        dE = f / fp
        E -= dE
        if np.all(np.abs(dE) < tol):
            break
    return E

def canonicalize_elements(inc_deg, Omega_deg, omega_deg):
    """
    Canonicalize (inc, Omega, omega) so that equivalent representations
    produce the same sky projection.

    Rule used:
    - Map inclinations into the [0,90] representative by:
        if inc in (90, 270]: replace inc -> 180 - inc and apply Omega += 180, omega += 180
      This handles the i=180 -> i=0 case and keeps the projection invariant.
    """
    inc = inc_deg % 360.0
    Omega = Omega_deg % 360.0
    omega = omega_deg % 360.0

    # If the inclination indicates the plane is flipped (i in (90,270]),
    # transform to the equivalent "canonical" branch: inc' = 180 - inc, and
    # shift node & periastron by +180° to preserve projection.
    if 90.0 < inc <= 270.0:
        inc = 180.0 - inc   # brings inc into [-90, 90], but symmetric; result will be <= 90 in magnitude
        # normalize to positive in [0,180)
        if inc < 0:
            inc = -inc
        Omega = (Omega + 180.0) % 360.0
        omega = (omega + 180.0) % 360.0

    # final safety: ensure inc in [0,180)
    inc = inc % 360.0
    if inc >= 180.0:
        inc -= 180.0

    return inc, Omega, omega

def orbit_secondary_relative_to_primary(times, a, e, P, T0,
                                        inc_deg, Omega_deg, omega_deg,
                                        return_arrays=True):
    """
    Compute secondary star positions projected on sky relative to primary at (0,0).
    This version canonicalizes elements so that i=0 and i=180 (and equivalent representations)
    yield identical projected ellipses.
    Returns dict {'x','y','z'} (arrays).
    """
    # canonicalize elements first
    inc_c, Omega_c, omega_c = canonicalize_elements(inc_deg, Omega_deg, omega_deg)

    t = np.asarray(times, dtype=float)
    inc = np.deg2rad(inc_c)
    Omega = np.deg2rad(Omega_c)
    omega = np.deg2rad(omega_c)

    # Mean anomaly
    M = 2.0 * np.pi * (t - T0) / P

    # Eccentric anomaly
    E = kepler_E(M, e)

    # True anomaly nu
    nu = 2.0 * np.arctan2(np.sqrt(1.0 + e) * np.sin(E/2.0),
                          np.sqrt(1.0 - e) * np.cos(E/2.0))

    # radius r (same units as a)
    r = a * (1.0 - e * np.cos(E))

    # perifocal coordinates
    xp = r * np.cos(nu)
    yp = r * np.sin(nu)
    zp = np.zeros_like(xp)

    # Combined rotation R = Rz(Omega) @ Rx(inc) @ Rz(omega)
    cosO = np.cos(Omega); sinO = np.sin(Omega)
    cosi = np.cos(inc);   sini = np.sin(inc)
    cosw = np.cos(omega); sinw = np.sin(omega)

    R11 =  cosO * cosw - sinO * sinw * cosi
    R12 = -cosO * sinw - sinO * cosw * cosi
    R13 =  sinO * sini

    R21 =  sinO * cosw + cosO * sinw * cosi
    R22 = -sinO * sinw + cosO * cosw * cosi
    R23 = -cosO * sini

    R31 =  sinw * sini
    R32 =  cosw * sini
    R33 =  cosi

    X = R11 * xp + R12 * yp + R13 * zp
    Y = R21 * xp + R22 * yp + R23 * zp
    Z = R31 * xp + R32 * yp + R33 * zp

    out = {'x': X, 'y': Y, 'z': Z}

    if not return_arrays:
        out = {k: v.tolist() for k, v in out.items()}

    return out

def plot_binaries(info, fibre=1.4, seeing=0.0, inc=None, traj='new', source='COMPOSITE', t_eval=2026, print_source=True):
    t = np.linspace(0,2*np.pi,1000)
    sb = info.copy()
    ref = sb.loc[sb['origin']=='COMPOSITE'].reset_index(drop=True).loc[0].copy()
    if source!='COMPOSITE':
        ref1 = sb.loc[sb['origin']==source].reset_index(drop=True).loc[0].copy()
        for kw in ref1.keys():
            if ref1[kw]==ref1[kw]:
                ref[kw] = ref1[kw]
    print(ref)

    j = ref['ID']
    starname = ref['PRIMARY']
    hip  = ref['HIP']
    wds = ref['WDS']
    vmag = ref['vmag']
    obtp = ref['OBTP']
    plx = ref['plx']
    dist = 1000/plx
    M1 = ref['Ms']
    M2 = ref['Ms2']
    period = ref['period']
    ecc = ref['ecc']
    node = ref['node']
    omega = ref['omega']
    T0 = ref['T0']
    mv1 = ref['mv1']
    mv2 = ref['mv2']
    if inc is None:
        inc = ref['inc']
    
    if inc!=inc:
        inc=52
    if ecc!=ecc:
        ecc = 0.0
    if node!=node:
        node = 180
    if omega!=omega:
        omega = 90
    if T0!=T0:
        T0 = 0.0

    flux_ratio = np.round(10**(-0.4*(mv2-mv1))*100,2)
    a_au = np.array((M1+M2)**(1/3)*period**(2/3)).astype('float')
    a_arcsec = np.array(a_au*plx/1000).astype('float')
    amin_au = np.round(a_au*(1-ecc),2)
    amin_years = np.round(amin_au**(3/2)/(M1+M2)**(1/3),2)
    acrit_au,acrit_years = holman_stability(M1,M2,ecc,period)

    if traj=='new':
        times = np.linspace(T0,T0+period,1000)

        pos_THE = orbit_secondary_relative_to_primary(np.linspace(t_eval,t_eval+10,100), a_arcsec, ecc, period, T0, inc, node, omega)
        the_x_arcsec = pos_THE['x'] ; the_y_arcsec = pos_THE['y']        
        the_amin_arcsec = np.min(np.sqrt(the_x_arcsec**2+the_y_arcsec**2))

        pos = orbit_secondary_relative_to_primary(times, a_arcsec, ecc, period, T0, inc, node, omega)
        x_arcsec = pos['x'] ; y_arcsec = pos['y']
        amin_arcsec = np.min(np.sqrt(x_arcsec**2+y_arcsec**2))
        the_amin_arcsec = np.min(np.sqrt(the_x_arcsec**2+the_y_arcsec**2))

        pos = orbit_secondary_relative_to_primary(times, a_arcsec, ecc, period, T0, 0, node, omega)
        x_arcsec0 = pos['x'] ; y_arcsec0 = pos['y']

        pos = orbit_secondary_relative_to_primary(np.linspace(t_eval,t_eval+10,100), a_au, ecc, period, T0, 0, node, omega)
        the_x_au = pos['x'] ; the_y_au = pos['y']

        pos = orbit_secondary_relative_to_primary(times, a_au, ecc, period, T0, 0, node, omega)
        x_au = pos['x'] ; y_au = pos['y']
    else:
        c_arcsec = a_arcsec*ecc
        b_arcsec = a_arcsec*np.sqrt(1-ecc**2)

        c_au = a_au*ecc
        b_au = a_au*np.sqrt(1-ecc**2)

        x_arcsec = b_arcsec*np.cos(t)
        y_arcsec = (a_arcsec*np.sin(t)+c_arcsec)*abs(np.cos(inc*np.pi/180))
        the_x_arcsec = x_arcsec
        the_y_arcsec = y_arcsec
        x_au = b_au*np.cos(t)
        y_au = a_au*np.sin(t)+c_au
        the_x_au = x_au
        the_y_au = y_au
        x_arcsec0 = x_arcsec
        y_arcsec0 = y_arcsec/abs(np.cos(inc*np.pi/180))

        amin_arcsec = np.min(np.sqrt(x_arcsec**2+y_arcsec**2))
        the_amin_arcsec = amin_arcsec

    thickness = 1 + abs(np.sin(inc*np.pi/180))*4*(1 - (y_arcsec - y_arcsec.min())/(y_arcsec.max() - y_arcsec.min()))  # 1→5

    seeing_x = np.ravel(x_arcsec+0.5*seeing/np.sqrt(2)*np.random.randn(1000)[:,np.newaxis])
    seeing_y = np.ravel(y_arcsec+0.5*seeing/np.sqrt(2)*np.random.randn(1000)[:,np.newaxis])
    selected = np.random.choice(np.arange(len(seeing_x)),10000,replace=False)
    seeing_x = seeing_x[selected]
    seeing_y = seeing_y[selected]

    points = np.array([x_arcsec, y_arcsec]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    fig = plt.figure(figsize=(13,9))
    fig.suptitle('ID = %.0f | %s | %s | %s | plx = %.1f mas (%.1f pc) | $m_{v}$ = %.2f | %s'%(j,wds,starname,hip,plx,dist,vmag,obtp),fontsize=12)
    plt.subplot(1,2,1)
    plt.plot(np.cos(np.linspace(0,2*np.pi,9)+np.pi/8)*fibre*0.5,np.sin(np.linspace(0,2*np.pi,9)+np.pi/8)*fibre*0.5,color='k',ls='-',zorder=100)
    plt.plot(x_arcsec0,y_arcsec0,color='k',ls='-.',alpha=0.8,lw=1.5)
    plt.scatter(seeing_x,seeing_y,color='k',alpha=0.1,s=3)

    plt.axhline(y=0,color='k',lw=1)
    plt.axvline(x=0,color='k',lw=1)
    ax = plt.gca()
    lc = LineCollection(segments, linewidths=thickness, color='black')
    ax.add_collection(lc)
    ax.set_aspect('equal')
    plt.scatter([0],[0],color='k',marker='o',s=100,label='$M_1$=%.2f'%(M1),zorder=10)
    if (T0!=0)&(traj=='new'):
        plt.scatter(the_x_arcsec,the_y_arcsec,color='gold',marker='o',s=30,zorder=10)        
        plt.text(the_x_arcsec[0]+0.2,the_y_arcsec[0]+0.2,'%.0f'%(t_eval),ha='left',zorder=101,fontweight='bold')
    plt.scatter(x_arcsec[0],y_arcsec[0],color='r',marker='o',ec='k',s=100,zorder=10)
    plt.scatter(the_x_arcsec[0],the_y_arcsec[0],color='gray',marker='o',ec='k',s=100,label='$M_2$=%.2f ($F_1$/$F_2$=%.1f%%)'%(M2,flux_ratio),zorder=10)

    plt.text(0.6,0.4,'%.1f" fibre'%(fibre),ha='left')
    danger1 = int(amin_arcsec<((fibre+seeing)*0.5))
    if danger1==1:
        contam = np.sqrt(x_arcsec**2+y_arcsec**2)<(0.5*(fibre+seeing))
        plt.scatter(x_arcsec[contam],y_arcsec[contam],color='r',marker='x',zorder=7,label='CONTAM DANGER!')

    plt.legend()
    plt.xlabel('X ["]',fontsize=14)
    plt.ylabel('Y ["]',fontsize=14)
    plt.title('i = %.0f° | P = %.1f yrs | e = %.2f | T0=%.0f yrs | $\omega$=%.0f° | $\Omega$=%.0f° \n a = %.2f" | $a_{min}$ = %.2f" | $a_{THE}$ = %.2f"'%(inc, period, ecc, T0, omega, node, a_arcsec, amin_arcsec, the_amin_arcsec))

    plt.subplot(1,2,2)
    plt.axis('equal')
    plt.plot(x_au,y_au,color='k',ls='-')
    plt.plot(np.cos(t)*amin_au,np.sin(t)*amin_au,color='r',ls=':')
    plt.plot(np.cos(t)*acrit_au,np.sin(t)*acrit_au,color='r',ls='-')
    plt.text(0,amin_au,'%.1f yrs'%(amin_years),color='r',ha='center',va='bottom')
    plt.text(0,acrit_au,'$a_{crit}$ = %.1f yrs'%(acrit_years),color='r',ha='center',va='bottom')
    plt.plot(np.cos(t)*1,np.sin(t)*1,color='k',ls='-')
    plt.scatter([0],[0],color='k',marker='o',s=100,zorder=10)
    plt.scatter(x_au[0],y_au[0],color='r',marker='o',ec='k',s=100,zorder=10)
    plt.scatter(the_x_au[0],the_y_au[0],color='gray',marker='o',ec='k',s=100,label='Ms=%.2f (%.1f%%)'%(M2,flux_ratio),zorder=10)
    plt.grid()
    plt.xlabel('X [AU]',fontsize=14)
    plt.ylabel('Y [AU]',fontsize=14)
    plt.title('Orbital plan\na = %.1f AU | $a_{min}$ = %.1f AU | $a_{crit}$ = %.1f AU'%(a_au, amin_au, acrit_au))

    plt.subplots_adjust(bottom=0.10+0.20*int(print_source),top=0.86)

    if print_source:
        summary = ''
        sb = sb.sort_values(by='origin').reset_index(drop=True)
        for source in sb.index:
            text = ''
            ss = sb.loc[source][['sep_arcsec','mv1','mv2','Ms','Ms2','omega','node','inc','ecc','period','origin']]
            for kw in ss.keys():
                value = ss[kw]
                if value==value:
                    value = str(value).replace('.','․')
                    value = value + ' '*(5-len(value))
                else:
                    value = ' '*4+' '
                text = text + '%s = %s | '%(kw,value) 
            summary = summary+text+'\n'
        plt.axes([0,0,1,0.22])
        plt.axis('off')
        plt.xlim(0,1) ; plt.ylim(0,1)
        plt.text(0.05,1,summary,ha='left',va='top',fontsize=9)

    return danger1