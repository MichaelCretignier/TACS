# TaCS (Targets Characterisation and Selection) for the Terra Hunting Experiment (v2.10)

<p align="center">
  <img src="logo.png" alt="Project logo" width="400">
</p>

## ⓵ Contact Me

If you have any problem, please contact me at:

michael.cretignier@physics.ox.ac.uk

## ⓶ Installation

Download the directory and try to run `THE_TCS_main.py` with your own Python installation.
If it crashes, install a Python environment:

### [Option 1] Conda install (Python 3.12.5)

#### [Mac M4 Chip] (Python 3.12.5)

```bash
conda create -n tcs -c conda-forge python=3.12.5 numpy=1.26.4 pandas=2.3.2 scipy=1.16.2 astropy=7.1.0 matplotlib=3.10.6 ipython=9.5.0 colorama=0.4.4 pyqt=5.15 -y 
```

#### [Mac M2 Chip] (Python 3.10)

```bash
conda create -n tcs -c conda-forge python=3.10.0 numpy=1.23.5 pandas=1.4.1 scipy=1.8 astropy=5.2.2 matplotlib=3.5 ipython=8.11.0 colorama=0.4.4 pyqt=5.15 -y
```

#### [Mac Intel Chip] (Python 3.8.8)

```bash
conda create -n tcs python=3.8.8 
conda activate tcs 
pip install -r requirements_3.8.8.txt
```

*Check if the snaky environment exists and is active:*

```bash
conda env list
```

### [Option 2] Venv install (Python 3.8.8)

```bash
python3 -m venv tcs 
source tcs/bin/activate 
pip install --upgrade pip 
pip install -r requirements_3.8.8.txt
```

## ⓷ Test file

Let's run the main test file `THE_TCS_main.py` containing all the useful features in a iPython shell (with the local python environment if needed). Move inside the `TACS` directory:

```bash
conda activate tcs 
cd .../GitHub/TACS/
```

Now launch iPython:

```bash
ipython 
```

And run the test file.

```python
run THE_TCS_main.py
```

A lot of figures will pop but just check you don't get any error message.

## ⓸ Tutorial

### Getting access to the target selection tables

*Let's enter the TACS Git Clone directory*

```bash
cd .../GitHub/TACS/
```

*Launch an IPython shell:*

```bash
ipython
```

*Let's initiate a `.tcs()` object:*

```python
import tacs

presurvey = tacs.tcs(version='5.2') #version of the catalogue 
```

*A lot have already be done from here! \
In `tacs`, all the relevant information are stored in attributes that all started with the `.info_` prefix:*

1) `.info_SC_` (scalars values)
2) `.info_XY_` (time-series)
3) `.info_IM_` (images)
4) `.info_TA_` (tables)

*Let's check the tables for the Terra Hunting Experiment available already:*
```python
# presurvey.info_ #and then press "tab" to see all the products available
print(presurvey.info_TA_stars_selected)
```
*This is a dictionary that contains several tables initiated by `tacs` when calling `.tcs()`*

```python
for kw in presurvey.info_TA_stars_selected.keys():
  print(kw)
```
*A lot of TACS tables already exists! Those are defined in `THE_TCS_variables.py`. The most important are:*

1) `.info_TA_stars_selected['GR8']`        (the initial sample of 1418 stars)
2) `.info_TA_stars_selected['solartwins']` (the solar twins sample)
3) `.info_TA_stars_selected['RVopti']`     (the RV optimised sample)
4) `.info_TA_stars_selected['presurvey']`  (the Union of solartwins and RVopti)

*To access the table, just get the `.data` attribute:*

```python
table = presurvey.info_TA_stars_selected['presurvey'].data
print(table)

#why not plotting a selection in sky?
presurvey.info_TA_stars_selected['presurvey'].plot_space_mission(newfig=False)

```

### Accessing peculiar star-by-star information


*Let's start with the most basic information that is the stellar names known by TACS:*

```python
starnames = tacs.get_info_starname('51Peg') #aka HD217014
```

*`tacs` can work with a large variety of names conventions. It's however recommand to use the HD one for simplicity.*

*We can ask the code if a given star still belongs to a given selection of stars:*

```python
cutoff_suntwins = presurvey.info_TA_cutoff['solartwins'].copy()
tacs.which_cutoff('HD217014', cutoff_suntwins)

#Nice!
#What about this star?

tacs.which_cutoff('HD16160', cutoff_suntwins)

#Indeed, it is not a solar twin and it has other issues as well.
```

*Remark: you can't use the `.which_cutoff()` method on the 'presurvey' since the 'presurvey' is the union of two samples without a proper cutoff list.*

*Let's plot all the information collected by the TACS for a given target:*

```python
selection = presurvey.info_TA_stars_selected['presurvey'].data.copy()
cutoff = presurvey.info_TA_cutoff['RVopti']

tacs.plot_summary('HD16160',selection=selection,cutoff=cutoff)
```

### Playing with star visibility

*While `tacs` was primarly developed for THE, it contains useful functions going well beyond THE mission. One of them is the stellar visibility.*

*Let's compute the visibility of a star for different spectrograph:*

```python
plt.figure(figsize=(18,5))
for n,ins in enumerate(['HARPS3','HARPS','NEID','ESPRESSO','KPF','SOPHIE']):
    star = tacs.tcs(sun_elevation=-12, instrument=ins)
    star.set_star(ra=8,dec=5) # specify the RA and DEC in degree
    plt.subplot(1,6,n+1) ; star.compute_nights(airmass_max=1.5, weather=False, plot=True) ; plt.title(ins)
plt.subplots_adjust(left=0.05,right=0.96)
```

*Interested to know the stars with the longest night duration in April for HARPS3?:*

```python
star3 = tacs.tcs(sun_elevation=-6, instrument='HARPS3')
star3.compute_SG_calendar(sun_elevation=-6, airmass_max=1.75, alpha_step=5, dec_step=5)

star3.compute_SG_month(month=4,plot=True) # April
```

*Want to know the standard stars of HARPS3 along the year?*

```python
standards1 = tacs.THE_standards 
star1 = tacs.tcs(sun_elevation=-6) 
plt.figure(figsize=(16,8))
s1 = plt.subplot(1,1,1)
for n,s in enumerate(standards1): 
    star1.set_star(starname=s,verbose=False)
    star1.plot_night_length(figure=s1,legend=False,airmass_max=[1.5],sun_elevation=[-12, -18],color='C%.0f'%(n),showname=True) #peak in April
    plt.ylim(-1,10)
plt.subplots_adjust(hspace=0.45,top=0.95,bottom=0.10)
```

## References

GR8 table is coming from [Freckelton et al. +25](https://ui.adsabs.harvard.edu/abs/2025yCat..75401786F/abstract).

The computation of the RV budget is made using:
 
1) ARVE ([Al Moulla + 25](https://ui.adsabs.harvard.edu/abs/2025A%26A...701A.266A/abstract))
2) GP ([O'Sullivan et al. in prep.]())
3) ExTEMPO ([Rajukar et al. in prep.](https://github.com/BRajkumar041992/ExTEMPO))

## Uninstall

```bash
[TERMINAL] 
conda remove --name tcs --all
```