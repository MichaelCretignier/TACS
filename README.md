# TaCS (Targets Characterisation and Selection) for the Terra Hunting Experiment (v2.02)

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

### Step-by-step

```bash
cd .../GitHub/TACS/
```

*Launch an IPython shell:*

```bash
ipython
```

*Let's initiate the a tacs object:*

```python
import tacs

presurvey = tacs.tcs(version='5.2') #version of the catalogue 
```

*A lot have already be done from here! \
In `tacs`, all the relevant information are stored in attributes that all started with `.info_`*

1) `.info_SC_` (scalars values)
2) `.info_XY_` (time-series)
3) `.info_IM_` (images)
4) `.info_TA_` (tables)

*Let's check the tables for the Terra Hunting available already:*
```python
# presurvey.info_ #and then press "tab" to see all the products available
print(presurvey.info_TA_stars_selected)
```
*This is a dictionary that contains several tables initiated by `tacs` when calling `.tcs()`*

```python
for kw in presurvey.info_TA_stars_selected.keys():
  print(kw)
```
*A lot of TACS tables already exists! The most important are:*

1) `.info_TA_stars_selected['GR8']`        # the initial sample of 1418 stars
2) `.info_TA_stars_selected['solartwins']` # the solar twins sample
3) `.info_TA_stars_selected['RVopti']`     # RV optimised sample
4) `.info_TA_stars_selected['presurvey']`  # Union of solartwins and RVopti 

*To access the table, just get the .data attribute:*

```python
table = presurvey.info_TA_stars_selected['presurvey'].data
print(table)
```

## References

GR8 table is coming from Freckelton et al. +25 (2025yCat..75401786F) \
https://ui.adsabs.harvard.edu/abs/2025yCat..75401786F/abstract

The computation of the RV budget is made using:
 
1) ARVE (Al Moulla + 25, 2025A&A...701A.266A)  \
(https://ui.adsabs.harvard.edu/abs/2025A%26A...701A.266A/abstract)
2) GP (O'Sullivan et al. in prep.)
3) ExTEMPO (https://github.com/BRajkumar041992/ExTEMPO)

## Uninstall

```bash
[TERMINAL] 
conda remove --name tcs --all
```