# TaCS (Targets Characterisation and Selection) for the Terra Hunting Experiment (v2.02)

<p align="center">
  <img src="logo.png" alt="Project logo" width="400">
</p>

## Contact Me

If you have any problem, please contact me at:

michael.cretignier@physics.ox.ac.uk

## Installation

Download the directory and try to run `THE_TCS_main.py` with your own Python installation.
If it crashes, install a Python environment:

 [Mac M4 Chip] Python environment (Conda install) (Python 3.12.5)

```bash
conda create -n tcs -c conda-forge python=3.12.5 numpy=1.26.4 pandas=2.3.2 scipy=1.16.2 astropy=7.1.0 matplotlib=3.10.6 ipython=9.5.0 colorama=0.4.4 pyqt=5.15 -y 
```

 [Mac M2 Chip] Python environment (Conda forge install) (Python 3.10)

```bash
conda create -n tcs -c conda-forge python=3.10.0 numpy=1.23.5 pandas=1.4.1 scipy=1.8 astropy=5.2.2 matplotlib=3.5 ipython=8.11.0 colorama=0.4.4 pyqt=5.15 -y
```

 [Mac Intel Chip] Python environment (Conda install) (Python 3.8.8)

```bash
conda create -n tcs python=3.8.8 
conda activate tcs 
pip install -r requirements_3.8.8.txt
```

[Alternative to conda] Python environment (Venv install)

```bash
python3 -m venv tcs 
source tcs/bin/activate 
pip install --upgrade pip 
pip install -r requirements_3.8.8.txt
```

## Test file

Let's run the main test file `THE_TCS_main.py` containing all the useful features in a iPython shell (with the local python environment if needed). Move inside the `TACS` directory:

```bash
conda activate tcs 
cd TACS/
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

## ⑤ Tutorial

### Step-by-step

```bash
cd .../GitHub/TACS/
```

*Launch an IPython shell:*

```bash
ipython
```

*This information is specified by:*

```python
import tacs

#let's initiate a tacs object
presurvey = tacs.tcs(version='5.2') 
```

*A lot have already be done from here! In `tacs`, all the relevant information are stored in attributes that all started with `.info_`*

1) .info_SC (scalars values)
2) .info_XY (time-series)
3) .info_IM (images)
4) .info_TA (tables)

*Let's check the tables for the Terra Hunting available*
```python
# presurvey.info_ #and then press "tab"
print(presurvey.info_TA_stars_selected)
```
*This is a dictionary that contains several tables initiated by tacs when calling `.tcs()`*
```python
# presurvey.info_ #and then press "tab"
print(presurvey.info_TA_stars_selected)
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