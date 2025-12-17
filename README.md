# Contact Me

If you have any problem, please contact me at:

michael.cretignier@physics.ox.ac.uk

# Installation

Download the directory and try to run THE_TCS_main.py with your own Python installation.
If bug, install a Python environment:

# Python environment (Conda install) <----- Best option

[TERMINAL] conda create -n tcs python=3.8.8 \
[TERMINAL]conda activate tcs \
[TERMINAL]pip install -r requirements.txt

# Python environment (Venv install)

[TERMINAL] python3 -m venv tcs \
[TERMINAL] source tcs/bin/activate \
[TERMINAL] pip install --upgrade pip \
[TERMINAL] pip install -r requirements.txt

# Tutorial (run the main.py in a iPython shell)

[TERMINAL] cd ../TACS/Python \
[TERMINAL] ipython
[IPYTHON] run THE_TCS_main.py

# References

GR8 table is coming from Freckelton et al. +25 (2025yCat..75401786F)

The computation of the RV budget is made using:
 
1) ARVE (Al Moulla + 25, 2025A&A...701A.266A)
2) GP (O'Sullivan et al. in prep.)
3) ExTEMPO (https://github.com/BRajkumar041992/ExTEMPO)

