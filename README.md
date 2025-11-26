Download the directory and try to run THE_TCS_main.py with your own Python installation.
If bug, install a Python environment:

[INSTALL by CONDA] <----- Best option

conda create -n tcs python=3.8.8 \
conda activate tcs \
pip install -r requirements.txt


[INSTALL by VENV]

python3 -m venv tcs \
source tcs/bin/activate \
pip install --upgrade pip \
pip install -r requirements.txt

[REFERENCES]

GR8 table is coming from Freckelton et al. +25 (2025yCat..75401786F)

The computation of the RV budget is made using:
 
1) ARVE (Al Moulla + 25, 2025A&A...701A.266A)
2) GP (O'Sullivan et al. in prep.)

If you have any problem, please contact me at:

michael.cretignier@physics.ox.ac.uk