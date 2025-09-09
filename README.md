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

[IMPORTANT]

THE members should download the Master_table_v2.0.csv on the Google drive

[REFERENCES]

The computation of the RV budget is made using:
 
1) ARVE (Al Moulla, K et al. in prep.)
2) GP (O'Sullivan et al. in prep.)

If you have any problem, please contact me at:

michael.cretignier@physics.ox.ac.uk