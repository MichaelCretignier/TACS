[INSTALL by VENV]

python3 -m venv tcs
source tcs/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

[INSTALL by CONDA]

conda create -n tcs python=3.8.8
conda activate tcs
pip install -r requirements.txt

[REFERENCES]

GR8 Master table containing atmospheric parameters were obtained as described in: 

[ADS] https://ui.adsabs.harvard.edu/abs/2025MNRAS.540.1786F/abstract

If you're a member of THE
Download in the TACS Google Drive repository the full tables:
	1) THE_Master_table.csv
	2) THE_DACE.csv
	3) THE_YARARA.csv

The computation of the RV budget is made using:
 
1) ARVE (Al Moulla, K et al. in prep.)
2) GP (O'Sullivan et al. in prep.)

If you have any problem, please contact me at:

michael.cretignier@physics.ox.ac.uk