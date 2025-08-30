Scripts and data for angular momentum paper

# Info on SEVN and SEVNpy version for reproducibility
Branch:             'SEVN' (main)
Last commit:         99c2ca17
Conda env name:      SEVN_25ago

This branch is the main SEVN version with small modifications to 
- allow for rlo_gamma_angmom > 1
- enforce CE and then merger for tight systems where He core radius > Semimajor (that could occur if collisions are disabled)
These modifications are included in the 25th August commits (https://gitlab.com/sevncodes/sevn/-/commit/794cb2623bc14beb089e5aae53011a109f97cff2)
+
OFFLINE, I added also the dadt outputs for SEVNpy, manually cherrypicking the commit then reverted of the 6th of August
https://gitlab.com/sevncodes/sevn/-/commit/c9d46ac9d73d8b2cde33724d10a42ef2948b7dd7

NB: To run the simulations, recall to also add manually the EvolveFunctor for MSBH (EvolveBHMSandBBH)



# Pipeline for analysis and plots

## Pipeline to reduce output dataframes
- need auxiliary libraries lib_out_extract.py, lib_log_extract.py, lib_mt_extract.py

1) out_extract.py	to extract from all output and evolved files the MSBH and BBH binaries (+eventually collect the remnants manually in a single folder)
2) log_extract.py	to extract from all logfiles the mt_history_full.csv file for all simulated binaries
3) mt_extract.py	to filter the mt information only of the binaries of interest and w.r.t. first SN (run for both label_names eventually)
4) write_summary.py	to write summary file with number of binaries undergoing a sepcific evolution, in each simulated set
5) mass_simulated.py     to get the simulated mass, then corrected and used for efficiency estimations

## Pipeline to plot large simulations
1) Correction_factor_for_efficiency.ipynb	to calculate the correction factor for each simulated set
2) plot_efficiencies.ipynb			to plot formation efficiencies, accounting to the correction factor just calculated
3) plot_hist.ipynb				to plot 1D histogram distributions
4) plot_hist2D.ipynb				to plot 2D histogram distributions in the MS-BH configuration

## Pipeline to plot subsect of binaries with SEVNpy
- need for auxiliary libraries calclib.py and plotlib.py
1) plot_ZAMSBH_grid.ipynb			to plot ZAMS-BH grids through SEVNpy
2) plot_single.ipynb				to plot single binary evolution with SEVNpy

## Other plotting scripts
1) plot_gamma.ipynb	to plot Sobermann and Pols trends
