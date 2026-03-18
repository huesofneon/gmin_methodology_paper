====================================
INSTRUCTIONS FOR RUNNING CODE
====================================
Legend: 
# = run number (3=JK, 4=DK)
$ = Sitename/Datasetname

Note: code has already been run to produce outputs/figures available on the repository.

=========R data files=========
Setup (done already):
1: "#_LoggerName.txt" is downloaded from logger to '\data\raw\'
	-Edit txt to be read by R such that the first line follows this format:
	LoggerName,Time,Celsius(C),Humidity(%rh),Dew Point(C),Serial Number

2: "v#_Min Conductance Datasheet_$.xlsx" is downloaded from google sheets to '\data\raw\'.

Running code:
3: "\massloss_run#.R" outputs '\data\mass_data_run4.rds' by pulling data from '\data\raw\' and combining:
	-drydown data (v?_Min_Conductance_Datasheet_$.xlsx), 
	-LA correction data (#_LA_Corrections_$), 
	-LMA ratio data (#_LMA_data_$.csv)
	-Leaf dry weights of samples (#_LEAF_dryweights_$.xlsx)
	-Physiological data [run 3, i.e. JK spp + ELYRHI only] (Arens_RWC_limits_JK+run4ELYRHI.xlsx)
NB: "massloss_run4.R" must be run first to obtain ELYRHI data

4: "\#_phytotron_temp.R" outputs '#_metdata_summary.rds' by pulling data from "\data\raw\#_LoggerName.txt"

5: "fitmodel_drydown_run#.R" fits models to drydown curves, outputting "modelcoeffs_gmin_run#.rds" and "df_comb_run#.rds"

6: "gmin_modelled_run#.R" calculates gmin & gmin.se across specified intervals, outputting all data in "gmin_modelled_summ_run#.rds"

=========R stats & figure files=========
7: "Fig3AB_Drydown_curves.R" plots Fig3A_ & Fig3B_ in '\figures'.

8: "Fig4_gmin_interval_comp.R" plots and compares gmin+se across all specified intervals for run3+run4ELYRHI in Fig4_ in '\figures'.

9: "Fig5_SMRI_calculation+plot.R", outputs the following in '\data\':
-"Table1_physiological_thresholds_run3+run4ELYRHI.rds"
-calculates SMRI and outputs "Table2_gmin+SMRI88_data_summ_run3+run4ELYRHI.rds"
(Note: gmin data reported for run4 species in Table3 can be found in "gmin_modelled_summ_run4.rds")

10: "Fig5_SMRI_calculation+plot.R" also performs stats analysis and plots SMRI+se for run3+run4ELYRHI as Fig5_ in '\figures'.

11: "Fig6_gmin_comp_run4.R" performs stats analysis and plots gmin_RWC80-50 for run4 species as Fig6 in '\figures'.
