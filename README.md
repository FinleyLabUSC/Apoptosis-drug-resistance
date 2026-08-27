ODE model for simulating caspase-mediated signaling, calibrated to published experimental data from the prostate cancer cell line, PC3. The model is implemented in MATLAB.

Accompanies manuscript entitled "Mechanistic Modeling of Intrinsic Drug Resistance in Prostate Cancer Apoptosis Signaling".

Fitted parameter set that gives the best fit to the experimental data is given in "gp_best.mat". The error for this best parameter set is given in "gf_best.mat". 

All of the 61 fitted parameter sets in which none of the parameters hit the search bounds are given in "gp_best_All.mat". The errors for these best-fit parameter sets are given in "gf_best_All.mat".

These fitted parameters are used to produce the main results of the manuscript: 
- Figure 6: "run_model_best_params_withTreatment.m"
- Figures 7 and 8: "run_model_best_params_VaryDrugStrength.m"
- Figures 9 and 10: run_model_best_params_TreatmentHeterogeneity_AllDrugs.m"

The reaction rates and ordinary differential equations for the model are given in "coreFile_DISC.m".
