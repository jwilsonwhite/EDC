# EDC
Population models for organisms exposed to endocrine disrupting compounds

This code performs the analyses described in White, Cole, Cherr, Connon, and Brander - “Scaling up the individual-level effects of endocrine disruptors: how many males does a population need?” (submitted to Ecological Applications)

File descriptions:

beta_fert_curve.m creates Fig. 1 in the manuscript

BH.m returns the Beverton-Holt stock-recruit function

EDC_params.m contains the model parameters

EDC_popmodel_as.m runs the age-structured population model. It is designed to be called by EDC_runs.m

EDC_runs.m calls EDC_popmodel_as.m to run the various model analyses described in the manuscript, and also create the figures.

m_comb.m generates the various combinations of genotype-specific matings presented in manuscript Table 2.

Summary_data_from_BCole.xlsx is a summary of field sex ratio data. Also shown in Appendix S1.

Directory Spawning_trials:

spawning_fr.m Fits a Beta CDF to the laboratory spawning data

Data_Oct2013 contains datasheets with the spawning trial data used in spawning_fr.m.

Eggs_F2013.csv describes how many eggs were produced in each trial
FishIDs_F2013.csv lists the sex, tag ID, length, and mass of each fish used in the trials.
TrialIDs_F2013.csv lists which fish (by tag ID) were used in each trial replicate.
