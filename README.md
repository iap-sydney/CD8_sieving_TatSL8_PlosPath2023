# CD8_sieving_TatSL8_PlosPath2023

To generate the figures in the manuscript, run the following codes in order. Further details of codes and what they produce are included below.
- Data_importation/save_UPenn_FTY_bcode_data_210602.m
- SL8_WT_timecourse_analysis_R_export_220318.m
- WT_timecourse_plots.R
- Barcode_1ntMutation_breakdown_all_times_code_221110.m
- SL8_Variant_timecourse_220812.R
- VL_and_WT_dynamics_at_CD8depletion_230120.m
- CD8_depletion_barcode_analysis_210831.R
- Summary_statistics_220222.m
 
Data:
- Penn SIV barcode data summary.SD_20230208.xlsx: contains high throughput sequencing results for barcode and Tat-SL8 regions
- Copy of 1st Betts FTY study 1 VL_day106 post MT807R1_SSD_necropsy_dpi_added_210426.xlsx: contains the viral load data for each animal.
- 3-9-22 TL8-PBMC.csv: contains data on total count and percent of CD8s specific for Tat-TL8 (same response as SL8).
- 3-9-22 TL8-PBMC_SSD_names_corrected.xlsx: same as "3-9-33 TL8-PBMC.csv" but also with a sheet for LNMC. xlsx can be read by MATLAB better than R.
- Immonen_etal_2020_SL8_evolution/Epitopes.xlsx: Data from Immonen et al. 2020 on SL8 epitope fractions over time.
- Allen_etal_2000/Allen_etal_2000_SL8_decay_data.csv: Data from Allen et al. 2000 on WT SL8 epitope fractions over time.
- d190-d373_Sheet_d190.csv and d190-d373_Sheet_d373.csv: contains phenotype breakdown of TL8/SL8-specific CD8s from samples taken at approximately day 190 and 373, respectively.

Note: if the cleaned data file does not exist (UPenn_FTY_barcode_SL8_Mutation_data_210607.mat), the code save_UPenn_FTY_bcode_data_210602.m (described below in "Data_importation Codes") must be run first.

Fig. 1 & 2 Functions:
- SL8_WT_timecourse_analysis_R_export_220318.m: This code loads the cleaned data saved in UPenn_FTY_barcode_SL8_Mutation_data_210607.mat and generates .csv files of WT percentages, VL, etc. for use in R (specifically WT_timecourse_plots.R).
- WT_timecourse_plots.R: visualize the time course data of WT SL8, VL, and total TL8 specific CD8s
- Time to detection panel: “Survival proportions: Survival of ATI1 vs ATI2 (FTY) percent” in “Results_Time To Detection FTY v2_SSD_edits.pzfx”
- Reactivation rate panel: “RR  ATI1,2 All monkeys” in “Results_RR FTY v3_with_FTY_control_comparison.pzfx”

Fig. 3 Functions:
- Barcode_1ntMutation_breakdown_all_times_code_221110.m: This code generates the .csv files for plotting dynamics of various Tat-SL8 variants. This code calls functions in the Poisson_model_variant_reactivation folder
- SL8_Variant_timecourse_220812.R: This code takes the output from Barcode_1ntMutation_breakdown_all_times_code_221110.m and generates the plots of percent each Tat-SL8 variant over time within an animal.
- Poisson_model_tests_221109.m: Code to fit Poisson models for detection of variants in reactivation following ATI-2

Fig. 4 Functions:
- CD8_depletion_barcode_analysis_210831.R: This code examines the barcodes detected before and after CD8 depletion. This code uses the data organised in CD8_analysis/VL_and_WT_dynamics_at_CD8depletion_230120.m Note, panels A and D are generated by WT_timecourse_plots.R and a prism file, respectively.
- Reactivation rate panel: “RR  ATI1,2 CD8 Deal. All monkeys” in “Results_RR FTY v3_with_FTY_control_comparison.pzfx”

Table S1 Function:
- CD8_TL8-spec_PBMC_d190d373_analysis_211216.R: This code compares CD8 phenotype levels during ART-1 and late in ART-2.

S2 Fig:
- “Barcode number ATI1,2 CD8 Deal. All monkeys” in “Results_RR FTY v3_with_FTY_control_comparison.pzfx”

- Summary_statistics_220222.m: This code calculates most of the statistics quoted in the manuscript. It also gathers the list of all SL8 variants ever greater than 20% of the viral load. Output is saved in Summary_statistics.mat

General_Functions Codes:
- SL8_variant_count_210607.m: code generates arrays of counts for barcode, SL8 amino acid variant, and barcode-SL8 variant combos
- background_subtraction.m: code does background subtraction based on counts and LOD provided.
- AUC_exp_between_measure_210507.m: calculates the integral of viral load (or some other data) when we assume the variable changes exponentially between time points
- AIC_calc_generic.m: Code to calculate AIC for a generic model.
- find_react_dpi.m: Code to identify dpi of reactivation (react_dpi), index of the corresponding time point in the VL data array (react_ind), and indicator of if reactivation was detected or if the animal was censored at the given time point (0 = censored, 1 = reactivation detected; react_det). Definition for reactivation is the earliest of 1) the first of at least two viral load quantifications above limit of detection or 2) viral load quantification of above 1000 copies/ml on a single day. In the case that ART was resumed, CD8s were depleted, or animal was necropsied before one reactivation was detected, animals were right censored at the day of said event.

Poisson_model_variant_reactivation Codes:
- SL8_variants_of_interest_lists_221014.m: Identifies detected SL8 amino acid variants that result from a single nucleotide mutation. Results are saved in Single_aa_muations.mat.
- SL8_single_mut_variant_count_221107.m: similar to SL8_variant_count_210607.m, but only including SL8 variants created by a single nucleotide mutation.
- loglikelihood_Poisson_restricted.m: Function to calculate Log-Likelihood of various Poisson models of reactivation based on parameters of model and model used is set by the definition vector. This function calls loglikelihood_Poisson_full.m and uses '0' for parameters not used in the current model being evaluated.
- loglikelihood_Poisson_full.m: Function to calculate Log-likelihood of Poisson model of reactivation with all terms included. Calls react_Poisson_Model_full_21118.m to calculate the probability of particular variants reactivating given parameters and variant-specific variables.
- react_Poisson_Model_full_211118.m: Calculate the probability of particular variants reactivating given parameters and variant-specific variables.

Data_importation Codes:
- save_UPenn_FTY_bcode_data_210602.m: Code to import sequencing and viral load data into MATLAB
- bcode_SL8ID_list_identification_210604.m: Code to identify unique barcodes and SL8 variants seen in sequencing. Each SL8 variant is assigned a 17-dimensional vector ID as defined above.
- bcode_SL8ID_total_count_210604.m: Code to record the count of barcodes and SL8 sequences seen in sequencing
- bcode_SL8_LOD_210604.m: generates the LODs for all barcodes, SL8 sequences, and barcode-SL8 combinations at current time point
- Best_seq_run_indices_210604.m: Code to collect indices for best sequencing runs (most templates) of short and SL8 sequencing runs separately
- Global_params_seq_coding_SL8WTaa.m: generates Global_params_seq_coding_SL8WTaa.mat, which is the look up arrays for determining which amino acid a codon corresponds to.
- aa_translation_210603.m: translates amino acid IDs between letter name and numeric name
- nt_translation_210603.m: takes a codon sequence (in either numeric or alpha form) and translates to the other form
- nt_aa_translation_210604.m: translates the codon sequence in T, C, A, G notation into the amino acid the codon codes for in numeric notation
- sort_bcode_210607.m: reorder barcodes based on prevelance over time (barcodes seen at first time point are ordered by size at that time point, then barcodes seen later, then those seen in LNMC but not Plasma, and then those seen only in PBMC
- Sequence_data_used_report_210604.m: Code to report what sequences we have used.

CD8_analysis Codes:
- VL_and_WT_dynamics_at_CD8depletion_230120.m: code gathers information about VL WT composition around CD8 depletion and the barcodes present before and after.

Notes/Definitions:
- SL8 ID vectors: 17-dimensional vectors such that i) the first coordinate is a unique identifier of the amino acid sequence (WT defined as 0), ii) the next 8 coordinates indicate the amino acid at the respective locations (the WT amino acids are indicated as 0 and mutated amino acids are numerically coded by the definition listed in aa_translation_210603.m), iii) the last 8 coordinates indicate the nucleotide sequence of the respective codons (WT codon sequences are indicated as 0, and mutated codon sequences are numerically coded by the mapping indicated in nt_translation_210603.m).
