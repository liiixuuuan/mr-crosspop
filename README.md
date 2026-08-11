# mr-crosspop

### Project Title
Bidirectional Cross-Population Mendelian Randomization of Cardiovascular-Kidney-Metabolic Syndrome

### Project Description
This project is a MR study that uses GWAS summary statistics to evaluate the directionality and potential causal effects of SNP-associated traits on CKM disease outcomes. Two-sample MR was applied across all 182 directed pairs among CKM nodes (PREVENT risk factors plus MASLD, eGFR/CKD, CAD, HF). IVW, MR-Egger, weighted median, MR-PRESSO, CAUSE, Steiger directionality were tested for robustness.

### Data Source
'Data List.xlsx' listed out the data used in this project.

### Example MR Analysis
'EUR Exposure _ BMI.R' and 'EAS Exposure _ BMI.R' are example scripts demonstrating how the Mendelian Randomization (MR) analysis pipeline is performed for European (EUR) and East Asian (EAS) populations, respectively.

### Results of MR pipeline
'Phase 1 Results.xlsx' listed out the pairs tested. 
'Tier1_robust_pairs.xlsx' included those that passed IVW Bonferroni significant, CAUSE analysis z > 1.96 AND 95% CI of γ excludes 0, Steiger directionality, and MR-Egger intercept p-value > 0.05.
'Tier2_discovery_pairs.xlsx' included those that only passed IVW FDR q-value < 0.05, CAUSE analysis z > 1.96 AND 95% CI of γ excludes 0, but not Tier 1 requirement.

### Data Visualisation
'UVMR_forest_EUR_IVW_all_pairs.png' and 'UVMR_forest_EAS_IVW_all_pairs.png' showed IVW results of all the tested pairs. (Function script: IVW Overall Forest Plot.R)
'Tier1_Tier2_forest_EUR_IVW_Egger_WM.png' and Tier1_Tier2_forest_EAS_IVW_Egger_WM.png' showed results of sensitivity test in Tier 1 and Tier 2 pairs. (Function script: MR Forest Plot.R)
'MR_network_EUR.png' and 'MR_network_EAS.png' showed causal CKM network in Tier 1 and Tier 2 pairs. (Function script: Phase 2 Script.R)
