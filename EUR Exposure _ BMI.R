# BMI -> CAD --------------------------------------------------------------

setwd('~/MR_data/ieu_opengwas/instruments')
eur_bmi_exposure_data = vcf2mr('ukb-b-19953.vcf.gz',clump=T,population='EUR')

eur_CAD_outcome_data = gwas_outcome_data('ieu-a-7',eur_bmi_exposure_data)

eur_bmiCAD_harm = harmonise_data(eur_bmi_exposure_data,eur_CAD_outcome_data, action = 3)

eur_bmiCAD_fstat <- Fstats(eur_bmiCAD_harm[eur_bmiCAD_harm$mr_keep == TRUE, ])
eur_bmiCAD_fstat

eur_bmiCAD_res = mr(eur_bmiCAD_harm)
eur_bmiCAD_res

eur_bmiCAD_egger = mr_pleiotropy_test(eur_bmiCAD_harm) # to get egger intercept
eur_bmiCAD_egger

library(MRPRESSO)
set.seed(123)
eur_bmiCAD_presso = mr_presso(BetaOutcome   = "beta.outcome",
                              BetaExposure  = "beta.exposure",
                              SdOutcome     = "se.outcome",
                              SdExposure    = "se.exposure",
                              OUTLIERtest   = TRUE,          # outlier detection
                              DISTORTIONtest = TRUE,         # distortion test
                              data           = eur_bmiCAD_harm[eur_bmiCAD_harm$mr_keep == TRUE, ],
                              NbDistribution = 1000,         # number of simulations
                              SignifThreshold = 0.05
                              )

eur_bmiCAD_presso$`Main MR results`          # global + outlier-corrected estimates
eur_bmiCAD_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eur_bmiCAD_het = mr_heterogeneity(eur_bmiCAD_harm)
eur_bmiCAD_het

eur_bmiCAD_harm <- add_metadata(eur_bmiCAD_harm)

# Data for Steiger
eur_bmiCAD_harm$samplesize.exposure = 461460 # add losing data
eur_bmiCAD_harm$r.exposure <- get_r_from_bsen(
  b  = eur_bmiCAD_harm$beta.exposure,
  se = eur_bmiCAD_harm$se.exposure,
  n  = eur_bmiCAD_harm$samplesize.exposure
)
eur_bmiCAD_harm$r.outcome <- get_r_from_lor(
  lor     = eur_bmiCAD_harm$beta.outcome,
  af      = eur_bmiCAD_harm$eaf.outcome,
  ncase   = eur_bmiCAD_harm$ncase.outcome,
  ncontrol = eur_bmiCAD_harm$ncontrol.outcome,
  prevalence = 0.06
)
eur_bmiCAD_harm$prevalence.outcome <- 0.06 # Change in different outcome

# Steiger test for each SNPs
eur_bmiCAD_steiger = steiger_filtering(eur_bmiCAD_harm[eur_bmiCAD_harm$mr_keep == TRUE, ])
nrow(eur_bmiCAD_steiger[eur_bmiCAD_steiger$steiger_dir==T,])
SNP_non_steiger = eur_bmiCAD_steiger$SNP[eur_bmiCAD_steiger$steiger_dir==F]

# Overall Steiger stats for all SNPs
eur_bmiCAD_direction <- directionality_test(eur_bmiCAD_harm[eur_bmiCAD_harm$mr_keep == TRUE, ])
eur_bmiCAD_direction

# Cause
setwd('~/MR_data/ieu_opengwas/instruments')
eur_bmi_cause = read_vcf_for_cause('ukb-b-19953.vcf.gz')
setwd('~/MR_data/EUR_new_outcome')
eur_CAD_cause = read_vcf_for_cause('ieu-a-7.vcf.gz')
eur_bmiCAD_merge = CAUSE_merge(eur_bmi_cause, eur_CAD_cause)
eur_bmiCAD_params = CAUSE_params(eur_bmiCAD_merge)
eur_bmiCAD_clump = X_clump(eur_bmiCAD_merge, 'EUR')

eur_bmiCAD_top_vars <- eur_bmiCAD_clump$snp
set.seed(123)
eur_bmiCAD_res <- cause(X=eur_bmiCAD_clump, variants = eur_bmiCAD_top_vars, param_ests = eur_bmiCAD_params)

summary(eur_bmiCAD_res, ci_size = 0.95)
class(eur_bmiCAD_res)
eur_bmiCAD_res$elpd # Zscore



# BMI -> Stroke -----------------------------------------------------------

eur_stroke_outcome_data = gwas_outcome_data('ebi-a-GCST005838', eur_bmi_exposure_data)
eur_bmistroke_harm = harmonise_data(eur_bmi_exposure_data, eur_stroke_outcome_data)
eur_bmistroke_fstat = Fstats(eur_bmistroke_harm[eur_bmistroke_harm$mr_keep == T,])
eur_bmistroke_fstat
eur_bmistroke_res = mr(eur_bmistroke_harm)
eur_bmistroke_res
eur_bmistroke_egger = mr_pleiotropy_test(eur_bmistroke_harm)
eur_bmistroke_egger

library(MRPRESSO)
set.seed(123)
eur_bmistroke_presso = mr_presso(BetaOutcome   = "beta.outcome",
                              BetaExposure  = "beta.exposure",
                              SdOutcome     = "se.outcome",
                              SdExposure    = "se.exposure",
                              OUTLIERtest   = TRUE,         
                              DISTORTIONtest = TRUE,       
                              data           = eur_bmistroke_harm[eur_bmistroke_harm$mr_keep == TRUE, ],
                              NbDistribution = 1000,        
                              SignifThreshold = 0.05)
eur_bmistroke_presso$`Main MR results`          # global + outlier-corrected estimates
eur_bmistroke_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eur_bmistroke_het = mr_heterogeneity(eur_bmistroke_harm)
eur_bmistroke_het
eur_bmistroke_harm <- add_metadata(eur_bmistroke_harm)
eur_bmistroke_harm$samplesize.exposure = 461460
eur_bmistroke_harm$r.exposure <- get_r_from_bsen(
  b  = eur_bmistroke_harm$beta.exposure,
  se = eur_bmistroke_harm$se.exposure,
  n  = eur_bmistroke_harm$samplesize.exposure
)
eur_bmistroke_harm$r.outcome <- get_r_from_lor(
  lor     = eur_bmistroke_harm$beta.outcome,
  af      = eur_bmistroke_harm$eaf.outcome,
  ncase   = eur_bmistroke_harm$ncase.outcome,
  ncontrol = eur_bmistroke_harm$ncontrol.outcome,
  prevalence = 0.09
)
eur_bmistroke_harm$prevalence.outcome <- 0.09 
eur_bmistroke_steiger = steiger_filtering(eur_bmistroke_harm[eur_bmistroke_harm$mr_keep == TRUE, ])
nrow(eur_bmistroke_steiger[eur_bmistroke_steiger$steiger_dir==T,])
SNP_non_steiger = eur_bmistroke_steiger$SNP[eur_bmistroke_steiger$steiger_dir==F]
eur_bmistroke_direction <- directionality_test(eur_bmistroke_harm[eur_bmistroke_harm$mr_keep == TRUE, ])
eur_bmistroke_direction

setwd('~/MR_data/ieu_opengwas/outcomes (full)')
eur_stroke_cause = read_vcf_for_cause('ebi-a-GCST005838.vcf.gz')
eur_bmistroke_merge = CAUSE_merge(eur_bmi_cause, eur_stroke_cause)
eur_bmistroke_params = CAUSE_params(eur_bmistroke_merge)
eur_bmistroke_clump = X_clump(eur_bmistroke_merge, 'EUR')

eur_bmistroke_top_vars <- eur_bmistroke_clump$snp
set.seed(123)
eur_bmistroke_res <- cause(X=eur_bmistroke_clump, variants = eur_bmistroke_top_vars, param_ests = eur_bmistroke_params)

summary(eur_bmistroke_res, ci_size = 0.95)
class(eur_bmistroke_res)
eur_bmistroke_res$elpd # Zscore



# BMI -> HF ---------------------------------------------------------------

setwd('~/MR_data/finngen/finngen_liftover')
eur_HF_outcome_data = finngen2mr_outcome(filepath= "summary_stats_release_finngen_R12_I9_CHD_hg19.gz",
                                         snps = eur_bmi_exposure_data$SNP,
                                         maf_threshold = 0.01, 
                                         outcome_name  = "CHD_FinnGen_R12",
                                         proxy_r2      = 0.8,
                                         proxy_kb      = 5000,
                                         proxy_pop     = "EUR")
eur_bmiHF_harm = harmonise_data(eur_bmi_exposure_data, eur_HF_outcome_data)
eur_bmiHF_fstat = Fstats(eur_bmiHF_harm[eur_bmiHF_harm$mr_keep == T,])
eur_bmiHF_fstat
eur_bmiHF_res = mr(eur_bmiHF_harm)
eur_bmiHF_res
eur_bmiHF_egger = mr_pleiotropy_test(eur_bmiHF_harm)
eur_bmiHF_egger

library(MRPRESSO)
set.seed(123)
eur_bmiHF_presso = mr_presso(BetaOutcome   = "beta.outcome",
                                 BetaExposure  = "beta.exposure",
                                 SdOutcome     = "se.outcome",
                                 SdExposure    = "se.exposure",
                                 OUTLIERtest   = TRUE,         
                                 DISTORTIONtest = TRUE,       
                                 data           = eur_bmiHF_harm[eur_bmiHF_harm$mr_keep == TRUE, ],
                                 NbDistribution = 1000,        
                                 SignifThreshold = 0.05)
eur_bmiHF_presso$`Main MR results`          # global + outlier-corrected estimates
eur_bmiHF_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eur_bmiHF_het = mr_heterogeneity(eur_bmiHF_harm)
eur_bmiHF_het
eur_bmiHF_harm <- add_metadata(eur_bmiHF_harm)
eur_bmiHF_harm$samplesize.exposure = 461460
eur_bmiHF_harm$r.exposure <- get_r_from_bsen(
  b  = eur_bmiHF_harm$beta.exposure,
  se = eur_bmiHF_harm$se.exposure,
  n  = eur_bmiHF_harm$samplesize.exposure
)
eur_bmiHF_harm$ncase.outcome = 37653
eur_bmiHF_harm$ncontrol.outcome = 462695
eur_bmiHF_harm$samplesize.outcome = 462695 + 37653
eur_bmiHF_harm$r.outcome <- get_r_from_lor(
  lor     = eur_bmiHF_harm$beta.outcome,
  af      = eur_bmiHF_harm$eaf.outcome,     
  ncase   = eur_bmiHF_harm$ncase.outcome, 
  ncontrol = eur_bmiHF_harm$ncontrol.outcome,
  prevalence = 0.02
)
eur_bmiHF_harm$prevalence.outcome <- 0.02
eur_bmiHF_steiger = steiger_filtering(eur_bmiHF_harm[eur_bmiHF_harm$mr_keep == TRUE, ])
nrow(eur_bmiHF_steiger[eur_bmiHF_steiger$steiger_dir==T,])
SNP_non_steiger = eur_bmiHF_steiger$SNP[eur_bmiHF_steiger$steiger_dir==F]
eur_bmiHF_direction <- directionality_test(eur_bmiHF_harm[eur_bmiHF_harm$mr_keep == TRUE, ])
eur_bmiHF_direction

setwd('~/MR_data/finngen/finngen_liftover')
eur_HF_cause = read_finngen_for_cause("summary_stats_release_finngen_R12_I9_CHD_hg19.gz")
eur_bmiHF_merge = CAUSE_merge(eur_bmi_cause, eur_HF_cause)
eur_bmiHF_params = CAUSE_params(eur_bmiHF_merge)
eur_bmiHF_clump = X_clump(eur_bmiHF_merge, 'EUR')

eur_bmiHF_top_vars <- eur_bmiHF_clump$snp
set.seed(123)
eur_bmiHF_res <- cause(X=eur_bmiHF_clump, variants = eur_bmiHF_top_vars, param_ests = eur_bmiHF_params)

summary(eur_bmiHF_res, ci_size = 0.95)
class(eur_bmiHF_res)
eur_bmiHF_res$elpd # Zscore



# BMI -> T2D --------------------------------------------------------------

setwd('~/MR_data/EUR_new_outcome/Mahajan_T2D')
eur_T2D_outcome_data = mahajan2mr_outcome('Mahajan_T2D_rsid.txt',
                                          snps          = eur_bmi_exposure_data$SNP,
                                          maf_threshold = 0.01,
                                          outcome_name  = "T2D_Mahajan",
                                          proxy_r2      = 0.8,
                                          proxy_kb      = 5000,
                                          proxy_pop     = "EUR")
eur_bmiT2D_harm = harmonise_data(eur_bmi_exposure_data, eur_T2D_outcome_data)
eur_bmiT2D_fstat = Fstats(eur_bmiT2D_harm[eur_bmiT2D_harm$mr_keep == T,])
eur_bmiT2D_fstat
eur_bmiT2D_res = mr(eur_bmiT2D_harm)
eur_bmiT2D_res
eur_bmiT2D_egger = mr_pleiotropy_test(eur_bmiT2D_harm)
eur_bmiT2D_egger

library(MRPRESSO)
set.seed(123)
eur_bmiT2D_presso = mr_presso(BetaOutcome   = "beta.outcome",
                             BetaExposure  = "beta.exposure",
                             SdOutcome     = "se.outcome",
                             SdExposure    = "se.exposure",
                             OUTLIERtest   = TRUE,         
                             DISTORTIONtest = TRUE,       
                             data           = eur_bmiT2D_harm[eur_bmiT2D_harm$mr_keep == TRUE, ],
                             NbDistribution = 1000,        
                             SignifThreshold = 0.05)
eur_bmiT2D_presso$`Main MR results`          # global + outlier-corrected estimates
eur_bmiT2D_presso$`MR-PRESSO results`$`Global Test`

eur_bmiT2D_het = mr_heterogeneity(eur_bmiT2D_harm)
eur_bmiT2D_het
eur_bmiT2D_harm <- add_metadata(eur_bmiT2D_harm)
eur_bmiT2D_harm$samplesize.exposure = 461460
eur_bmiT2D_harm$r.exposure <- get_r_from_bsen(
  b  = eur_bmiT2D_harm$beta.exposure,
  se = eur_bmiT2D_harm$se.exposure,
  n  = eur_bmiT2D_harm$samplesize.exposure
)
eur_bmiT2D_harm$ncase.outcome = 74124 
eur_bmiT2D_harm$ncontrol.outcome = 824006
eur_bmiT2D_harm$samplesize.outcome = 898130
eur_bmiT2D_harm$r.outcome <- get_r_from_lor(
  lor     = eur_bmiT2D_harm$beta.outcome,
  af      = eur_bmiT2D_harm$eaf.outcome,     
  ncase   = eur_bmiT2D_harm$ncase.outcome,  ######ask where to get ncase and ncontrol
  ncontrol = eur_bmiT2D_harm$ncontrol.outcome,
  prevalence = 0.10
)
eur_bmiT2D_harm$prevalence.outcome <- 0.10
eur_bmiT2D_direction <- directionality_test(eur_bmiT2D_harm[eur_bmiT2D_harm$mr_keep == TRUE, ])
eur_bmiT2D_direction

eur_T2D_cause = read_mahajan_for_cause('Mahajan_T2D_rsid.txt')
eur_bmiT2D_merge = CAUSE_merge(eur_bmi_cause, eur_T2D_cause)
eur_bmiT2D_params = CAUSE_params(eur_bmiT2D_merge)
eur_bmiT2D_clump = X_clump(eur_bmiT2D_merge, 'EUR')

eur_bmiT2D_top_vars <- eur_bmiT2D_clump$snp
set.seed(123)
eur_bmiT2D_res <- cause(X=eur_bmiT2D_clump, variants = eur_bmiT2D_top_vars, param_ests = eur_bmiT2D_params)

summary(eur_bmiT2D_res, ci_size = 0.95)
class(eur_bmiT2D_res)
eur_bmiT2D_res$elpd # Zscore


# BMI -> GFR --------------------------------------------------------------

eur_GFR_outcome_data = gwas_outcome_data('ebi-a-GCST003372',eur_bmi_exposure_data)
eur_bmiGFR_harm = harmonise_data(eur_bmi_exposure_data,eur_GFR_outcome_data, action = 3)
eur_bmiGFR_fstat <- Fstats(eur_bmiGFR_harm[eur_bmiGFR_harm$mr_keep == TRUE, ])
eur_bmiGFR_fstat
eur_bmiGFR_res = mr(eur_bmiGFR_harm)
eur_bmiGFR_res
eur_bmiGFR_egger = mr_pleiotropy_test(eur_bmiGFR_harm) # to get egger intercept
eur_bmiGFR_egger

library(MRPRESSO)
set.seed(123)
eur_bmiGFR_presso = mr_presso(BetaOutcome   = "beta.outcome",
                              BetaExposure  = "beta.exposure",
                              SdOutcome     = "se.outcome",
                              SdExposure    = "se.exposure",
                              OUTLIERtest   = TRUE,          # outlier detection
                              DISTORTIONtest = TRUE,         # distortion test
                              data           = eur_bmiGFR_harm[eur_bmiGFR_harm$mr_keep == TRUE, ],
                              NbDistribution = 1000,         # number of simulations
                              SignifThreshold = 0.05
)

eur_bmiGFR_presso$`Main MR results`          # global + outlier-corrected estimates
eur_bmiGFR_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eur_bmiGFR_het = mr_heterogeneity(eur_bmiGFR_harm)
eur_bmiGFR_het

eur_bmiGFR_harm <- add_metadata(eur_bmiGFR_harm)
eur_bmiGFR_harm$samplesize.exposure = 461460 # add losing data
eur_bmiGFR_harm$r.exposure <- get_r_from_bsen(
  b  = eur_bmiGFR_harm$beta.exposure,
  se = eur_bmiGFR_harm$se.exposure,
  n  = eur_bmiGFR_harm$samplesize.exposure
)
eur_bmiGFR_harm$r.outcome <- get_r_from_bsen(
  b  = eur_bmiGFR_harm$beta.outcome,
  se = eur_bmiGFR_harm$se.outcome,
  n  = eur_bmiGFR_harm$samplesize.outcome
)
eur_bmiGFR_direction <- directionality_test(eur_bmiGFR_harm[eur_bmiGFR_harm$mr_keep == TRUE, ])
eur_bmiGFR_direction

setwd('~/MR_data/ieu_opengwas/outcomes (full)')
eur_GFR_cause = read_vcf_for_cause('ebi-a-GCST003372.vcf.gz')
eur_bmiGFR_merge = CAUSE_merge(eur_bmi_cause, eur_GFR_cause)
eur_bmiGFR_params = CAUSE_params(eur_bmiGFR_merge)
eur_bmiGFR_clump = X_clump(eur_bmiGFR_merge, 'EUR')

eur_bmiGFR_top_vars <- eur_bmiGFR_clump$snp
set.seed(123)
eur_bmiGFR_res <- cause(X=eur_bmiGFR_clump, variants = eur_bmiGFR_top_vars, param_ests = eur_bmiGFR_params)

summary(eur_bmiGFR_res, ci_size = 0.95)
class(eur_bmiGFR_res)
eur_bmiGFR_res$elpd # Zscore


# BMI -> AF ---------------------------------------------------------------

setwd('~/MR_data/gwas_catalog')
eur_AF_outcome_data = gwastsv2mr_outcome('GCST90624413.tsv.gz',
                                          snps          = eur_bmi_exposure_data$SNP,
                                          maf_threshold = 0.01,
                                          outcome_name  = "AF_GCST90624413",
                                          proxy_r2      = 0.8,
                                          proxy_kb      = 5000,
                                          proxy_pop     = "EUR")
eur_bmiAF_harm = harmonise_data(eur_bmi_exposure_data, eur_AF_outcome_data)
eur_bmiAF_fstat = Fstats(eur_bmiAF_harm[eur_bmiAF_harm$mr_keep == T,])
eur_bmiAF_fstat
eur_bmiAF_res = mr(eur_bmiAF_harm)
eur_bmiAF_res
eur_bmiAF_egger = mr_pleiotropy_test(eur_bmiAF_harm)
eur_bmiAF_egger

library(MRPRESSO)
set.seed(123)
eur_bmiAF_presso = mr_presso(BetaOutcome   = "beta.outcome",
                              BetaExposure  = "beta.exposure",
                              SdOutcome     = "se.outcome",
                              SdExposure    = "se.exposure",
                              OUTLIERtest   = TRUE,         
                              DISTORTIONtest = TRUE,       
                              data           = eur_bmiAF_harm[eur_bmiAF_harm$mr_keep == TRUE, ],
                              NbDistribution = 1000,        
                              SignifThreshold = 0.05)
eur_bmiAF_presso$`Main MR results`          # global + outlier-corrected estimates
eur_bmiAF_presso$`MR-PRESSO results`$`Global Test`

eur_bmiAF_het = mr_heterogeneity(eur_bmiAF_harm)
eur_bmiAF_het
eur_bmiAF_harm <- add_metadata(eur_bmiAF_harm)
eur_bmiAF_harm$samplesize.exposure = 461460
eur_bmiAF_harm$r.exposure <- get_r_from_bsen(
  b  = eur_bmiAF_harm$beta.exposure,
  se = eur_bmiAF_harm$se.exposure,
  n  = eur_bmiAF_harm$samplesize.exposure
)
eur_bmiAF_harm$ncase.outcome = 192851
eur_bmiAF_harm$ncontrol.outcome = 1239541
eur_bmiAF_harm$r.outcome <- get_r_from_lor(
  lor     = eur_bmiAF_harm$beta.outcome,
  af      = eur_bmiAF_harm$eaf.outcome,     
  ncase   = eur_bmiAF_harm$ncase.outcome,  
  ncontrol = eur_bmiAF_harm$ncontrol.outcome,
  prevalence = 0.02
)
eur_bmiAF_harm$prevalence.outcome <- 0.02
eur_bmiAF_direction <- directionality_test(eur_bmiAF_harm[eur_bmiAF_harm$mr_keep == TRUE, ])
eur_bmiAF_direction

eur_AF_cause = read_gwastsv_for_cause('GCST90624413.tsv.gz')
eur_bmiAF_merge = CAUSE_merge(eur_bmi_cause, eur_AF_cause)
eur_bmiAF_params = CAUSE_params(eur_bmiAF_merge)
eur_bmiAF_clump = X_clump(eur_bmiAF_merge, 'EUR')

eur_bmiAF_top_vars <- eur_bmiAF_clump$snp
set.seed(123)
eur_bmiAF_res <- cause(X=eur_bmiAF_clump, variants = eur_bmiAF_top_vars, param_ests = eur_bmiAF_params)

summary(eur_bmiAF_res, ci_size = 0.95)
class(eur_bmiAF_res)
eur_bmiAF_res$elpd # Zscore



# BMI -> NAFLDa -----------------------------------------------------------

setwd('~/MR_data/finngen/finngen_liftover')
eur_NAFLDa_outcome_data = finngen2mr_outcome('summary_stats_release_finngen_R12_NAFLD_hg19.gz',
                                             snps          = eur_bmi_exposure_data$SNP,
                                             outcome_name  = "FinnGen_NAFLD",
                                             proxy_pop     = "EUR")
eur_bmiNAFLDa_harm = harmonise_data(eur_bmi_exposure_data, eur_NAFLDa_outcome_data)
eur_bmiNAFLDa_fstat = Fstats(eur_bmiNAFLDa_harm[eur_bmiNAFLDa_harm$mr_keep == T,])
eur_bmiNAFLDa_fstat
eur_bmiNAFLDa_res = mr(eur_bmiNAFLDa_harm)
eur_bmiNAFLDa_res
eur_bmiNAFLDa_egger = mr_pleiotropy_test(eur_bmiNAFLDa_harm)
eur_bmiNAFLDa_egger

library(MRPRESSO)
set.seed(123)
eur_bmiNAFLDa_presso = mr_presso(BetaOutcome   = "beta.outcome",
                                 BetaExposure  = "beta.exposure",
                                 SdOutcome     = "se.outcome",
                                 SdExposure    = "se.exposure",
                                 OUTLIERtest   = TRUE,         
                                 DISTORTIONtest = TRUE,       
                                 data           = eur_bmiNAFLDa_harm[eur_bmiNAFLDa_harm$mr_keep == TRUE, ],
                                 NbDistribution = 1000,        
                                 SignifThreshold = 0.05)
eur_bmiNAFLDa_presso$`Main MR results`          # global + outlier-corrected estimates
eur_bmiNAFLDa_presso$`MR-PRESSO results`$`Global Test`

eur_bmiNAFLDa_het = mr_heterogeneity(eur_bmiNAFLDa_harm)
eur_bmiNAFLDa_het
eur_bmiNAFLDa_harm <- add_metadata(eur_bmiNAFLDa_harm)
eur_bmiNAFLDa_harm$samplesize.exposure = 461460
eur_bmiNAFLDa_harm$r.exposure <- get_r_from_bsen(
  b  = eur_bmiNAFLDa_harm$beta.exposure,
  se = eur_bmiNAFLDa_harm$se.exposure,
  n  = eur_bmiNAFLDa_harm$samplesize.exposure
)
eur_bmiNAFLDa_harm$samplesize.outcome = 218792
eur_bmiNAFLDa_harm$ncase.outcome = 894
eur_bmiNAFLDa_harm$ncontrol.outcome = 217898
eur_bmiNAFLDa_harm$r.outcome <- get_r_from_lor(
  lor     = eur_bmiNAFLDa_harm$beta.outcome,
  af      = eur_bmiNAFLDa_harm$eaf.outcome,     
  ncase   = eur_bmiNAFLDa_harm$ncase.outcome,  
  ncontrol = eur_bmiNAFLDa_harm$ncontrol.outcome,
  prevalence = 0.27
)
eur_bmiNAFLDa_harm$prevalence.outcome <- 0.27
eur_bmiNAFLDa_direction <- directionality_test(eur_bmiNAFLDa_harm[eur_bmiNAFLDa_harm$mr_keep == TRUE, ])
eur_bmiNAFLDa_direction

eur_bmiNAFLDa_merge = CAUSE_merge(eur_bmi_cause, eur_NAFLDa_cause)
eur_bmiNAFLDa_params = CAUSE_params(eur_bmiNAFLDa_merge)
eur_bmiNAFLDa_clump = X_clump(eur_bmiNAFLDa_merge, 'EUR')

eur_bmiNAFLDa_top_vars <- eur_bmiNAFLDa_clump$snp
set.seed(123)
eur_bmiNAFLDa_res <- cause(X=eur_bmiNAFLDa_clump, variants = eur_bmiNAFLDa_top_vars, param_ests = eur_bmiNAFLDa_params)

summary(eur_bmiNAFLDa_res, ci_size = 0.95)
class(eur_bmiNAFLDa_res)
eur_bmiNAFLDa_res$elpd # Zscore



# BMI -> NAFLDb -----------------------------------------------------------

setwd('~/MR_data/finngen/finngen_liftover')
eur_NAFLDb_outcome_data = meta2mr_outcome('meta_analysis_ukbb_summary_stats_finngen_R12_NAFLD_meta_out_hg19.gz',
                                          snps = eur_bmi_exposure_data$SNP,
                                          outcome_name = 'FinnGen_NAFLD_meta',
                                          pop = 'EUR')
eur_bmiNAFLDb_harm = harmonise_data(eur_bmi_exposure_data, eur_NAFLDb_outcome_data)
eur_bmiNAFLDb_fstat = Fstats(eur_bmiNAFLDb_harm[eur_bmiNAFLDb_harm$mr_keep == T,])
eur_bmiNAFLDb_res = mr(eur_bmiNAFLDb_harm)
eur_bmiNAFLDb_res
eur_bmiNAFLDb_egger = mr_pleiotropy_test(eur_bmiNAFLDb_harm)
eur_bmiNAFLDb_egger

library(MRPRESSO)
set.seed(123)
eur_bmiNAFLDb_presso = mr_presso(BetaOutcome   = "beta.outcome",
                                 BetaExposure  = "beta.exposure",
                                 SdOutcome     = "se.outcome",
                                 SdExposure    = "se.exposure",
                                 OUTLIERtest   = TRUE,         
                                 DISTORTIONtest = TRUE,       
                                 data           = eur_bmiNAFLDb_harm[eur_bmiNAFLDb_harm$mr_keep == TRUE, ],
                                 NbDistribution = 1000,        
                                 SignifThreshold = 0.05)
eur_bmiNAFLDb_presso$`Main MR results`          # global + outlier-corrected estimates
eur_bmiNAFLDb_presso$`MR-PRESSO results`$`Global Test`

eur_bmiNAFLDb_het = mr_heterogeneity(eur_bmiNAFLDb_harm)
eur_bmiNAFLDb_het

# eur_bmiNAFLDb_harm <- add_metadata(eur_bmiNAFLDb_harm)
# eur_bmiNAFLDb_harm$samplesize.exposure = 461460
# eur_bmiNAFLDb_harm$r.exposure <- get_r_from_bsen(
#   b  = eur_bmiNAFLDb_harm$beta.exposure,
#   se = eur_bmiNAFLDb_harm$se.exposure,
#   n  = eur_bmiNAFLDb_harm$samplesize.exposure
# )
# eur_bmiNAFLDb_harm$samplesize.outcome = 218792
# eur_bmiNAFLDb_harm$ncase.outcome = 894
# eur_bmiNAFLDb_harm$ncontrol.outcome = 217898
# eur_bmiNAFLDb_harm$r.outcome <- get_r_from_lor(
#   lor     = eur_bmiNAFLDb_harm$beta.outcome,
#   af      = eur_bmiNAFLDb_harm$eaf.outcome,     
#   ncase   = eur_bmiNAFLDb_harm$ncase.outcome,  
#   ncontrol = eur_bmiNAFLDb_harm$ncontrol.outcome,
#   prevalence = 0.27
# )
# eur_bmiNAFLDb_harm$prevalence.outcome <- 0.27
# eur_bmiNAFLDb_direction <- directionality_test(eur_bmiNAFLDb_harm[eur_bmiNAFLDb_harm$mr_keep == TRUE, ])
# eur_bmiNAFLDb_direction

setwd('~/MR_data/finngen/finngen_liftover')
eur_NAFLDb_cause = read_meta_for_cause('meta_analysis_ukbb_summary_stats_finngen_R12_NAFLD_meta_out_hg19.gz')
eur_bmiNAFLDb_merge = CAUSE_merge(eur_bmi_cause, eur_NAFLDb_cause)
eur_bmiNAFLDb_params = CAUSE_params(eur_bmiNAFLDb_merge)
eur_bmiNAFLDb_clump = X_clump(eur_bmiNAFLDb_merge, 'EUR')

eur_bmiNAFLDb_top_vars <- eur_bmiNAFLDb_clump$snp
set.seed(123)
eur_bmiNAFLDb_res <- cause(X=eur_bmiNAFLDb_clump, variants = eur_bmiNAFLDb_top_vars, param_ests = eur_bmiNAFLDb_params)

summary(eur_bmiNAFLDb_res, ci_size = 0.95)
class(eur_bmiNAFLDb_res)
eur_bmiNAFLDb_res$elpd # Zscore



# BMI -> NAFLDc -----------------------------------------------------------

setwd('~/MR_data/gwas_catalog')
eur_NAFLDc_outcome_data = gwastsv2mr_outcome('GCST90091033_buildGRCh37.tsv.gz',
                                             snps = eur_bmi_exposure_data$SNP,
                                             outcome_name = 'NAFLD_GOLD',
                                             proxy_pop = 'EUR')
eur_bmiNAFLDc_harm = harmonise_data(eur_bmi_exposure_data, eur_NAFLDc_outcome_data)
eur_bmiNAFLDc_fstat = Fstats(eur_bmiNAFLDc_harm[eur_bmiNAFLDc_harm$mr_keep == T,])
eur_bmiNAFLDc_res = mr(eur_bmiNAFLDc_harm)
eur_bmiNAFLDc_res
eur_bmiNAFLDc_egger = mr_pleiotropy_test(eur_bmiNAFLDc_harm)
eur_bmiNAFLDc_egger

library(MRPRESSO)
set.seed(123)
eur_bmiNAFLDc_presso = mr_presso(BetaOutcome   = "beta.outcome",
                                 BetaExposure  = "beta.exposure",
                                 SdOutcome     = "se.outcome",
                                 SdExposure    = "se.exposure",
                                 OUTLIERtest   = TRUE,         
                                 DISTORTIONtest = TRUE,       
                                 data           = eur_bmiNAFLDc_harm[eur_bmiNAFLDc_harm$mr_keep == TRUE, ],
                                 NbDistribution = 1000,        
                                 SignifThreshold = 0.05)
eur_bmiNAFLDc_presso$`Main MR results`          # global + outlier-corrected estimates
eur_bmiNAFLDc_presso$`MR-PRESSO results`$`Global Test`

eur_bmiNAFLDc_het = mr_heterogeneity(eur_bmiNAFLDc_harm)
eur_bmiNAFLDc_het
eur_bmiNAFLDc_harm <- add_metadata(eur_bmiNAFLDc_harm)
eur_bmiNAFLDc_harm$samplesize.exposure = 461460
eur_bmiNAFLDc_harm$r.exposure <- get_r_from_bsen(
  b  = eur_bmiNAFLDc_harm$beta.exposure,
  se = eur_bmiNAFLDc_harm$se.exposure,
  n  = eur_bmiNAFLDc_harm$samplesize.exposure
)
eur_bmiNAFLDc_harm$samplesize.outcome = 778614
eur_bmiNAFLDc_harm$ncase.outcome = 8434 
eur_bmiNAFLDc_harm$ncontrol.outcome= 770180 
eur_bmiNAFLDc_harm$r.outcome <- get_r_from_lor(
  lor     = eur_bmiNAFLDc_harm$beta.outcome,
  af      = eur_bmiNAFLDc_harm$eaf.outcome,     
  ncase   = eur_bmiNAFLDc_harm$ncase.outcome,  
  ncontrol = eur_bmiNAFLDc_harm$ncontrol.outcome,
  prevalence = 0.27
)
eur_bmiNAFLDc_harm$prevalence.outcome <- 0.27
eur_bmiNAFLDc_direction <- directionality_test(eur_bmiNAFLDc_harm[eur_bmiNAFLDc_harm$mr_keep == TRUE, ])
eur_bmiNAFLDc_direction

setwd('~/MR_data/gwas_catalog')
eur_NAFLDc_cause = read_gwastsv_for_cause_2.0('GCST90091033_buildGRCh37.tsv.gz')
eur_bmiNAFLDc_merge = CAUSE_merge(eur_bmi_cause, eur_NAFLDc_cause)
eur_bmiNAFLDc_params = CAUSE_params(eur_bmiNAFLDc_merge)
eur_bmiNAFLDc_clump = X_clump(eur_bmiNAFLDc_merge, 'EUR')

eur_bmiNAFLDc_top_vars <- eur_bmiNAFLDc_clump$snp
set.seed(123)
eur_bmiNAFLDc_res <- cause(X=eur_bmiNAFLDc_clump, variants = eur_bmiNAFLDc_top_vars, param_ests = eur_bmiNAFLDc_params)

summary(eur_bmiNAFLDc_res, ci_size = 0.95)
class(eur_bmiNAFLDc_res)
eur_bmiNAFLDc_res$elpd # Zscore