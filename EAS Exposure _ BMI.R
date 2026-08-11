# BMI -> CAD --------------------------------------------------------------

setwd('~/MR_data/bbj/hum0197.v3.BBJ.BMI.v1/hum0197.v3.BBJ.BMI.v1')
eas_bmi_exposure_data = sakaue2mr('GWASsummary_BMI_Japanese_SakaueKanai2020.auto.txt.gz',
                                  clump = T,
                                  population = 'EAS')

setwd('~/MR_data/bbj')
eas_CAD_outcome_data = hum2mr_outcome('hum0197.v5.gwas.CAD.v1.txt.gz',
                                      snps       = eas_bmi_exposure_data$SNP,
                                      trait_name = "CAD_Japanese")

eas_bmiCAD_harm = harmonise_data(eas_bmi_exposure_data,eas_CAD_outcome_data, action = 3)

eas_bmiCAD_fstat <- Fstats(eas_bmiCAD_harm[eas_bmiCAD_harm$mr_keep == TRUE, ])

eas_bmiCAD_res = mr(eas_bmiCAD_harm)
eas_bmiCAD_res

eas_bmiCAD_egger = mr_pleiotropy_test(eas_bmiCAD_harm) # to get egger intercept
eas_bmiCAD_egger

library(MRPRESSO)
set.seed(123)
eas_bmiCAD_presso = mr_presso(BetaOutcome   = "beta.outcome",
                              BetaExposure  = "beta.exposure",
                              SdOutcome     = "se.outcome",
                              SdExposure    = "se.exposure",
                              OUTLIERtest   = TRUE,          # outlier detection
                              DISTORTIONtest = TRUE,         # distortion test
                              data           = eas_bmiCAD_harm[eas_bmiCAD_harm$mr_keep == TRUE, ],
                              NbDistribution = 1000,         # number of simulations
                              SignifThreshold = 0.05
)

eas_bmiCAD_presso$`Main MR results`          # global + outlier-corrected estimates
eas_bmiCAD_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eas_bmiCAD_het = mr_heterogeneity(eas_bmiCAD_harm)
eas_bmiCAD_het

eas_bmiCAD_harm <- add_metadata(eas_bmiCAD_harm)

# Data for Steiger
eas_bmiCAD_harm$samplesize.exposure =163836 # add losing data
eas_bmiCAD_harm$ncase.outcome = 32512
eas_bmiCAD_harm$ncontrol.outcome = 146214
eas_bmiCAD_harm$r.exposure <- get_r_from_bsen(
  b  = eas_bmiCAD_harm$beta.exposure,
  se = eas_bmiCAD_harm$se.exposure,
  n  = eas_bmiCAD_harm$samplesize.exposure
)
eas_bmiCAD_harm$r.outcome <- get_r_from_lor(
  lor     = eas_bmiCAD_harm$beta.outcome,
  af      = eas_bmiCAD_harm$eaf.outcome,
  ncase   = eas_bmiCAD_harm$ncase.outcome,
  ncontrol = eas_bmiCAD_harm$ncontrol.outcome,
  prevalence = 0.06
)
eas_bmiCAD_harm$prevalence.outcome <- 0.06 # Change in different outcome

eas_bmiCAD_direction <- directionality_test(eas_bmiCAD_harm[eas_bmiCAD_harm$mr_keep == TRUE, ])
eas_bmiCAD_direction

setwd('~/MR_data/bbj/hum0197.v3.BBJ.BMI.v1/hum0197.v3.BBJ.BMI.v1')
eas_bmi_cause = sakaue2cause('GWASsummary_BMI_Japanese_SakaueKanai2020.auto.txt.gz')
setwd('~/MR_data/bbj')
eas_CAD_cause = hum2cause('hum0197.v5.gwas.CAD.v1.txt.gz')
eas_bmiCAD_merge = CAUSE_merge(eas_bmi_cause, eas_CAD_cause)
eas_bmiCAD_params = CAUSE_params(eas_bmiCAD_merge)
eas_bmiCAD_clump = X_clump(eas_bmiCAD_merge, 'EAS')

eas_bmiCAD_top_vars <- eas_bmiCAD_clump$snp
set.seed(123)
eas_bmiCAD_cause_res <- cause(X=eas_bmiCAD_clump, variants = eas_bmiCAD_top_vars, param_ests = eas_bmiCAD_params)

summary(eas_bmiCAD_cause_res, ci_size = 0.95)
class(eas_bmiCAD_cause_res)
eas_bmiCAD_cause_res$elpd # Zscore



# BMI -> Stroke -----------------------------------------------------------

setwd('~/MR_data/bbj/hum0197.v3.BBJ.IS.v1/hum0197.v3.BBJ.IS.v1')
eas_stroke_outcome_data = hum2mr_outcome('GWASsummary_IS_Japanese_SakaueKanai2020.auto.txt.gz',
                                            snps       = eas_bmi_exposure_data$SNP,
                                            trait_name = "Stroke_IS_Japanese")
eas_bmistroke_harm = harmonise_data(eas_bmi_exposure_data, eas_stroke_outcome_data)
eas_bmistroke_fstat = Fstats(eas_bmistroke_harm[eas_bmistroke_harm$mr_keep == T,])
eas_bmistroke_res = mr(eas_bmistroke_harm)
eas_bmistroke_res
eas_bmistroke_egger = mr_pleiotropy_test(eas_bmistroke_harm)
eas_bmistroke_egger

library(MRPRESSO)
set.seed(123)
eas_bmistroke_presso = mr_presso(BetaOutcome   = "beta.outcome",
                                 BetaExposure  = "beta.exposure",
                                 SdOutcome     = "se.outcome",
                                 SdExposure    = "se.exposure",
                                 OUTLIERtest   = TRUE,         
                                 DISTORTIONtest = TRUE,       
                                 data           = eas_bmistroke_harm[eas_bmistroke_harm$mr_keep == TRUE, ],
                                 NbDistribution = 1000,        
                                 SignifThreshold = 0.05)
eas_bmistroke_presso$`Main MR results`          # global + outlier-corrected estimates
eas_bmistroke_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eas_bmistroke_het = mr_heterogeneity(eas_bmistroke_harm)
eas_bmistroke_het
eas_bmistroke_harm <- add_metadata(eas_bmistroke_harm)
eas_bmistroke_harm$samplesize.exposure =163836 # add losing data
eas_bmistroke_harm$ncase.outcome = 22664
eas_bmistroke_harm$ncontrol.outcome = 152022
eas_bmistroke_harm$r.exposure <- get_r_from_bsen(
  b  = eas_bmistroke_harm$beta.exposure,
  se = eas_bmistroke_harm$se.exposure,
  n  = eas_bmistroke_harm$samplesize.exposure
)
eas_bmistroke_harm$r.outcome <- get_r_from_lor(
  lor     = eas_bmistroke_harm$beta.outcome,
  af      = eas_bmistroke_harm$eaf.outcome,
  ncase   = eas_bmistroke_harm$ncase.outcome,
  ncontrol = eas_bmistroke_harm$ncontrol.outcome,
  prevalence = 0.002
)
eas_bmistroke_harm$prevalence.outcome <- 0.002
eas_bmistroke_direction <- directionality_test(eas_bmistroke_harm[eas_bmistroke_harm$mr_keep == TRUE, ])
eas_bmistroke_direction

eas_stroke_cause = hum2cause('GWASsummary_IS_Japanese_SakaueKanai2020.auto.txt.gz')
eas_bmistroke_merge = CAUSE_merge(eas_bmi_cause, eas_stroke_cause)
eas_bmistroke_params = CAUSE_params(eas_bmistroke_merge)
eas_bmistroke_clump = X_clump(eas_bmistroke_merge, 'EAS')

eas_bmistroke_top_vars <- eas_bmistroke_clump$snp
set.seed(123)
eas_bmistroke_cause_res <- cause(X=eas_bmistroke_clump, variants = eas_bmistroke_top_vars, param_ests = eas_bmistroke_params)

summary(eas_bmistroke_cause_res, ci_size = 0.95)
class(eas_bmistroke_cause_res)
eas_bmistroke_cause_res$elpd # Zscore



# BMI -> HF ---------------------------------------------------------------

setwd('~/MR_data/bbj/hum0197.v3.BBJ.CHF.v1/hum0197.v3.BBJ.CHF.v1')
eas_HF_outcome_data = hum2mr_outcome('GWASsummary_CHF_Japanese_SakaueKanai2020.auto.txt.gz',
                                      snps       = eas_BMI_exposure_data$SNP,
                                      trait_name = "HF_Japanese")
eas_bmiHF_harm = harmonise_data(eas_bmi_exposure_data, eas_HF_outcome_data)
eas_bmiHF_fstat = Fstats(eas_bmiHF_harm[eas_bmiHF_harm$mr_keep == T,])
eas_bmiHF_res = mr(eas_bmiHF_harm)
eas_bmiHF_res
eas_bmiHF_egger = mr_pleiotropy_test(eas_bmiHF_harm)
eas_bmiHF_egger

library(MRPRESSO)
set.seed(123)
eas_bmiHF_presso = mr_presso(BetaOutcome   = "beta.outcome",
                                 BetaExposure  = "beta.exposure",
                                 SdOutcome     = "se.outcome",
                                 SdExposure    = "se.exposure",
                                 OUTLIERtest   = TRUE,         
                                 DISTORTIONtest = TRUE,       
                                 data           = eas_bmiHF_harm[eas_bmiHF_harm$mr_keep == TRUE, ],
                                 NbDistribution = 1000,        
                                 SignifThreshold = 0.05)
eas_bmiHF_presso$`Main MR results`          # global + outlier-corrected estimates
eas_bmiHF_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eas_bmiHF_het = mr_heterogeneity(eas_bmiHF_harm)
eas_bmiHF_het
eas_bmiHF_harm <- add_metadata(eas_bmiHF_harm)
eas_bmiHF_harm$samplesize.exposure = 163836 # add losing data
eas_bmiHF_harm$ncase.outcome = 10540 
eas_bmiHF_harm$ncontrol.outcome = 168186
eas_bmiHF_harm$r.exposure <- get_r_from_bsen(
  b  = eas_bmiHF_harm$beta.exposure,
  se = eas_bmiHF_harm$se.exposure,
  n  = eas_bmiHF_harm$samplesize.exposure
)
eas_bmiHF_harm$r.outcome <- get_r_from_lor(
  lor     = eas_bmiHF_harm$beta.outcome,
  af      = eas_bmiHF_harm$eaf.outcome,
  ncase   = eas_bmiHF_harm$ncase.outcome,
  ncontrol = eas_bmiHF_harm$ncontrol.outcome,
  prevalence = 0.02
)
eas_bmiHF_harm$prevalence.outcome <- 0.02
eas_bmiHF_direction <- directionality_test(eas_bmiHF_harm[eas_bmiHF_harm$mr_keep == TRUE, ])
eas_bmiHF_direction

eas_HF_cause = hum2cause('GWASsummary_CHF_Japanese_SakaueKanai2020.auto.txt.gz')
eas_bmiHF_merge = CAUSE_merge(eas_bmi_cause, eas_HF_cause)
eas_bmiHF_params = CAUSE_params(eas_bmiHF_merge)
eas_bmiHF_clump = X_clump(eas_bmiHF_merge, 'EAS')

eas_bmiHF_top_vars <- eas_bmiHF_clump$snp
set.seed(123)
eas_bmiHF_cause_res <- cause(X=eas_bmiHF_clump, variants = eas_bmiHF_top_vars, param_ests = eas_bmiHF_params)

summary(eas_bmiHF_cause_res, ci_size = 0.95)
class(eas_bmiHF_cause_res)
eas_bmiHF_cause_res$elpd # Zscore



# BMI -> T2D --------------------------------------------------------------

setwd('~/MR_data/bbj/hum0197.v3.BBJ.T2D.v1/hum0197.v3.BBJ.T2D.v1')
eas_T2D_outcome_data = hum2mr_outcome('GWASsummary_T2D_Japanese_SakaueKanai2020.auto.txt.gz',
                                     snps       = eas_bmi_exposure_data$SNP,
                                     trait_name = "T2D_Japanese")
eas_bmiT2D_harm = harmonise_data(eas_bmi_exposure_data, eas_T2D_outcome_data)
eas_bmiT2D_fstat = Fstats(eas_bmiT2D_harm[eas_bmiT2D_harm$mr_keep == T,])
eas_bmiT2D_res = mr(eas_bmiT2D_harm)
eas_bmiT2D_res
eas_bmiT2D_egger = mr_pleiotropy_test(eas_bmiT2D_harm)
eas_bmiT2D_egger

library(MRPRESSO)
set.seed(123)
eas_bmiT2D_presso = mr_presso(BetaOutcome   = "beta.outcome",
                             BetaExposure  = "beta.exposure",
                             SdOutcome     = "se.outcome",
                             SdExposure    = "se.exposure",
                             OUTLIERtest   = TRUE,         
                             DISTORTIONtest = TRUE,       
                             data           = eas_bmiT2D_harm[eas_bmiT2D_harm$mr_keep == TRUE, ],
                             NbDistribution = 1000,        
                             SignifThreshold = 0.05)
eas_bmiT2D_presso$`Main MR results`          # global + outlier-corrected estimates
eas_bmiT2D_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eas_bmiT2D_het = mr_heterogeneity(eas_bmiT2D_harm)
eas_bmiT2D_het
eas_bmiT2D_harm <- add_metadata(eas_bmiT2D_harm)
eas_bmiT2D_harm$samplesize.exposure = 163836 # add losing data
eas_bmiT2D_harm$ncase.outcome = 45383
eas_bmiT2D_harm$ncontrol.outcome = 132032
eas_bmiT2D_harm$r.exposure <- get_r_from_bsen(
  b  = eas_bmiT2D_harm$beta.exposure,
  se = eas_bmiT2D_harm$se.exposure,
  n  = eas_bmiT2D_harm$samplesize.exposure
)
eas_bmiT2D_harm$r.outcome <- get_r_from_lor(
  lor     = eas_bmiT2D_harm$beta.outcome,
  af      = eas_bmiT2D_harm$eaf.outcome,
  ncase   = eas_bmiT2D_harm$ncase.outcome,
  ncontrol = eas_bmiT2D_harm$ncontrol.outcome,
  prevalence = 0.15
)
eas_bmiT2D_harm$prevalence.outcome <- 0.15
eas_bmiT2D_direction <- directionality_test(eas_bmiT2D_harm[eas_bmiT2D_harm$mr_keep == TRUE, ])
eas_bmiT2D_direction

eas_T2D_cause = hum2cause('GWASsummary_T2D_Japanese_SakaueKanai2020.auto.txt.gz')
eas_bmiT2D_merge = CAUSE_merge(eas_bmi_cause, eas_T2D_cause)
eas_bmiT2D_params = CAUSE_params(eas_bmiT2D_merge)
eas_bmiT2D_clump = X_clump(eas_bmiT2D_merge, 'EAS')

eas_bmiT2D_top_vars <- eas_bmiT2D_clump$snp
set.seed(123)
eas_bmiT2D_cause_res <- cause(X=eas_bmiT2D_clump, variants = eas_bmiT2D_top_vars, param_ests = eas_bmiT2D_params)

summary(eas_bmiT2D_cause_res, ci_size = 0.95)
class(eas_bmiT2D_cause_res)
eas_bmiT2D_cause_res$elpd # Zscore



# BMI -> eGFR -------------------------------------------------------------

setwd('~/MR_data/bbj')
eas_GFR_outcome_data = sakaue2mr_outcome('hum0197.v5.gwas.eGFR.v1.txt.gz',
                                         snps       = eas_bmi_exposure_data$SNP,
                                         trait_name = 'GFR_Japanese')
eas_bmiGFR_harm = harmonise_data(eas_bmi_exposure_data, eas_GFR_outcome_data)
eas_bmiGFR_fstat = Fstats(eas_bmiGFR_harm[eas_bmiGFR_harm$mr_keep == T,])
eas_bmiGFR_res = mr(eas_bmiGFR_harm)
eas_bmiGFR_res
eas_bmiGFR_egger = mr_pleiotropy_test(eas_bmiGFR_harm)
eas_bmiGFR_egger

library(MRPRESSO)
set.seed(123)
eas_bmiGFR_presso = mr_presso(BetaOutcome   = "beta.outcome",
                              BetaExposure  = "beta.exposure",
                              SdOutcome     = "se.outcome",
                              SdExposure    = "se.exposure",
                              OUTLIERtest   = TRUE,         
                              DISTORTIONtest = TRUE,       
                              data           = eas_bmiGFR_harm[eas_bmiGFR_harm$mr_keep == TRUE, ],
                              NbDistribution = 1000,        
                              SignifThreshold = 0.05)
eas_bmiGFR_presso$`Main MR results`          # global + outlier-corrected estimates
eas_bmiGFR_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eas_bmiGFR_het = mr_heterogeneity(eas_bmiGFR_harm)
eas_bmiGFR_het
eas_bmiGFR_harm <- add_metadata(eas_bmiGFR_harm)
eas_bmiGFR_harm$samplesize.exposure = 163836 # add losing data
eas_bmiGFR_harm$samplesize.outcome = 1154633
eas_bmiGFR_harm$r.exposure <- get_r_from_bsen(
  b  = eas_bmiGFR_harm$beta.exposure,
  se = eas_bmiGFR_harm$se.exposure,
  n  = eas_bmiGFR_harm$samplesize.exposure
)
eas_bmiGFR_harm$r.outcome <- get_r_from_bsen(
  b  = eas_bmiGFR_harm$beta.outcome,
  se = eas_bmiGFR_harm$se.outcome,
  n  = eas_bmiGFR_harm$samplesize.outcome
)
eas_bmiGFR_direction <- directionality_test(eas_bmiGFR_harm[eas_bmiGFR_harm$mr_keep == TRUE, ])
eas_bmiGFR_direction

eas_GFR_cause = sakaue2cause('hum0197.v5.gwas.eGFR.v1.txt.gz')
eas_bmiGFR_merge = CAUSE_merge(eas_bmi_cause, eas_GFR_cause)
eas_bmiGFR_params = CAUSE_params(eas_bmiGFR_merge)
eas_bmiGFR_clump = X_clump(eas_bmiGFR_merge, 'EAS')

eas_bmiGFR_top_vars <- eas_bmiGFR_clump$snp
set.seed(123)
eas_bmiGFR_cause_res <- cause(X=eas_bmiGFR_clump, variants = eas_bmiGFR_top_vars, param_ests = eas_bmiGFR_params)

summary(eas_bmiGFR_cause_res, ci_size = 0.95)
class(eas_bmiGFR_cause_res)
eas_bmiGFR_cause_res$elpd # Zscore



# BMI -> AF ---------------------------------------------------------------

setwd('~/MR_data/bbj/hum0197.v3.BBJ.AF.v1/hum0197.v3.BBJ.AF.v1')
eas_AF_outcome_data = hum2mr_outcome('GWASsummary_Atrial_Flutter_Japanese_SakaueKanai2020.auto.txt.gz',
                                      snps = eas_bmi_exposure_data$SNP,
                                      trait_name = 'AF Japanese')
eas_bmiAF_harm = harmonise_data(eas_bmi_exposure_data, eas_AF_outcome_data)
eas_bmiAF_fstat = Fstats(eas_bmiAF_harm[eas_bmiAF_harm$mr_keep == T,])
eas_bmiAF_res = mr(eas_bmiAF_harm)
eas_bmiAF_res
eas_bmiAF_egger = mr_pleiotropy_test(eas_bmiAF_harm)
eas_bmiAF_egger

library(MRPRESSO)
set.seed(123)
eas_bmiAF_presso = mr_presso(BetaOutcome   = "beta.outcome",
                             BetaExposure  = "beta.exposure",
                             SdOutcome     = "se.outcome",
                             SdExposure    = "se.exposure",
                             OUTLIERtest   = TRUE,         
                             DISTORTIONtest = TRUE,       
                             data           = eas_bmiAF_harm[eas_bmiAF_harm$mr_keep == TRUE, ],
                             NbDistribution = 1000,        
                             SignifThreshold = 0.05)
eas_bmiAF_presso$`Main MR results`          # global + outlier-corrected estimates
eas_bmiAF_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eas_bmiAF_het = mr_heterogeneity(eas_bmiAF_harm)
eas_bmiAF_het
eas_bmiAF_harm <- add_metadata(eas_bmiAF_harm)
eas_bmiAF_harm$samplesize.exposure = 163836 # add losing data
eas_bmiAF_harm$ncase.outcome = 45383
eas_bmiAF_harm$ncontrol.outcome = 132032
eas_bmiAF_harm$r.exposure <- get_r_from_bsen(
  b  = eas_bmiAF_harm$beta.exposure,
  se = eas_bmiAF_harm$se.exposure,
  n  = eas_bmiAF_harm$samplesize.exposure
)
eas_bmiAF_harm$r.outcome <- get_r_from_lor(
  lor     = eas_bmiAF_harm$beta.outcome,
  af      = eas_bmiAF_harm$eaf.outcome,
  ncase   = eas_bmiAF_harm$ncase.outcome,
  ncontrol = eas_bmiAF_harm$ncontrol.outcome,
  prevalence = 0.15
)
eas_bmiAF_harm$prevalence.outcome <- 0.15
eas_bmiAF_direction <- directionality_test(eas_bmiAF_harm[eas_bmiAF_harm$mr_keep == TRUE, ])
eas_bmiAF_direction

eas_AF_cause = hum2cause('GWASsummary_Atrial_Flutter_Japanese_SakaueKanai2020.auto.txt.gz')
eas_bmiAF_merge = CAUSE_merge(eas_bmi_cause, eas_AF_cause)
eas_bmiAF_params = CAUSE_params(eas_bmiAF_merge)
eas_bmiAF_clump = X_clump(eas_bmiAF_merge, 'EAS')

eas_bmiAF_top_vars <- eas_bmiAF_clump$snp
set.seed(123)
eas_bmiAF_cause_res <- cause(X=eas_bmiAF_clump, variants = eas_bmiAF_top_vars, param_ests = eas_bmiAF_params)

summary(eas_bmiAF_cause_res, ci_size = 0.95)
class(eas_bmiAF_cause_res)
eas_bmiAF_cause_res$elpd # Zscore



# BMI -> MASLDa -----------------------------------------------------------

setwd('~/MR_data/bbj/hum0197.v3.BBJ.BMI.v1/hum0197.v3.BBJ.BMI.v1')
eas_bmi_exposure_data = sakaue2mr('GWASsummary_BMI_Japanese_SakaueKanai2020.auto.txt.gz',
                                  clump = T,
                                  population = 'EAS')

setwd('~/MR_data/bbj/hum0197.v3.BBJ.ALT.v1/hum0197.v3.BBJ.ALT.v1')
eas_MASLDa_outcome_data = sakaue2mr_outcome('GWASsummary_ALT_Japanese_SakaueKanai2020.auto.txt.gz',
                                      snps       = eas_bmi_exposure_data$SNP,
                                      trait_name = "ALT_Japanese")

eas_bmiMASLDa_harm = harmonise_data(eas_bmi_exposure_data,eas_MASLDa_outcome_data, action = 3)
eas_bmiMASLDa_fstat <- Fstats(eas_bmiMASLDa_harm[eas_bmiMASLDa_harm$mr_keep == TRUE, ])
eas_bmiMASLDa_res = mr(eas_bmiMASLDa_harm)
eas_bmiMASLDa_res
eas_bmiMASLDa_egger = mr_pleiotropy_test(eas_bmiMASLDa_harm) # to get egger intercept
eas_bmiMASLDa_egger

library(MRPRESSO)
set.seed(123)
eas_bmiMASLDa_presso = mr_presso(BetaOutcome   = "beta.outcome",
                                 BetaExposure  = "beta.exposure",
                                 SdOutcome     = "se.outcome",
                                 SdExposure    = "se.exposure",
                                 OUTLIERtest   = TRUE,          # outlier detection
                                 DISTORTIONtest = TRUE,         # distortion test
                                 data           = eas_bmiMASLDa_harm[eas_bmiMASLDa_harm$mr_keep == TRUE, ],
                                 NbDistribution = 1000,         # number of simulations
                                 SignifThreshold = 0.05
)
eas_bmiMASLDa_presso$`Main MR results`          # global + outlier-corrected estimates
eas_bmiMASLDa_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eas_bmiMASLDa_het = mr_heterogeneity(eas_bmiMASLDa_harm)
eas_bmiMASLDa_het

eas_bmiMASLDa_harm <- add_metadata(eas_bmiMASLDa_harm)
eas_bmiMASLDa_harm$samplesize.exposure = 163836 # add losing data
eas_bmiMASLDa_harm$r.exposure <- get_r_from_bsen(
  b  = eas_bmiMASLDa_harm$beta.exposure,
  se = eas_bmiMASLDa_harm$se.exposure,
  n  = eas_bmiMASLDa_harm$samplesize.exposure
)
eas_bmiMASLDa_harm$samplesize.outcome = 150545 # add losing data
eas_bmiMASLDa_harm$r.outcome <- get_r_from_bsen(
  b  = eas_bmiMASLDa_harm$beta.outcome,
  se = eas_bmiMASLDa_harm$se.outcome,
  n  = eas_bmiMASLDa_harm$samplesize.outcome
)
eas_bmiMASLDa_direction <- directionality_test(eas_bmiMASLDa_harm[eas_bmiMASLDa_harm$mr_keep == TRUE, ])
eas_bmiMASLDa_direction

setwd('~/MR_data/bbj/hum0197.v3.BBJ.BMI.v1/hum0197.v3.BBJ.BMI.v1')
eas_bmi_cause = sakaue2cause('GWASsummary_BMI_Japanese_SakaueKanai2020.auto.txt.gz')
setwd('~/MR_data/bbj/hum0197.v3.BBJ.ALT.v1/hum0197.v3.BBJ.ALT.v1')
eas_MASLDa_cause = sakaue2cause('GWASsummary_ALT_Japanese_SakaueKanai2020.auto.txt.gz')
eas_bmiMASLDa_merge = CAUSE_merge(eas_bmi_cause, eas_MASLDa_cause)
eas_bmiMASLDa_params = CAUSE_params(eas_bmiMASLDa_merge)
eas_bmiMASLDa_clump = X_clump(eas_bmiMASLDa_merge, 'EAS')

eas_bmiMASLDa_top_vars <- eas_bmiMASLDa_clump$snp
set.seed(123)
eas_bmiMASLDa_cause_res <- cause(X=eas_bmiMASLDa_clump, variants = eas_bmiMASLDa_top_vars, param_ests = eas_bmiMASLDa_params)

summary(eas_bmiMASLDa_cause_res, ci_size = 0.95)
class(eas_bmiMASLDa_cause_res)
eas_bmiMASLDa_cause_res$elpd # Zscore



# BMI -> MASLDb -----------------------------------------------------------

eas_MASLDb_outcome_data = load_tpmi_saige(logistic_file = '~/MR_data/twb/TWB TPMI PheWeb summary stats PLINK2/571.5.glm.logistic.gz',
                                           saige_file = '~/MR_data/twb/saige/571.5.saige.txt.gz',
                                           type = 'outcome')
eas_bmiMASLDb_harm = harmonise_data(eas_bmi_exposure_data,eas_MASLDb_outcome_data, action = 3)
eas_bmiMASLDb_fstat <- Fstats(eas_bmiMASLDb_harm[eas_bmiMASLDb_harm$mr_keep == TRUE, ])
eas_bmiMASLDb_res = mr(eas_bmiMASLDb_harm)
eas_bmiMASLDb_res
eas_bmiMASLDb_egger = mr_pleiotropy_test(eas_bmiMASLDb_harm) # to get egger intercept
eas_bmiMASLDb_egger

library(MRPRESSO)
set.seed(123)
eas_bmiMASLDb_presso = mr_presso(BetaOutcome   = "beta.outcome",
                                 BetaExposure  = "beta.exposure",
                                 SdOutcome     = "se.outcome",
                                 SdExposure    = "se.exposure",
                                 OUTLIERtest   = TRUE,          # outlier detection
                                 DISTORTIONtest = TRUE,         # distortion test
                                 data           = eas_bmiMASLDb_harm[eas_bmiMASLDb_harm$mr_keep == TRUE, ],
                                 NbDistribution = 1000,         # number of simulations
                                 SignifThreshold = 0.05
)
eas_bmiMASLDb_presso$`Main MR results`          # global + outlier-corrected estimates
eas_bmiMASLDb_presso$`MR-PRESSO results`$`Global Test`   # global pleiotropy test (i din include this in my report)

eas_bmiMASLDb_het = mr_heterogeneity(eas_bmiMASLDb_harm)
eas_bmiMASLDb_het

eas_bmiMASLDb_harm <- add_metadata(eas_bmiMASLDb_harm)
eas_bmiMASLDb_harm$samplesize.exposure = 163836 # add losing data
eas_bmiMASLDb_harm$r.exposure <- get_r_from_bsen(
  b  = eas_bmiMASLDb_harm$beta.exposure,
  se = eas_bmiMASLDb_harm$se.exposure,
  n  = eas_bmiMASLDb_harm$samplesize.exposure
)
eas_bmiMASLDb_harm$ncase.outcome = 18927 
eas_bmiMASLDb_harm$ncontrol.outcome = 289822
eas_bmiMASLDb_harm$samplesize.outcome = 289822 +18927
eas_bmiMASLDb_harm$r.outcome <- get_r_from_lor(
  lor     = eas_bmiMASLDb_harm$beta.outcome,
  af      = eas_bmiMASLDb_harm$eaf.outcome,
  ncase   = eas_bmiMASLDb_harm$ncase.outcome,
  ncontrol = eas_bmiMASLDb_harm$ncontrol.outcome,
  prevalence = 0.35
)
eas_bmiMASLDb_harm$prevalence.outcome <- 0.35 # Change in different outcome

eas_bmiMASLDb_direction <- directionality_test(eas_bmiMASLDb_harm[eas_bmiMASLDb_harm$mr_keep == TRUE, ])
eas_bmiMASLDb_direction

eas_MASLDb_cause = load_tpmi_saige_cause(logistic_file = '~/MR_data/twb/TWB TPMI PheWeb summary stats PLINK2/571.5.glm.logistic.gz',
                                         saige_file = '~/MR_data/twb/saige/571.5.saige.txt.gz')
eas_bmiMASLDb_merge = CAUSE_merge(eas_bmi_cause, eas_MASLDb_cause)
eas_bmiMASLDb_params = CAUSE_params(eas_bmiMASLDb_merge)
eas_bmiMASLDb_clump = X_clump(eas_bmiMASLDb_merge, 'EAS')

eas_bmiMASLDb_top_vars <- eas_bmiMASLDb_clump$snp
set.seed(123)
eas_bmiMASLDb_cause_res <- cause(X=eas_bmiMASLDb_clump, variants = eas_bmiMASLDb_top_vars, param_ests = eas_bmiMASLDb_params)

summary(eas_bmiMASLDb_cause_res, ci_size = 0.95)
class(eas_bmiMASLDb_cause_res)
eas_bmiMASLDb_cause_res$elpd # Zscore