library(vcfR)
library(TwoSampleMR)
library(data.table)


# EUR File ----------------------------------------------------------------

vcf2mr <- function(vcffile, 
                   trait_name = NULL,
                   min_pval   = 5e-8,
                   min_maf    = 0.01,
                   clump      = TRUE,
                   clump_r2   = 0.001,
                   clump_kb   = 10000,
                   population = "EUR",
                   verbose    = TRUE) {
  
  if (!file.exists(vcffile)) stop("File not found: ", vcffile)
  if (verbose) message("Reading VCF: ", vcffile)
  
  # --- 1. Read VCF ---
  vcf <- read.vcfR(vcffile, verbose = FALSE)
  
  # --- 2. Fix duplicate/missing IDs ---
  ids <- vcf@fix[, "ID"]
  missing_id <- is.na(ids) | ids == "."
  vcf@fix[missing_id, "ID"] <- paste0(vcf@fix[missing_id, "CHROM"], ":",
                                      vcf@fix[missing_id, "POS"])
  vcf@fix[, "ID"] <- make.unique(vcf@fix[, "ID"], sep = "_")
  if (verbose) message("Duplicates fixed: ", sum(duplicated(ids)))
  
  # --- 3. Extract FORMAT fields ---
  fix_df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  
  es <- extract.gt(vcf, element = "ES", as.numeric = TRUE)
  se <- extract.gt(vcf, element = "SE", as.numeric = TRUE)
  lp <- extract.gt(vcf, element = "LP", as.numeric = TRUE)
  
  # --- 4. Extract AF --- check FORMAT first, fall back to INFO if missing ---
  af_format_available <- "AF" %in% vcf@gt[1, ]
  
  if (af_format_available) {
    af <- extract.gt(vcf, element = "AF", as.numeric = TRUE)
    eaf_vals <- as.numeric(af[, 1])
  } else {
    if (verbose) message("AF not in FORMAT — extracting from INFO field")
    af_info <- extract.info(vcf, element = "AF", as.numeric = TRUE)
    eaf_vals <- as.numeric(af_info)
  }
  
  # --- 5. Extract SS if available ---
  ss_available <- "SS" %in% strsplit(vcf@gt[1, 1], ":")[[1]]
  if (ss_available) {
    ss <- extract.gt(vcf, element = "SS", as.numeric = TRUE)
    ss_vals <- as.numeric(ss[, 1])
  } else {
    if (verbose) message("SS not in FORMAT — samplesize set to NA")
    ss_vals <- rep(NA, nrow(fix_df))
  }
  
  # --- 6. Build data frame ---
  mr_data <- data.frame(
    SNP           = fix_df$ID,
    CHR           = fix_df$CHROM,
    POS           = as.numeric(fix_df$POS),
    effect_allele = fix_df$ALT,
    other_allele  = fix_df$REF,
    beta          = as.numeric(es[, 1]),
    se            = as.numeric(se[, 1]),
    pval          = 10^(-as.numeric(lp[, 1])),
    eaf           = eaf_vals,
    samplesize    = ss_vals,
    stringsAsFactors = FALSE
  )
  
  mr_data$phenotype <- if (!is.null(trait_name)) trait_name else
    tools::file_path_sans_ext(tools::file_path_sans_ext(basename(vcffile)))
  
  # --- 7. QC filters ---
  n_before <- nrow(mr_data)
  mr_data <- mr_data[!is.na(mr_data$beta) & !is.na(mr_data$se) &
                       !is.na(mr_data$pval), ]
  mr_data <- mr_data[mr_data$pval < min_pval, ]
  
  # MAF filter only if EAF available
  if (!all(is.na(mr_data$eaf))) {
    mr_data$maf <- ifelse(mr_data$eaf < 0.5, mr_data$eaf, 1 - mr_data$eaf)
    mr_data <- mr_data[is.na(mr_data$maf) | mr_data$maf >= min_maf, ]
  } else {
    if (verbose) message("EAF not available — skipping MAF filter")
  }
  
  if (verbose) message("Variants after QC: ", nrow(mr_data),
                       " (removed ", n_before - nrow(mr_data), ")")
  
  # --- 8. Format for TwoSampleMR ---
  exposure_dat <- format_data(
    mr_data,
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    samplesize_col    = "samplesize",
    chr_col           = "CHR",
    pos_col           = "POS",
    phenotype_col     = "phenotype"
  )
  
  # --- 9. Clumping ---
  if (clump) {
    if (verbose) message("Clumping with ", population, " panel...")
    exposure_dat <- clump_data(exposure_dat, clump_r2 = clump_r2,
                               clump_kb = clump_kb, pop = population)
    if (verbose) message("Variants after clumping: ", nrow(exposure_dat))
  }
  
  return(exposure_dat)
}


gz2mr <- function(filepath, exposure_name, pop,
                           pval_threshold = 5e-8,
                           clump_kb = 10000,
                           clump_r2 = 0.001) {
  
  dt <- fread(filepath, select = c("rsID", "REF", "ALT", "POOLED_ALT_AF",
                                   "EFFECT_SIZE", "SE", "pvalue", "N"))
  
  setnames(dt,
           old = c("rsID", "REF", "ALT", "POOLED_ALT_AF", "EFFECT_SIZE", "SE", "pvalue", "N"),
           new = c("SNP",  "other_allele", "effect_allele", "eaf",
                   "beta", "se", "pval", "samplesize")
  )
  
  # pvalue is character in this file — convert
  dt[, pval := as.numeric(pval)]
  
  # Step 1: filter to GWS variants
  dt <- dt[pval < pval_threshold]
  message(sprintf("%d variants passing P < %s before clumping", nrow(dt), pval_threshold))
  
  out <- format_data(
    as.data.frame(dt),
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    samplesize_col    = "samplesize"
  )
  
  out$exposure <- exposure_name
  
  # Step 1: LD clumping — r2 < 0.001, 10,000 kb, ancestry-matched panel
  out <- clump_data(
    out,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p1 = pval_threshold,
    pop      = pop        # "EUR" or "EAS" passed by caller
  )
  message(sprintf("%d instruments after clumping (r2 < %s, %s kb, %s panel)",
                  nrow(out), clump_r2, clump_kb, pop))
  
  return(out)
}



finngen2mr <- function(filepath,
                       exposure_name  = "CHD_FinnGen_R12",
                       maf_threshold  = 0.01,
                       pval_threshold = 5e-8,
                       clump          = TRUE,
                       clump_kb       = 10000,
                       clump_r2       = 0.001,
                       pop            = "EUR")
{
  dt <- fread(filepath, sep = "\t", showProgress = TRUE)
  
  setnames(dt,
           old = c("#chrom", "pos",  "ref",           "alt",            "rsids", "af_alt", "beta", "sebeta", "pval"),
           new = c("CHR",    "BP",   "other_allele",  "effect_allele",  "SNP",   "eaf",    "beta", "se",     "pval")
  )
  
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(SNP)]
  dt <- dt[se > 0]
  dt <- dt[maf >= maf_threshold]
  dt <- dt[!grepl("^\\.", SNP)]   # remove variants without rsID
  
  dt <- dt[pval < pval_threshold]
  message(sprintf("%d variants passing P < %s before clumping", nrow(dt), pval_threshold))
  
  if (nrow(dt) == 0) stop("No variants passing GWS threshold — check file or threshold.")
  
  dt$exposure <- exposure_name
  
  out <- format_data(
    as.data.frame(dt),
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    chr_col           = "CHR",
    pos_col           = "BP",
    phenotype_col     = "exposure"
  )
  
  # --- Step 1: LD clumping — r2 < 0.001, 10,000 kb, ancestry-matched panel ---
  if (clump) {
    out <- clump_data(
      out,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p1 = pval_threshold,
      pop      = pop
    )
    message(sprintf("%d instruments after clumping (r2 < %s, %s kb, %s panel)",
                    nrow(out), clump_r2, clump_kb, pop))
  }
  
  return(out)
}

meta2mr <- function(filepath,
                    exposure_name,
                    pop,
                    pval_threshold = 5e-8,
                    clump          = TRUE,
                    clump_kb       = 10000,
                    clump_r2       = 0.001,
                    maf_threshold  = 0.01) {
  
  
  library(data.table)
  library(TwoSampleMR)
  
  dt <- fread(filepath, select = c("rsid", "REF", "ALT",
                                   "all_inv_var_meta_beta",
                                   "all_inv_var_meta_sebeta",
                                   "all_inv_var_meta_p",
                                   "FINNGEN_af_alt",
                                   "UKBB_af_alt",
                                   "all_meta_N",
                                   "#CHR", "POS"))
  
  setnames(dt,
           old = c("rsid", "REF",         "ALT",            
                   "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "all_inv_var_meta_p",
                   "#CHR", "POS"),
           new = c("SNP",  "other_allele", "effect_allele",  
                   "beta", "se",           "pval",
                   "CHR",  "BP")
  )
  
  dt[, eaf := ifelse(!is.na(FINNGEN_af_alt), FINNGEN_af_alt, UKBB_af_alt)]
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval)]
  dt <- dt[!is.na(SNP) & SNP != ""]
  dt <- dt[se > 0]
  dt <- dt[maf >= maf_threshold]
  
  dt <- dt[pval < pval_threshold]
  message(sprintf("%d variants passing P < %s before clumping", nrow(dt), pval_threshold))
  
  if (nrow(dt) == 0) stop("No variants passing GWS threshold — check file or threshold.")
  
  dt$exposure <- exposure_name
  
  out <- format_data(
    as.data.frame(dt),
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    samplesize_col    = "all_meta_N",
    chr_col           = "CHR",
    pos_col           = "BP",
    phenotype_col     = "exposure"
  )
  
  # --- Step 1: LD clumping — r2 < 0.001, 10,000 kb, ancestry-matched panel ---
  if (clump) {
    out <- clump_data(
      out,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p1 = pval_threshold,
      pop      = pop
    )
    message(sprintf("%d instruments after clumping (r2 < %s, %s kb, %s panel)",
                    nrow(out), clump_r2, clump_kb, pop))
  }
  
  return(out)
}


gwastsv2mr <- function(filepath,
                       exposure_name,
                       pop,
                       pval_threshold = 5e-8,
                       clump          = TRUE,
                       clump_kb       = 10000,
                       clump_r2       = 0.001,
                       maf_threshold  = 0.01) {
  
  library(data.table)
  library(TwoSampleMR)
  
  dt <- fread(filepath, sep = "\t", header = TRUE, showProgress = TRUE)
  
  if ("variant_id" %in% colnames(dt) & !"rs_id" %in% colnames(dt)) {
    setnames(dt, "variant_id", "rs_id")
  }
  
  setnames(dt,
           old = c("rs_id", "chromosome", "base_pair_location",
                   "effect_allele", "other_allele",
                   "effect_allele_frequency", "beta", "standard_error", "p_value"),
           new = c("SNP",   "CHR",        "BP",
                   "effect_allele", "other_allele",
                   "eaf",           "beta", "se",   "pval")
  )
  
  dt[, effect_allele := toupper(effect_allele)]
  dt[, other_allele  := toupper(other_allele)]
  
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval)]
  dt <- dt[se > 0]
  dt <- dt[maf >= maf_threshold]
  dt <- dt[!is.na(SNP) & SNP != ""]
  
  dt <- dt[pval < pval_threshold]
  message(sprintf("%d variants passing P < %s before clumping", nrow(dt), pval_threshold))
  
  if (nrow(dt) == 0) stop("No variants passing GWS threshold — check file or threshold.")
  
  dt$exposure <- exposure_name
  
  out <- format_data(
    as.data.frame(dt),
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    chr_col           = "CHR",
    pos_col           = "BP",
    phenotype_col     = "exposure"
  )
  
  # --- Step 1: LD clumping — r2 < 0.001, 10,000 kb, ancestry-matched panel ---
  if (clump) {
    out <- clump_data(
      out,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p1 = pval_threshold,
      pop      = pop
    )
    message(sprintf("%d instruments after clumping (r2 < %s, %s kb, %s panel)",
                    nrow(out), clump_r2, clump_kb, pop))
  }
  
  return(out)
}



mahajan2mr <- function(filepath,
                       exposure_name  = "T2D_Mahajan",
                       pop,
                       pval_threshold = 5e-8,
                       clump          = TRUE,
                       clump_kb       = 10000,
                       clump_r2       = 0.001,
                       maf_threshold  = 0.01) {
  
  library(data.table)
  library(TwoSampleMR)
  
  dt <- fread(filepath, sep = "\t", header = TRUE, showProgress = TRUE)
  
  setnames(dt,
           old = c("rsID", "CHR", "BP", "EA",            "NEA",          "EAF", "Beta", "SE", "Pvalue"),
           new = c("SNP",  "CHR", "BP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval")
  )
  
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(SNP)]
  dt <- dt[se > 0]
  dt <- dt[maf >= maf_threshold]
  dt <- dt[!grepl("^\\.", SNP)]
  
  dt <- dt[pval < pval_threshold]
  message(sprintf("%d variants passing P < %s before clumping", nrow(dt), pval_threshold))
  
  if (nrow(dt) == 0) stop("No variants passing GWS threshold — check file or threshold.")
  
  dt$exposure <- exposure_name
  
  out <- format_data(
    as.data.frame(dt),
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    chr_col           = "CHR",
    pos_col           = "BP",
    phenotype_col     = "exposure"
  )
  
  # --- Step 1: LD clumping — r2 < 0.001, 10,000 kb, ancestry-matched panel ---
  if (clump) {
    out <- clump_data(
      out,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p1 = pval_threshold,
      pop      = pop
    )
    message(sprintf("%d instruments after clumping (r2 < %s, %s kb, %s panel)",
                    nrow(out), clump_r2, clump_kb, pop))
  }
  
  return(out)
}


gwas_outcome_data <- function(outcome_id, exposure_dat, proxies = TRUE) {
  
  snps <- exposure_dat$SNP
  message(sprintf("Querying %d SNPs against %s...", length(snps), outcome_id))
  
  out <- extract_outcome_data(
    snps     = snps,
    outcomes = outcome_id,
    proxies  = proxies,
    rsq      = 0.8
  )
  
  message(sprintf("Returned %d outcome SNPs.", nrow(out)))
  return(out)
}



finngen2mr_outcome <- function(filepath,
                               snps          = NULL,
                               maf_threshold = 0.01,
                               outcome_name  = "CHD_FinnGen_R12",
                               proxy_r2      = 0.8,
                               proxy_kb      = 5000,
                               proxy_pop     = "EUR") 
{
  
  library(data.table)
  library(TwoSampleMR)
  
  dt <- fread(filepath, sep = "\t", showProgress = TRUE)
  
  setnames(dt,
           old = c("#chrom", "pos", "ref", "alt", "rsids", "af_alt", "beta", "sebeta", "pval"),
           new = c("CHR",    "BP",  "other_allele", "effect_allele", "SNP",  "eaf",    "beta", "se", "pval")
  )
  
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(SNP)]
  dt <- dt[se > 0]
  dt <- dt[maf >= maf_threshold]
  dt <- dt[!grepl("^\\.", SNP)]
  
  if (!is.null(snps)) {
    found   <- dt[SNP %in% snps]
    missing <- setdiff(snps, dt$SNP)
    
    message(sprintf("SNPs found directly : %d / %d", nrow(found), length(snps)))
    message(sprintf("SNPs missing        : %d — searching for proxies (r2 >= %.1f)...", 
                    length(missing), proxy_r2))
  } else {
    found   <- dt
    missing <- character(0)
  }
  
  found$outcome <- outcome_name                  # add label as column first
  
  found_formatted <- format_data(
    as.data.frame(found),
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    chr_col           = "CHR",
    pos_col           = "BP",
    phenotype_col     = "outcome"                # tell format_data where the label lives
  )
  
  if (length(missing) > 0) {
    message("Querying IEU OpenGWAS for proxy SNPs...")
    
    proxy_formatted <- tryCatch({
      extract_outcome_data(
        snps            = missing,
        outcomes        = outcome_name,
        proxies         = TRUE,
        rsq             = proxy_r2,
        kb              = proxy_kb,
        maf_threshold   = maf_threshold,
        opengwas_jwt    = ieugwasr::get_opengwas_jwt()
      )
    }, error = function(e) {
      message("Proxy lookup failed: ", e$message)
      message("Returning directly found SNPs only.")
      NULL
    })
    
    if (!is.null(proxy_formatted)) {
      proxy_formatted$is_proxy <- !(proxy_formatted$SNP %in% snps)
      
      n_proxies       <- sum(proxy_formatted$is_proxy)
      n_still_missing <- length(setdiff(missing, proxy_formatted$SNP))
      
      message(sprintf("Proxies found       : %d", n_proxies))
      message(sprintf("SNPs still missing  : %d", n_still_missing))
      
      found_formatted$is_proxy <- FALSE
      outcome_dat <- rbind(found_formatted, proxy_formatted)
    } else {
      found_formatted$is_proxy <- FALSE
      outcome_dat <- found_formatted
    }
    
  } else {
    found_formatted$is_proxy <- FALSE
    outcome_dat <- found_formatted
  }
  
  message(sprintf("Total outcome SNPs returned: %d", nrow(outcome_dat)))
  return(outcome_dat)
}



mahajan2mr_outcome <- function(filepath,
                               snps          = NULL,
                               maf_threshold = 0.01,
                               outcome_name  = "T2D_Mahajan",
                               proxy_r2      = 0.8,
                               proxy_kb      = 5000,
                               proxy_pop     = "EUR") {
  
  library(data.table)
  library(TwoSampleMR)
  
  dt <- fread(filepath, sep = "\t", header = TRUE, showProgress = TRUE)
  
  setnames(dt,
           old = c("rsID", "CHR", "BP",  "EA",            "NEA",           "EAF", "Beta", "SE", "Pvalue"),
           new = c("SNP",  "CHR", "BP",  "effect_allele", "other_allele",  "eaf", "beta", "se", "pval")
  )
  
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(SNP)]
  dt <- dt[se > 0]
  dt <- dt[maf >= maf_threshold]
  dt <- dt[!grepl("^\\.", SNP)]
  
  if (!is.null(snps)) {
    found   <- dt[SNP %in% snps]
    missing <- setdiff(snps, dt$SNP)
    
    message(sprintf("SNPs found directly : %d / %d", nrow(found), length(snps)))
    message(sprintf("SNPs missing        : %d — searching for proxies (r2 >= %.1f)...",
                    length(missing), proxy_r2))
  } else {
    found   <- dt
    missing <- character(0)
  }
  
  found$outcome <- outcome_name
  
  found_formatted <- format_data(
    as.data.frame(found),
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    chr_col           = "CHR",
    pos_col           = "BP",
    phenotype_col     = "outcome"
  )
  
  if (length(missing) > 0) {
    message("Querying IEU OpenGWAS for proxy SNPs...")
    
    proxy_formatted <- tryCatch({
      extract_outcome_data(
        snps         = missing,
        outcomes     = outcome_name,
        proxies      = TRUE,
        rsq          = proxy_r2,
        kb           = proxy_kb,
        maf_threshold = maf_threshold,
        opengwas_jwt = ieugwasr::get_opengwas_jwt()
      )
    }, error = function(e) {
      message("Proxy lookup failed: ", e$message)
      message("Returning directly found SNPs only.")
      NULL
    })
    
    if (!is.null(proxy_formatted)) {
      proxy_formatted$is_proxy <- !(proxy_formatted$SNP %in% snps)
      
      n_proxies       <- sum(proxy_formatted$is_proxy)
      n_still_missing <- length(setdiff(missing, proxy_formatted$SNP))
      
      message(sprintf("Proxies found       : %d", n_proxies))
      message(sprintf("SNPs still missing  : %d", n_still_missing))
      
      found_formatted$is_proxy <- FALSE
      outcome_dat <- rbind(found_formatted, proxy_formatted)
    } else {
      found_formatted$is_proxy <- FALSE
      outcome_dat <- found_formatted
    }
    
  } else {
    found_formatted$is_proxy <- FALSE
    outcome_dat <- found_formatted
  }
  
  message(sprintf("Total outcome SNPs returned: %d", nrow(outcome_dat)))
  return(outcome_dat)
}



gwastsv2mr_outcome <- function(filepath,
                               snps          = NULL,
                               maf_threshold = 0.01,
                               outcome_name  = "eGFR_GCST90624413",
                               proxy_r2      = 0.8,
                               proxy_kb      = 5000,
                               proxy_pop     = "EUR") 
{  
  library(data.table)
  library(TwoSampleMR)
  
  dt <- fread(filepath, sep = "\t", header = TRUE, showProgress = TRUE)
  
  setnames(dt,
           old = c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value"),
           new = c("SNP",        "CHR",        "BP",                 "effect_allele", "other_allele", "eaf",                    "beta", "se",             "pval")
  )
  
  dt[, effect_allele := toupper(effect_allele)]
  dt[, other_allele  := toupper(other_allele)]
  
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval)]
  dt <- dt[se > 0]
  dt <- dt[maf >= maf_threshold]
  dt <- dt[!is.na(SNP) & SNP != ""]          # drop variants with no rsID
  
  if (!is.null(snps)) {
    found   <- dt[SNP %in% snps]
    missing <- setdiff(snps, dt$SNP)
    
    message(sprintf("SNPs found directly : %d / %d", nrow(found), length(snps)))
    message(sprintf("SNPs missing        : %d — searching for proxies (r2 >= %.1f)...",
                    length(missing), proxy_r2))
  } else {
    found   <- dt
    missing <- character(0)
  }
  
  found$outcome <- outcome_name
  
  found_formatted <- format_data(
    as.data.frame(found),
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    chr_col           = "CHR",
    pos_col           = "BP",
    samplesize_col    = "n",               # per-SNP sample size available
    phenotype_col     = "outcome"
  )
  
  if (length(missing) > 0) {
    message("Querying IEU OpenGWAS for proxy SNPs...")
    
    proxy_formatted <- tryCatch({
      extract_outcome_data(
        snps          = missing,
        outcomes      = outcome_name,
        proxies       = TRUE,
        rsq           = proxy_r2,
        kb            = proxy_kb,
        maf_threshold = maf_threshold,
        opengwas_jwt  = ieugwasr::get_opengwas_jwt()
      )
    }, error = function(e) {
      message("Proxy lookup failed: ", e$message)
      message("Returning directly found SNPs only.")
      NULL
    })
    
    if (!is.null(proxy_formatted)) {
      proxy_formatted$is_proxy <- !(proxy_formatted$SNP %in% snps)
      
      n_proxies       <- sum(proxy_formatted$is_proxy)
      n_still_missing <- length(setdiff(missing, proxy_formatted$SNP))
      
      message(sprintf("Proxies found       : %d", n_proxies))
      message(sprintf("SNPs still missing  : %d", n_still_missing))
      
      found_formatted$is_proxy <- FALSE
      outcome_dat <- rbind(found_formatted, proxy_formatted)
    } else {
      found_formatted$is_proxy <- FALSE
      outcome_dat <- found_formatted
    }
    
  } else {
    found_formatted$is_proxy <- FALSE
    outcome_dat <- found_formatted
  }
  
  message(sprintf("Total outcome SNPs returned: %d", nrow(outcome_dat)))
  return(outcome_dat)
}



meta2mr_outcome <- function(filepath,
                            outcome_name,
                            snps,
                            pop,
                            proxy_r2      = 0.8,
                            maf_threshold = 0.01) {
  
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
  
  dt <- fread(filepath, select = c("rsid", "REF", "ALT",
                                   "all_inv_var_meta_beta",
                                   "all_inv_var_meta_sebeta",
                                   "all_inv_var_meta_p",
                                   "FINNGEN_af_alt",
                                   "UKBB_af_alt",
                                   "all_meta_N",
                                   "#CHR", "POS"))
  
  setnames(dt,
           old = c("rsid", "REF",         "ALT",            
                   "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "all_inv_var_meta_p",
                   "#CHR", "POS"),
           new = c("SNP",  "other_allele", "effect_allele",  
                   "beta", "se",           "pval",
                   "CHR",  "BP")
  )
  
  dt[, eaf := ifelse(!is.na(FINNGEN_af_alt), FINNGEN_af_alt, UKBB_af_alt)]
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval)]
  dt <- dt[!is.na(SNP) & SNP != ""]
  dt <- dt[se > 0]
  dt <- dt[maf >= maf_threshold]
  
  # --- Step 1: direct match of exposure instruments in outcome file ---
  present    <- intersect(snps, dt$SNP)
  missing    <- setdiff(snps, dt$SNP)
  message(sprintf("%d / %d instruments found directly in outcome file",
                  length(present), length(snps)))
  
  dt_matched <- dt[SNP %in% present]
  
  # --- Step 2: proxy search for missing SNPs (r2 >= proxy_r2) ---
  if (length(missing) > 0) {
    message(sprintf("Searching proxies for %d missing SNP(s) (r2 >= %s, %s panel)...",
                    length(missing), proxy_r2, pop))
    
    proxy_hits <- tryCatch(
      ieugwasr::ld_proxies(
        rsid    = missing,
        pop     = pop,
        r2      = proxy_r2,
        searchspace = dt$SNP  # restrict candidate proxies to SNPs present in outcome file
      ),
      error = function(e) {
        message("Proxy lookup failed: ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if (!is.null(proxy_hits) && nrow(proxy_hits) > 0) {
      # keep best proxy (highest r2) per missing SNP
      setDT(proxy_hits)
      proxy_hits <- proxy_hits[order(-R2)]
      proxy_hits <- proxy_hits[!duplicated(rsid)]
      
      dt_proxy <- dt[SNP %in% proxy_hits$proxy_snp]
      # map proxy SNP back to the original (target) exposure SNP name
      map <- setNames(proxy_hits$rsid, proxy_hits$proxy_snp)
      dt_proxy[, target_SNP := map[SNP]]
      
      message(sprintf("%d / %d missing instruments recovered via proxy",
                      nrow(dt_proxy), length(missing)))
      
      still_missing <- setdiff(missing, proxy_hits$rsid)
      if (length(still_missing) > 0) {
        message(sprintf("%d instrument(s) with no direct match or proxy: %s",
                        length(still_missing), paste(still_missing, collapse = ", ")))
      }
      
      # NOTE: allele harmonisation for proxy SNPs (effect/other allele phase
      # relative to the target SNP) is NOT resolved here — ld_proxies() output
      # includes proxy allele correlation fields (e.g. proxy_a1/proxy_a2 vs
      # target_a1/target_a2) that should be checked before trusting sign/beta
      # direction on proxied SNPs. harmonise_data() with action=2/3 will catch
      # simple strand issues but not r2<1 allele-phase flips — verify manually
      # for your proxy set if precision on a small number of SNPs matters.
      
      dt_proxy[, SNP := target_SNP]
      dt_proxy[, target_SNP := NULL]
      
      dt_matched <- rbind(dt_matched, dt_proxy, fill = TRUE)
    } else {
      message("No proxies found for missing instruments.")
    }
  }
  
  if (nrow(dt_matched) == 0) stop("No exposure instruments (direct or proxy) found in outcome file.")
  
  dt_matched$outcome <- outcome_name
  
  out <- format_data(
    as.data.frame(dt_matched),
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    samplesize_col    = "all_meta_N",
    chr_col           = "CHR",
    pos_col           = "BP",
    phenotype_col     = "outcome"
  )
  
  return(out)
}



##### CAUSE #####

library(cause)
library(dplyr)

read_vcf_for_cause <- function(vcffile) {
  
  vcf <- read.vcfR(vcffile, verbose = FALSE)
  
  # Fix duplicate/missing IDs
  ids <- vcf@fix[, "ID"]
  missing_id <- is.na(ids) | ids == "."
  vcf@fix[missing_id, "ID"] <- paste0(vcf@fix[missing_id, "CHROM"], ":",
                                      vcf@fix[missing_id, "POS"])
  vcf@fix[, "ID"] <- make.unique(vcf@fix[, "ID"], sep = "_")
  
  # Extract fields manually
  fix_df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  es <- extract.gt(vcf, element = "ES", as.numeric = TRUE)
  se <- extract.gt(vcf, element = "SE", as.numeric = TRUE)
  lp <- extract.gt(vcf, element = "LP", as.numeric = TRUE)
  af <- extract.gt(vcf, element = "AF", as.numeric = TRUE)
  
  # Build data frame in CAUSE-compatible format
  dat <- data.frame(
    ID  = fix_df$ID,
    ALT = fix_df$ALT,
    REF = fix_df$REF,
    ES  = as.numeric(es[, 1]),
    SE  = as.numeric(se[, 1]),
    p   = 10^(-as.numeric(lp[, 1])),
    AF  = as.numeric(af[, 1]),
    stringsAsFactors = FALSE
  )
  
  return(dat)
}



read_finngen_for_cause <- function(filepath) {
  
  library(data.table)
  
  # --- Read -------------------------------------------------------------------
  dt <- fread(filepath, sep = "\t", showProgress = TRUE)
  
  # --- Rename -----------------------------------------------------------------
  setnames(dt,
           old = c("#chrom", "pos", "ref", "alt", "rsids", "af_alt", "beta", "sebeta", "pval"),
           new = c("CHR",    "BP",  "REF", "ALT", "ID",    "AF",     "ES",   "SE",     "p")
  )
  
  # --- Fix duplicate/missing IDs ----------------------------------------------
  missing_id <- is.na(dt$ID) | dt$ID == "."
  dt$ID[missing_id] <- paste0(dt$CHR[missing_id], ":", dt$BP[missing_id])
  dt$ID <- make.unique(dt$ID, sep = "_")
  
  # --- Minimal QC (same implicit standards as VCF function) -------------------
  dt <- dt[!is.na(ES) & !is.na(SE) & !is.na(p)]
  dt <- dt[SE > 0]
  
  # --- Build CAUSE-compatible data frame --------------------------------------
  dat <- data.frame(
    ID  = dt$ID,
    ALT = dt$ALT,
    REF = dt$REF,
    ES  = dt$ES,
    SE  = dt$SE,
    p   = dt$p,
    AF  = dt$AF,
    stringsAsFactors = FALSE
  )
  
  message(sprintf("Variants returned: %s", format(nrow(dat), big.mark = ",")))
  
  return(dat)
}



read_mahajan_for_cause <- function(filepath) {
  
  library(data.table)
  
  # --- Read -------------------------------------------------------------------
  dt <- fread(filepath, sep = "\t", header = TRUE, showProgress = TRUE)
  
  # --- Rename to CAUSE conventions --------------------------------------------
  setnames(dt,
           old = c("rsID", "EA",  "NEA", "EAF", "Beta", "SE", "Pvalue"),
           new = c("ID",   "ALT", "REF", "AF",  "ES",   "SE", "p")
  )
  
  # --- Fix duplicate/missing IDs ----------------------------------------------
  missing_id <- is.na(dt$ID) | dt$ID == "."
  dt$ID[missing_id] <- paste0(dt$CHR[missing_id], ":", dt$BP[missing_id])
  dt$ID <- make.unique(dt$ID, sep = "_")
  
  # --- Minimal QC (matching read_vcf_for_cause and read_finngen_for_cause) ----
  dt <- dt[!is.na(ES) & !is.na(SE) & !is.na(p)]
  dt <- dt[SE > 0]
  
  # --- Build CAUSE-compatible data frame --------------------------------------
  dat <- data.frame(
    ID  = dt$ID,
    ALT = dt$ALT,
    REF = dt$REF,
    ES  = dt$ES,
    SE  = dt$SE,
    p   = dt$p,
    AF  = dt$AF,
    stringsAsFactors = FALSE
  )
  
  message(sprintf("Variants returned: %s", format(nrow(dat), big.mark = ",")))
  
  return(dat)
}



read_gwastsv_for_cause <- function(filepath) {
  
  library(data.table)
  
  # --- Read -------------------------------------------------------------------
  dt <- fread(filepath, sep = "\t", header = TRUE, showProgress = TRUE)
  
  # --- Rename to CAUSE conventions --------------------------------------------
  setnames(dt,
           old = c("rs_id", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value"),
           new = c("ID",    "ALT",           "REF",           "AF",                     "ES",   "SE",             "p")
  )
  
  # --- Uppercase alleles (GWAS Catalog uses lowercase) ------------------------
  dt[, ALT := toupper(ALT)]
  dt[, REF := toupper(REF)]
  
  # --- Fix duplicate/missing IDs ----------------------------------------------
  missing_id <- is.na(dt$ID) | dt$ID == ""
  dt$ID[missing_id] <- paste0(dt$chromosome[missing_id], ":", dt$base_pair_location[missing_id])
  dt$ID <- make.unique(dt$ID, sep = "_")
  
  # --- Minimal QC (matching all other read_*_for_cause functions) -------------
  dt <- dt[!is.na(ES) & !is.na(SE) & !is.na(p)]
  dt <- dt[SE > 0]
  
  # --- Build CAUSE-compatible data frame --------------------------------------
  dat <- data.frame(
    ID  = dt$ID,
    ALT = dt$ALT,
    REF = dt$REF,
    ES  = dt$ES,
    SE  = dt$SE,
    p   = dt$p,
    AF  = dt$AF,
    stringsAsFactors = FALSE
  )
  
  message(sprintf("Variants returned: %s", format(nrow(dat), big.mark = ",")))
  
  return(dat)
}



read_gwastsv_for_cause_2.0 <- function(filepath) {
  
  library(data.table)
  
  # --- Read -------------------------------------------------------------------
  dt <- fread(filepath, sep = "\t", header = TRUE, showProgress = TRUE,
              select = c("variant_id", "other_allele", "effect_allele", 
                         "beta", "standard_error"))
  
  # --- Flexible SNP column: handle both "rs_id" and "variant_id" -------------
  if ("variant_id" %in% colnames(dt) & !"rs_id" %in% colnames(dt)) {
    setnames(dt, "variant_id", "rs_id")
  }
  
  # --- Rename -----------------------------------------------------------------
  setnames(dt,
           old = c("rs_id", "other_allele", "effect_allele", "beta", "standard_error"),
           new = c("snp",   "A2",           "A1",            "beta_hat", "seb")
           # other_allele = reference = A2, effect_allele = A1
  )
  
  # --- QC filters -------------------------------------------------------------
  dt <- dt[!is.na(beta_hat) & !is.na(seb) & seb > 0]
  dt <- dt[!is.na(snp) & snp != ""]
  
  return(dt)
}



gz2cause <- function(filepath) {
  
  dt <- fread(filepath, select = c("rsID", "REF", "ALT", "EFFECT_SIZE", "SE"))
  
  setnames(dt,
           old = c("rsID", "REF", "ALT", "EFFECT_SIZE", "SE"),
           new = c("snp",  "A2",  "A1",  "beta_hat",    "seb")
           # REF = reference (other) = A2, ALT = effect = A1
  )
  
  return(dt)
}



read_meta_for_cause <- function(filepath) {
  
  dt <- fread(filepath, select = c("rsid", "REF", "ALT",
                                   "all_inv_var_meta_beta",
                                   "all_inv_var_meta_sebeta"))
  
  setnames(dt,
           old = c("rsid", "REF", "ALT", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta"),
           new = c("snp",  "A2",  "A1",  "beta_hat",              "seb")
           # REF = reference (other) = A2, ALT = effect = A1
  )
  
  # Remove variants without rsID
  dt <- dt[!is.na(snp) & !grepl("^\\.", snp) & snp != ""]
  
  # Remove rows with missing beta or se
  dt <- dt[!is.na(beta_hat) & !is.na(seb) & seb > 0]
  
  return(dt)
}



# EAS File ----------------------------------------------------------------

sakaue2mr <- function(filepath,
                      trait_name = "BMI_Japanese",
                      min_pval   = 5e-8,
                      min_maf    = 0.01,
                      min_info   = 0.8,
                      clump      = TRUE,
                      clump_r2   = 0.001,
                      clump_kb   = 10000,
                      population = "EAS",        # Japanese GWAS → EAS panel
                      verbose    = TRUE) {
  
  library(data.table)
  library(TwoSampleMR)
  
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  if (verbose) message("Reading: ", filepath)
  
  dt <- fread(filepath, sep = "\t", showProgress = verbose)
  
  if (verbose) message("Total variants: ", nrow(dt))
  
  # --- Rename columns ---
  dt <- dt[, .(
    SNP           = SNP,
    CHR           = CHR,
    POS           = BP,
    effect_allele = ALLELE1,      # ALLELE1 = effect allele (A1)
    other_allele  = ALLELE0,      # ALLELE0 = other allele
    eaf           = A1FREQ,       # effect allele frequency
    beta          = BETA,
    se            = SE,
    pval          = P_BOLT_LMM_INF,  # use BOLT-LMM p-value (more accurate than LINREG)
    info          = INFO
  )]
  
  # --- QC filters ---
  n_before <- nrow(dt)
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(SNP)]
  dt <- dt[se > 0]
  dt <- dt[info >= min_info]                          # imputation quality filter
  dt <- dt[pval < min_pval]                           # GWS threshold
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  dt <- dt[maf >= min_maf]
  
  if (verbose) message("Variants after QC: ", nrow(dt),
                       " (removed ", n_before - nrow(dt), ")")
  
  dt$phenotype <- trait_name
  
  # --- Format for TwoSampleMR ---
  exposure_dat <- format_data(
    as.data.frame(dt),
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    chr_col           = "CHR",
    pos_col           = "POS",
    phenotype_col     = "phenotype"
  )
  
  # --- Clumping ---
  if (clump) {
    if (verbose) message("Clumping with ", population, " panel...")
    exposure_dat <- clump_data(exposure_dat,
                               clump_r2 = clump_r2,
                               clump_kb = clump_kb,
                               pop      = population)
    if (verbose) message("Variants after clumping: ", nrow(exposure_dat))
  }
  
  return(exposure_dat)
}



sakaue2mr2.0 <- function(
    gwas_file,
    clump         = TRUE,
    pval_threshold = 5e-8,
    min_maf       = 0.001,
    imp_filter    = 0.3,
    clump_r2      = 0.001,
    clump_kb      = 10000,
    pop           = "EAS"
) {
  library(data.table)
  library(TwoSampleMR)
  
  # ── 1. Load ──────────────────────────────────────────────────────────────────
  message("Loading Sakaue BOLT-LMM file...")
  df <- fread(gwas_file)
  
  # ── 2. QC filters ────────────────────────────────────────────────────────────
  df <- df[INFO >= imp_filter]
  df <- df[A1FREQ >= min_maf & A1FREQ <= (1 - min_maf)]
  
  # P_BOLT_LMM_INF is stored as character ("5.6E-01") — convert
  df[, pval := as.numeric(P_BOLT_LMM_INF)]
  df <- df[!is.na(pval)]
  
  # ── 3. Filter to genome-wide significant SNPs ─────────────────────────────
  df <- df[pval < pval_threshold]
  message(sprintf("%d SNPs pass p < %g", nrow(df), pval_threshold))
  
  if (nrow(df) == 0) stop("No SNPs pass the p-value threshold.")
  
  # ── 4. Format for TwoSampleMR ─────────────────────────────────────────────
  # BOLT-LMM: ALLELE1 = effect allele, ALLELE0 = other allele
  exposure_dat <- format_data(
    dat                  = as.data.frame(df),
    type                 = "exposure",
    snp_col              = "SNP",
    beta_col             = "BETA",
    se_col               = "SE",
    eaf_col              = "A1FREQ",
    effect_allele_col    = "ALLELE1",
    other_allele_col     = "ALLELE0",
    pval_col             = "pval",
    chr_col              = "CHR",
    pos_col              = "BP",       # hg19
    min_pval             = 1e-300
  )
  
  exposure_dat$exposure <- tools::file_path_sans_ext(
    tools::file_path_sans_ext(basename(gwas_file))  # strips .txt.gz
  )
  
  # ── 5. Clumping ───────────────────────────────────────────────────────────
  if (clump) {
    message("Clumping...")
    exposure_dat <- clump_data(
      exposure_dat,
      clump_r2  = clump_r2,
      clump_kb  = clump_kb,
      clump_p1  = pval_threshold,
      pop       = pop
    )
    message(sprintf("%d SNPs retained after clumping", nrow(exposure_dat)))
  }
  
  return(exposure_dat)
}



hum2mr <- function(filepath,
                   trait_name    = "CAD_Japanese",
                   min_maf       = 0.01,
                   min_info      = 0.8,
                   clump         = TRUE,
                   pop           = "EAS",          # 1000G ancestry for LD clumping
                   clump_r2      = 0.001,
                   clump_kb      = 10000,
                   pval_threshold = 5e-8,
                   verbose       = TRUE) {
  
  library(data.table)
  library(TwoSampleMR)
  
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  if (verbose) message("Reading: ", filepath)
  
  dt <- fread(filepath, sep = "\t", showProgress = verbose)
  
  if (verbose) message("Total variants: ", nrow(dt))
  
  # --- Rename columns (mirrors hum2mr_outcome) ---
  dt <- dt[, .(
    SNP           = SNPID,
    CHR           = CHR,
    POS           = POS,
    effect_allele = Allele2,       # Allele2 = effect allele (BETA refers to this)
    other_allele  = Allele1,       # Allele1 = reference allele
    eaf           = AF_Allele2,
    beta          = BETA,
    se            = SE,
    pval          = p.value,
    info          = imputationInfo,
    samplesize    = N
  )]
  
  # --- QC filters ---
  n_before <- nrow(dt)
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(SNP)]
  dt <- dt[se > 0]
  dt <- dt[info >= min_info]
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  dt <- dt[maf >= min_maf]
  
  valid_alleles <- c("A", "T", "C", "G")
  dt <- dt[effect_allele %in% valid_alleles & other_allele %in% valid_alleles]
  
  if (verbose) message("Variants after QC: ", nrow(dt),
                       " (removed ", n_before - nrow(dt), ")")
  
  # --- Genome-wide significance filter ---
  dt <- dt[pval < pval_threshold]
  if (verbose) message("Variants at P < ", pval_threshold, ": ", nrow(dt))
  
  if (nrow(dt) == 0) {
    warning("No variants survive the p-value threshold.")
    return(NULL)
  }
  
  dt$phenotype <- trait_name
  
  # --- Format as EXPOSURE for TwoSampleMR ---
  exposure_dat <- format_data(
    as.data.frame(dt),
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    samplesize_col    = "samplesize",
    chr_col           = "CHR",
    pos_col           = "POS",
    phenotype_col     = "phenotype"
  )
  
  if (verbose) message("Variants before clumping: ", nrow(exposure_dat))
  
  # --- LD clumping (protocol: r² < 0.001, 10,000 kb, ancestry-matched) ---
  if (clump) {
    exposure_dat <- clump_data(
      exposure_dat,
      clump_r2  = clump_r2,
      clump_kb  = clump_kb,
      pop       = pop
    )
    if (verbose) message("Instruments after clumping: ", nrow(exposure_dat))
  }
  
  return(exposure_dat)
}



hum2mr_outcome <- function(filepath,
                           snps          = NULL,
                           trait_name    = "CAD_Japanese",
                           min_maf       = 0.01,
                           min_info      = 0.8,
                           verbose       = TRUE) {
  
  library(data.table)
  library(TwoSampleMR)
  
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  if (verbose) message("Reading: ", filepath)
  
  dt <- fread(filepath, sep = "\t", showProgress = verbose)
  
  if (verbose) message("Total variants: ", nrow(dt))
  
  # --- Rename columns ---
  dt <- dt[, .(
    SNP           = SNPID,
    CHR           = CHR,
    POS           = POS,
    effect_allele = Allele2,        # Allele2 = effect allele (BETA refers to this)
    other_allele  = Allele1,        # Allele1 = reference allele
    eaf           = AF_Allele2,     # effect allele frequency
    beta          = BETA,
    se            = SE,
    pval          = p.value,
    info          = imputationInfo,
    samplesize    = N
  )]
  
  # --- QC filters ---
  n_before <- nrow(dt)
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(SNP)]
  dt <- dt[se > 0]
  dt <- dt[info >= min_info]
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  dt <- dt[maf >= min_maf]
  
  # Remove non-standard alleles
  valid_alleles <- c("A", "T", "C", "G")
  dt <- dt[effect_allele %in% valid_alleles & other_allele %in% valid_alleles]
  
  if (verbose) message("Variants after QC: ", nrow(dt),
                       " (removed ", n_before - nrow(dt), ")")
  
  # --- Filter to exposure SNPs only (no p-value threshold for outcome) ---
  if (!is.null(snps)) {
    dt_found   <- dt[SNP %in% snps]
    missing    <- setdiff(snps, dt$SNP)
    
    message(sprintf("SNPs found directly : %d / %d", nrow(dt_found), length(snps)))
    message(sprintf("SNPs missing        : %d", length(missing)))
    
    dt <- dt_found
  }
  
  dt$phenotype <- trait_name
  
  # --- Format as OUTCOME for TwoSampleMR ---
  outcome_dat <- format_data(
    as.data.frame(dt),
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    samplesize_col    = "samplesize",
    chr_col           = "CHR",
    pos_col           = "POS",
    phenotype_col     = "phenotype"
  )
  
  if (verbose) message("Outcome SNPs returned: ", nrow(outcome_dat))
  return(outcome_dat)
}



sakaue2mr_outcome <- function(filepath,
                              snps       = NULL,
                              trait_name = "IS_Japanese",
                              min_maf    = 0.01,
                              min_info   = 0.8,
                              verbose    = TRUE) {
  
  library(data.table)
  library(TwoSampleMR)
  
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  if (verbose) message("Reading: ", filepath)
  
  dt <- fread(filepath, sep = "\t", showProgress = verbose)
  if (verbose) message("Total variants: ", nrow(dt))
  
  # --- Rename ---
  dt <- dt[, .(
    SNP           = SNP,
    CHR           = CHR,
    POS           = BP,
    effect_allele = ALLELE1,
    other_allele  = ALLELE0,
    eaf           = A1FREQ,
    beta          = BETA,
    se            = SE,
    pval          = P_BOLT_LMM_INF,
    info          = INFO
  )]
  
  # --- QC ---
  n_before <- nrow(dt)
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(SNP)]
  dt <- dt[se > 0]
  dt <- dt[info >= min_info]
  dt[, maf := ifelse(eaf < 0.5, eaf, 1 - eaf)]
  dt <- dt[maf >= min_maf]
  
  valid_alleles <- c("A", "T", "C", "G")
  dt <- dt[effect_allele %in% valid_alleles & other_allele %in% valid_alleles]
  
  if (verbose) message("Variants after QC: ", nrow(dt),
                       " (removed ", n_before - nrow(dt), ")")
  
  # --- Filter to exposure SNPs only ---
  if (!is.null(snps)) {
    dt_found <- dt[SNP %in% snps]
    missing  <- setdiff(snps, dt$SNP)
    message(sprintf("SNPs found directly : %d / %d", nrow(dt_found), length(snps)))
    message(sprintf("SNPs missing        : %d", length(missing)))
    dt <- dt_found
  }
  
  dt$phenotype <- trait_name
  
  # --- Format as OUTCOME ---
  outcome_dat <- format_data(
    as.data.frame(dt),
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    chr_col           = "CHR",
    pos_col           = "POS",
    phenotype_col     = "phenotype"
  )
  
  if (verbose) message("Outcome SNPs returned: ", nrow(outcome_dat))
  return(outcome_dat)
}



egfr2mr_outcome <- function(filepath, snps = NULL, proxies = TRUE) {
  
  dt <- fread(filepath, select = c("RSID", "Allele1", "Allele2", "Freq1", 
                                   "Effect", "StdErr", "P-value", "n_total_sum"))
  
  setnames(dt, old = c("RSID", "Allele1", "Allele2", "Freq1", 
                       "Effect", "StdErr", "P-value", "n_total_sum"),
           new = c("SNP",  "effect_allele", "other_allele", "eaf",
                   "beta", "se", "pval", "samplesize"))
  
  if (!is.null(snps)) dt <- dt[SNP %in% snps]
  
  dt[, outcome := "eGFR"]
  
  out <- format_data(
    as.data.frame(dt),
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    samplesize_col    = "samplesize"
  )
  
  # --- Proxy search for missing SNPs ---
  if (!is.null(snps) && proxies) {
    missing_snps <- setdiff(snps, out$SNP)
    
    if (length(missing_snps) > 0) {
      message(sprintf("Searching proxies for %d missing SNPs (r2 >= 0.8)...", 
                      length(missing_snps)))
      
      proxy_out <- extract_outcome_data(
        snps             = missing_snps,
        outcomes         = "ebi-a-GCST90104571",  # placeholder, not used for local file
        proxies          = TRUE,
        rsq              = 0.8,
        align_alleles    = 1,
        palindromes      = 1,
        maf_threshold    = 0.3
      )
      
      if (!is.null(proxy_out) && nrow(proxy_out) > 0) {
        out <- bind_rows(out, proxy_out)
        message(sprintf("Recovered %d SNPs via proxies.", nrow(proxy_out)))
      }
    } else {
      message("All SNPs found directly — no proxy search needed.")
    }
  }
  
  return(out)
}



# SAIGE_TPMI_FILE ---------------------------------------------------------

# ── Core helper: build CHR:POS → rsID map from logistic file ─────────────────
build_tpmi_rsid_map <- function(logistic_file) {
  message("Building rsID map from logistic file...")
  log_df <- fread(logistic_file, select = c("#CHROM", "POS", "rsID"))
  setnames(log_df, "#CHROM", "CHR")
  log_df[, chrpos := paste0(CHR, ":", POS)]
  
  # Drop any missing rsIDs and duplicated positions
  log_df <- log_df[!is.na(rsID) & rsID != "" & rsID != "."]
  log_df <- log_df[!duplicated(chrpos)]
  
  log_df[, .(chrpos, rsID)]
}

# ── annotate_rsids: attach rsID to SAIGE data.table by CHR:POS ───────────────
annotate_rsids <- function(saige_dt, rsid_map) {
  saige_dt[, chrpos := paste0(CHR, ":", POS)]
  saige_dt <- merge(saige_dt, rsid_map, by = "chrpos", all.x = FALSE)
  # all.x = FALSE: drops rows with no rsID match (can't use in MR anyway)
  
  n_matched <- nrow(saige_dt)
  message(sprintf("  rsID annotation: %d SNPs matched", n_matched))
  saige_dt
}

# ── Main loader ───────────────────────────────────────────────────────────────
load_tpmi_saige <- function(
    saige_file,
    logistic_file,
    type           = c("exposure", "outcome"),  # which MR role
    pval_threshold = 5e-8,                      # only used for exposure
    min_maf        = 0.001,
    imp_filter     = 0.3,
    qc_pass_only   = TRUE
) {
  type <- match.arg(type)
  library(data.table)
  library(TwoSampleMR)
  
  # 1. Load SAIGE
  message("Loading SAIGE file...")
  d <- fread(saige_file)
  message(sprintf("  %d rows loaded", nrow(d)))
  
  # 2. QC filters
  if (qc_pass_only) {
    d <- d[QC == "PASS"]
    message(sprintf("  %d rows after QC == PASS", nrow(d)))
  }
  d <- d[imputationInfo >= imp_filter]
  d <- d[AF_Allele2 >= min_maf & AF_Allele2 <= (1 - min_maf)]
  
  # 3. For exposure: pre-filter to p < threshold before rsID merge (faster)
  if (type == "exposure") {
    d <- d[p.value < pval_threshold]
    message(sprintf("  %d rows pass p < %g", nrow(d), pval_threshold))
    if (nrow(d) == 0) stop("No SNPs pass the p-value threshold.")
  }
  
  # 4. Annotate rsIDs via logistic file
  rsid_map <- build_tpmi_rsid_map(logistic_file)
  d <- annotate_rsids(d, rsid_map)
  if (nrow(d) == 0) stop("No SNPs remained after rsID annotation.")
  
  # 5. Format for TwoSampleMR
  # SAIGE: Allele2 = effect allele, Allele1 = other allele
  out <- format_data(
    dat                  = as.data.frame(d),
    type                 = type,
    snp_col              = "rsID",
    beta_col             = "BETA",
    se_col               = "SE",
    eaf_col              = "AF_Allele2",
    effect_allele_col    = "Allele2",
    other_allele_col     = "Allele1",
    pval_col             = "p.value",
    chr_col              = "CHR",
    pos_col              = "POS",
    min_pval             = 1e-300
  )
  
  out$phenotype <- tools::file_path_sans_ext(
    tools::file_path_sans_ext(basename(saige_file))
  )
  
  message(sprintf("Done. %d SNPs in %s object.", nrow(out), type))
  out
}



##### CAUSE #####

sakaue2cause <- function(filepath, verbose = TRUE) {
  
  library(data.table)
  
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  if (verbose) message("Reading: ", filepath)
  
  dt <- fread(filepath, sep = "\t", showProgress = verbose)
  
  # --- Rename to CAUSE-compatible format ---
  dat <- data.frame(
    ID  = dt$SNP,
    ALT = dt$ALLELE1,       # effect allele
    REF = dt$ALLELE0,       # other allele
    ES  = dt$BETA,
    SE  = dt$SE,
    p   = dt$P_BOLT_LMM_INF,
    AF  = dt$A1FREQ,
    stringsAsFactors = FALSE
  )
  
  # --- Basic QC ---
  dat <- dat[!is.na(dat$ES) & !is.na(dat$SE) & !is.na(dat$p), ]
  dat <- dat[dat$SE > 0, ]
  dat <- dat[dt$INFO >= 0.8, ]   # imputation quality filter
  
  # Remove non-standard alleles
  valid <- c("A", "T", "C", "G")
  dat <- dat[dat$ALT %in% valid & dat$REF %in% valid, ]
  
  if (verbose) message("Variants returned: ", nrow(dat))
  return(dat)
}



hum2cause <- function(filepath, verbose = TRUE) {
  
  library(data.table)
  
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  if (verbose) message("Reading: ", filepath)
  
  dt <- fread(filepath, sep = "\t", showProgress = verbose)
  
  # --- Rename to CAUSE-compatible format ---
  dat <- data.frame(
    ID  = dt$SNPID,
    ALT = dt$Allele2,       # effect allele
    REF = dt$Allele1,       # other allele
    ES  = dt$BETA,
    SE  = dt$SE,
    p   = dt$p.value,
    AF  = dt$AF_Allele2,
    stringsAsFactors = FALSE
  )
  
  # --- Basic QC ---
  dat <- dat[!is.na(dat$ES) & !is.na(dat$SE) & !is.na(dat$p), ]
  dat <- dat[dat$SE > 0, ]
  dat <- dat[dt$imputationInfo >= 0.8, ]
  
  # Remove non-standard alleles
  valid <- c("A", "T", "C", "G")
  dat <- dat[dat$ALT %in% valid & dat$REF %in% valid, ]
  
  if (verbose) message("Variants returned: ", nrow(dat))
  return(dat)
}



egfr2cause <- function(filepath) {
  
  dt <- fread(filepath, select = c("RSID", "Allele1", "Allele2", 
                                   "Effect", "StdErr"))
  
  setnames(dt,
           old = c("RSID", "Allele1", "Allele2", "Effect", "StdErr"),
           new = c("snp",  "A1",      "A2",      "beta_hat", "seb")
  )
  
  return(dt)
}



load_tpmi_saige_cause <- function(
    saige_file,
    logistic_file,
    min_maf    = 0.001,
    imp_filter = 0.3
) {
  library(data.table)
  library(cause)
  
  # ── 1. Load SAIGE ────────────────────────────────────────────────────────────
  message("Loading SAIGE file for CAUSE...")
  d <- fread(saige_file)
  
  # ── 2. QC filters (no QC=="PASS" filter for CAUSE — needs genome-wide SNPs) ─
  d <- d[imputationInfo >= imp_filter]
  d <- d[AF_Allele2 >= min_maf & AF_Allele2 <= (1 - min_maf)]
  message(sprintf("  %d SNPs after QC filters", nrow(d)))
  
  # ── 3. Annotate rsIDs via logistic file ──────────────────────────────────────
  rsid_map <- build_tpmi_rsid_map(logistic_file)  # reuse helper from before
  d <- annotate_rsids(d, rsid_map)
  message(sprintf("  %d SNPs after rsID annotation", nrow(d)))
  
  if (nrow(d) == 0) stop("No SNPs remained after rsID annotation.")
  
  # ── 4. Remove duplicated rsIDs (CAUSE can't handle them) ─────────────────────
  n_before <- nrow(d)
  d <- d[!duplicated(rsID)]
  message(sprintf("  %d SNPs after removing duplicates (%d removed)",
                  nrow(d), n_before - nrow(d)))
  
  # ── 5. Rename to CAUSE-required column names ─────────────────────────────────
  # CAUSE gwas_merge() expects: snp, beta, se, A1, A2, p
  # SAIGE: Allele2 = effect allele (A1 in CAUSE convention)
  out <- d[, .(
    snp      = rsID,
    beta_hat = BETA,    
    seb      = SE,     
    A1       = Allele2,
    A2       = Allele1,
    p        = p.value
  )]
  
  message(sprintf("Done. %d SNPs ready for CAUSE gwas_merge().", nrow(out)))
  return(out)
}



# Usage Function ----------------------------------------------------------

##### F stats #####

Fstats <- function(harm_dat)
{
  harm_dat$Fstat = (harm_dat$beta.exposure/harm_dat$se.exposure)^2
  
  mean_F = mean(harm_dat$Fstat, na.rm =T)
  
  message("Number of instruments: ", nrow(harm_dat))
  message("Mean F-statistic: ", round(mean_F, 2))
  
  if (mean_F > 10) {
    message("PASS: Mean F > 10, relevance assumption satisfied")
  } else {
    message("FAIL: Mean F < 10, weak instrument bias likely")
  }
  
  return(harm_dat)  # returns data with Fstat column added
}



##### CAUSE #####

CAUSE_merge <- function(exposure_dat, outcome_dat) {
  
  .standardise <- function(df) {
    if ("ID" %in% colnames(df)) {
      df <- df %>% rename(snp = ID, beta_hat = ES, seb = SE, A1 = ALT, A2 = REF)
    }
    return(df)
  }
  
  exposure_dat <- .standardise(exposure_dat)
  outcome_dat  <- .standardise(outcome_dat)
  
  X <- gwas_merge(
    exposure_dat, outcome_dat,
    snp_name_cols = c("snp",      "snp"),
    beta_hat_cols = c("beta_hat", "beta_hat"),
    se_cols       = c("seb",      "seb"),
    A1_cols       = c("A1",       "A1"),
    A2_cols       = c("A2",       "A2")
  )
  
  return(X)
}



CAUSE_params <- function(X)
{
  set.seed(123)
  varlist <- with(X, sample(snp, size = min(1000000, nrow(X)), replace = FALSE))
  params <- est_cause_params(X, varlist)
  
  return(params)
}


library(dplyr)
library(ieugwasr)
library(genetics.binaRies)

X_clump <- function(X, pop) {
  
  X_filtered <- X %>%
    mutate(p1 = 2 * pnorm(-abs(beta_hat_1 / seb1))) %>%
    filter(!is.na(p1), p1 < 5e-8) %>%
    distinct(snp, .keep_all = TRUE)
  
  cat("Variants before clumping:", nrow(X_filtered), "\n")
  
  clump_input <- X_filtered %>%
    transmute(
      rsid = snp,
      pval = p1,
      id = "exposure"
    )
  
  clumped_ids <- ieugwasr::ld_clump(
    dat = clump_input,
    clump_kb = 10000,
    clump_r2 = 0.01,
    clump_p = 5e-8,
    pop = pop
  )
  
  X_clumped <- X_filtered %>%
    semi_join(clumped_ids, by = c("snp" = "rsid"))
  
  return(X_clumped)
}
