tidy_inputGWAS <- function (GWAS, verbose = FALSE) 
{
  if (verbose) 
    cat("# Preparation of the data... \n")
  GWASnames = list(SNPID = c("rsid", "snpid", "snp", "rnpid", 
                             "rs"), CHR = c("chr"), POS = c("pos"), ALT = c("a1", 
                                                                            "alt", "alts"), REF = c("ref", "a0", "a2"), BETA = c("beta", 
                                                                                                                                 "b", "beta1", "or"), SE = c("se", "std"), Z = c("z", 
                                                                                                                                                                                 "zscore"), N = c("n", "neff"), Ncases = c("n_cases", 
                                                                                                                                                                                                                           "ncases", "n_case", "ncase"), Ncontrols = c("n_controls", 
                                                                                                                                                                                                                                                                       "ncontrols", "n_control", "ncontrol"))
  if (is.character(GWAS)) {
    if (!file.exists(GWAS)) 
      stop("the file does not exist", call. = FALSE)
    if (verbose) 
      cat(paste0("The file used as input is: \"", strsplit(GWAS, 
                                                           "/")[[1]][length(strsplit(GWAS, "/")[[1]])], 
                 "\".  \n"))
    HeaderGWAS = colnames(data.table::fread(GWAS, nrows = 1, 
                                            showProgress = FALSE, data.table = F))
  }
  else if (is.data.frame(GWAS)) {
    if (data.table::is.data.table(GWAS)) 
      GWAS = as.data.frame(GWAS)
    if (verbose) 
      cat(paste0("The data.frame used as input is: \"", 
                 attributes(GWAS)$GName, "\".  \n"))
    HeaderGWAS = colnames(GWAS)
  }
  HeaderGWAS = tolower(HeaderGWAS)
  if (all(!HeaderGWAS %in% GWASnames[["SNPID"]])) 
    stop("no SNPID column", call. = FALSE)
  if (sum(HeaderGWAS %in% GWASnames[["SNPID"]]) > 1) 
    stop("multiple SNPID columns, please provide only one", 
         call. = FALSE)
  tmp = paste0("   SNPID column, ok")
  if (all(!HeaderGWAS %in% GWASnames[["CHR"]])) 
    stop("no CHR column", call. = FALSE)
  tmp = c(tmp, "CHR column, ok")
  if (all(!HeaderGWAS %in% GWASnames[["POS"]])) 
    stop("no POS column", call. = FALSE)
  tmp = c(tmp, "POS column, ok")
  if (all(!HeaderGWAS %in% GWASnames[["ALT"]])) 
    stop("no ALT column", call. = FALSE)
  if (sum(HeaderGWAS %in% GWASnames[["ALT"]]) > 1) 
    stop("multiple ALT columns, please provide only one", 
         call. = FALSE)
  tmp = c(tmp, "ALT column, ok")
  if (all(!HeaderGWAS %in% GWASnames[["REF"]])) 
    stop("no REF column", call. = FALSE)
  if (sum(HeaderGWAS %in% GWASnames[["REF"]]) > 1) 
    stop("multiple REF columns, please provide only one", 
         call. = FALSE)
  tmp = c(tmp, "REF column, ok")
  if (all(!HeaderGWAS %in% GWASnames[["Z"]])) {
    if (!all(!HeaderGWAS %in% GWASnames[["BETA"]]) & !all(!HeaderGWAS %in% 
                                                          GWASnames[["SE"]])) {
      tmp = c(tmp, "BETA column, ok")
      tmp = c(tmp, "SE column, ok")
      getZ = TRUE
    }
    else {
      stop("no effect (BETA/SE or Z) column(s)", call. = FALSE)
    }
  }
  else if (!all(!HeaderGWAS %in% GWASnames[["Z"]])) {
    tmp = c(tmp, "Z column, ok")
    getZ = FALSE
  }
  else {
    stop("no effect (BETA/SE or Z) column(s)", call. = FALSE)
  }
  if (sum(HeaderGWAS %in% GWASnames[["BETA"]]) > 1) 
    stop("multiple BETA columns, please provide only one", 
         call. = FALSE)
  if (sum(HeaderGWAS %in% GWASnames[["SE"]]) > 1) 
    stop("multiple SE columns, please provide only one", 
         call. = FALSE)
  if (sum(HeaderGWAS %in% GWASnames[["Z"]]) > 1) 
    stop("multiple Z columns, please provide only one", 
         call. = FALSE)
  if (all(!HeaderGWAS %in% GWASnames[["N"]])) 
    stop("no N column", call. = FALSE)
  tmp = c(tmp, paste0("N column, ok \n"))
  if (sum(HeaderGWAS %in% GWASnames[["N"]]) > 1) 
    stop("multiple N columns, please provide only one", 
         call. = FALSE)
  if (verbose) 
    cat(paste(tmp, collapse = " - "))
  if (is.character(GWAS)) {
    GWASData = tibble::as_tibble(data.table::fread(GWAS, 
                                                   showProgress = FALSE, data.table = F))
    attributes(GWASData)$GName = basename(GWAS)
  }
  else if (is.data.frame(GWAS)) {
    GWASData = tibble::as_tibble(GWAS)
    rm(GWAS)
  }
  SNPID = match(HeaderGWAS, GWASnames[["SNPID"]])
  SNPID = which(!is.na(SNPID))[1]
  CHR = match(HeaderGWAS, GWASnames[["CHR"]])
  CHR = which(!is.na(CHR))[1]
  POS = match(HeaderGWAS, GWASnames[["POS"]])
  POS = which(!is.na(POS))[1]
  ALT = match(HeaderGWAS, GWASnames[["ALT"]])
  ALT = which(!is.na(ALT))[1]
  REF = match(HeaderGWAS, GWASnames[["REF"]])
  REF = which(!is.na(REF))[1]
  BETA = match(HeaderGWAS, GWASnames[["BETA"]])
  BETA = which(!is.na(BETA))[1]
  SE = match(HeaderGWAS, GWASnames[["SE"]])
  SE = which(!is.na(SE))[1]
  ZSTAT = match(HeaderGWAS, GWASnames[["Z"]])
  ZSTAT = which(!is.na(ZSTAT))[1]
  N = match(HeaderGWAS, GWASnames[["N"]])
  N = which(!is.na(N))[1]
  colNumbers = c(SNPID, CHR, POS, ALT, REF, BETA, SE, ZSTAT, 
                 N)
  colNames = c("rsid", "chr", "pos", "alt", "ref", "beta", 
               "se", "Z", "N")
  colNames = colNames[!is.na(colNumbers)]
  colNumbers = colNumbers[!is.na(colNumbers)]
  GWASData_clean <- GWASData %>% dplyr::select(tidyselect::all_of(colNumbers)) %>% 
    stats::setNames(colNames)
  if (getZ) {
    if ("or" %in% HeaderGWAS) {
      GWASData_clean <- dplyr::mutate(Z = log(.data$beta)/.data$se, 
                                      beta = NULL, se = NULL)
    }
    else {
      GWASData_clean <- GWASData_clean %>% dplyr::mutate(Z = .data$beta/.data$se, 
                                                         beta = NULL, se = NULL)
    }
  }
  else {
    GWASData_clean <- GWASData_clean %>% dplyr::mutate(beta = NULL, 
                                                       se = NULL)
  }
  GWASData_clean <- GWASData_clean %>% dplyr::mutate(std_beta = .data$Z/sqrt(.data$N), 
                                                     std_SE = 1/sqrt(.data$N), p = 2 * stats::pnorm(-abs(.data$Z)))
  res = GWASData_clean
  return(res)
}


run_LDSC <- function (exposure_data, exposure_name, outcome_data, outcome_name, 
                      ld, hm3, save_logfiles, verbose) 
{
  utils::write.table(exposure_data, paste0(exposure_name, 
                                           ".tsv"), sep = "\t", quote = F, row.names = F)
  utils::write.table(outcome_data, paste0(outcome_name, ".tsv"), 
                     sep = "\t", quote = F, row.names = F)
  if (verbose) 
    cat("> Munging exposure data... \n")
  invisible(utils::capture.output(GenomicSEM::munge(paste0(exposure_name, 
                                                           ".tsv"), hm3, exposure_name)))
  if (verbose & save_logfiles) 
    cat("  Please check the log file", paste0(exposure_name, 
                                              "_munge.log"), "to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files\n")
  if (verbose) 
    cat("> Munging outcome data... \n")
  invisible(utils::capture.output(GenomicSEM::munge(paste0(outcome_name, 
                                                           ".tsv"), hm3, outcome_name)))
  if (verbose & save_logfiles) 
    cat("  Please check the log file", paste0(outcome_name, 
                                              "_munge.log"), "to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files\n")
  traits <- c(paste0(exposure_name, ".sumstats.gz"), paste0(outcome_name, 
                                                            ".sumstats.gz"))
  sample.prev <- c(NA, NA)
  population.prev <- c(NA, NA)
  trait.names <- c(exposure_name, outcome_name)
  if (verbose) 
    cat("> Running cross-trait LDSC... \n")
  invisible(utils::capture.output(LDSCoutput <- GenomicSEM::ldsc(traits = traits, 
                                                                 sample.prev = sample.prev, population.prev = population.prev, 
                                                                 ld = ld, wld = ld, trait.names = trait.names)))
  if (verbose & save_logfiles) 
    cat("  Please check the log file", paste0(c(traits, 
                                                "ldsc.log"), collapse = "_"), "for detailed results of the cross-trait LDSC analysis\n")
  if (verbose & !save_logfiles) 
    cat("  Please consider saving the log files and checking them to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files\n")
  h2_exp = as.numeric(LDSCoutput$S[1, 1])
  h2_exp_SE = sqrt(LDSCoutput$V[1, 1])
  gcov_int = as.numeric(LDSCoutput$I[1, 2])
  output = readLines(paste0(c(traits, "ldsc.log"), collapse = "_"))
  gcov_int_SE = as.numeric(stringr::str_replace(stringr::str_split(output[stringr::str_detect(output, 
                                                                                              "Cross trait Intercept")], "\\(")[[1]][2], "\\)", ""))
  int_exp = LDSCoutput$I[1, 1]
  int_out = LDSCoutput$I[2, 2]
  h2_out = as.numeric(LDSCoutput$S[2, 2])
  h2_out_SE = sqrt(LDSCoutput$V[3, 3])
  rgcov = as.numeric(LDSCoutput$S[1, 2])
  rg = as.numeric(LDSCoutput$S[1, 2]/sqrt(LDSCoutput$S[1, 
                                                       1] * LDSCoutput$S[2, 2]))
  rgcov_SE = sqrt(LDSCoutput$V[2, 2])
  if (verbose) 
    cat("> Cleaning temporary files... \n")
  if (save_logfiles) {
    file.remove(paste0(exposure_name, c(".tsv", ".sumstats.gz")))
  }
  else {
    file.remove(paste0(exposure_name, c(".tsv", "_munge.log", 
                                        ".sumstats.gz")))
  }
  if (save_logfiles) {
    file.remove(paste0(outcome_name, c(".tsv", ".sumstats.gz")))
  }
  else {
    file.remove(paste0(outcome_name, c(".tsv", "_munge.log", 
                                       ".sumstats.gz")))
  }
  if (!save_logfiles) 
    file.remove(paste0(c(traits, "ldsc.log"), collapse = "_"))
  return(list(h2_LDSC = h2_exp, h2_LDSC_se = h2_exp_SE, lambda = gcov_int, 
              lambda_se = gcov_int_SE, int_exp = int_exp, int_out = int_out, 
              h2_out = h2_out, h2_out_se = h2_out_SE, rgcov = rgcov, 
              rgcov_se = rgcov_SE, rg = rg))
}

run_MR <- function (exposure_data, outcome_data, MR_threshold = 5e-08, 
          MR_pruning_dist = 500, MR_pruning_LD = 0, MR_reverse = NULL, 
          verbose = TRUE) 
{
  data <- dplyr::inner_join(exposure_data, outcome_data, by = c("rsid", 
                                                                "chr", "pos"), suffix = c(".exp", ".out"))
  data <- data %>% dplyr::mutate(std_beta.out = dplyr::case_when(alt.exp == 
                                                                   alt.out & ref.exp == ref.out ~ std_beta.out, ref.exp == 
                                                                   alt.out & alt.exp == ref.out ~ -std_beta.out, TRUE ~ 
                                                                   NA_real_)) %>% dplyr::filter(!is.na(std_beta.out))
  if (verbose) 
    cat("> Identifying IVs... \n")
  data_thresholded <- data %>% dplyr::filter(.data$p.exp < 
                                               MR_threshold)
  if (verbose) 
    cat("   ", format(nrow(data_thresholded), big.mark = ","), 
        "IVs with p <", format(MR_threshold, scientific = T), 
        "\n")
  if (nrow(data_thresholded) == 0) 
    stop("no IV left after thresholding")
  if (!is.null(MR_reverse)) {
    reverse_t_threshold = stats::qnorm(MR_reverse)
    data_thresholded_filtered <- data_thresholded %>% dplyr::filter((abs(.data$std_beta.exp) - 
                                                                       abs(.data$std_beta.out))/sqrt(.data$std_SE.exp^2 + 
                                                                                                       .data$std_SE.out^2) > reverse_t_threshold)
    if (verbose) 
      cat("   ", nrow(data_thresholded) - nrow(data_thresholded_filtered), 
          "IVs excluded - more strongly associated with the outcome than with the exposure, p <", 
          format(MR_reverse, scientific = T), "\n")
    data_thresholded = data_thresholded_filtered
    rm(data_thresholded_filtered)
  }
  if (nrow(data_thresholded) == 0) 
    stop("no IV left after excluding IVs more strongly associated with the outcome than with the exposure")
  ToPrune <- data_thresholded %>% dplyr::transmute(SNP = .data$rsid, 
                                                   chr_name = .data$chr, chr_start = .data$pos, pval.exposure = .data$p.exp)
  if (MR_pruning_LD > 0) {
    if (verbose) 
      cat("   Pruning : distance : ", MR_pruning_dist, 
          "Kb", " - LD threshold : ", MR_pruning_LD, "\n")
    SNPsToKeep = c()
    for (chr in unique(ToPrune$chr_name)) {
      SNPsToKeep = c(SNPsToKeep, suppressMessages(TwoSampleMR::clump_data(ToPrune[ToPrune$chr_name == 
                                                                                    chr, ], clump_kb = MR_pruning_dist, clump_r2 = MR_pruning_LD, pop='EAS')$SNP))
    }
  }
  else {
    prune_byDistance <- function(data, prune.dist = 100, 
                                 byP = T) {
      if (byP) {
        SNP_order = order(data %>% dplyr::pull(4))
      }
      else {
        SNP_order = order(data %>% dplyr::pull(4), decreasing = T)
      }
      data = data[SNP_order, ]
      snp = 0
      while (T) {
        snp = snp + 1
        ToRemove = which(data$chr_name == data$chr_name[snp] & 
                           abs(data$chr_start - data$chr_start[snp]) < 
                           prune.dist * 1000)
        if (length(ToRemove) > 1) {
          ToRemove = ToRemove[-1]
          data = data[-ToRemove, ]
        }
        if (snp == nrow(data)) 
          break
      }
      return(unlist(data[, 1]))
    }
    if (verbose) 
      cat("   Pruning : distance : ", MR_pruning_dist, 
          "Kb \n")
    SNPsToKeep = prune_byDistance(ToPrune, prune.dist = MR_pruning_dist, 
                                  byP = T)
  }
  data_pruned <- data_thresholded %>% dplyr::filter(.data$rsid %in% 
                                                      SNPsToKeep)
  if (verbose) 
    cat("   ", format(nrow(data_pruned), big.mark = ","), 
        "IVs left after pruning \n")
  if (verbose) 
    cat("> Performing MR \n")
  res_MR_TwoSampleMR <- TwoSampleMR::mr_ivw(data_pruned$std_beta.exp, 
                                            data_pruned$std_beta.out, data_pruned$std_SE.exp, data_pruned$std_SE.out)
  if (verbose) 
    cat("   ", "IVW-MR observed effect:", format(res_MR_TwoSampleMR$b, 
                                                 digits = 3), "(", format(res_MR_TwoSampleMR$se, 
                                                                          digits = 3), ")\n")
  return(list(alpha_obs = res_MR_TwoSampleMR$b, alpha_obs_se = res_MR_TwoSampleMR$se, 
              n_exp = mean(data_pruned$N.exp), n_out = mean(data_pruned$N.out), 
              IVs = data_pruned %>% dplyr::select(.data$std_beta.exp, 
                                                  .data$std_SE.exp), IVs_rs = data_pruned$rsid))
}

get_correction <- function (IVs, lambda, lambda_se, h2_LDSC, h2_LDSC_se, alpha_obs, 
          alpha_obs_se, n_exp, n_out, MR_threshold, verbose) 
{
  M = 1150000
  Tr = -stats::qnorm(MR_threshold/2)
  lambdaPrime = lambda/sqrt(n_exp * n_out)
  if (verbose) 
    cat("> Estimating genetic architecture parameters... \n")
  get_pi <- function(my_pi, sumbeta2, Tr, n_exp, h2_LDSC, 
                     M) {
    if (0 >= my_pi) 
      return(1e+06)
    sigma = sqrt(h2_LDSC/(my_pi * M))
    denominator = (my_pi * (2 * (sigma^2 + 1/n_exp) * stats::pnorm(-Tr/sqrt(1 + 
                                                                              n_exp * sigma^2)) + 2 * Tr * (n_exp * sigma^4 + 
                                                                                                              2 * sigma^2 + 1/n_exp) * exp(-Tr^2/(2 * (n_exp * 
                                                                                                                                                         sigma^2 + 1)))/(sqrt(2 * pi) * (1 + n_exp * sigma^2)^(3/2))) + 
                     (1 - my_pi) * 1/n_exp * (2 * stats::pnorm(-Tr) + 
                                                2 * Tr * stats::dnorm(Tr))) * M
    return(abs(denominator - sumbeta2))
  }
  get_geneticArchitecture <- function(theta, Nexp, M, Tr) {
    h2_LDSC = theta[length(theta)]
    effects = theta[-length(theta)]
    sumBeta2 = sum(effects^2)
    nSP = 5
    Res_SP = data.frame(SP = 1:nSP, SP_pi = NA_real_, diff = NA_real_, 
                        pi = NA_real_)
    for (i in 1:nSP) {
      theta = 3 * 10^(stats::runif(1, -7, -1))
      res_optim = stats::optimise(get_pi, interval = c(1e-07, 
                                                       0.3), sumbeta2 = sumBeta2, Tr = Tr, n_exp = n_exp, 
                                  h2_LDSC = h2_LDSC, M = M, tol = 1e-06)
      Res_SP[i, 2:4] = c(theta, res_optim$objective, res_optim$minimum)
    }
    pi_x <- Res_SP %>% dplyr::arrange(diff) %>% dplyr::slice(1) %>% 
      dplyr::pull(pi)
    sigma = sqrt(h2_LDSC/(M * pi_x))
    return(c(as.numeric(pi_x), as.numeric(sigma)))
  }
  theta = c(IVs$std_beta.exp, h2_LDSC)
  res_genA = get_geneticArchitecture(theta, n_exp, M, Tr)
  get_alpha <- function(n_exp, lambdaPrime, pi_x, sigma, alpha_obs, 
                        Tr) {
    A = stats::pnorm(-Tr/sqrt(1 + n_exp * sigma^2))
    B = 2 * Tr * exp(-Tr^2/(2 * (n_exp * sigma^2 + 1)))/(sqrt(2 * 
                                                                pi) * (1 + n_exp * sigma^2)^(3/2))
    C = (stats::pnorm(-Tr) + Tr * stats::dnorm(Tr))
    a = (pi_x * (2 * (sigma^2 + 1/n_exp) * A + B * (n_exp * 
                                                      sigma^4 + 2 * sigma^2 + 1/n_exp)))
    b = (1 - pi_x) * 2/n_exp * C
    d = a + b
    alpha = (alpha_obs * d - (lambdaPrime * pi_x * (2 * 
                                                      A + B + sigma^2 * n_exp * B) + lambdaPrime * (1 - 
                                                                                                      pi_x) * 2 * C))/(pi_x * sigma^2 * (2 * A + B * n_exp * 
                                                                                                                                           sigma^2 + B))
    return(alpha)
  }
  if (verbose) 
    cat("> Estimating corrected effect... \n")
  alpha_corrected = get_alpha(n_exp, lambdaPrime, res_genA[1], 
                              res_genA[2], alpha_obs, Tr)
  get_correctedSE <- function(IVs, lambda, lambda_se, h2_LDSC, 
                              h2_LDSC_se, alpha_obs, alpha_obs_se, n_exp, n_out, M, 
                              Tr, s = 1000, sthreshold = 0.05, extracheck = T) {
    get_s <- function(s) {
      effects = IVs$std_beta.exp
      effects_se = IVs$std_SE.exp
      L = matrix(stats::rnorm(s, lambda, lambda_se), ncol = s)/sqrt(n_exp * 
                                                                      n_out)
      E = matrix(stats::rnorm(nrow(IVs) * s, effects, 
                              effects_se), ncol = s)
      negativeH2 = FALSE
      H = matrix(stats::rnorm(s, h2_LDSC, h2_LDSC_se), 
                 ncol = s)
      if (any(H < 0)) {
        negativeH2 = TRUE
        while (!all(H > 0)) {
          H[H < 0] = stats::rnorm(length(H[H < 0]), 
                                  h2_LDSC, h2_LDSC_se)
        }
      }
      D = rbind(E, H)
      pis = apply(D, 2, function(x) get_geneticArchitecture(x, 
                                                            n_exp, M, Tr))
      B = matrix(stats::rnorm(s, alpha_obs, alpha_obs_se), 
                 ncol = s)
      all_params = data.frame(pi = pis[1, ], sigma = pis[2, 
      ], alpha = B[1, ], lambda = L[1, ])
      all_params$corrected = apply(all_params, 1, function(x) get_alpha(n_exp, 
                                                                        x[4], x[1], x[2], x[3], Tr))
      all_params$warning_negH2 = negativeH2
      return(all_params)
    }
    res = get_s(s)
    num_groups = 10
    res_subsets <- res %>% dplyr::group_by((dplyr::row_number() - 
                                              1)%/%(dplyr::n()/num_groups)) %>% tidyr::nest() %>% 
      dplyr::pull(data)
    subsets_var = unlist(lapply(res_subsets, function(x) stats::var(x$corrected)))
    subsets_cov = unlist(lapply(res_subsets, function(x) stats::cov(x$corrected, 
                                                                    x$alpha)))
    needmore = F
    if (stats::sd(subsets_var)/base::mean(subsets_var) > 
        sthreshold) 
      needmore = T
    if (stats::sd(subsets_cov)/base::mean(subsets_cov) > 
        sthreshold) 
      needmore = T
    if (extracheck & (alpha_obs_se^2 + stats::sd(res$corrected)^2 - 
                      2 * stats::cov(res$alpha, res$corrected)) < 0) 
      needmore = T
    tmp_sd_corrected = stats::sd(res$corrected)
    while (needmore) {
      res = rbind(res, get_s(s))
      res_subsets <- res %>% dplyr::filter(.data$corrected < 
                                             alpha_corrected + 10 * tmp_sd_corrected, .data$corrected > 
                                             alpha_corrected - 10 * tmp_sd_corrected) %>% 
        dplyr::group_by((dplyr::row_number() - 1)%/%(dplyr::n()/num_groups)) %>% 
        tidyr::nest() %>% dplyr::pull(data)
      subsets_var = unlist(lapply(res_subsets, function(x) stats::var(x$corrected)))
      subsets_cov = unlist(lapply(res_subsets, function(x) stats::cov(x$corrected, 
                                                                      x$alpha)))
      needmore = F
      if (stats::sd(subsets_var)/base::mean(subsets_var) > 
          sthreshold) 
        needmore = T
      if (stats::sd(subsets_cov)/base::mean(subsets_cov) > 
          sthreshold) 
        needmore = T
      if (extracheck & (alpha_obs_se^2 + stats::sd(res$corrected)^2 - 
                        2 * stats::cov(res$alpha, res$corrected)) < 
          0) 
        needmore = T
    }
    all_res = c(stats::sd(res$corrected), stats::cov(res$alpha, 
                                                     res$corrected), nrow(res), any(res$warning_negH2))
    return(all_res)
  }
  se_cov = get_correctedSE(IVs, lambda, lambda_se, h2_LDSC, 
                           h2_LDSC_se, alpha_obs, alpha_obs_se, n_exp, n_out, M, 
                           Tr)
  if (verbose) 
    cat("   ", "corrected effect:", format(alpha_corrected, 
                                           digits = 3), "(", format(se_cov[1], digits = 3), 
        ")\n")
  if (verbose) 
    cat("   ", "covariance between observed and corrected effect:", 
        format(se_cov[2], digits = 3), "\n")
  if (verbose) 
    cat("           ", se_cov[3], "simulations were used to estimate the variance and the covariance.\n")
  if (se_cov[4]) 
    warning("some h2 were negative in the parametric bootstrap, please check that the heritability of the exposure is not too low as this might compromise the results.\n")
  if (verbose) 
    cat("> Testing difference between observed and corrected effect... \n")
  test_diff = (alpha_obs - alpha_corrected)/sqrt(alpha_obs_se^2 + 
                                                   se_cov[1]^2 - 2 * se_cov[2])
  p_diff = 2 * stats::pnorm(-abs(test_diff), lower.tail = T)
  return(list(alpha_corrected = alpha_corrected, alpha_corrected_se = se_cov[1], 
              cov_obs_corrected = se_cov[2], test_diff = test_diff, 
              p_diff = p_diff, pi_x = res_genA[1], sigma2_x = res_genA[2]^2))
}



MRlap <- function (exposure, exposure_name = NULL, outcome, outcome_name = NULL, 
          ld, hm3, MR_threshold = 5e-08, MR_pruning_dist = 500, MR_pruning_LD = 0, 
          MR_reverse = 0.001, save_logfiles = FALSE, verbose = TRUE) 
{
  InitPath = getwd()
  on.exit(setwd(InitPath))
  StartTime = proc.time()
  platform = .Platform$OS.type
  if (!is.logical(verbose)) 
    stop("verbose : should be logical", call. = FALSE)
  if (!is.logical(save_logfiles)) 
    stop("save_logfiles : should be logical", call. = FALSE)
  if (verbose) 
    cat("<<< Preparation of analysis >>> \n")
  if (verbose) 
    cat("> Checking parameters \n")
  if (is.data.frame(exposure)) {
    attributes(exposure)$GName = deparse(substitute(exposure))
  }
  if (is.data.frame(outcome)) {
    attributes(outcome)$GName = deparse(substitute(outcome))
  }
  if (is.character(ld)) {
    if (!dir.exists(ld)) 
      stop("ld : the folder does not exist", call. = FALSE)
    ld = normalizePath(ld)
  }
  else stop("ld : wrong format, should be character", call. = FALSE)
  if (is.character(hm3)) {
    if (!file.exists(hm3)) 
      stop("hm3 : the file does not exist", call. = FALSE)
    hm3 = normalizePath(hm3)
  }
  else stop("hm3 : wrong format, should be character", call. = FALSE)
  if (!is.numeric(MR_threshold)) 
    stop("MR_threshold : non-numeric argument", call. = FALSE)
  if (MR_threshold > 10^-5) 
    stop("MR_threshold : superior to the threshold limit", 
         call. = FALSE)
  if (verbose) 
    cat("The p-value threshold used for selecting MR instruments is:", 
        format(MR_threshold, scientific = T), "\n")
  if (!is.numeric(MR_pruning_dist)) 
    stop("MR_pruning_dist : non-numeric argument", call. = FALSE)
  if (MR_pruning_dist < 10) 
    stop("MR_pruning_dist : should be higher than 10Kb", 
         call. = FALSE)
  if (MR_pruning_dist > 50000) 
    stop("MR_pruning_dist : should be lower than 50Mb", 
         call. = FALSE)
  if (verbose) 
    cat("The distance used for pruning MR instruments is: ", 
        MR_pruning_dist, "Kb \n")
  if (!is.numeric(MR_pruning_LD)) 
    stop("MR_pruning_LD : non-numeric argument", call. = FALSE)
  if (MR_pruning_LD < 0) 
    stop("MR_pruning_LD : should be positive", call. = FALSE)
  if (MR_pruning_LD > 1) 
    stop("MR_pruning_LD : should not be larger than 1", 
         call. = FALSE)
  if (MR_pruning_LD > 0) {
    if (verbose) 
      cat("The LD threshold used for pruning MR instruments is:", 
          MR_pruning_LD, "\n")
  }
  else {
    if (verbose) 
      cat("Distance-based pruning will be used for MR instruments \n")
  }
  if (!is.null(exposure_name) & !is.character(exposure_name)) 
    stop("exposure_name : non-character argument", call. = FALSE)
  if (!is.null(outcome_name) & !is.character(outcome_name)) 
    stop("outcome_name : non-character argument", call. = FALSE)
  if (verbose) 
    cat(paste0("> Processing exposure ", ifelse(is.null(exposure_name), 
                                                "", paste0("(", exposure_name, ") ")), "summary statistics... \n"))
  if (is.null(exposure_name)) 
    exposure_name = "exposure"
  exposure_data = tidy_inputGWAS(exposure, verbose)
  if (verbose) 
    cat(paste0("> Processing outcome ", ifelse(is.null(outcome_name), 
                                               "", paste0("(", outcome_name, ") ")), "summary statistics... \n"))
  if (is.null(outcome_name)) 
    outcome_name = "outcome"
  outcome_data = tidy_inputGWAS(outcome, verbose)
  if (verbose) 
    cat("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n")
  if (verbose) 
    cat("<<< Performing cross-trait LDSC >>>  \n")
  LDSC_results = run_LDSC(exposure_data, exposure_name, outcome_data, 
                          outcome_name, ld, hm3, save_logfiles, verbose)
  if (verbose) 
    cat("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n")
  if (verbose) 
    cat("<<< Running IVW-MR >>>  \n")
  MR_results = run_MR(exposure_data, outcome_data, MR_threshold, 
                      MR_pruning_dist, MR_pruning_LD, MR_reverse, verbose)
  if (verbose) 
    cat("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n")
  if (verbose) 
    cat("<<< Estimating corrected effect >>>  \n")
  correction_results = with(c(MR_results, LDSC_results), get_correction(IVs, 
                                                                        lambda, lambda_se, h2_LDSC, h2_LDSC_se, alpha_obs, alpha_obs_se, 
                                                                        n_exp, n_out, MR_threshold, verbose))
  tmp = "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"
  if (is.na(correction_results$alpha_corrected_se)) {
    if (verbose) 
      cat("WARNING: the sampling strategy failed in the estimation of the standard error.\n")
    if (verbose) 
      cat("Please try to increase the number of simulations s. \n")
  }
  Time = as.integer((proc.time() - StartTime)[3])
  minutes <- as.integer(trunc(Time/60))
  seconds <- Time - minutes * 60
  if (verbose) 
    cat("Runtime of the analysis: ", minutes, " minute(s) and ", 
        seconds, " second(s).  \n")
  results_MR = with(c(MR_results, correction_results), list(observed_effect = alpha_obs, 
                                                            observed_effect_se = alpha_obs_se, m_IVs = nrow(IVs), 
                                                            IVs = IVs_rs, observed_effect_p = 2 * stats::pnorm(-abs(alpha_obs/alpha_obs_se)), 
                                                            corrected_effect = alpha_corrected, corrected_effect_se = alpha_corrected_se, 
                                                            corrected_effect_p = 2 * stats::pnorm(-abs(alpha_corrected/alpha_corrected_se)), 
                                                            test_difference = test_diff, p_difference = p_diff))
  results_LDSC = with(LDSC_results, list(h2_exp = h2_LDSC, 
                                         h2_exp_se = h2_LDSC_se, int_exp = int_exp, h2_out = h2_out, 
                                         h2_out_se = h2_out_se, int_out = int_out, gcov = rgcov, 
                                         gcov_se = rgcov_se, rg = rg, int_crosstrait = lambda, 
                                         int_crosstrait_se = lambda_se))
  results_GeneticArchitecture = with(correction_results, list(polygenicity = pi_x, 
                                                              perSNP_heritability = sigma2_x))
  results = list(MRcorrection = results_MR, LDSC = results_LDSC, 
                 GeneticArchitecture = results_GeneticArchitecture)
  return(results)
}