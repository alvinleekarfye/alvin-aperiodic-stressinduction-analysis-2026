rm(list = ls())

library(pbapply)
library(parallel)
library(mgcv)
library(emmeans) 

df_combined <- readRDS("C:/df_combined_delta.rds")
df_combined$Subject <- as.factor(df_combined$Subject)
df_combined$Gender <- as.factor(df_combined$Gender)
df_combined$ROI <- as.factor(df_combined$ROI)
df_combined$Time <- as.factor(df_combined$Time)

gam_model_delta <- gam(
  Residual_Delta_bc ~
    ROI * Time + Age + Gender + Cortisol_Baseline +
    s(Subject, bs = "re"),
  data = df_combined,
  family = scat,
  method = "REML"
)

gam_model_delta_2 <- gam(
  Residual_Delta_bc ~ 
    ROI + Time + Age + Gender + Cortisol_Baseline +
    s(Subject, bs = "re"),
  data = df_combined,
  family = scat,
  method = "REML"
)

wild_boot_pTerms_table_parallel <- function(model, null_model, data,
                                            subject_var = "Subject",
                                            nboot = 1000,
                                            ncores = 4,
                                            weight_type = c("rademacher", "mammen", "webb"),
                                            seed = 123) {
  
  weight_type <- match.arg(weight_type)
  
  anova_orig <- anova(model)
  ptab0 <- anova_orig$pTerms.table
  term_names <- rownames(ptab0)
  Chi_orig <- ptab0[, "Chi.sq"]
  
  ptab <- summary(model)$p.table
  coef_names <- rownames(ptab)
  coef_orig <- ptab[, "Estimate"]
  
  resp_name <- as.character(formula(model)[[2]])
  form_orig <- formula(model)
  orig_sp   <- model$sp
  
  yhat  <- fitted(null_model)
  resid <- residuals(null_model, type = "response")
  
  yhat_full  <- fitted(model)
  resid_full <- residuals(model, type = "response")
  
  subjects <- unique(data[[subject_var]])
  n_sub <- length(subjects)
  
  em_orig <- emmeans::emmeans(model, revpairwise ~ Time | ROI)
  posthoc_df <- as.data.frame(em_orig$contrasts)
  contrast_names <- paste(posthoc_df$ROI, posthoc_df$contrast, sep = " | ")
  posthoc_orig <- posthoc_df$estimate
  
  bootstrap_one <- function(b) {
    
    if (weight_type == "rademacher") {
      w <- sample(c(-1, 1), size = n_sub, replace = TRUE)
    } else if (weight_type == "mammen") {
      w_vals <- c((1 - sqrt(5))/2, (1 + sqrt(5))/2)
      probs <- c((sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)))
      w <- sample(w_vals, size = n_sub, replace = TRUE, prob = probs)
    } else if (weight_type == "webb") {
      w_vals <- c(-0.857, -0.571, -0.429, 0.429, 0.571, 0.857)
      w <- sample(w_vals, size = n_sub, replace = TRUE)
    }
    
    names(w) <- subjects
    w_i <- unname(w[as.character(data[[subject_var]])])
    
    if(length(w_i) != nrow(data)) stop("Weight mapping failed")
    
    data_boot <- data
    data_boot[[resp_name]] <- yhat + resid * w_i
    
    m_boot <- try(
      mgcv::gam(formula = form_orig,
                data = data_boot,
                family = model$family,
                method = "REML",
                sp = orig_sp),
      silent = TRUE
    )
    
    if (inherits(m_boot, "try-error")) {
      return(list(
        Chi = rep(NA, length(term_names)),
        coef = rep(NA, length(coef_names)),
        posthoc = rep(NA, length(contrast_names)),
        coef_full = rep(NA, length(coef_names)),
        posthoc_full = rep(NA, length(contrast_names))
      ))
    }
    
    a_boot <- try(anova(m_boot)$pTerms.table, silent = TRUE)
    if (inherits(a_boot, "try-error") || !all(term_names %in% rownames(a_boot))) {
      Chi_boot <- rep(NA, length(term_names))
    } else {
      Chi_boot <- a_boot[term_names, "Chi.sq"]
    }
    
    s <- summary(m_boot)$p.table
    if (!all(coef_names %in% rownames(s))) {
      coef_boot <- rep(NA, length(coef_names))
    } else {
      coef_boot <- s[coef_names, "Estimate"]
    }
    
    em_boot <- try(emmeans::emmeans(m_boot, revpairwise ~ Time | ROI), silent = TRUE)
    
    if (inherits(em_boot, "try-error")) {
      posthoc_boot <- rep(NA, length(contrast_names))
    } else {
      tmp <- as.data.frame(em_boot$contrasts)
      tmp_names <- paste(tmp$ROI, tmp$contrast, sep = " | ")
      
      if (!all(contrast_names %in% tmp_names)) {
        posthoc_boot <- rep(NA, length(contrast_names))
      } else {
        tmp <- tmp[match(contrast_names, tmp_names), ]
        posthoc_boot <- tmp$estimate
      }
    }

    data_boot_full <- data
    data_boot_full[[resp_name]] <- yhat_full + resid_full * w_i
    
    m_full <- try(
      mgcv::gam(formula = form_orig,
                data = data_boot_full,
                family = model$family,
                method = "REML",
                sp = orig_sp),
      silent = TRUE
    )
    
    if (inherits(m_full, "try-error")) {
      coef_full <- rep(NA, length(coef_names))
      posthoc_full <- rep(NA, length(contrast_names))
    } else {
      
      s2 <- summary(m_full)$p.table
      if (!all(coef_names %in% rownames(s2))) {
        coef_full <- rep(NA, length(coef_names))
      } else {
        coef_full <- s2[coef_names, "Estimate"]
      }
      
      em_full <- try(emmeans::emmeans(m_full, revpairwise ~ Time | ROI), silent = TRUE)
      
      if (inherits(em_full, "try-error")) {
        posthoc_full <- rep(NA, length(contrast_names))
      } else {
        tmp2 <- as.data.frame(em_full$contrasts)
        tmp2_names <- paste(tmp2$ROI, tmp2$contrast, sep = " | ")
        
        if (!all(contrast_names %in% tmp2_names)) {
          posthoc_full <- rep(NA, length(contrast_names))
        } else {
          tmp2 <- tmp2[match(contrast_names, tmp2_names), ]
          posthoc_full <- tmp2$estimate
        }
      }
    }
    
    list(
      Chi = Chi_boot,
      coef = coef_boot,
      posthoc = posthoc_boot,
      
      coef_full = coef_full,
      posthoc_full = posthoc_full
    )
  }
  
  set.seed(seed)
  cl <- makeCluster(ncores)
  clusterSetRNGStream(cl, seed)
  
  clusterExport(
    cl,
    varlist = c("data", "yhat", "resid",
                "yhat_full", "resid_full",
                "subjects", "n_sub",
                "resp_name", "form_orig",
                "term_names", "coef_names",
                "model", "weight_type",
                "orig_sp", "contrast_names"),
    envir = environment()
  )
  
  clusterEvalQ(cl, { library(mgcv); library(emmeans) })
  
  boot_list <- pbapply::pblapply(seq_len(nboot), bootstrap_one, cl = cl)
  stopCluster(cl)
  
  Chi_boot <- do.call(rbind, lapply(boot_list, `[[`, "Chi"))
  coef_boot <- do.call(rbind, lapply(boot_list, `[[`, "coef"))
  posthoc_boot <- do.call(rbind, lapply(boot_list, `[[`, "posthoc"))
  
  coef_boot_full <- do.call(rbind, lapply(boot_list, `[[`, "coef_full"))
  posthoc_boot_full <- do.call(rbind, lapply(boot_list, `[[`, "posthoc_full"))
  
  colnames(Chi_boot) <- term_names
  colnames(coef_boot) <- coef_names
  colnames(posthoc_boot) <- contrast_names
  colnames(coef_boot_full) <- coef_names
  colnames(posthoc_boot_full) <- contrast_names
  
  p_terms_boot <- sapply(seq_along(Chi_orig), function(j) {
    mean(Chi_boot[, j] >= Chi_orig[j], na.rm = TRUE)
  })
  names(p_terms_boot) <- term_names
  
  bootstrap_p_coef <- sapply(seq_along(coef_names), function(j) {
    mean(abs(coef_boot[, j]) >= abs(coef_orig[j]), na.rm = TRUE)
  })
  names(bootstrap_p_coef) <- coef_names
  
  posthoc_p <- sapply(seq_along(posthoc_orig), function(j) {
    mean(abs(posthoc_boot[, j]) >= abs(posthoc_orig[j]), na.rm = TRUE)
  })
  names(posthoc_p) <- contrast_names
  
  boot_se <- apply(coef_boot_full, 2, sd, na.rm = TRUE)
  boot_ci <- t(apply(coef_boot_full, 2, quantile,
                     probs = c(0.005, 0.995), na.rm = TRUE))
  
  posthoc_se <- apply(posthoc_boot_full, 2, sd, na.rm = TRUE)
  posthoc_ci <- t(apply(posthoc_boot_full, 2, quantile,
                        probs = c(0.005, 0.995), na.rm = TRUE))

  posthoc_p_holm <- p.adjust(posthoc_p, method = "holm")
  
  return(list(
    original_pTerms_table = ptab0,
    bootstrap_p_terms = p_terms_boot,
    bootstrap_p_coef = bootstrap_p_coef,
    boot_se = boot_se,
    boot_ci = boot_ci,
    posthoc_p = posthoc_p,
    posthoc_p_holm = posthoc_p_holm,
    posthoc_se = posthoc_se,
    posthoc_ci = posthoc_ci
  ))
}

boot_parallel_delta <- wild_boot_pTerms_table_parallel(
  gam_model_delta,
  gam_model_delta_2,
  df_combined,
  subject_var = "Subject",
  nboot = 10000,
  ncores = 16,
  weight_type = "rademacher",
  seed = 888
)
saveRDS(boot_parallel_delta, "C:/bootstrap_delta_results.rds")
