rm(list = ls())

library(mgcv)
library(emmeans)
library(eegUtils)
library(ggplot2)
library(dplyr)
library(patchwork)
library(e1071)
library(DHARMa)

df_combined <- readRDS("C:/df_combined_covariation.rds")
df_combined$Subject <- as.factor(df_combined$Subject)
df_combined$Gender <- as.factor(df_combined$Gender)
df_combined$ROI <- as.factor(df_combined$ROI)
df_combined$Time <- as.factor(df_combined$Time)

df_diff <- df_combined %>%
  group_by(Subject, ROI) %>%
  summarise(
    Exponent_diff = Exponent[Time == "T3"] - Exponent[Time == "T2"],
    Theta_diff    = Residual_Theta[Time == "T3"] - Residual_Theta[Time == "T2"],
    Delta_diff    = Residual_Delta[Time == "T3"] - Residual_Delta[Time == "T2"],
    Beta_diff    = Residual_Beta[Time == "T3"] - Residual_Beta[Time == "T2"],
    
    Age           = first(Age),
    Gender        = first(Gender),
    Cortisol_Baseline   = first(Cortisol_Baseline),
    .groups = "drop"
  )

gam_model_delta_coupling <- gam(
  Delta_diff ~
    Exponent_diff * ROI + Age + Gender + Cortisol_Baseline +
    s(Subject, bs = "re"),
  data = df_diff,
  family = scat,
  method = "REML"
)
summary(gam_model_delta_coupling)
anova(gam_model_delta_coupling)
gam.check(gam_model_delta_coupling)

print(gam_model_delta_coupling$df.residual)
skewness(resid(gam_model_delta_coupling))
kurtosis(resid(gam_model_delta_coupling))

sim <- simulateResiduals(fittedModel = gam_model_delta_coupling, re.form = NULL, n = 10000)
plot(sim)  

roi_slopes <- emtrends(gam_model_delta_coupling, ~ ROI, var = "Exponent_diff")
roi_slopes_summary <- summary(roi_slopes, infer = c(TRUE, TRUE), level = 0.99)
roi_slopes_df <- as.data.frame(roi_slopes_summary)
roi_slopes_df$p.adjusted <- p.adjust(roi_slopes_df$p.value, method = "holm")
significant_rois <- roi_slopes_df %>%
  filter(p.adjusted < 0.05)
significant_rois

posthoc_df <- roi_slopes_df %>%
  rename(estimate = Exponent_diff.trend) %>%  
  mutate(
    signif = case_when(
      p.adjusted < 0.001 ~ "***",
      p.adjusted < 0.01  ~ "**",
      p.adjusted < 0.05  ~ "*",
      TRUE ~ " "
    )
  )

roi_order <- c("O2", "O1", "P8", "P4", "Pz", "P3", "P7", "T4", "T3", "C4", "Cz", "C3", "F8", "F4", "Fz", "F3", "F7", "Fp2", "Fp1")

posthoc_df <- posthoc_df %>%
  mutate(ROI = factor(ROI, levels = roi_order))

p_coupling_forest <- ggplot(posthoc_df, aes(x = estimate, y = ROI)) +
  geom_point(size = 3, aes(color = signif)) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.5) +
  geom_text(aes(x = upper.CL, label = signif), hjust = -0.3, size = 4) +
  scale_color_manual(values = c("*" = "green", "**" = "green", "***" = "green")) +
  labs(x = "Aperiodic-Oscillatory β", y = "Site", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none")
p_coupling_forest <- p_coupling_forest + coord_cartesian(xlim = c(-0.50, 0.50))

plot_emm <- emtrends(gam_model_delta_coupling, ~ ROI, var = "Exponent_diff")
emm_df <- as.data.frame(plot_emm)

df_topo <- emm_df %>%
  rename(electrode = ROI, quantity = Exponent_diff.trend) %>%
  select(electrode, quantity) %>%
  mutate(
    electrode = as.factor(electrode),
  )

p_coupling <- topoplot(df_topo,
                       quantity = "quantity",
                       chan_marker  = "name",
                       method = "Biharmonic",
                       head = TRUE,
                       palette = "YlOrRd",
                       fill_title="β")
p_coupling <- p_coupling +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.50,0.50),
    oob = scales::squish,
    name = "β"
  )
p_coupling <- p_coupling + labs(title = "Exponent-Delta") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p_coupling

p_coupling2 <- topoplot(df_topo,
                        quantity = "quantity",
                        chan_marker  = "name",
                        method = "Biharmonic",
                        head = TRUE,
                        palette = "YlOrRd",
                        fill_title="β")
p_coupling2 <- p_coupling2 +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.50,0.50),
    oob = scales::squish,
    name = "β"
  )
p_coupling2 <- p_coupling2 + labs(title = "Exponent-Delta") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p_coupling2 <- p_coupling2 + theme(legend.position = "none")

combined <- p_coupling_forest + p_coupling + p_coupling2 +
  plot_layout(widths = c(0.75, 1, 1))
combined

gam_model_theta_coupling <- gam(
  Theta_diff ~ 
    Exponent_diff * ROI + Age + Gender + Cortisol_Baseline +
    s(Subject, bs = "re"),    
  data = df_diff,
  family = scat,
  method = "REML"
)
summary(gam_model_theta_coupling)
anova(gam_model_theta_coupling)
gam.check(gam_model_theta_coupling)

print(gam_model_theta_coupling$df.residual)
skewness(resid(gam_model_theta_coupling))
kurtosis(resid(gam_model_theta_coupling))

sim <- simulateResiduals(fittedModel = gam_model_theta_coupling, re.form = NULL, n = 10000)
plot(sim)  

roi_slopes <- emtrends(gam_model_theta_coupling, ~ ROI, var = "Exponent_diff")
roi_slopes_summary <- summary(roi_slopes, infer = c(TRUE, TRUE), level = 0.99)
roi_slopes_df <- as.data.frame(roi_slopes_summary)
roi_slopes_df$p.adjusted <- p.adjust(roi_slopes_df$p.value, method = "holm")
significant_rois <- roi_slopes_df %>%
  filter(p.adjusted < 0.05)
significant_rois

posthoc_df <- roi_slopes_df %>%
  rename(estimate = Exponent_diff.trend) %>%
  mutate(
    signif = case_when(
      p.adjusted < 0.001 ~ "***",
      p.adjusted < 0.01  ~ "**",
      p.adjusted < 0.05  ~ "*",
      TRUE ~ " "
    )
  )

roi_order <- c("O2", "O1", "P8", "P4", "Pz", "P3", "P7", "T4", "T3", "C4", "Cz", "C3", "F8", "F4", "Fz", "F3", "F7", "Fp2", "Fp1")

posthoc_df <- posthoc_df %>%
  mutate(ROI = factor(ROI, levels = roi_order))

p_coupling_forest <- ggplot(posthoc_df, aes(x = estimate, y = ROI)) +
  geom_point(size = 3, aes(color = signif)) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.5) +
  geom_text(aes(x = upper.CL, label = signif), hjust = -0.3, size = 4) +
  scale_color_manual(values = c("*" = "green", "**" = "green", "***" = "green")) +
  labs(x = "Aperiodic-Oscillatory β", y = "Site", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none") 
p_coupling_forest <- p_coupling_forest + coord_cartesian(xlim = c(-0.50, 0.50))

plot_emm <- emtrends(gam_model_theta_coupling, ~ ROI, var = "Exponent_diff")
emm_df <- as.data.frame(plot_emm)

df_topo <- emm_df %>%
  rename(electrode = ROI, quantity = Exponent_diff.trend) %>%
  select(electrode, quantity) %>%
  mutate(
    electrode = as.factor(electrode),
  )

p_coupling <- topoplot(df_topo, 
                       quantity = "quantity", 
                       chan_marker  = "name",
                       method = "Biharmonic",
                       head = TRUE,
                       palette = "YlOrRd",
                       fill_title="β")
p_coupling <- p_coupling +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.50,0.50),
    oob = scales::squish,
    name = "β"
  )
p_coupling <- p_coupling + labs(title = "Exponent-Theta") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p_coupling 

p_coupling2 <- topoplot(df_topo, 
                        quantity = "quantity", 
                        chan_marker  = "name",
                        method = "Biharmonic",
                        head = TRUE,
                        palette = "YlOrRd",
                        fill_title="β")
p_coupling2 <- p_coupling2 +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.50,0.50),
    oob = scales::squish,
    name = "β"
  )
p_coupling2 <- p_coupling2 + labs(title = "Exponent-Theta") +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) 
p_coupling2 <- p_coupling2 + theme(legend.position = "none")

combined <- p_coupling_forest + p_coupling + p_coupling2 +
  plot_layout(widths = c(0.75, 1, 1)) 
combined

gam_model_beta_coupling <- gam(
  Beta_diff ~ 
    Exponent_diff * ROI + Age + Gender + Cortisol_Baseline +
    s(Subject, bs = "re"),    
  data = df_diff,
  family = scat,
  method = "REML"
)
summary(gam_model_beta_coupling)
anova(gam_model_beta_coupling)
gam.check(gam_model_beta_coupling)

print(gam_model_beta_coupling$df.residual)
skewness(resid(gam_model_beta_coupling))
kurtosis(resid(gam_model_beta_coupling))

sim <- simulateResiduals(fittedModel = gam_model_beta_coupling, re.form = NULL, n = 10000)
plot(sim)  

roi_slopes <- emtrends(gam_model_beta_coupling, ~ ROI, var = "Exponent_diff")
roi_slopes_summary <- summary(roi_slopes, infer = c(TRUE, TRUE), level = 0.99)
roi_slopes_df <- as.data.frame(roi_slopes_summary)
roi_slopes_df$p.adjusted <- p.adjust(roi_slopes_df$p.value, method = "holm")
significant_rois <- roi_slopes_df %>%
  filter(p.adjusted < 0.05)
significant_rois

posthoc_df <- roi_slopes_df %>%
  rename(estimate = Exponent_diff.trend) %>%  
  mutate(
    signif = case_when(
      p.adjusted < 0.001 ~ "***",
      p.adjusted < 0.01  ~ "**",
      p.adjusted < 0.05  ~ "*",
      TRUE ~ " "
    )
  )

roi_order <- c("O2", "O1", "P8", "P4", "Pz", "P3", "P7", "T4", "T3", "C4", "Cz", "C3", "F8", "F4", "Fz", "F3", "F7", "Fp2", "Fp1")

posthoc_df <- posthoc_df %>%
  mutate(ROI = factor(ROI, levels = roi_order))

p_coupling_forest <- ggplot(posthoc_df, aes(x = estimate, y = ROI)) +
  geom_point(size = 3, aes(color = signif)) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.5) +
  geom_text(aes(x = upper.CL, label = signif), hjust = -0.3, size = 4) +
  scale_color_manual(values = c("*" = "green", "**" = "green", "***" = "green")) +
  labs(x = "Aperiodic-Oscillatory β", y = "Site", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none") 
p_coupling_forest <- p_coupling_forest + coord_cartesian(xlim = c(-0.50, 0.50))

plot_emm <- emtrends(gam_model_beta_coupling, ~ ROI, var = "Exponent_diff")
emm_df <- as.data.frame(plot_emm)

df_topo <- emm_df %>%
  rename(electrode = ROI, quantity = Exponent_diff.trend) %>%
  select(electrode, quantity) %>%
  mutate(
    electrode = as.factor(electrode),
  )

p_coupling <- topoplot(df_topo, 
                       quantity = "quantity", 
                       chan_marker  = "name",
                       method = "Biharmonic",
                       head = TRUE,
                       palette = "YlOrRd",
                       fill_title="β")
p_coupling <- p_coupling +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.50,0.50),
    oob = scales::squish,
    name = "β"
  )
p_coupling <- p_coupling + labs(title = "Exponent-Beta") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p_coupling 

p_coupling2 <- topoplot(df_topo, 
                        quantity = "quantity", 
                        chan_marker  = "name",
                        method = "Biharmonic",
                        head = TRUE,
                        palette = "YlOrRd",
                        fill_title="β")
p_coupling2 <- p_coupling2 +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.50,0.50),
    oob = scales::squish,
    name = "β"
  )
p_coupling2 <- p_coupling2 + labs(title = "Exponent-Beta") +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) 
p_coupling2 <- p_coupling2 + theme(legend.position = "none")

combined <- p_coupling_forest + p_coupling + p_coupling2 +
  plot_layout(widths = c(0.75, 1, 1)) 
combined