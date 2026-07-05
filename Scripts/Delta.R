library(mgcv)
library(emmeans)
library(eegUtils)
library(ggplot2)
library(dplyr)
library(patchwork)
library(e1071)
library(DHARMa)

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
summary(gam_model_delta)
anova(gam_model_delta)
gam.check(gam_model_delta)

print(gam_model_delta$df.residual)
skewness(resid(gam_model_delta))
kurtosis(resid(gam_model_delta))

sim <- simulateResiduals(fittedModel = gam_model_delta, re.form = NULL, n = 10000)
plot(sim)      

boot_parallel_delta <- readRDS("C:/bootstrap_delta_results.rds")
ci_df <- as.data.frame(boot_parallel_delta$posthoc_ci)
ci_df$label <- rownames(ci_df)
ci_df$ROI <- sub(" \\| .*", "", ci_df$label)
ci_df$contrast <- sub(".*\\| ", "", ci_df$label)
ci_df <- ci_df[, c("ROI", "contrast", "0.5%", "99.5%")]
colnames(ci_df)[3:4] <- c("wild_ci_low", "wild_ci_high")

p_holm_df <- as.data.frame(boot_parallel_delta$posthoc_p_holm)
p_holm_df$label <- rownames(p_holm_df)
p_holm_df$ROI <- sub(" \\| .*", "", p_holm_df$label)
p_holm_df$contrast <- sub(".*\\| ", "", p_holm_df$label)
p_holm_df <- p_holm_df[, c("ROI", "contrast", "boot_parallel_delta$posthoc_p_holm")]
colnames(p_holm_df)[3] <- c("posthoc_p_holm")

wild_df <- merge(
  p_holm_df,
  ci_df,
  by = c("ROI", "contrast"),
  all.x = TRUE
)

posthoc <- emmeans(gam_model_delta, revpairwise ~ Time | ROI, adjust = "none") 
posthoc_df <- as.data.frame(summary(posthoc$contrasts, infer = c(TRUE, TRUE), level = 0.99))
pvals <- posthoc_df$p.value
pvals_adj <- p.adjust(pvals, method = "holm")
posthoc_df$p.value_holm <- pvals_adj
posthoc_sig <- posthoc_df %>%
  filter(p.value_holm < 0.05)
posthoc_sig

posthoc_df <- posthoc_df %>%
  mutate(
    signif = case_when(
      p.value_holm < 0.001 ~ "***",
      p.value_holm < 0.01 ~ "**",
      p.value_holm < 0.05 ~ "*",
      TRUE ~ " "
    )
  )

posthoc_df <- merge(
  posthoc_df,
  wild_df,
  by = c("ROI", "contrast"),
  all.x = TRUE
)

posthoc_df <- posthoc_df %>%
  mutate(
    wildsignif = case_when(
      posthoc_p_holm < 0.001 ~ "***",
      posthoc_p_holm < 0.01 ~ "**",
      posthoc_p_holm < 0.05 ~ "*",
      TRUE ~ " "
    )
  )

roi_order <- c("O2", "O1", "P8", "P4", "Pz", "P3", "P7", "T4", "T3", "C4", "Cz", "C3", "F8", "F4", "Fz", "F3", "F7", "Fp2", "Fp1")

posthoc_df <- posthoc_df %>%
  mutate(ROI = factor(ROI, levels = roi_order))

p7 <- ggplot(posthoc_df, aes(x = estimate, y = ROI)) +
  geom_point(size = 3, aes(color = wildsignif)) +
  geom_errorbarh(aes(xmin = wild_ci_low, xmax = wild_ci_high), height = 0.5) +
  geom_text(aes(x = upper.CL, label = wildsignif), hjust = -0.3, size = 4) +
  scale_color_manual(values = c("*" = "green", "**" = "green", "***" = "green")) +
  labs(x = "Estimated Difference (Stress - Training)", y = "Site", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none") 
p7 <- p7 + coord_cartesian(xlim = c(-0.2, 0.2))

plot_emm <- emmeans(gam_model_delta, ~ Time | ROI)
emm_df <- as.data.frame(plot_emm)
head(emm_df)
df_topo <- emm_df %>%
  rename(segment = Time, electrode = ROI, quantity = emmean) %>%
  select(segment, electrode, quantity) %>%
  mutate(
    electrode = as.factor(electrode),
    segment = recode(segment,
                     "T2" = "training",
                     "T3" = "stress")
  )

df_training <- subset(df_topo, segment == "training")
p8 <- topoplot(df_training, 
               quantity = "quantity", 
               chan_marker  = "name",
               method = "Biharmonic",
               head = TRUE,
               palette = "YlOrRd",
               fill_title="Residual Delta\nBaseline-corrected")
p8 <- p8 +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.1,0.1),
    oob = scales::squish,
    name = "β"
  )
p8 <- p8 + theme(legend.position = "none")
p8 <- p8 + labs(title = "Training") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p8

df_stress <- subset(df_topo, segment == "stress")
p9 <- topoplot(df_stress, 
               quantity = "quantity", 
               chan_marker  = "name",
               method = "Biharmonic",
               head = TRUE,
               palette = "YlOrRd",
               fill_title="Residual Delta\nBaseline-corrected")
p9 <- p9 +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-0.1,0.1),
    oob = scales::squish,
    name = "β"
  )
p9 <- p9 + labs(title = "Stress") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p9

combined <- p7 + p8 + p9 +
  plot_layout(widths = c(0.75, 1, 1)) 
combined