library(mgcv)
library(emmeans)
library(eegUtils)
library(ggplot2)
library(dplyr)
library(patchwork)
library(e1071)
library(DHARMa)

df_combined <- readRDS("C:/df_combined_exponent.rds")
df_combined$Subject <- as.factor(df_combined$Subject)
df_combined$Gender <- as.factor(df_combined$Gender)
df_combined$ROI <- as.factor(df_combined$ROI)
df_combined$Time <- as.factor(df_combined$Time)

gam_model_R2 <- gam(
  R2 ~ 
    ROI * Time + 
    s(Subject, bs = "re"),    
  data = df_combined,
  family = scat,
  method = "REML"
)
summary(gam_model_R2)
anova(gam_model_R2)
gam.check(gam_model_R2)

print(gam_model_R2$df.residual)
skewness(resid(gam_model_R2))
kurtosis(resid(gam_model_R2))

sim <- simulateResiduals(fittedModel = gam_model_R2, re.form = NULL, n = 10000)
plot(sim)  

posthoc <- emmeans(gam_model_R2, revpairwise ~ Time | ROI, adjust = "none") 
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

roi_order <- c("O2", "O1", "P8", "P4", "Pz", "P3", "P7", "T4", "T3", "C4", "Cz", "C3", "F8", "F4", "Fz", "F3", "F7", "Fp2", "Fp1")

posthoc_df <- posthoc_df %>%
  mutate(ROI = factor(ROI, levels = roi_order))

p19 <- ggplot(posthoc_df, aes(x = estimate, y = ROI)) +
  geom_point(size = 3, aes(color = signif)) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.5) +
  geom_text(aes(x = upper.CL, label = signif), hjust = -0.3, size = 4) +
  scale_color_manual(values = c("*" = "green", "**" = "green", "***" = "green")) +
  labs(x = "Estimated Difference (Stress - Training)", y = "Site", color = "Significance") +
  theme_minimal() +
  theme(legend.position = "none") 
p19 <- p19 + coord_cartesian(xlim = c(-0.2, 0.2))

plot_emm <- emmeans(gam_model_R2, ~ Time | ROI)
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
p20 <- topoplot(df_training, 
               quantity = "quantity", 
               chan_marker  = "name",
               method = "Biharmonic",
               head = TRUE,
               palette = "YlOrRd",
               fill_title="R2")
p20 <- p20 +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0.90,
    limits = c(0.985,1),
    oob = scales::squish,
    name = "β"
  )
p20 <- p20 + theme(legend.position = "none")
p20 <- p20 + labs(title = "Training") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p20

df_stress <- subset(df_topo, segment == "stress")
p21 <- topoplot(df_stress, 
               quantity = "quantity", 
               chan_marker  = "name",
               method = "Biharmonic",
               head = TRUE,
               palette = "YlOrRd",
               fill_title="R2")
p21 <- p21 +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0.90,
    limits = c(0.985,1),
    oob = scales::squish,
    name = "β"
  )
p21 <- p21 + labs(title = "Stress") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
p21

combined <- p19 + p20 + p21 +
  plot_layout(widths = c(0.75, 1, 1)) 
combined