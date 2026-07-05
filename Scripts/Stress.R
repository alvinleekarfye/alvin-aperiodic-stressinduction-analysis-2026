rm(list=ls())

library(tidyverse)
library(mgcv)
library(emmeans)

stress_df_long <- readRDS("C:/stress_df_long.rds")
stress_df_long$ID <- as.factor(stress_df_long$ID)
stress_df_long$Time <- as.factor(stress_df_long$Time)

model1 <- gam(VAS ~ Time + s(ID, bs="re"), data = stress_df_long, method ='REML', family=gaussian(link=identity))
summary <- summary(model1)
summary
plot.gam(model1, seWithMean=TRUE)

anova_results <- anova(model1)
print(anova_results$pTerms.table)
print(model1$df.residual)

emm <- emmeans(model1, ~ Time)
pairwise_comparisons <- contrast(emm, method = "pairwise", adjust = "holm")
pairwise_df <- as.data.frame(pairwise_comparisons)
pairwise_df

df_summary <- stress_df_long %>%
  group_by(Time) %>%
  summarise(
    mean_VAS = mean(VAS, na.rm = TRUE),
    sd_VAS = sd(VAS, na.rm = TRUE),
    n = sum(!is.na(VAS)),
    se_VAS = sd_VAS / sqrt(n),
    ci_lower = mean_VAS - qt(0.995, df = n-1) * se_VAS,
    ci_upper = mean_VAS + qt(0.995, df = n-1) * se_VAS
  )

ggplot(df_summary, aes(x = Time, y = mean_VAS)) +
  geom_col(fill = "grey") +               
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.2) +                
  labs(x = "Timepoint", y = "Perceived Stress") +
  theme_minimal() +
  ylim(c=0,100) +
  scale_x_discrete(labels = c("Baseline", "MIST Training\nCondition",
                              "MIST Stress\nCondition 1", "MIST Stress\nCondition 2", "MIST Stress\nCondition 3")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



