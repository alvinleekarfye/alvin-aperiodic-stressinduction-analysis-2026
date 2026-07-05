rm(list=ls())

library(tidyverse)
library(mgcv)
library(emmeans)

mist_df <- read.csv("C:/MIST.csv")

round(mean(mist_df$overall_acc_1, na.rm = TRUE), 2)
round(mean(mist_df$overall_acc_2, na.rm = TRUE), 2)
round(mean(mist_df$overall_acc_3, na.rm = TRUE), 2)
round(mean(mist_df$overall_acc_4, na.rm = TRUE), 2)

round(mean(mist_df$overall_RT_1, na.rm = TRUE), 2)
round(mean(mist_df$overall_RT_2, na.rm = TRUE), 2)
round(mean(mist_df$overall_RT_3, na.rm = TRUE), 2)
round(mean(mist_df$overall_RT_4, na.rm = TRUE), 2)

round(sd(mist_df$overall_acc_1, na.rm = TRUE), 2)
round(sd(mist_df$overall_acc_2, na.rm = TRUE), 2)
round(sd(mist_df$overall_acc_3, na.rm = TRUE), 2)
round(sd(mist_df$overall_acc_4, na.rm = TRUE), 2)

round(sd(mist_df$overall_RT_1, na.rm = TRUE), 2)
round(sd(mist_df$overall_RT_2, na.rm = TRUE), 2)
round(sd(mist_df$overall_RT_3, na.rm = TRUE), 2)
round(sd(mist_df$overall_RT_4, na.rm = TRUE), 2)

df_acc_long <- mist_df %>%
  select(id, overall_acc_1, overall_acc_2, overall_acc_3, overall_acc_4) %>%
  pivot_longer(cols = starts_with("overall_acc"), names_to = "time", values_to = "ACC")

df_acc_long <- df_acc_long %>%
  mutate(time = case_when(
    time == "overall_acc_1" ~ 5,
    time == "overall_acc_2" ~ 10,
    time == "overall_acc_3" ~ 15,
    time == "overall_acc_4" ~ 20
  ))

df_acc_long$id <- factor(df_acc_long$id)
df_acc_long$time <- factor(as.character(df_acc_long$time), levels = c("5", "10", "15", "20"))

model1 <- gam(ACC ~ time + s(id, bs='re'), data = df_acc_long, method="REML")
results <- summary(model1)
print(results$pTerms.table)
anova(model1)
print(model1$df.residual)

posthoc <- emmeans(model1, pairwise~time, at = list(time = c(5, 10, 15, 20)), data=df_acc_long, adjust = "holm")
summary(posthoc$contrasts)
contrasts_df <- as.data.frame(posthoc$contrasts)
contrasts_df

df_RT_long <- mist_df %>%
  select(id, overall_RT_1, overall_RT_2, overall_RT_3, overall_RT_4) %>%
  pivot_longer(cols = starts_with("overall_RT"), names_to = "time", values_to = "RT")

df_RT_long <- df_RT_long %>%
  mutate(time = case_when(
    time == "overall_RT_1" ~ 5,
    time == "overall_RT_2" ~ 10,
    time == "overall_RT_3" ~ 15,
    time == "overall_RT_4" ~ 20
  ))

df_RT_long$id <- factor(df_RT_long$id)
df_RT_long$time <- factor(as.character(df_RT_long$time), levels = c("5", "10", "15", "20"))

model2 <- gam(RT ~ time + s(id, bs='re'), data = df_RT_long, family=inverse.gaussian(link="inverse"), method="REML")
results<- summary(model2)
print(results$pTerms.table)
anova(model2)
print(model2$df.residual)

posthoc <- emmeans(model2, pairwise~time, at = list(time = c(5, 10, 15, 20)), data=df_RT_long, adjust = "holm")
summary(posthoc$contrasts)
contrasts_df <- as.data.frame(posthoc$contrasts)
contrasts_df
