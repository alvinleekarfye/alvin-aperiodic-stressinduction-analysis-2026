library(dplyr)
library(tidyr)
library(ggplot2)

folders <- c(
  baseline = "D:/aperiod/export_eeg_psd/baseline",
  stress   = "D:/aperiod/export_eeg_psd/stress",
  training = "D:/aperiod/export_eeg_psd/training"
)

selected_participant = 9999

read_psd <- function(folder_path, condition_name) {
  
  file <- list.files(
    folder_path,
    pattern = paste0("^", selected_participant, "_T3.*\\.csv$"),
    full.names = TRUE
  )
  
  df <- read.csv(file[1])
  colnames(df) <- trimws(colnames(df))
  
  freq_col <- grep("freq", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  psd_col  <- grep("psd", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  
  df %>%
    rename(
      Frequency_Hz = all_of(freq_col),
      PSD = all_of(psd_col)
    ) %>%
    mutate(condition = condition_name)
}

psd_data <- bind_rows(
  lapply(names(folders), function(nm) {
    read_psd(folders[[nm]], nm)
  })
) %>%
  filter(Frequency_Hz >= 1, Frequency_Hz <= 30)


df_aperiod <- readRDS("C:/df_aperiod.rds")

df_aperiod <- df_aperiod %>%
  filter(Subject == selected_participant & ROI == "T3")

df_aperiod_cond <- df_aperiod %>%
  group_by(segment) %>%
  summarise(
    Offset = mean(Offset, na.rm = TRUE),
    Exponent = mean(Exponent, na.rm = TRUE),
    R2 = mean(R2, na.rm = TRUE),
    .groups = "drop"
  )

freq <- seq(1, 30, length.out = 200)

df_fit <- df_aperiod_cond %>%
  crossing(Frequency_Hz = freq) %>%
  mutate(
    log10_f = log10(Frequency_Hz),
    log_PSD = Offset - Exponent * log10_f,
    PSD = 10^log_PSD
  )

r2_lookup <- df_fit %>%
  dplyr::group_by(segment) %>%
  dplyr::summarise(R2 = mean(R2, na.rm = TRUE)) %>%
  dplyr::mutate(label = paste0(segment, " (R² = ", round(R2, 3), ")"))
label_vec <- setNames(r2_lookup$label, r2_lookup$segment)

df_fit$segment <- factor(
  df_fit$segment,
  levels = c("baseline", "training", "stress")
)
psd_data$condition <- factor(
  psd_data$condition,
  levels = c("baseline", "training", "stress")
)

ggplot() +
  geom_line(data = psd_data, aes(Frequency_Hz, PSD, color = condition), size = 0.8, alpha = 0.8) +
  geom_line(data = df_fit, aes(Frequency_Hz, PSD, color = segment), linetype = 5, size = 0.8, alpha = 0.8) +
  scale_x_log10(breaks = c(1, 4, 8, 12, 30)) +
  scale_y_log10(limits = c(0.2, 35), breaks = c(0.2, 1, 3, 10, 30)) +
  labs(x = "Frequency (Hz)", y = "μV²/Hz") +
  geom_vline(xintercept = c(1, 4, 8, 12, 30), linetype = "dotted", alpha = 0.8) + 
  annotate("text", x = sqrt(1*4),   y = 0.2, label = "Delta", size = 3, hjust = 0.5) +
  annotate("text", x = sqrt(4*8),   y = 0.2, label = "Theta", size = 3, hjust = 0.5) +
  annotate("text", x = sqrt(8*12),  y = 0.2, label = "Alpha", size = 3, hjust = 0.5) +
  annotate("text", x = sqrt(12*30), y = 0.2, label = "Beta",  size = 3, hjust = 0.5) +
  scale_color_manual(
    values = c(
      baseline = "#1f77b4",
      stress   = "#d62728",
      training = "#2ca02c"),
    labels = c(label_vec)
  ) +
  labs(color = "Conditions") +
  theme_minimal() +
  theme(legend.position = c(0.8, 0.75),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))
