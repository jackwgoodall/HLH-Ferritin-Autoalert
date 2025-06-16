# HLH Diagnosis Rate Analysis
# Date: 14-06-2025
# Author: Jack Goodall
# Contact: j.goodall@sheffield.ac.uk

# Load packages and set global options
if (!require("pacman", quietly = TRUE)) { install.packages("pacman"); library("pacman") }
pacman::p_load(ggplot2, tidyverse, lmtest, AER)
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
options(scipen = 999, digits = 12)  # Disable scientific notation

# Set global ggplot theme
theme_set(
  theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    )
)

# -------------------------
# Emergency Admissions Data
# -------------------------

# Load and prepare STH ED attendance and admission data
STH_attendance <- read_csv("data/STH_attendance.csv",
                           col_types = cols(Total_Emergency_Admissions = col_number(),
                                            Total_Attendances = col_number(),
                                            Month = col_number(), Year = col_number(),
                                            Date = col_date(format = "%d/%m/%Y")))

# Pivot to long format for plotting
STH_attendance_long <- STH_attendance %>%
  mutate("Total Emergency Admissions" = Total_Emergency_Admissions,
         "Total ED Attendances" = Total_Attendances) %>%
  plotly::select(Date, "Total Emergency Admissions") %>%
  pivot_longer(cols = c("Total Emergency Admissions"),
               names_to = "name",
               values_to = "count")

# Plot Emergency Admissions
STH_attendance_long %>%
  ggplot(aes(x = Date, y = count, colour = name)) +
  geom_segment(aes(x = Date, xend = Date, y = 0, yend = count), color = "black", size = 0.8) +
  geom_point(size = 3) +
  labs(title = "Monthly STH Emergency Admissions", x = "Date", y = "Count", colour = "Category")

# ------------------------------------
# Interrupted Time Series: Main Model
# ------------------------------------

# Create quarter-year variable and summarise emergency admissions
STH_attendance_quarter <- STH_attendance %>%
  mutate(
    quarter_year = case_when(
      Month %in% c(1, 2, 3) ~ paste0(Year, "-Q1"),
      Month %in% c(4, 5, 6) ~ paste0(Year, "-Q2"),
      Month %in% c(7, 8, 9) ~ paste0(Year, "-Q3"),
      Month %in% c(10, 11, 12) ~ paste0(Year, "-Q4")
    )
  ) %>%
  group_by(quarter_year) %>%
  summarise(emergency_admissions = sum(Total_Emergency_Admissions)) %>%
  ungroup()

# Merge with main outcome data
main_analysis <- read.csv("data/main_analysis.csv")
main_analysis <- STH_attendance_quarter %>%
  inner_join(main_analysis, join_by(quarter_year))

# Poisson ITS model
poisson_model_quarter_ITS <- glm(
  count ~ time + intervention + time_from + offset(log(emergency_admissions)),
  family = poisson,
  data = main_analysis
)

# Model diagnostics
summary(poisson_model_quarter_ITS)
dispersiontest(poisson_model_quarter_ITS)
dwtest(poisson_model_quarter_ITS)

# -------------------------------
# Sensitivity: No Respiratory Virus
# -------------------------------

HLH_noresp <- read.csv("data/HLH_noresp.csv")
HLH_noresp <- STH_attendance_quarter %>%
  full_join(HLH_noresp)

poisson_model_quarter_ITS_noresp <- glm(
  count ~ time + intervention + time_from + offset(log(emergency_admissions)),
  family = poisson,
  data = HLH_noresp
)

summary(poisson_model_quarter_ITS_noresp)
dispersiontest(poisson_model_quarter_ITS_noresp)
dwtest(poisson_model_quarter_ITS_noresp)

# ------------------------
# Sensitivity: Only Coded
# ------------------------

HLH_coded <- read.csv("data/HLH_coded.csv")
HLH_coded <- STH_attendance_quarter %>%
  full_join(HLH_coded)

poisson_model_quarter_ITS_coded <- glm(
  count ~ time + intervention + time_from + offset(log(emergency_admissions)),
  family = poisson,
  data = HLH_coded
)

summary(poisson_model_quarter_ITS_coded)
dispersiontest(poisson_model_quarter_ITS_coded)
dwtest(poisson_model_quarter_ITS_coded)

# -------------------
# Visualisation Setup
# -------------------

main_analysis <- main_analysis %>%
  mutate(intervention_2 = case_when(intervention == 0 ~ "Pre",
                                    intervention == 1 ~ "Post"))

# Get predicted counts and confidence intervals
preds <- predict(poisson_model_quarter_ITS, newdata = main_analysis, type = "link", se.fit = TRUE)

main_analysis$poisson_pred <- exp(preds$fit)
main_analysis$poisson_pred_lower <- exp(preds$fit - 1.96 * preds$se.fit)
main_analysis$poisson_pred_upper <- exp(preds$fit + 1.96 * preds$se.fit)

# Counterfactual predictions (no intervention)
main_analysis_no_intervention <- main_analysis
main_analysis_no_intervention$intervention <- 0

preds_no_intervention <- predict(poisson_model_quarter_ITS,
                                 newdata = main_analysis_no_intervention,
                                 type = "link", se.fit = TRUE)

main_analysis$poisson_pred_no_intervention <- exp(preds_no_intervention$fit)

# Plot modelled HLH cases
ggplot(main_analysis, aes(x = quarter_year, y = count)) +
  geom_line(aes(y = poisson_pred_no_intervention, group = 1), linetype = "dashed", size = 1, color = "orange") +
  geom_point(aes(colour = intervention_2)) +
  geom_line(aes(y = poisson_pred, group = intervention_2, color = intervention_2), size = 1) +
  geom_ribbon(aes(ymin = poisson_pred_lower, ymax = poisson_pred_upper, group = intervention_2, fill = intervention_2), alpha = 0.2) +
  geom_vline(xintercept = "2022-Q1", linetype = "dashed", color = "black") +
  labs(x = "Quarter Year", y = "Count of Diagnosed HLH Cases") +
  scale_fill_manual(values = c("cadetblue", "orange")) +
  scale_colour_manual(values = c("cadetblue", "orange")) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 5), labels = scales::label_number(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(linewidth = 0.5),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'solid', colour = "grey"))

# ------------------------------------
# Estimate Missed Cases Without Alert
# ------------------------------------

acutally_diagnosed <- sum(main_analysis$count[main_analysis$intervention_2 == "Post"])
predicted_diagnosed <- sum(main_analysis$poisson_pred_no_intervention[main_analysis$intervention_2 == "Post"])
missed_cases_annualised <- (acutally_diagnosed - predicted_diagnosed) / 21 * 12
missed_cases_annualised

# ----------------------------------
# HLH Cases During COVID-19 Pandemic
# ----------------------------------

print("Median cases before the intervention:")
median(main_analysis$count[main_analysis$intervention == 0])

print("Median cases during 2020:")
median(main_analysis$count[grepl("2020", main_analysis$quarter_year)])

print("Mean cases before the intervention:")
mean(main_analysis$count[main_analysis$intervention == 0])

print("Mean cases during 2020:")
mean(main_analysis$count[grepl("2020", main_analysis$quarter_year)])
