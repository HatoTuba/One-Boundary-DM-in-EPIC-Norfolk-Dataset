# Run logistic regression and survival analyses here
setwd("/srv/projects/2024_10_TH_ENDR021_2024/analysis/r_analyses")
library(dplyr); library(ggplot2); library(ggsurvfit)
library(tidyr); library(forcats); library(here)
library(ggeffects); library(sjPlot); library(ggcorrplot)
library(caret); library(brms); library(pscl)
library(yardstick); library(pROC)
# require(lme4); require(survminer) # not existent
require(gee); require(broom); library(rsample)
require(survival); require(performance); require(rms); require(survRM2) 
source("Functions.R")

####-----HC3 Analysis-----####
# data_original <- read.csv("data/demog_HC3_full.csv")
data <- read.csv("data/demog_HC3_full_py.csv")

# prepare the data for survival analyses
data <- data |> mutate(EOF_Y10TH = case_when(EOF_Y10TH == 0.0 ~ 2023.2,
                                         TRUE ~ EOF_Y10TH)) |>
  mutate(time_to_event = case_when(INC_DEMENTIA_ALL == "No Dementia" & EOF_Y10TH < DOC_HC3 ~ 0,
                          TRUE ~ EOF_Y10TH - DOC_HC3)) |>
  filter(!if_any(c(SHORT_EMSE, meanRT_comp, EDLEVEL09), ~is.na(.)))

# summarize time to event data of dementia diagnosis 
data |> filter(INC_DEMENTIA_ALL == "Dementia") |> 
  summarise(tm      = mean(time_to_event), 
            tmedian = median(time_to_event), 
            mintm   = min(time_to_event),
            maxtm   = max(time_to_event),
            stm     = sd(time_to_event))


# uses function to prepare the data for survival analysis
data_surv <- surv_prepare(data)

# standardize the numeric data
data_st <- data_surv |> mutate(
  across(
    .cols = where(is.numeric) & !c(event, time_to_event),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )
)



#####-----STANDARD SURVIVAL ANALYSIS FOR COMPLEX VST-----#####
# data <- read.csv("data/demog_HC3.csv")
data <- data_st |> select(-X, -X_z) # standardized data
### survival analysis with complex task ###
# create surv object
surv_obj <- with(data, 
                 Surv(time=time_to_event, 
                             event=event))

# cox model with mean RT 
cox_mod_meanrt <- coxph(
  surv_obj ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    meanRT_comp_z + DOC_AGE_HC3_z:meanRT_comp_z,
  data=data
)

# cox model with RT quantiles
cox_mod_qRT <- coxph(
   surv_obj ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    rt_q1_z + rt_q2_z + rt_q3_z + DOC_AGE_HC3_z:rt_q1_z + DOC_AGE_HC3_z:rt_q2_z  + DOC_AGE_HC3_z:rt_q3_z,
    data=data
)

# cox model with cognitive parameters
cox_mod_cogpar <- coxph(
  surv_obj ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    dv_comp_z + a_comp_z + t0_comp_z +
    DOC_AGE_HC3_z:dv_comp_z+ DOC_AGE_HC3_z:a_comp_z + DOC_AGE_HC3_z:t0_comp_z, 
  data=data
)

# prepare tables
tsurv_qRT <- tidy(cox_mod_qRT, 
               exponentiate = TRUE,
              conf.int=T, 
              conf.level=0.95) |> 
  mutate(
    HR = round(estimate, 2),
    CI= paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)

tsurv_cogpar <- tidy(cox_mod_cogpar, 
              exponentiate = T,
              conf.int=T, 
              conf.level=0.95
              ) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)

tsurv_meanrt <- tidy(cox_mod_meanrt, 
               exponentiate = T,
               conf.int=T, 
               conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)


comparison_surv <- full_join(
  tsurv_meanrt |> rename_with(~paste0(., "_mmeanRT"), -term),
  tsurv_qRT |> rename_with(~paste0(., "_mqRT"), -term),
  join_by("term")) |> 
  full_join(tsurv_cogpar |> rename_with(~paste0(., "_mcog"), -term), join_by("term"))

metrics_surv <- tibble(
  term = c("AIC", "BIC", "Concordance"),
  HR_mmeanRT = c(AIC(cox_mod_meanrt), BIC(cox_mod_meanrt), summary(cox_mod_meanrt)$concordance[1]),
  CI_mmeanRT = NA_character_,
  p_mmeanRT = NA_real_,
  HR_mqRT = c(AIC(cox_mod_qRT), BIC(cox_mod_qRT), summary(cox_mod_qRT)$concordance[1]),
  CI_mqRT = NA_character_,
  p_mqRT = NA_real_,
  HR_mcog = c(AIC(cox_mod_cogpar), BIC(cox_mod_cogpar), summary(cox_mod_cogpar)$concordance[1]),
  CI_mcog = NA_character_,
  p_mcog = NA_real_,
)


comparison_tab_surv <- bind_rows(comparison_surv, metrics_surv)

write.table(comparison_tab_surv, 'data/results/results_survival_complex.txt', sep=",")

#######----SIMPLE SURVIVAL----######
data <- data_st |> select(-X, -X_z)
# create surv object
surv_obj_simp <- with(data, 
                 Surv(time=time_to_event, 
                      event=event))
# cox model with mean RT 
cox_mod_meanrt_simple <-coxph(
  surv_obj_simp ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    meanRT_simple_HC3_z + DOC_AGE_HC3_z:meanRT_simple_HC3_z,
  data=data
)

# cox model with quantiles
cox_mod_qRT_simple <- coxph(
  surv_obj_simp ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    rt_q1_simp_z + rt_q2_simp_z + rt_q3_simp_z + 
    DOC_AGE_HC3_z:rt_q1_simp_z + DOC_AGE_HC3_z:rt_q2_simp_z  + 
    DOC_AGE_HC3_z:rt_q3_simp_z,
  data=data
)

# cox model with parameters
cox_mod_cogpar_simple <- coxph(
  surv_obj_simp ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    v_simple_z + a_simple_z + t0_simple_z +
    DOC_AGE_HC3_z:v_simple_z + DOC_AGE_HC3_z:a_simple_z + DOC_AGE_HC3_z:t0_simple_z, 
  data=data
)

## prepare tables
tsurv_meanrt_simple <- tidy(cox_mod_meanrt_simple, 
               exponentiate = TRUE,
               conf.int=T, 
               conf.level=0.95) |> 
  mutate(
    HR = round(estimate, 2),
    CI= paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)

tsurv_qRT_simple <- tidy(cox_mod_qRT_simple, 
               exponentiate = T,
               conf.int=T, 
               conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)

tsurv_cogpar_simple <- tidy(cox_mod_cogpar_simple, 
               exponentiate = T,
               conf.int=T, 
               conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)


comparison_surv_simple <- full_join(
  tsurv_meanrt_simple |> rename_with(~paste0(., "_meanrt"), -term),
  tsurv_qRT_simple |> rename_with(~paste0(., "_mqRT"), -term),
  join_by("term")) |> 
  full_join(tsurv_cogpar_simple |> rename_with(~paste0(., "_mcogpar"), -term), join_by("term"))

metrics_surv_simple <- tibble(
  term = c("AIC", "BIC", "Concordance"),
  HR_meanrt = c(AIC(cox_mod_meanrt_simple), BIC(cox_mod_meanrt_simple), summary(cox_mod_meanrt_simple)$concordance[1]),
  CI_meanrt = NA_character_,
  `p value_meanrt` = NA_real_,
  HR_mqRT = c(AIC(cox_mod_qRT_simple), BIC(cox_mod_qRT_simple), summary(cox_mod_qRT_simple)$concordance[1]),
  CI_mqRT = NA_character_,
  `p value_mqRT` = NA_real_,
  HR_mcogpar = c(AIC(cox_mod_cogpar_simple), BIC(cox_mod_cogpar_simple), summary(cox_mod_cogpar_simple)$concordance[1]),
  CI_mcogpar = NA_character_,
  `p value_mcogpar` = NA_real_,
)

comparison_tab_surv_simple <- bind_rows(comparison_surv_simple, metrics_surv_simple)
write.table(comparison_tab_surv_simple, "data/results/results_survival_simple.txt", sep=",")



