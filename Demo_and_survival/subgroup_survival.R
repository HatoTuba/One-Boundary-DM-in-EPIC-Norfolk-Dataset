## This script divides the dataset into two groups, one with participants younger than 75,
## one with participants older than 75 (because mean age of dementia group was 75)
####-----Subgroup 75+-----####
# define older group using standardized data from the full survival analysis
older_group <- data_st |> filter(DOC_AGE_HC3 >= 75) |> group_by(INC_DEMENTIA_ALL) 

older_group |> group_by(INC_DEMENTIA_ALL, SEX) |> count()

older_group |> filter(INC_DEMENTIA_ALL == 1) |> 
  summarise(avg_surv= mean(time_to_event), 
            sd_surv = sd(time_to_event), 
            med_surv = median(time_to_event))

#####-----STANDARD SURVIVAL ANALYSIS subgroup 75+-----#####
# data <- read.csv("data/demog_HC3.csv")
data_old <- older_group |> select(-X, -X_z) # standardized data

### survival analysis with complex task ###
# create surv object
surv_obj_old <- with(data_old, 
                 Surv(time=time_to_event, 
                      event=event))

# cox model with mean RT 
cox_mod_meanrt_old <- coxph(
  surv_obj_old ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    meanRT_comp_z + DOC_AGE_HC3_z:meanRT_comp_z,
  data=data_old
)

cox_mod_qRT_old <- coxph(
  surv_obj_old ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    rt_q1_z + rt_q2_z + rt_q3_z + DOC_AGE_HC3_z:rt_q1_z + 
    DOC_AGE_HC3_z:rt_q2_z  + DOC_AGE_HC3_z:rt_q3_z,
  data=data_old
)

# cox model with parameters
cox_mod_cogpar_old <- coxph(
  surv_obj_old ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    dv_comp_z + a_comp_z + t0_comp_z +
    DOC_AGE_HC3_z:dv_comp_z+ DOC_AGE_HC3_z:a_comp_z + DOC_AGE_HC3_z:t0_comp_z, 
  data=data_old
)

# prepare tables
tsurv_meanrt_old <- tidy(cox_mod_meanrt_old, 
                         exponentiate = T,
                         conf.int=T, 
                         conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)


tsurv_qRT_old <- tidy(cox_mod_qRT_old, 
                  exponentiate = TRUE,
                  conf.int=T, 
                  conf.level=0.95) |> 
  mutate(
    HR = round(estimate, 2),
    CI= paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)

tsurv_cogpar_old <- tidy(cox_mod_cogpar_old, 
                     exponentiate = T,
                     conf.int=T, 
                     conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)


comparison_surv_old <- full_join(
  tsurv_meanrt_old |> rename_with(~paste0(., "_mmeanRT"), -term),
  tsurv_qRT_old |> rename_with(~paste0(., "_mqRT"), -term),
  join_by("term")) |> 
  full_join(tsurv_cogpar_old |> rename_with(~paste0(., "_mcog"), -term), join_by("term"))

metrics_surv_old <- tibble(
  term = c("AIC", "BIC", "Concordance"),
  HR_mmeanRT = c(AIC(cox_mod_meanrt_old), BIC(cox_mod_meanrt_old), summary(cox_mod_meanrt_old)$concordance[1]),
  CI_mmeanRT = NA_character_,
  p_mmeanRT = NA_real_,
  HR_mqRT = c(AIC(cox_mod_qRT_old), BIC(cox_mod_qRT_old), summary(cox_mod_qRT_old)$concordance[1]),
  CI_mqRT = NA_character_,
  p_mqRT = NA_real_,
  HR_mcog = c(AIC(cox_mod_cogpar_old), BIC(cox_mod_cogpar_old), summary(cox_mod_cogpar_old)$concordance[1]),
  CI_mcog = NA_character_,
  p_mcog = NA_real_,
)

comparison_tab_surv_old <- bind_rows(comparison_surv_old, metrics_surv_old)

write.table(comparison_tab_surv_old, 'data/results/results_survival_complex_old.txt', sep=",")

####-----Subgroup <75-----####
young <- data_st |> filter(DOC_AGE_HC3 < 75) 

young |> group_by(INC_DEMENTIA_ALL) |> count()
young |> filter(INC_DEMENTIA_ALL == 1) |> 
  summarise(avg_surv= mean(time_to_event), 
            sd_surv = sd(time_to_event), 
            med_surv = median(time_to_event))

data_young <- young |> select(-X, -X_z) # standardized data; clean up unnecessary cols

### survival analysis with complex task ###
# create surv object
surv_obj_young <- with(data_young, 
                     Surv(time=time_to_event, 
                          event=event))

# cox model with mean RT 
cox_mod_meanrt_young <- coxph(
  surv_obj_young ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    meanRT_comp_z + DOC_AGE_HC3_z:meanRT_comp_z,
  data=data_young
)

cox_mod_qRT_young <- coxph(
  surv_obj_young ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    rt_q1_z + rt_q2_z + rt_q3_z + DOC_AGE_HC3_z:rt_q1_z + 
    DOC_AGE_HC3_z:rt_q2_z  + DOC_AGE_HC3_z:rt_q3_z,
  data=data_young
)

# cox model with parameters
cox_mod_cogpar_young <- coxph(
  surv_obj_young ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    dv_comp_z + a_comp_z + t0_comp_z +
    DOC_AGE_HC3_z:dv_comp_z+ DOC_AGE_HC3_z:a_comp_z + DOC_AGE_HC3_z:t0_comp_z, 
  data=data_young
)


# prepare tables
tsurv_meanrt_young <- tidy(cox_mod_meanrt_young, 
                         exponentiate = T,
                         conf.int=T, 
                         conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)


tsurv_qRT_young <- tidy(cox_mod_qRT_young, 
                      exponentiate = TRUE,
                      conf.int=T, 
                      conf.level=0.95) |> 
  mutate(
    HR = round(estimate, 2),
    CI= paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)

tsurv_cogpar_young <- tidy(cox_mod_cogpar_young, 
                         exponentiate = T,
                         conf.int=T, 
                         conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)



comparison_surv_young <- full_join(
  tsurv_meanrt_young |> rename_with(~paste0(., "_mmeanRT"), -term),
  tsurv_qRT_young |> rename_with(~paste0(., "_mqRT"), -term),
  join_by("term")) |> 
  full_join(tsurv_cogpar_young |> rename_with(~paste0(., "_mcog"), -term), join_by("term"))

metrics_surv_young <- tibble(
  term = c("AIC", "BIC", "Concordance"),
  HR_mmeanRT = c(AIC(cox_mod_meanrt_young), BIC(cox_mod_meanrt_young), summary(cox_mod_meanrt_young)$concordance[1]),
  CI_mmeanRT = NA_character_,
  p_mmeanRT = NA_real_,
  HR_mqRT = c(AIC(cox_mod_qRT_young), BIC(cox_mod_qRT_young), summary(cox_mod_qRT_young)$concordance[1]),
  CI_mqRT = NA_character_,
  p_mqRT = NA_real_,
  HR_mcog = c(AIC(cox_mod_cogpar_young), BIC(cox_mod_cogpar_young), summary(cox_mod_cogpar_young)$concordance[1]),
  CI_mcog = NA_character_,
  p_mcog = NA_real_,
)

comparison_tab_surv_young <- bind_rows(comparison_surv_young, metrics_surv_young)

write.table(comparison_tab_surv_young, 'data/results/results_survival_complex_young.txt', sep=",")


####################################################################
####-----Subgroup Analysis Simple 75+-----####
# data <- data_st |> select(-X, -X_z)
# create surv object
surv_obj_simp_old <- with(data_old, 
                      Surv(time=time_to_event, 
                           event=event))
# cox model with mean RT 
cox_mod_meanrt_simple_old <-coxph(
  surv_obj_simp_old ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    meanRT_simple_HC3_z + DOC_AGE_HC3_z:meanRT_simple_HC3_z,
  data=data_old
)

# cox model with quantiles
cox_mod_qRT_simple_old <- coxph(
  surv_obj_simp_old ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    rt_q1_simp_z + rt_q2_simp_z + rt_q3_simp_z + 
    DOC_AGE_HC3_z:rt_q1_simp_z + DOC_AGE_HC3_z:rt_q2_simp_z  + 
    DOC_AGE_HC3_z:rt_q3_simp_z,
  data=data_old
)

# cox model with parameters
cox_mod_cogpar_simple_old <- coxph(
  surv_obj_simp_old ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    v_simple_z + a_simple_z + t0_simple_z +
    DOC_AGE_HC3_z:v_simple_z + DOC_AGE_HC3_z:a_simple_z + 
    DOC_AGE_HC3_z:t0_simple_z, 
  data=data_old
)


tsurv_meanrt_simple_old <- tidy(cox_mod_meanrt_simple_old, 
                            exponentiate = TRUE,
                            conf.int=T, 
                            conf.level=0.95) |> 
  mutate(
    HR = round(estimate, 2),
    CI= paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)

tsurv_qRT_simple_old <- tidy(cox_mod_qRT_simple_old, 
                         exponentiate = T,
                         conf.int=T, 
                         conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)

tsurv_cogpar_simple_old <- tidy(cox_mod_cogpar_simple_old, 
                            exponentiate = T,
                            conf.int=T, 
                            conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)


comparison_surv_simple_old <- full_join(
  tsurv_meanrt_simple_old |> rename_with(~paste0(., "_meanrt"), -term),
  tsurv_qRT_simple_old |> rename_with(~paste0(., "_mqRT"), -term),
  join_by("term")) |> 
  full_join(tsurv_cogpar_simple_old |> rename_with(~paste0(., "_mcogpar"), -term), join_by("term"))

metrics_surv_simple_old <- tibble(
  term = c("AIC", "BIC", "Concordance"),
  HR_meanrt = c(AIC(cox_mod_meanrt_simple_old), BIC(cox_mod_meanrt_simple_old), 
                summary(cox_mod_meanrt_simple_old)$concordance[1]),
  CI_meanrt = NA_character_,
  `p value_meanrt` = NA_real_,
  HR_mqRT = c(AIC(cox_mod_qRT_simple_old), BIC(cox_mod_qRT_simple_old), 
              summary(cox_mod_qRT_simple_old)$concordance[1]),
  CI_mqRT = NA_character_,
  `p value_mqRT` = NA_real_,
  HR_mcogpar = c(AIC(cox_mod_cogpar_simple_old), BIC(cox_mod_cogpar_simple_old), 
                 summary(cox_mod_cogpar_simple_old)$concordance[1]),
  CI_mcogpar = NA_character_,
  `p value_mcogpar` = NA_real_,
)
# 
# cox.zph(cox_mod1_simple)
# cox.zph(cox_mod2_simple)
# cox.zph(cox_mod_meanrt_simple)

comparison_tab_surv_simple_old <- bind_rows(comparison_surv_simple, metrics_surv_simple_old)
write.table(comparison_tab_surv_simple_old, 
            "data/results/results_survival_simple_subgroupold.txt", sep=",")
View(comparison_tab_surv_simple_old)

####-----Subgroup Analysis Simple 75------####
# create surv object
surv_obj_simp_young <- with(data_young, 
                      Surv(time=time_to_event, 
                           event=event))
# cox model with mean RT 
cox_mod_meanrt_simple_young <-coxph(
  surv_obj_simp_young ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    meanRT_simple_HC3_z + DOC_AGE_HC3_z:meanRT_simple_HC3_z,
  data=data_young
)

# cox model with quantiles
cox_mod_qRT_simple_young <- coxph(
  surv_obj_simp_young ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    rt_q1_simp_z + rt_q2_simp_z + rt_q3_simp_z + 
    DOC_AGE_HC3_z:rt_q1_simp_z + DOC_AGE_HC3_z:rt_q2_simp_z  + 
    DOC_AGE_HC3_z:rt_q3_simp_z,
  data=data_young
)

# cox model with parameters
cox_mod_cogpar_simple_young <- coxph(
  surv_obj_simp_young ~ DOC_AGE_HC3_z + SEX + EDLEVEL09 + SHORT_EMSE_z +
    v_simple_z + a_simple_z + t0_simple_z +
    DOC_AGE_HC3_z:v_simple_z + DOC_AGE_HC3_z:a_simple_z + DOC_AGE_HC3_z:t0_simple_z, 
  data=data_young
)


tsurv_meanrt_simple_young <- tidy(cox_mod_meanrt_simple_young, 
                            exponentiate = TRUE,
                            conf.int=T, 
                            conf.level=0.95) |> 
  mutate(
    HR = round(estimate, 2),
    CI= paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)

tsurv_qRT_simple_young <- tidy(cox_mod_qRT_simple_young, 
                         exponentiate = T,
                         conf.int=T, 
                         conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)

tsurv_cogpar_simple_young <- tidy(cox_mod_cogpar_simple_young, 
                            exponentiate = T,
                            conf.int=T, 
                            conf.level=0.95
) |> 
  mutate(
    HR = round(estimate, 2),
    CI = paste0(round(conf.low, 2),"-", round(conf.high, 2)),
    `p value` = round(p.value, 3)) |>
  select(term, HR, CI, `p value`)


comparison_surv_simple_young <- full_join(
  tsurv_meanrt_simple_young |> rename_with(~paste0(., "_meanrt"), -term),
  tsurv_qRT_simple_young |> rename_with(~paste0(., "_mqRT"), -term),
  join_by("term")) |> 
  full_join(tsurv_cogpar_simple_young |> rename_with(~paste0(., "_mcogpar"), -term), join_by("term"))

metrics_surv_simple_young <- tibble(
  term = c("AIC", "BIC", "Concordance"),
  HR_meanrt = c(AIC(cox_mod_meanrt_simple_young), BIC(cox_mod_meanrt_simple_young), 
                summary(cox_mod_meanrt_simple_young)$concordance[1]),
  CI_meanrt = NA_character_,
  `p value_meanrt` = NA_real_,
  HR_mqRT = c(AIC(cox_mod_qRT_simple_young), BIC(cox_mod_qRT_simple_young), 
              summary(cox_mod_qRT_simple_young)$concordance[1]),
  CI_mqRT = NA_character_,
  `p value_mqRT` = NA_real_,
  HR_mcogpar = c(AIC(cox_mod_cogpar_simple_young), BIC(cox_mod_cogpar_simple_young), 
                 summary(cox_mod_cogpar_simple_young)$concordance[1]),
  CI_mcogpar = NA_character_,
  `p value_mcogpar` = NA_real_,
)

comparison_tab_surv_simple_young <- bind_rows(comparison_surv_simple_young, metrics_surv_simple_young)
View(comparison_tab_surv_simple_young)
write.table(comparison_tab_surv_simple_young, 
            "data/results/results_survival_simple_subgroupyoung.txt", sep=",")
