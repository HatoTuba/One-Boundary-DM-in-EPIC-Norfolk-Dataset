setwd("/srv/projects/2024_10_TH_ENDR021_2024/analysis/r_analyses")
require(ggplot2); require(here); library(tidyr); require(dplyr)

# This script illustrates data cleaning and demographic summary preparations

# Read and clean the data, get necessary variables.
# whole data files
desc1 <- read.csv(file_path_desc1)
desc2 <- read.csv(file_path_desc2)

#---------------------
## participants who participated in vst testing
subj_HC3 <- desc1 |> filter(FLAG_VST ==1)

# select the necessary variables for later analysis
demog_HC3 <- subj_HC3 %>% select(RELEASEID, INC_DEMENTIA_ALL, DOC_AGE_HC3, SEX, MARITAL_STATUS3, 
                                 EDLEVEL09, TOWNINDX, SHORT_EMSE, PAL_FTMS, RUPE3, TOTALACT_3HC,
                                 VST_SIMPLE_MTIME, VST_SIMPLE_SDTIME, VST_COMPLEX_MTIME, VST_COMPLEX_SDTIME, 
                                 VST_GLAUC, VST_CACT, VST_MACDEG, VST_RTPATH, VST_GLASSES, HEARING_FU4, # this is 3rd HC
                                 # diseases: prevalent heart attack, CA3:prevalent cancer, prevalent stroke, self reported diabetes
                                 CVA3, CA3, MI3, DM3, DEPRESSION, TOTALACT_3HC, SF36PF_FOLLOW4,
                                 DOC_HC3, INC_DEMENTIA_ALL_DATE_Y10TH, EOF_Y10TH)

#####----OUTLIER ELIMINATION COMPLEX VST----#####

#-----ALL PARTICIPANTS-----
## there are 7171 participants in total
# nrow(dat_complex_hc3) 
# length(unique(dat_complex_hc3$RELEASEID))
# gives the number of participants
complex_HC3_ids <- unique(dat_complex_hc3$RELEASEID)
# eliminate trials with hard limits
all_data_clean <- eliminate_outliers(dat_complex_hc3)
mean(all_data_clean$TIME)


## ELIMINATE BASED ON VALID N_OBS ##
## HC3 ##
# gives out how many valid trials each participant has
total_nobs = all_data_clean %>% summarise(n_obs = n())

total_nobs %>% group_by(n_obs) %>% summarise(count = n())

# these are the people who were excluded from the analysis
see <- total_nobs[total_nobs$n_obs<40,]
write.csv(see, "data/excluded_HC3complex.csv")


# participants who has n_obs < 40 are eliminated
clean_data <- all_data_clean %>% filter(!RELEASEID %in% see$RELEASEID) %>% 
  select(-is_outlier)

clean_data$TIME = clean_data$TIME/1000

hist(clean_data$TIME, breaks=100)
# summary stats for each participants after eliminating the outliers
sum_data_clean <- clean_data %>% summarise(across(TIME, list(
  mean = mean,
  median = median,
  min = min,
  max = max
)))

write.csv(clean_data, "HC3_complex.csv")



####----- OUTLIER DETECTION SIMPLE VST-----####
# these are not exactly quantiles but I loose too many data points otherwise
# Quantile based, does not work well for me
find_outlier1 <- function(x) {
  Q1 <-quantile(x, 0.11)
  Q3 <- quantile(x, 0.99)
  IQR <- Q3 - Q1
  
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  is_outlier = x < lower_bound | x > upper_bound
}


# add another column outlier=TRUE/FALSE
data_with_outliers_HC3 <- dat_simple_HC3 |>
  group_by(RELEASEID) |>
  mutate(is_outlier = find_outlier2_simp(TIME)) # outliers with cutoffs

# to see the number of outlier trials per participant
outlier_sum_HC3 <- data_with_outliers_HC3 |> 
  group_by(RELEASEID) |> 
  summarize(outlier_count = sum(is_outlier))

# summary stats for each participants after eliminating the outliers
sum_data_clean_HC3 <- data_clean_HC3 |> summarise(across(TIME, list(
  mean = mean,
  median = median,
  min = min,
  max = max
)))


#####----Data wrangling----#####

# this treats death as right-censoring meaning standard survival analysis
demog_HC3 <- demog_HC3 |> mutate(event=INC_DEMENTIA_ALL) |>
  mutate(#ifelse(!is.na(INC_DEMENTIA_ALL_DATE_Y10TH), 1, 0), 
    time_to_event = ifelse(event==1, INC_DEMENTIA_ALL_DATE_Y10TH - DOC_HC3,
                           EOF_Y10TH - DOC_HC3))


# combine all visual problems in one column
demog_HC3 <- demog_HC3 |> rowwise() |>
  mutate(visual_problems = {
    vals <- c_across(c(VST_GLAUC, VST_CACT, VST_MACDEG, VST_RTPATH))
    first(vals[vals !=0], default=0)
  }) |> ungroup()

# combine marital statuses
demog_HC3 <- demog_HC3 |> 
  mutate(marital_st = recode(MARITAL_STATUS3, `1`="Single", `2`="Married", `3`="Seperated/divorced/widowed", 
                             `4`="Seperated/divorced/widowed", `5`="Seperated/divorced/widowed"))

# recode categorical variables
demog_HC3 <- demog_HC3 |>  
  mutate(INC_DEMENTIA_ALL = recode(INC_DEMENTIA_ALL, `0`="No Dementia", `1`="Dementia")) |>
  mutate(EDLEVEL09 = recode(EDLEVEL09, `0`="No Ed", `1`="O-A Level", `2`="O-A Level", `3`="Degree")) |>
  mutate(age_group5 = cut(DOC_AGE_HC3, breaks=seq(1, 100, by=5), # add 5 year age groups
                          right=FALSE)) |>
  mutate(SEX = recode(SEX, `1`="Men", `2`="Women")) 

# write everything into csv file
write.csv(demog_HC3, "data/demog_HC3all_large.csv")


####----- COMPLEX VST -----####
# increasing drift
post_sumHC3_simple <- read.csv("data/postsumHC3_simple.csv")
post_sumHC3_inc <- read.csv("data/postsumHC3_complex_inc.csv")

HC3_rt <- read.csv("data/HC3_complex.csv")
HC3_rt_simple <- read.csv("data/HC3_simple.csv")

# mean RTs of clean data
meanRTHC3 <- HC3_rt |> 
  group_by(RELEASEID) |> 
  summarise(meanRT_comp = mean(TIME))

meanRTHC3_simple <- HC3_rt_simple |> 
  group_by(RELEASEID) |> 
  summarise(meanRT_simple_HC3 = mean(TIME))

# pivot wide for parameter values
post_sumHC3_wide <- post_sumHC3_inc |> 
  pivot_wider(id_cols=RELEASEID, names_from = param, values_from=mean) |>
  rename("a_comp"= "a", "dv_comp"="dv", "t0_comp"="t0")
post_sumHC3_simple_wide <- post_sumHC3_simple |> 
  pivot_wider(id_cols = RELEASEID, names_from = param, values_from = mean) |> 
  rename("a_simple"= "a", "v_simple"="v", "t0_simple"="t0")
# combine parameter values and mean RTs 
demog_HC3 <- left_join(demog_HC3, meanRTHC3_simple)
demog_HC3 <- left_join(demog_HC3, meanRTHC3)
# demog_HC3 <- demog_HC3 |> filter(!is.na(meanRT)) # remove outlier participants
demog_HC3 <- left_join(demog_HC3, post_sumHC3_wide)
demog_HC3 <- left_join(demog_HC3, post_sumHC3_simple_wide)
# write.csv(demog_HC3, "data/demog_HC3_full.csv", fileEncoding = "UTF-8")

# Prepare the data for cognitive modeling
demog_HC3_python <- demog_HC3 |> filter(!is.na(meanRT_comp))

demog_HC3_python <- demog_HC3_python |> filter(!is.na(EDLEVEL09))


# Prepare RT quantiles
# complex VST
rt_quant <- HC3_rt |> 
  group_by(RELEASEID) |> 
  summarise(rt_q1 = quantile(TIME, probs = 0.1, na.rm=T),
            rt_q2 = quantile(TIME, probs = 0.5, na.rm=T),
            rt_q3 = quantile(TIME, probs = 0.9, na.rm=T)) 
# simple VST
rt_quant_simple <- HC3_rt_simple |>
  group_by(RELEASEID) |> 
  dplyr::summarize(rt_q1_simp = quantile(TIME, probs = 0.1, na.rm=T),
                   rt_q2_simp = quantile(TIME, probs = 0.5, na.rm=T),
                   rt_q3_simp = quantile(TIME, probs = 0.9, na.rm=T)) 

write.csv(rt_quant_simple, "data/rt_quant_simple.csv")
demog_HC3_python <- left_join(demog_HC3_python, rt_quant)
demog_HC3_python <- left_join(demog_HC3_python, rt_quant_simple)

write.csv(demog_HC3_python, "data/demog_HC3_full_py.csv", fileEncoding = "UTF-8")
demog_HC3_python <- read.csv("data/demog_HC3_full_py.csv")
read.csv("data/rt_quant_simple.csv")


#######-----Demographics Table Preparation-----#######
data <- read.csv("data/demog_HC3_full_py.csv")
data_original <- read.csv("data/demog_HC3_full.csv")

age_by_sex <- data |>
  summarise(mean_age = mean(DOC_AGE_HC3, na.rm=TRUE),
            sd_age   = sd(DOC_AGE_HC3, na.rm=TRUE),
            .by      = SEX) 
age_by_dementia <- data |> 
  summarise(mean_age = mean(DOC_AGE_HC3, na.rm=TRUE),
            sd_age   = sd(DOC_AGE_HC3, na.rm=TRUE),
            .by      = INC_DEMENTIA_ALL) 

age <- bind_rows(age_by_dementia, age_by_sex)

age_sum <- age |> pivot_longer(cols = -c(SEX, INC_DEMENTIA_ALL), 
                       names_to = "condition", 
                       values_to = "count") |>  
  pivot_wider(names_from=c(SEX, INC_DEMENTIA_ALL),
              values_from = count)

education <- data |> 
  group_by(EDLEVEL09, INC_DEMENTIA_ALL) |> 
  count() |> pivot_wider(names_from = INC_DEMENTIA_ALL,
                         values_from = n)
education_sex <- data |> group_by(EDLEVEL09, SEX) |> 
  count() |> pivot_wider(names_from = SEX, values_from = n)

education_all <- left_join(education, education_sex)

write.table(education_all, "data/results/education.txt")

marriage <- data |> group_by(INC_DEMENTIA_ALL, marital_st) |> count()
marriage_sex <- data |> group_by(SEX, marital_st) |> count()
marriage_all <- bind_rows(marriage, marriage_sex)
write.table(marriage_all, "data/results/marriage_all.txt", sep=",")


#####----comorbidity----#####
diseases <- data |> summarise(across(
  c(CVA3, CA3, MI3, DM3, DEPRESSION), 
  ~sum(.x==1, na.rm=TRUE)),
  .by= INC_DEMENTIA_ALL)

diseases_sex <- data |>
  summarise(across(
    c(CVA3, CA3, MI3, DM3, DEPRESSION), 
    ~sum(.x==1, na.rm=TRUE)),
    .by= SEX)

diseases_all <- bind_rows(diseases, diseases_sex) |> 
  rename(c(
    Depression= DEPRESSION,
    Prevalent_cancer = CA3,  
    Prevalent_stroke = CVA3,  
    Diabetes = DM3,  
    Prevalent_heart_attack = MI3
  ))

diseases_all <- diseases_all |> pivot_longer(cols = -c(SEX, INC_DEMENTIA_ALL), 
             names_to = "condition", 
             values_to = "count") |>  
  pivot_wider(names_from=c(SEX, INC_DEMENTIA_ALL),
              values_from = count)

## Prepare hearing and visual problems for table
# combine all visual problems in one column
data <- data |> rowwise() |>
  mutate(visual_problems = {
    vals <- c_across(c(VST_GLAUC, VST_CACT, VST_MACDEG, VST_RTPATH))
    first(vals[vals !=0], default=0)
  }) |> ungroup()

hearing_problems = data |> 
  summarise(hearing_problems = sum(HEARING_FU4==1, na.rm=T),
            .by=INC_DEMENTIA_ALL)
hearing_problems_sex <- data |> 
  summarise(hearing_problems = sum(HEARING_FU4==1, na.rm=T),
            .by=SEX)

visual_problems = data |>
  summarise(visual_prob = sum(visual_problems != 0, na.rm = T),
            .by=INC_DEMENTIA_ALL)
visual_problems_sex = data |>
  summarise(visual_prob = sum(visual_problems != 0, na.rm = T),
            .by=SEX)

hearing <- bind_rows(hearing_problems, hearing_problems_sex)
visual <- bind_rows(visual_problems, visual_problems_sex)

visual_prob <- visual |> pivot_longer(cols = -c(SEX, INC_DEMENTIA_ALL), 
                        names_to = "condition", 
                        values_to = "count") |>  
  pivot_wider(names_from=c(SEX, INC_DEMENTIA_ALL),
              values_from = count)

hearing_prob <- hearing |> pivot_longer(cols = -c(SEX, INC_DEMENTIA_ALL), 
                                             names_to = "condition", 
                                             values_to = "count") |>  
  pivot_wider(names_from=c(SEX, INC_DEMENTIA_ALL),
              values_from = count)


hear_visual <- bind_rows(visual_prob, hearing_prob)

comorbidity <- bind_rows(age_sum, diseases_all, hear_visual)

write.table(comorbidity, "data/results/comorbidity.txt", sep=",")


#####----MEAN RT----#####
meanRTs <- data |> 
  summarise(mean_comp   = mean(meanRT_comp), 
            sd_comp     = sd(meanRT_comp),
            mean_simple = mean(meanRT_simple_HC3, na.rm =T), 
            sd_simple   = sd(meanRT_simple_HC3,   na.rm =T),
            .by         = INC_DEMENTIA_ALL) 

meanRTs_sex <- data |> 
  summarise(mean_comp   = mean(meanRT_comp), 
            sd_comp     = sd(meanRT_comp),
            mean_simple = mean(meanRT_simple_HC3, na.rm =T), 
            sd_simple   = sd(meanRT_simple_HC3,   na.rm =T),
            .by         = SEX) 

meanRTs_all <- data |> 
  summarise(mean_comp   = mean(meanRT_comp), 
            sd_comp     = sd(meanRT_comp),
            mean_simple = mean(meanRT_simple_HC3, na.rm =T), 
            sd_simple   = sd(meanRT_simple_HC3,   na.rm =T)
            )

# combine mean RTs
all_meanRT_sum <- bind_rows(meanRTs, meanRTs_sex, meanRTs_all)

all_rt <- all_meanRT_sum |> pivot_longer(cols = -c(SEX, INC_DEMENTIA_ALL), 
                                             names_to = "mean_rt", 
                                             values_to = "count") |>  
  pivot_wider(names_from=c(SEX, INC_DEMENTIA_ALL),
              values_from = count)

#####----PARAMETERS----#####

params_sum <- data |>
  summarise(
    drfit_simple    = mean(v_simple, na.rm=T),
    drift_simple_sd = sd(v_simple, na.rm=T),
    bound_simp      = mean(a_simple, na.rm=T),
    bound_simp_sd   = sd(a_simple, na.rm=T),
    ndt_simp        = mean(t0_simple, na.rm=T),
    nst_simp_sd     = sd(t0_simple, na.rm=T),
    inc_drift       = mean(dv_comp, na.rm=T),
    inc_drift_Sd    = sd(dv_comp, na.rm=T),
    bound_comp      = mean(a_comp, na.rm=T),
    bound_comp_sd   = sd(a_comp, na.rm=T),
    ndt_comp        = mean(t0_comp, na.rm=T),
    ndt_comp_sd     = sd(t0_comp, na.rm=T),
    .by             = c(INC_DEMENTIA_ALL, SEX)
  )


params_sum_dement <- data |>
  summarise(
    drfit_simple    = mean(v_simple, na.rm=T),
    drift_simple_sd = sd(v_simple, na.rm=T),
    bound_simp      = mean(a_simple, na.rm=T),
    bound_simp_sd   = sd(a_simple, na.rm=T),
    ndt_simp        = mean(t0_simple, na.rm=T),
    nst_simp_sd     = sd(t0_simple, na.rm=T),
    inc_drift       = mean(dv_comp, na.rm=TRUE),
    inc_drift_Sd    = sd(dv_comp, na.rm=T),
    bound_comp      = mean(a_comp, na.rm=T),
    bound_comp_sd   = sd(a_comp, na.rm=T),
    ndt_comp        = mean(t0_comp, na.rm=T),
    ndt_comp_sd     = sd(t0_comp, na.rm=T),
    .by             = INC_DEMENTIA_ALL
  )

# combine men and women + dementia outcomes
all_params_sum <- bind_rows(params_sum, params_sum_dement)

all_params <- all_params_sum |> pivot_longer(cols = -c(SEX, INC_DEMENTIA_ALL), 
                              names_to = "condition", 
                              values_to = "count") |>  
  pivot_wider(names_from=c(SEX, INC_DEMENTIA_ALL),
              values_from = count)

write.csv(all_params, 'data/results/all_params.csv')


