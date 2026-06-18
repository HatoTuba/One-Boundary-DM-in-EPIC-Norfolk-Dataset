### THIS SCRIPT CONTAINTS CUSTOM FUNCTIONS THAT HAS BEEN USED IN THE ANALYSES 

# these are not exactly quantiles but I lose too many data points otherwise
# Quantile based, does not work well for this case
find_outlier1 <- function(x) {
  Q1 <-quantile(x, 0.01)
  IQR <- Q3 - Q1
  Q3 <- quantile(x, 0.999)
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  is_outlier = x < lower_bound | x > upper_bound
}

# OUTLIER FOR SIMPLE VST
find_outlier2_simp <- function(x) {
  lower_bound = 250
  upper_bound = 2500
  is_outlier = x > upper_bound
  is_outlier = x < lower_bound | x > upper_bound
}

####-----Data Elimination Complex Task-----####
# hard limits
find_outlier2 <- function(x) {
  lower_bound = 300
  upper_bound = 10000 # quantile .999 was smaller than 10 secs
  is_outlier = x < lower_bound | x > upper_bound
}

# to eliminate outliers in the complex VST RTs 
eliminate_outliers <- function(data){
  data_outlier <- data %>% 
    group_by(RELEASEID) %>% 
    mutate(is_outlier = find_outlier2(TIME))
  outlier_sum <- data_outlier %>%
    group_by(RELEASEID) %>% 
    summarize(outlier_count = sum(is_outlier))
  clean_data <- data_outlier %>% 
    filter(!is_outlier) %>% select(-TEST, )
  return(clean_data)
}

####----------------------------------------------------------------------------------------
# To treat the data, give proper names to variables
rename <- function(data){
  data <- 
    data |> 
    mutate(INC_DEMENTIA_ALL = ifelse(INC_DEMENTIA_ALL==1, "Dementia", "No Dementia")) |>
    mutate(SEX = ifelse(SEX==1, "Men", "Women")) |>
    # combine O and A levels as in the other study
    mutate(EDLEVEL09 = recode(EDLEVEL09, `0`="No Ed", `1`="O-A Level", `2`="O-A Level", `3`="Degree")) |>
    mutate(age=round(DOC_AGE_HC3)) |> # age=round(DOC_AGE_HC5), 
    mutate(age_group10 = cut(DOC_AGE_HC3, 
                               breaks=seq(1, 100, by=10), right=FALSE)) |>
    rowwise() |>
      mutate(visual_problems = {
        vals <- c_across(c(VST_GLAUC, VST_CACT, VST_MACDEG, VST_RTPATH))
        first(vals[vals !=0], default=0)
      }) |> ungroup()
  # mutate(visual_problems = recode(visual_problems, `NA`="Unknown", `0`="No", `1`="Yes", `3`="Unknown"))

  return(data)
}

meanRT <- function(data){
  rt <- data |>
    group_by(RELEASEID) |>
    summarise(meanRT = mean(TIME))
  return(rt)
}


# prepare the data for survival analyses
surv_prepare <- function(data){
  # prepare the factors
  data <- data |>
    mutate(SEX = recode(SEX, `Men`=0, `Women`=1)) |>
    mutate(INC_DEMENTIA_ALL = recode(INC_DEMENTIA_ALL, `No Dementia`=0, `Dementia`=1)) |>
    mutate(EDLEVEL09 = recode(EDLEVEL09, `No Ed`=0, `O-A Level`=1, `Degree`=2)) |>
    mutate(SEX=factor(SEX) |> fct_recode("Men"="0", "Women"="1"),
           EDLEVEL09 = factor(EDLEVEL09) |>
             fct_recode("No degree" = "0",
                        "O-A  Level" = "1",
                        "Degree"="2"))
}


plot_RT <- function(data){
  # mean RT regression line to all people
  # add age group for each 10 years
  data <- data |> 
    mutate(age_group10 = cut(DOC_AGE_HC3, 
                             breaks=seq(1, 100, by=10), right=FALSE))
  # mean RT regression separately for each age group
  data |> #filter(INC_DEMENTIA_ALL=="No Dementia") |>
    ggplot(aes(x=age, y=meanRT, color=age_group10)) +
    geom_point(size=0.5)+
    # coord_cartesian(ylim=c(1, 5)) +
    geom_smooth(method="lm", se=T) +
    labs(title = "Mean RT as a function of age",
         x="Ages",
         y="Mean RT",
         color = "Age group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    # scale_color_brewer() +
    facet_wrap(~INC_DEMENTIA_ALL, scales="fixed") 
}

