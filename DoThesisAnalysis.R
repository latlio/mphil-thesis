library(tidyverse)
library(survival)
library(survminer)
library(compareGroups)
library(sjlabelled)
library(rlang)
library(naniar)
library(patchwork)
library(corrr)
library(cmprsk)
library(scales)

data_path <- "/Users/lathanliou/Desktop/Academic/Cambridge/Thesis/Data/FinalDataPrep"
setwd(data_path)
source("PrepareThesisData.R")
final_data_path <- "FinalDataPrep/final_thesis_data.csv"
analysis_output_path <- "/Users/lathanliou/Desktop/Academic/Cambridge/Thesis/Statistical Analyses/FinalAnalyses"
analysis_output_path_baseline <- "/Users/lathanliou/Desktop/Academic/Cambridge/Thesis/Statistical Analyses/FinalAnalyses/BaselineCharacteristics"
analysis_output_path_exploratory <- "/Users/lathanliou/Desktop/Academic/Cambridge/Thesis/Statistical Analyses/FinalAnalyses/Exploratory"
analysis_output_path_KM <- "/Users/lathanliou/Desktop/Academic/Cambridge/Thesis/Statistical Analyses/FinalAnalyses/KMCurves"
analysis_output_path_missing <- "/Users/lathanliou/Desktop/Academic/Cambridge/Thesis/Statistical Analyses/FinalAnalyses/Missingness"
analysis_output_path_model <- "/Users/lathanliou/Desktop/Academic/Cambridge/Thesis/Statistical Analyses/FinalAnalyses/Modelling"

PlotDistributions <- function(data) {
  var_names <- list(
    'AgeMenarche' = "Age at Menarche",
    'BMI' = "BMI",
    'DiagnosisAge' = "Age at Diagnosis",
    'prs_z' = "Standardized PRS"
  )
  
  var_labeller <- function(variable, value) {
    return(var_names[value])
  }
  
  plot <- data %>%
    select(DiagnosisAge, AgeMenarche, BMI, prs_z) %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    ggplot(aes(value)) + 
    facet_wrap(~ variable, scales = "free", labeller = var_labeller) +
    geom_histogram() + 
    theme_bw() + 
    labs(x = "Value",
         y = "Frequency")
  
  ggsave(paste0(analysis_output_path_exploratory, "/dist_covars.png"), 
         plot, device = "png", width = 5, height = 5)
  
  return(plot)
}

AssessCorrelation <- function(data) {
  corr <- data %>%
    select(prs_z, DiagnosisAge, BMI, AgeMenarche, Parity, IMD) %>%
    mutate(logbmi = log(BMI)) %>%
    select(-BMI) %>%
    corrr::correlate(x = ., 
                     method = "spearman", 
                     quiet = TRUE) %>% 
    corrr::rearrange(x = ., 
                     method = "MDS", 
                     absolute = FALSE)
  
  return(corr)
}

PrepareDescriptiveData <- function(data) {
  descriptive_data <- data %>%
    replace_na(., list(Radiotherapy = 0,
                       Chemotherapy = 0,
                       HormoneTherapy = 0)) %>%
    mutate(groupcause = case_when(
      completeI2025 == 1 ~ "Incident CAD",
      bc == 1 ~ "Died from Breast Cancer",
      completeI2025 == 0 & bc == 0 & VitalStatus == 1 ~ "Died from Other Causes",
      TRUE ~ "Survivors"
    ),
    grouptherapy = case_when(
      Chemotherapy == 1 & Radiotherapy == 0 & HormoneTherapy == 0 ~ "Only Chemotherapy",
      Radiotherapy == 1 & Chemotherapy == 0 & HormoneTherapy == 0 ~ "Only Radiotherapy",
      HormoneTherapy == 1 & Chemotherapy == 0 & Radiotherapy == 0 ~ "Only Hormone Therapy",
      Chemotherapy == 1 & Radiotherapy == 1 & HormoneTherapy == 0 ~ "Chemo + Radio",
      Chemotherapy == 1 & Radiotherapy == 0 & HormoneTherapy == 1 ~ "Chemo + Hormone",
      Chemotherapy == 0 & Radiotherapy == 1 & HormoneTherapy == 1 ~ "Radio + Hormone",
      Chemotherapy == 1 & Radiotherapy == 1 & HormoneTherapy == 1 ~ "All Therapy",
      Chemotherapy == 0 & Radiotherapy == 0 & HormoneTherapy == 0 ~ "No Therapy"
    ),
    groupcad = case_when(
      completeI2025 == 1 ~ "General CAD",
      completeI2025 == 1 ~ "Clinical CAD",
      completeI20 == 1 ~ "Angina",
      completeI2123 == 1 ~ "MI",
      completeI2425 == 1 ~ "Other CAD",
      TRUE ~ "Non-CAD"
    ),
    genotyped = case_when(
      !is.na(prs_z) ~ 1,
      TRUE ~ 0
    ),
    Chemotherapy = case_when(
      Chemotherapy == 1 ~ "Yes", 
      Chemotherapy == 0 ~ "No"
    ),
    Radiotherapy = case_when(
      Radiotherapy == 1 ~ "Yes",            
      Radiotherapy == 0 ~ "No"
    ),
    HormoneTherapy = case_when(
      HormoneTherapy == 1 ~ "Yes",
      HormoneTherapy == 0 ~ "No"
    ),
    TumourGrade = case_when(
      TumourGrade == 1 ~ "Well Differentiated",
      TumourGrade == 2 ~ "Moderately Differentiated",
      TumourGrade == 3 ~ "Poorly/Undifferentiated"
    ),
    DetectedByScreeningMammogram = case_when(
      DetectedByScreeningMammogram == 1 ~ "Yes",
      DetectedByScreeningMammogram == 0 ~ "No"
    ),
    eduCat = case_when(
      eduCat == 1 ~ "GCSE or similar",
      eduCat == 2 ~ "A-Level of similar",
      eduCat == 3 ~ "Graduate"
    ),
    HRTEver = case_when(
      HRTEver == 1 ~ "Yes",
      HRTEver == 0 ~ "No"
    ),
    smokingEver = case_when(
      smokingEver == 0 ~ "Never",
      smokingEver == 1 ~ "Past",
      smokingEver == 2 ~ "Current, in last year"
    ),
    AlcNow = case_when(
      AlcNow == 1 ~ "Yes",
      AlcNow == 0 ~ "No"
    ),
    EthnicityClass = case_when(
      EthnicityClass == 1 ~ "European",
      EthnicityClass == 2 ~ "Hispanic-American",
      EthnicityClass == 3 ~ "African",
      EthnicityClass == 4 ~ "Asian",
      EthnicityClass == 5 ~ "South-east Asian",
      EthnicityClass == 6 ~ "Other"
    ))
  
  labelled_data <- sjlabelled::set_label(descriptive_data, c("Unique individual ID",
                                                     "Years between diagnosis and study entry",
                                                     "Years between diagnosis and SEARCH follow-up",
                                                     "Years between diagnosis and I20 event",
                                                     "Years between diagnosis and I2125 event",
                                                     "Years between diagnosis and I2123 event",
                                                     "Years between diagnosis and I2425 event",
                                                     "Years between diagnosis and I2025 event",
                                                     "Age at diagnosis",
                                                     "Vital status at last follow-up",
                                                     "Received adjuvant chemotherapy",
                                                     "Received adjuvant radiotherapy",
                                                     "Received adjuvant hormone therapy",
                                                     "Histopathological grade",
                                                     "Tumour maximum diameter (mm)",
                                                     "Number of nodes excised",
                                                     "Number of nodes involved",
                                                     "Detected by screening mammogram",
                                                     "Primary Cause of Death 1a",
                                                     "Primary Cause of Death 1b",
                                                     "Primary Cause of Death 1c",
                                                     "Secondary Cause of Death", 
                                                     "Underlying Cause of Death",
                                                     "Age at diagnosis",
                                                     "Highest level of education received",
                                                     "Age at menarche",
                                                     "Number of full-term pregnancies",
                                                     "Height (cm)",
                                                     "Weight (kg)",
                                                     "BMI", 
                                                     "Use of hormonal replacement therapy",
                                                     "Smoking history", 
                                                     "Currently drinking alcohol",
                                                     "Ethnicity", 
                                                     "Menstrual cycle frequency",
                                                     "Thyroid disease",
                                                     "Index of Multiple Deprivation Score",
                                                     "Breast cancer-related death",
                                                     "Polygenic Risk Score",
                                                     "Standardized Polygenic Risk Score",
                                                     "Angina event",
                                                     "Clinical CAD event",
                                                     "MI event",
                                                     "Other CAD event",
                                                     "General CAD event",
                                                     "Standardized Polygenic Risk Score Quintile",
                                                     "Standardized BMI",
                                                     "BMI Category",
                                                     "Competing Risk"))
 return(labelled_data) 
}

MakeBaselineTable <- function(data, groupvar) {
  formula <- paste0(sym(groupvar), "~ DiagnosisAge + 
    Chemotherapy + 
    Radiotherapy + 
    HormoneTherapy + 
    TumourGrade + 
    TumourSize + 
    Nodes_excised + 
    Nodes_involved + 
    DetectedByScreeningMammogram + 
    eduCat + 
    IMD +
    AgeMenarche + 
    Parity + 
    height +
    weight + 
    BMI +
    HRTEver + 
    smokingEver + 
    AlcNow + 
    EthnicityClass + 
    RegularMenstrualCycle + 
    ThyroidDisease")
  res <- compareGroups::compareGroups(formula = as.formula(formula),
                                      data = data,
                                      method = c(DiagnosisAge = 1,
                                                 Chemotherapy = 3,
                                                 Radiotherapy = 3,
                                                 HormoneTherapy = 3,
                                                 TumourGrade = 3,
                                                 TumourSize = 2,
                                                 Nodes_excised = 2,
                                                 Nodes_involved = 2,
                                                 DetectedByScreeningMammogram = 3,
                                                 eduCat = 3,
                                                 Parity = 2,
                                                 HRTEver = 3,
                                                 smokingEver = 3,
                                                 AlcNow = 3,
                                                 EthnicityClass = 3,
                                                 RegularMenstrualCycle = 3,
                                                 ThyroidDisease = 3))
  
  restab <- compareGroups::createTable(res, hide.no = "No", show.all = TRUE)
  
  compareGroups::export2word(restab, 
                             file = paste0(analysis_output_path_baseline,
                                           "/", 
                                           groupvar,
                                           "-baseline_characteristics.docx"))
}

PlotPRSByMissing <- function(data, missingvar, label) {
  plot <- ggplot(data, aes(x = prs_z, 
                           fill = !!sym(missingvar))) + 
    geom_density(alpha = 0.5) + 
    scale_fill_discrete(name = paste0(missingvar),
                        labels = c("Not missing", "Missing")) +
    scale_y_continuous(limits = c(0, 0.5), expand = c(0,0)) +
    theme_bw() + 
    labs(x = "Standardized PRS",
         y = "Density",
         title = label) + 
    theme(legend.position = "none")
  # 
  # ggsave(paste0(analysis_output_path_missing, "/", missingvar, "-prsplot.png"), 
  #        plot, device = "png", width = 5, height = 5)
}

PerformMissignessAnalysis <- function(data) {
  variables_of_interest <- data %>%
    select(prs_z,
           AgeMenarche,
           eduCat,
           HRTEver,
           IMD,
           RegularMenstrualCycle,
           smokingEver, 
           ThyroidDisease)
  plot_missing_factors <- naniar::vis_miss(variables_of_interest)
  
  ggsave(paste0(analysis_output_path_missing, "/", "missingplot.png"), 
         plot_missing_factors, device = "png", width = 5, height = 5)
  
  voi_shadow <- naniar::bind_shadow(variables_of_interest)
  
  NAvars <- voi_shadow %>%
    select(ends_with("_NA")) %>%
    colnames()
  
  missingvarlabels <- c("PRS",
                        "Age at menarche",
                        "Education",
                        "Hormone replacement therapy",
                        "Index of Multiple Deprivation",
                        "Menstrual cycle frequency",
                        "Smoking",
                        "Thyroid disease")
  
  plot_list <- map2(NAvars, missingvarlabels,
                    PlotPRSByMissing,
                    data = voi_shadow)
  
  return(plot_list)
}

PlotACMKMCurves <- function(data, var, label, subgroup_labels) {
  #all-cause mortality
  #YearsToStatus-YearsToEntry
  surv_formula <- paste0("Surv(YearsToStatus, VitalStatus) ~ ", var) %>%
    as.formula()
  
  fit <- surv_fit(surv_formula, data = data)
  coxfit <- coxph(surv_formula, data = data)
  
  unadjcurves <- survminer::ggsurvplot(fit,
                                       conf.int = FALSE,
                                       censor = FALSE,
                                       legend = "right",
                                       legend.title = paste0(label),
                                       legend.labs = subgroup_labels,
                                       risk.table = TRUE,
                                       xlim = c(0,20),
                                       ylim = c(0.5,1),
                                       break.time.by = 5)
  
  ggsave(paste0(analysis_output_path_KM, "/", var, "-crude-ACM-curves.png"), 
         print(unadjcurves), device = "png", width = 8, height = 7)
}

PlotAgeAdjKMCurves <- function(data, var) {
  surv_formula <- paste0("Surv(DiagnosisAge + YearsToEntry,
                         DiagnosisAge + YearsToI2025, completeI2025) ~ ", var) %>%
    as.formula()

  fit <- surv_fit(surv_formula, data = data)
  coxfit <- coxph(surv_formula, data = data)

  unadjcurves <- survminer::ggsurvplot(fit,
                                       conf.int = TRUE,
                                       legend.title = paste0(var),
                                       risk.table = TRUE)

  adjcurves <- survminer::ggadjustedcurves(coxfit,
                                           data = as.data.frame(data),
                                           variable = paste0(var),
                                           legend = "right",
                                           legend.title = paste0(var),
                                           caption = "Adjusted Survival Curves (Age Time Scale)",
                                           font.caption = 8,
                                           xlab = "Time",
                                           ylab = "Survival Probability")

  ggsave(paste0(analysis_output_path, "/", var, "-ageunadjsurvcurves.png"),
         unadjcurves$plot, device = "png", width = 5, height = 5)

  ggsave(paste0(analysis_output_path, "/", var, "-ageadjsurvcurves.png"),
         adjcurves, device = "png", width = 5, height = 5)
}

ExploreKM <- function(data, var, label, subgroup_labels) {
  surv_formula <- paste0("Surv(YearsToI2025, completeI2025) ~ ", var) %>%
    as.formula()
  
  fit <- surv_fit(surv_formula, data = data)
  
  #plot KM
  km <- survminer::ggsurvplot(fit,
                              conf.int = FALSE,
                              censor = FALSE,
                              legend = "right",
                              legend.title = paste0(label),
                              ylab = "Survival Rate",
                              legend.labs = subgroup_labels,
                              risk.table = FALSE,
                              xlim = c(0,20),
                              ylim = c(0.75, 1),
                              break.time.by = 5)
  #plot cumhaz
  cumhaz <- survminer::ggsurvplot(fit,
                              conf.int = FALSE,
                              censor = FALSE,
                              legend = "right",
                              legend.title = paste0(label),
                              ylab = "Cumulative Probability of CAD",
                              legend.labs = subgroup_labels,
                              risk.table = TRUE,
                              xlim = c(0,20),
                              ylim = c(0, 0.25),
                              break.time.by = 5,
                              fun = "cumhaz")
  fin <- km$plot / cumhaz$plot
  
  ggsave(paste0(analysis_output_path_KM, "/", var, "-fullKM.png"), 
         fin, device = "png", width = 7, height = 10)
}

PlotPRSDensity <- function(data) {
  p <- data %>%
    mutate(completeI2025 = as.factor(completeI2025)) %>%
    ggplot(aes(x = prs_z, fill = completeI2025)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = mean(data %>% 
                                   filter(completeI2025 == 1) %>% 
                                   select(prs_z) %>% pull(),
                                 na.rm = TRUE), col = "#93DBDE") +
    geom_vline(xintercept = mean(data %>%
                                   filter(completeI2025 == 0) %>%
                                   select(prs_z) %>% pull(),
                                 na.rm = TRUE), col = "#F8B5B1") + 
    theme_bw() + 
    labs(x = "Polygenic Risk Score",
         y = "Density") + 
    scale_fill_discrete(name = "Status", labels = c("No Event", 
                                                    "Incident CAD Event")) + 
    scale_y_continuous(limits = c(0, 0.5), expand = c(0,0))
  ggsave(paste0(analysis_output_path_exploratory, "/PRS-density-plot.png"), p, 
         device = "png", width = 5, height = 5)
}

FitCoxSensitivity <- function(data, outcome) {
  data <- data %>%
    dplyr::rename(eventindicator = !!sym(paste0("complete", outcome)),
           YearsToEvent = !!sym(paste0("YearsTo", outcome)))
  
  fit <- coxph(Surv(DiagnosisAge + YearsToEntry,
                    DiagnosisAge + YearsToEvent,
                    eventindicator) ~ prs_z,
               data) %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = T)
  
  return(fit)
}

PlotCoxSensitivity <- function(coefs_df) {
  plot <- coefs_df %>%
    mutate(definition = case_when(
      outcome == "I2125" ~ "Clinical",
      outcome == "I2123" ~ "MI",
      outcome == "I2425" ~ "Other",
      outcome == "I2025" ~ "General"
    )) %>%
    ggplot(aes(x = definition, y = estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = 0.05,
                  position = position_dodge(width = 0.2)) +
    theme_bw() +
    labs(x = "CAD Definition",
         y = "Hazard Ratio") +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(1, 2)) + 
    coord_flip()
  
  # data <- structure(list(
  #   mean  = c(NA, 1.34, 1.27, 1.35, 1.36), 
  #   lower = c(NA, 1.23, 1.07, 1.24, 1.23),
  #   upper = c(NA, 1.46, 1.51, 1.48, 1.47)),
  #   .Names = c("mean", "lower", "upper"), 
  #   row.names = c(NA, -5L), 
  #   class = "data.frame")
  # 
  # tabletext <- cbind(
  #   c("Definition", "General (I20-25)", "MI (I21-23)", "Clinical (I21-25)", 
  #     "Other (I24-25)"),
  #   c("Events", "934", "170", "794", "763"),
  #   c("HR", "1.34", "1.27", "1.35", "1.36"),
  #   c("95% CI", "(1.23, 1.46)", "(1.07, 1.51)", "(1.24, 1.48)", "(1.23, 1.47)")
  # )
  # 
  # plot <- forestplot::forestplot(tabletext, 
  #                                data,
  #                                hrzl_lines = gpar(col="#444444"),
  #                                new_page = T,
  #                                is.summary=c(TRUE,rep(FALSE,4),TRUE),
  #                                clip=c(0.1,2.5),
  #                                col=fpColors(box="royalblue",
  #                                             line="darkblue", 
  #                                             summary="royalblue"))
  
  ggsave(paste0(analysis_output_path, "/", "sensitivityplot.png"), 
         plot, device = "png", width = 5, height = 5)
}

FitCoxControls <- function(data, outcome) {
formula <- paste0("Surv(YearsToEntry, YearsToStatus,",
                  outcome, ") ~ prs_z + DiagnosisAge")
  fit <- coxph(as.formula(formula),
               data) %>%
    broom::tidy(exponentiate = TRUE)
  
  return(fit)
}

AssessLTDep <- function(data) {
  
  fit <- coxph(Surv(YearsToI2025,
                    completeI2025) ~ YearsToEntry,
               data) %>%
    broom::tidy()

  return(fit)
}

ObtainSchoenfeldp <- function(fit) {
  output <- cox.zph(fit)$table[1,3]
  return(output)
}

AssessQuintileDoseResponse <- function(data, quintilevar, var, label, width) {
  data <- data %>%
    mutate(logbmi_quintile = cut(log(BMI), quantile(log(BMI), 
                                                    probs = 0:5/5,
                                                    na.rm = TRUE),
                                 include.lowest = TRUE),
           age_quintile = cut(DiagnosisAge, quantile(DiagnosisAge,
                                                     probs = 0:5/5,
                                                     na.rm = TRUE),
                              include.lowest = TRUE))
  
  formula <- paste0("Surv(YearsToEntry, YearsToI2025, completeI2025) ~ ",
                    quintilevar)
  fit <- coxph(as.formula(formula),
               data) %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = T)
  
  meanquintiles <- data %>%
    drop_na(!!sym(quintilevar)) %>%
    group_by(!!sym(quintilevar)) %>%
    summarize(mean = mean(eval(parse(text = var)), na.rm = T))
  
  doseplot <- fit %>%
    select(estimate, conf.low, conf.high) %>%
    add_row(estimate = 1, conf.low = 1, conf.high = 1, .before = 1) %>%
    mutate(Quintile = 1:nrow(.),
           mean = meanquintiles %>% pull()) %>%
    ggplot(aes(x = mean, y = estimate, col = as.factor(Quintile))) + 
    geom_point() +
    scale_y_continuous(trans = log_trans(),
                       labels = scales::number_format(accuracy = 1)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = width) + 
    theme_bw() + 
    labs(x = "Mean",
         y = "Hazard Ratio",
         col = "Quintile",
         title = label)
  
  return(doseplot)
}

AssessControlDoseResponse <- function(data, outcome, modelterms, quintilevar,
                                      var, label, width){
  data <- data %>%
    mutate(logbmi_quintile = cut(log(BMI), quantile(log(BMI), 
                                                    probs = 0:5/5,
                                                    na.rm = TRUE),
                                 include.lowest = TRUE),
           age_quintile = cut(DiagnosisAge, quantile(DiagnosisAge,
                                                     probs = 0:5/5,
                                                     na.rm = TRUE),
                              include.lowest = TRUE))
  if (outcome == "all") {
    formula <- paste0("Surv(YearsToEntry, YearsToStatus, VitalStatus) ~ ",
                      modelterms)
    fit <- coxph(as.formula(formula),
                 data) %>%
      broom::tidy(exponentiate = TRUE,
                  conf.int = T)
    
    meanquintiles <- data %>%
      drop_na(!!sym(quintilevar)) %>%
      group_by(!!sym(quintilevar)) %>%
      summarize(mean = mean(eval(parse(text = var)), na.rm = T))
    
    doseplot <- fit %>%
      select(estimate, conf.low, conf.high) %>%
      add_row(estimate = 1, conf.low = 1, conf.high = 1, .before = 1) %>%
      dplyr::slice(1:5) %>%
      mutate(Quintile = 1:nrow(.),
             mean = meanquintiles %>% pull()) %>%
      ggplot(aes(x = mean, y = estimate, col = as.factor(Quintile))) + 
      geom_point() +
      scale_y_continuous(trans = log_trans(),
                         labels = scales::number_format(accuracy = 0.1)) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                    width = width) + 
      theme_bw() + 
      labs(x = "Mean",
           y = "Hazard Ratio",
           col = "Quintile",
           title = label)
  }
  else {
    formula <- paste0("Surv(YearsToEntry, YearsToStatus, bc) ~ ",
                      modelterms)
    fit <- coxph(as.formula(formula),
                 data) %>%
      broom::tidy(exponentiate = TRUE,
                  conf.int = T)
    
    meanquintiles <- data %>%
      drop_na(!!sym(quintilevar)) %>%
      group_by(!!sym(quintilevar)) %>%
      summarize(mean = mean(eval(parse(text = var)), na.rm = T))
    
    doseplot <- fit %>%
      select(estimate, conf.low, conf.high) %>%
      add_row(estimate = 1, conf.low = 1, conf.high = 1, .before = 1) %>%
      dplyr::slice(1:5) %>%
      mutate(Quintile = 1:nrow(.),
             mean = meanquintiles %>% pull()) %>%
      ggplot(aes(x = mean, y = estimate, col = as.factor(Quintile))) + 
      geom_point() +
      scale_y_continuous(trans = log_trans(),
                         labels = scales::number_format(accuracy = 0.1)) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                    width = width) + 
      theme_bw() + 
      labs(x = "Mean",
           y = "Hazard Ratio",
           col = "Quintile",
           title = label)
  }
  return(doseplot)
}

AnalyzeSensitivityExp <- function(data) {
  fit_smoke <- coxph(Surv(YearsToEntry,
                          YearsToI2025,
                          completeI2025) ~ smokingEver,
                     data)
  fit_smoke2 <- coxph(Surv(YearsToEntry,
                           YearsToI2025,
                           completeI2025) ~ as.factor(smokingEver),
                      data)
  lrt <- anova(fit_smoke, fit_smoke2)
  
  res_smoke <- fit_smoke %>%
    broom::tidy(exponentiate = TRUE)
  res_smoke2 <- fit_smoke2 %>%
    broom::tidy(exponentiate = TRUE)
  
  return(list("lrt" = lrt,
              "smoking_ord" = res_smoke,
              "smoking_cat" = res_smoke2))
}

PerformUnivariateModelling <- function(data, var) {
  data <- data %>%
    drop_na(c("Chemotherapy",
              "Radiotherapy",
              "HormoneTherapy",
              "AgeMenarche",
              "BMI",
              "smokingEver",
              "AlcNow",
              "ThyroidDisease",
              "Parity",
              "HRTEver", 
              "IMD"))
  
  formula <- paste0("Surv(YearsToEntry, YearsToI2025, completeI2025) ~ ",
                    var)
  fit <- coxph(as.formula(formula),
               data) %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = T)
  
  return(fit)
}

PerformSequentialModelling <- function(data) {
  fit1 <- coxph(Surv(YearsToEntry,
                    YearsToI2025,
                    completeI2025) ~ prs_z + DiagnosisAge,
               data)
  
  fit2 <- coxph(Surv(YearsToEntry,
                      YearsToI2025,
                      completeI2025) ~ prs_z + DiagnosisAge +
                   Chemotherapy + Radiotherapy + HormoneTherapy ,
                 data %>%
                   drop_na(c("Chemotherapy",
                             "Radiotherapy",
                             "HormoneTherapy")))
  fit3a <- coxph(Surv(YearsToEntry,
                      YearsToI2025,
                      completeI2025) ~ prs_z + DiagnosisAge + 
                   Chemotherapy + Radiotherapy + HormoneTherapy,
                 data %>%
                   drop_na(c("Chemotherapy",
                             "Radiotherapy",
                             "HormoneTherapy",
                             "AgeMenarche",
                             "BMI",
                             "smokingEver",
                             "AlcNow",
                             "ThyroidDisease",
                             "Parity",
                             "HRTEver",
                             "IMD")))
  fit3b <- coxph(Surv(YearsToEntry,
                      YearsToI2025,
                      completeI2025) ~ prs_z + DiagnosisAge +
                   Chemotherapy + Radiotherapy + HormoneTherapy +
                   log(IMD) +
                   AgeMenarche + log(BMI) + as.factor(smokingEver) + AlcNow + 
                   ThyroidDisease + Parity + HRTEver,
                 data %>%
                   drop_na(c("Chemotherapy",
                             "Radiotherapy",
                             "HormoneTherapy",
                             "AgeMenarche",
                             "BMI",
                             "smokingEver",
                             "AlcNow",
                             "ThyroidDisease",
                             "Parity",
                             "HRTEver",
                             "IMD")))
  list_fits <- list("Model 1" = fit1, 
                    "Model 2" = fit2, 
                    "Model 3a" = fit3a, 
                    "Model 3b" = fit3b)
  primaryanalysis <-  list_fits %>% map_dfr(broom::tidy, exponentiate = TRUE,
                                            conf.int = T)
  counts <- list_fits %>% map_dfr(broom::glance)
  Schoenfeldres <- list_fits %>% map_dfr(ObtainSchoenfeldp)
  
  write_csv(primaryanalysis,
            paste0(analysis_output_path_model, "/primary_analysis.csv"))
  write_csv(counts,
            paste0(analysis_output_path_model, "/primary_analysis_counts.csv"))
  write_csv(Schoenfeldres, 
            paste0(analysis_output_path_model, "/primary_analysis_ph_assumption.csv"))
  
  return(primaryanalysis)
}

AnalyzeFinalModel <- function(data) {
  #statistically and nominally significant variables were kept
  finalfit <- coxph(Surv(YearsToEntry,
                         YearsToI2025,
                         completeI2025) ~ prs_z + DiagnosisAge + HormoneTherapy +
                      log(BMI) + as.factor(smokingEver) + AlcNow + log(IMD),
                    data)
  res <- finalfit %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = T)
  
  #prs quintile
  finalfit2 <- coxph(Surv(YearsToEntry,
                          YearsToI2025,
                          completeI2025) ~ prs_z_quintile + DiagnosisAge + HormoneTherapy +
                       log(BMI) + as.factor(smokingEver) + AlcNow + log(IMD),
                     data)
  res2 <- finalfit2 %>%
    broom::tidy(exponentiate = T,
                conf.int = T)
  
  #interact with anti-hormone therapy
  finalfit3 <- coxph(Surv(YearsToEntry,
                          YearsToI2025,
                          completeI2025) ~ prs_z*HormoneTherapy + DiagnosisAge +
                       log(BMI) + as.factor(smokingEver) + AlcNow + log(IMD),
                     data)
  res3 <- finalfit3 %>%
    broom::tidy(exponentiate = T,
                conf.int = T)
  
  #interaction with bmi
  finalfit4 <- coxph(Surv(YearsToEntry,
                          YearsToI2025,
                          completeI2025) ~ prs_z*log(BMI) + HormoneTherapy + DiagnosisAge +
                       as.factor(smokingEver) + AlcNow + log(IMD),
                     data)
  res4 <- finalfit4 %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = T)
  
  #interaction with bmi WHO categories
  finalfit5 <- coxph(Surv(YearsToEntry,
                          YearsToI2025,
                          completeI2025) ~ prs_z*as.factor(bmi_cat) + HormoneTherapy + 
                       DiagnosisAge + as.factor(smokingEver) + AlcNow + log(IMD),
                     data)
  res5 <- finalfit5 %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = T)
  
  #interaction with smoking
  finalfit6 <- coxph(Surv(YearsToEntry,
                          YearsToI2025,
                          completeI2025) ~ prs_z*as.factor(smokingEver) + HormoneTherapy +
                       DiagnosisAge + log(BMI) + AlcNow + log(IMD),
                     data)
  res6 <- finalfit6 %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = T)
  
  #interaction with time
  finalfit7 <- coxph(Surv(YearsToEntry,
                          YearsToI2025,
                          completeI2025) ~ prs_z*YearsToI2025 + DiagnosisAge + HormoneTherapy +
                       log(BMI) + as.factor(smokingEver) + AlcNow + log(IMD),
                     data)
  res7 <- finalfit7 %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = T)
  
  #interaction with alcohol
  finalfit8 <- coxph(Surv(YearsToEntry,
                          YearsToI2025,
                          completeI2025) ~ prs_z*AlcNow + DiagnosisAge +
                       HormoneTherapy + log(BMI) + as.factor(smokingEver) + 
                       log(IMD),
                     data)
  res8 <- finalfit8 %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = T)
  
  #adjusted curve
  curve <- survminer::ggadjustedcurves(finalfit3,
                                       data = as.data.frame(data),
                                       variable = "HormoneTherapy",
                                       legend = "right",
                                       legend.title = "Hormone Therapy",
                                       xlab = "Time",
                                       ylab = "Survival Probability",
                                       xlim = c(0, 20),
                                       ylim = c(0.7, 1))
  ggsave(paste0(analysis_output_path_KM, "/adj_horm_KM.png"), 
         curve, device = "png", width = 8, height = 7)
  
  # facetfit <- coxph(Surv(DiagnosisAge + YearsToEntry,
  #                        DiagnosisAge + YearsToI2025,
  #                        completeI2025) ~ prs_z +
  #                     BMI + smokingEver + AlcNow,
  #                   data)
  # 
  # #facet by hormone therapy looking at prs quintile
  # facet_curve1 <- survminer::ggadjustedcurves(facetfit,
  #                                            data = as.data.frame(data %>%
  #                                                                   filter(HormoneTherapy == 1)),
  #                                            variable = "prs_z_quintile",
  #                                            title = "Received Hormone Therapy",
  #                                            font.main = c(12, "bold"),
  #                                            legend = "right",
  #                                            legend.title = "PRS Quintile",
  #                                            caption = "Adjusted Survival Curves (Age Time Scale)",
  #                                            font.caption = 8,
  #                                            xlab = "Time",
  #                                            ylab = "Survival Probability")
  # facet_curve0 <- survminer::ggadjustedcurves(facetfit,
  #                                             data = as.data.frame(data %>%
  #                                                                    filter(HormoneTherapy == 0)),
  #                                             variable = "prs_z_quintile",
  #                                             title = "No Hormone Therapy",
  #                                             font.main = c(12, "bold"),
  #                                             legend = "right",
  #                                             legend.title = "PRS Quintile",
  #                                             caption = "Adjusted Survival Curves (Age Time Scale)",
  #                                             font.caption = 8,
  #                                             xlab = "Time",
  #                                             ylab = "Survival Probability")
  
  # ggsave(path = paste0(analysis_output_path, "/final_adjKMCurve.png"), curve,
  #        device = "png", width = 5, height = 5)
  
  # forest plot
  # forest <- survminer::ggforest(finalfit,
  #                               data = as.data.frame(data))
  output_fit <- list(res,
                     res2,
                     res3,
                     res4,
                     res5,
                     res6,
                     res7,
                     res8)
  output_files <- list("final",
                         "prs_quintile",
                         "horm",
                         "bmi",
                         "bmi_cat",
                         "smoking",
                         "time",
                         "alc")
  quick_save_csv <- function(obj, filename) {
    write_csv(obj, paste0(analysis_output_path_model, "/final_model_", 
                          filename, ".csv"))
  }
  purrr::map2(output_fit, output_files, quick_save_csv)
}

PerformInteractionModelling <- function(data) {
  fit_chemo <- coxph(Surv(YearsToEntry,
                          YearsToI2025,
                          completeI2025) ~ prs_z*Chemotherapy + DiagnosisAge,
                     data) %>%
    broom::tidy(exponentiate = TRUE)
  fit_radio <- coxph(Surv(YearsToEntry,
                          YearsToI2025,
                          completeI2025) ~ prs_z*Radiotherapy + DiagnosisAge,
                     data) %>%
    broom::tidy(exponentiate = TRUE)
  fit_hormone <- coxph(Surv(YearsToEntry,
                            YearsToI2025,
                            completeI2025) ~ prs_z*HormoneTherapy + DiagnosisAge,
                       data) %>%
    broom::tidy(exponentiate = TRUE)
  fit_hormone_adj <- coxph(Surv(YearsToEntry,
                                YearsToI2025,
                                completeI2025) ~ prs_z*HormoneTherapy + DiagnosisAge +
                             BMI + smokingEver + AlcNow,
                           data) %>%
    broom::tidy(exponentiate = TRUE)
  
  interactionanalysis <- bind_rows(fit_chemo,
                                   fit_radio,
                                   fit_hormone,
                                   fit_hormone_adj)

  write_csv(interactionanalysis, paste0(analysis_output_path,
                                        "/interaction_analysis.csv"))
  
  return(interactionanalysis)
}

PlotAdjKMPRS <- function(data) {
  fit <- coxph(Surv(YearsToEntry,
                    YearsToI2025,
                    completeI2025) ~ prs_z_quintile + DiagnosisAge +
                 HormoneTherapy + log(BMI) + as.factor(smokingEver) +
                 AlcNow + log(IMD),
               data)

  curve <- survminer::ggadjustedcurves(fit,
                                       data = as.data.frame(data),
                                       variable = "prs_z_quintile",
                                       legend = "right",
                                       legend.title = "PRS Quintile",
                                       xlab = "Time",
                                       ylab = "Survival Probability",
                                       xlim = c(0, 20),
                                       ylim = c(0.7, 1))
  
  ggsave(paste0(analysis_output_path_KM, "/final-adj-curve-prs.png"), 
         curve, device = "png", width = 8, height = 7)
}

AnalyzeTimeSplit <- function(data, episode) {
  #time split
  splitdata <- survSplit(Surv(YearsToEntry,
                              YearsToI2025,
                              completeI2025) ~ .,
                         data,
                         cut = c(5,10), 
                         episode = "timegroup")
  
  splitfit <- coxph(Surv(YearsToEntry,
                         YearsToI2025,
                         completeI2025) ~ prs_z + HormoneTherapy + 
                      BMI + smokingEver + AlcNow + DiagnosisAge,
                    splitdata %>%
                      filter(timegroup == episode)) %>%
    broom::tidy(exponentiate = TRUE)
  
  return(splitfit)
}

AnalyzeTimeSplitTherapy <- function(data, var, episode) {
  splitdata <- survSplit(Surv(YearsToEntry,
                              YearsToI2025,
                              completeI2025) ~ .,
                         data,
                         cut = c(5,10), 
                         episode = "timegroup")
  
  surv_formula <- paste0("Surv(YearsToEntry, YearsToI2025, completeI2025) ~ 
                          DiagnosisAge + log(BMI) + smokingEver +",
                         var) %>%
    as.formula()
  
  splitfit <- coxph(surv_formula,
                    splitdata %>%
                      filter(timegroup == episode)) %>%
    broom::tidy(exponentiate = TRUE)
  
  return(splitfit)
}

AnalyzeTimeSplitTherapy2 <- function(data, var, episode) {
  splitdata <- survSplit(Surv(YearsToEntry,
                              YearsToStatus,
                              VitalStatus) ~ .,
                         data,
                         cut = c(5,10), 
                         episode = "timegroup")
  
  surv_formula <- paste0("Surv(YearsToEntry, YearsToStatus, VitalStatus) ~ 
                         DiagnosisAge + log(BMI) + smokingEver +",
                         var) %>%
    as.formula()
  
  splitfit <- coxph(surv_formula,
                    splitdata %>%
                      filter(timegroup == episode)) %>%
    broom::tidy(exponentiate = TRUE)
  
  return(splitfit)
}

PlotTimeSplitTherapy <- function(data){
  #xaxis - time periods
  #yaxis - HR
  plot <- data %>%
    filter(term == "Chemotherapy" |
             term == "Radiotherapy" |
             term == "HormoneTherapy") %>%
    mutate(timeperiod = rep(c("0-5", "5-10", "10+"), 3),
           timeperiod = factor(timeperiod, 
                                  levels = c("0-5", "5-10", "10+"))) %>%
    ggplot(aes(x = timeperiod, y = estimate, col = term)) + 
    geom_point(position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = 0.3),
                  width = 0.1) + 
    theme_bw() + 
    coord_flip() +
    facet_wrap(~ term) +
    labs(x = "Time Period (years)",
         y = "Hazard Ratio",
         col = "Oncotherapy")
  
  ggsave(paste0(analysis_output_path, "/timesplit-therapy-plot.png"), 
         plot, device = "png", width = 8, height = 7)
}

AnalyzeTimeSplitInteract <- function(data) {
  splitdata <- survSplit(Surv(YearsToEntry,
                              YearsToI2025,
                              completeI2025) ~ .,
                         data,
                         cut = c(5,10), 
                         episode = "timegroup")
  
  splitfitinteract <- coxph(Surv(YearsToEntry,
                                 YearsToI2025,
                                 completeI2025) ~ prs_z*strata(timegroup) + HormoneTherapy +
                              log(BMI) + smokingEver + AlcNow + DiagnosisAge +
                              ThyroidDisease + HRTEver,
                            splitdata) %>%
    broom::tidy(exponentiate = TRUE)
  
  write_csv(splitfitinteract, paste0(analysis_output_path, 
                                     "/split_interact_finalmodel.csv"))
  
  return(splitfitinteract)
}

statistics.raplot <- function(x1, x2, y,s1,s2,t) {   
  require(Hmisc)
  require(pROC)
  require(pracma)
  require(caret)
  
  s <- is.na(x1 + x2 + y)  #Remove rows with missing data
  if (any(s)) {
    s <- !s
    x1 <- x1[s]
    x2 <- x2[s]
    y <- y[s]
  }
  n <- length(y)
  y <- as.numeric(y)
  u <- sort(unique(y))
  if (length(u) != 2 || u[1] != 0 || u[2] != 1)
    stop("y must have two values: 0 and 1")
  r <- range(x1, x2)
  if (r[1] < 0 || r[2] > 1)
    stop("x1 and x2 must be in [0,1]")
  if(length(x1)!=length(x2))
    stop("Reference (Null) and New (Alt) model vectors must be the same length")
  incidence<-sum(y)/n
  if (missing(t)) {t<-c(0, incidence,1) }
  a <- y == 1
  b <- y == 0
  na <- sum(a)
  nb <- sum(b)
  d <- x2 - x1
  pev<-na/n
  pne<-nb/n
  # NRI 
  n.thresh<-length(t)-1
  risk.class.x1.ev<-cut2(x1[a],t)
  risk.class.x2.ev<-cut2(x2[a],t)
  thresh<-c()
  for (i in 1:(length(t)-1)){
    ifelse(i==(length(t)-1),thresh[i]<-paste("[",toString(t[i]),",",toString(t[i+1]),"]"), 
           thresh[i]<-paste("[",toString(t[i]),",",toString(t[i+1]),")"))
  }
  levels(risk.class.x1.ev)<-thresh
  levels(risk.class.x2.ev)<-thresh
  cM.ev<-confusionMatrix(risk.class.x2.ev,risk.class.x1.ev)
  pup.ev<-0  # P(up|event)
  pdown.ev<-0 # P(down|event)
  for (i in 1:(n.thresh-1)){pup.ev<-pup.ev+sum(cM.ev$table[(i+1):n.thresh,i])}
  for (i in 2:n.thresh){pdown.ev<-pdown.ev+sum(cM.ev$table[1:(i-1),i])}
  npup.ev<-pup.ev
  pup.ev<-pup.ev/na 
  ndown.ev<-pup.ev
  pdown.ev<-pdown.ev/na
  risk.class.x1.ne<-cut2(x1[b],t)
  risk.class.x2.ne<-cut2(x2[b],t)  
  levels(risk.class.x1.ne)<-thresh
  levels(risk.class.x2.ne)<-thresh
  cM.ne<-confusionMatrix(risk.class.x2.ne,risk.class.x1.ne)
  pup.ne<-0 # P(up|nonevent)
  pdown.ne<-0 # P(down|nonevent)
  for (i in 1:(n.thresh-1)){pup.ne<-pup.ev+sum(cM.ne$table[(i+1):n.thresh,i])}
  for (i in 2:n.thresh){pdown.ne<-pdown.ne+sum(cM.ne$table[1:(i-1),i])} 
  ndown.ne<-pdown.ne
  pdown.ne<-pdown.ne/nb
  npup.ne<-pup.ne
  pup.ne<-pup.ne/nb
  nri <- pup.ev - pdown.ev - (pup.ne - pdown.ne)
  se.nri <- sqrt((pup.ev + pdown.ev)/na + (pup.ne + pdown.ne)/nb)
  z.nri <- nri/se.nri
  nri.ev <- pup.ev - pdown.ev
  se.nri.ev <- sqrt((pup.ev + pdown.ev)/na)
  z.nri.ev <- nri.ev/se.nri.ev
  nri.ne <- pdown.ne - pup.ne
  se.nri.ne <- sqrt((pdown.ne + pup.ne)/nb)
  z.nri.ne <- nri.ne/se.nri.ne 
  # weighted NRI
  pup<-(npup.ev+npup.ne)/n
  pdown<-(ndown.ev+ndown.ne)/n
  pevent.up<- pup.ev * pev / pup        # P(event|up)= P(up|event).P(event)/P(up)
  pevent.dn<- pdown.ev * pev / pdown      # P(event|dn)= P(dn|event).P(event)/P(dn)
  pnonevent.up<- pup.ne * pne / pup      
  pnonevent.dn<- pdown.ne * pne / pdown
  wnri<-s1*(pevent.up*pup - pevent.dn * pdown)+s2*(pnonevent.dn*pdown-pnonevent.up*pup)
  se.wnri<-NA
  z.wnri<-NA
  # Category Free NRI calculations
  cfpup.ev <- mean(d[a] > 0)
  cfpup.ne <- mean(d[b] > 0)
  cfpdown.ev <- mean(d[a] < 0)
  cfpdown.ne <- mean(d[b] < 0)
  cfnri <- cfpup.ev - cfpdown.ev - (cfpup.ne - cfpdown.ne)
  se.cfnri <- sqrt((cfpup.ev + cfpdown.ev)/na + (cfpup.ne + cfpdown.ne)/nb)
  z.cfnri <- cfnri/se.cfnri
  cfnri.ev <- cfpup.ev - cfpdown.ev
  se.cfnri.ev <- sqrt((cfpup.ev + cfpdown.ev)/na)
  z.cfnri.ev <- cfnri.ev/se.cfnri.ev
  cfnri.ne <- cfpdown.ne - cfpup.ne
  se.cfnri.ne <- sqrt((cfpdown.ne + cfpup.ne)/nb)
  z.cfnri.ne <- cfnri.ne/se.cfnri.ne
  # IDI calculations
  improveSens <- sum(d[a])/na
  improveSpec <- -sum(d[b])/nb
  idi.ev <- mean(improveSens)
  idi.ne <- mean(improveSpec)
  idi <- idi.ev + idi.ne
  relidi <- 100*((sum(x2[a])/na - sum(x2[b])/nb)/(sum(x1[a])/na-sum(x1[b])/nb)-1) #relative IDI expressed as a percentage
  se.relidi <-NA
  z.relidi <- NA
  var.ev <- var(d[a])/na
  se.idi.ev <- sqrt(var.ev)
  z.idi.ev <- idi.ev/se.idi.ev
  var.ne <- var(d[b])/nb
  se.idi.ne <- sqrt(var.ne)
  z.idi.ne <- idi.ne/se.idi.ne
  se.idi <- sqrt(var.ev + var.ne)
  z.idi <- idi/se.idi
  # AUC calculations
  roc.x1 <- roc(y, x1)
  auc.x1 <- auc(roc.x1)
  ci.auc.x1 <- ci.auc(roc.x1)
  se.auc.x1 <- (ci.auc.x1[3] - auc.x1)/qnorm(0.975)
  roc.x2 <- roc(y, x2)
  auc.x2 <- auc(roc.x2)
  ci.auc.x2 <- ci.auc(roc.x2)
  se.auc.x2 <- (ci.auc.x2[3] - auc.x2)/qnorm(0.975)
  roc.test.x1.x2 <- roc.test(roc.x1, roc.x2)  #Uses the default Delong method
  sens.x1 <- roc.x1$sensitivities
  spec.x1 <- 1 - roc.x1$specificities
  n.x1 <- length(sens.x1)
  x1 <- roc.x1$thresholds
  x1 <- x1[c(-1,-n.x1)]
  sens.x1 <- sens.x1[c(-1,-n.x1)]
  spec.x1 <- spec.x1[c(-1,-n.x1)]
  sens.x2 <- roc.x2$sensitivities
  spec.x2 <- 1 - roc.x2$specificities
  n.x2 <- length(sens.x2)
  x2 <- roc.x2$thresholds
  x2 <- x2[c(-1,-n.x2)]
  sens.x2 <- sens.x2[c(-1,-n.x2)]
  spec.x2 <- spec.x2[c(-1,-n.x2)]
  # Integrated sensitivity & 1-specificity calculations: Note a 1 and 0 are added to the beginning and end of the sens and spec, and a 0 & 1 to the risks, so that it is the area from 0 to 1 and not just partial.
  is.x1 <- trapz(x = c(0,x1,1), y = c(1,sens.x1,0))  # area under curves (relates to integrated sens, 1-spec)
  is.x2 <- trapz(x = c(0,x2,1), y = c(1,sens.x2,0))
  ip.x1 <- trapz(x = c(0,x1,1), y = c(1,spec.x1,0))
  ip.x2 <- trapz(x = c(0,x2,1), y = c(1,spec.x2,0))
  
  # Output
  output <- c(n, na, nb, pup.ev, pup.ne, pdown.ev, pdown.ne, nri, se.nri, z.nri,
              nri.ev, se.nri.ev, z.nri.ev, nri.ne, se.nri.ne, z.nri.ne, 
              cfpup.ev, cfpup.ne, cfpdown.ev, cfpdown.ne, cfnri, se.cfnri, z.cfnri,
              cfnri.ev, se.cfnri.ev, z.cfnri.ev, cfnri.ne, se.cfnri.ne, z.cfnri.ne, 
              improveSens, improveSpec, 
              s1,s2,wnri,se.wnri,z.wnri,
              idi.ev, se.idi.ev, z.idi.ev, idi.ne, 
              se.idi.ne, z.idi.ne, idi, se.idi, z.idi,relidi, se.relidi, z.relidi, is.x1, NA, is.x2, NA, 
              ip.x1, NA, ip.x2, NA, auc.x1, se.auc.x1, auc.x2, se.auc.x2, 
              roc.test.x1.x2$p.value,incidence)
  names(output) <- c("n", "na", "nb", "pup.ev", "pup.ne", "pdown.ev", "pdown.ne", 
                     "nri", "se.nri", "z.nri", "nri.ev", "se.nri.ev", "z.nri.ev",
                     "nri.ne", "se.nri.ne", "z.nri.ne",
                     "cfpup.ev", "cfpup.ne", "cfpdown.ev", "cfpdown.ne", 
                     "cfnri", "se.cfnri", "z.cfnri", "cfnri.ev", "se.cfnri.ev", "z.cfnri.ev",
                     "cfnri.ne", "se.cfnri.ne", "z.cfnri.ne", "improveSens", "improveSpec",
                     "s1","s2","wnri","se.wnri","z.wnri",
                     "idi.ev", "se.idi.ev", "z.idi.ev", "idi.ne", "se.idi.ne", 
                     "z.idi.ne", "idi", "se.idi", "z.idi","relidi", "se.relidi", "z.relidi", 
                     "is.x1", "se.is.x1",
                     "is.x2", "se.is.x2", "ip.x1", "se.ip.x1", "ip.x2", "se.ip.x2", 
                     "auc.x1", "se.auc.x1", "auc.x2", "se.auc.x2", 
                     "roc.test.x1.x2.pvalue","incidence")
  return(output)
}

CI.raplot <- function(x1, x2, y, s1,s2, t, cis = c("asymptotic", "boot"), conf.level = 0.95, n.boot = 2000, dp = 4) {
  if(length(x1)!=length(x2))
    stop("Reference (Null) and New (Alt) model vectors must be the same length")
  if (missing(s1)){
    s1<-0
    s2<-0
  }
  ifelse(missing(t),
         results <- statistics.raplot(x1, x2,y,s1,s2),  
         results <- statistics.raplot(x1, x2, y, s1,s2,t)
  )
  
  if (cis == "boot") {
    
    results.boot <- matrix(NA, n.boot, length(results))
    
    colnames(results.boot) <- names(results)
    
    for (i in 1:n.boot) {
      boot.index <- sample(length(y), replace = TRUE)
      risk.model1.boot <- x1[boot.index]
      risk.model2.boot <- x2[boot.index]
      cc.status.boot <- y[boot.index]           
      results.boot[i, ] <- statistics.raplot(x1 = risk.model1.boot, 
                                             x2 = risk.model2.boot, 
                                             y = cc.status.boot,
                                             s1,
                                             s2,
                                             t)
    }
    
    results.se.boot <- apply(results.boot, 2, sd) 
    results[grep("se", names(results))] <- results.se.boot[grep("se", names(results)) - 1]
  }
  
  # calculate cis and return 
  z <- abs(qnorm((1 - conf.level)/2))
  
  results.matrix <- matrix(NA, 27, 2)
  
  results.matrix[1, ] <- c("Total (n)", results["n"])
  results.matrix[2, ] <- c("Events (n)", results["na"])
  results.matrix[3, ] <- c("Non-events (n)", results["nb"])
  results.matrix[4, ] <- c("cfNRI and summary statistics","-------------------------")
  results.matrix[5, ] <- c("cfNRI events (%)", 
                           paste(round(100*results["cfnri.ev"], dp-2), " (", 
                                 round(100*results["cfnri.ev"] - z * 100*results["se.cfnri.ev"], dp-2),
                                 " to ", round(100*results["cfnri.ev"] + 
                                                 z * 100*results["se.cfnri.ev"], dp-2), ")", sep = ""))
  results.matrix[6, ] <- c("cfNRI non-events (%)", 
                           paste(round(100*results["cfnri.ne"], dp-2), " (",
                                 round(100*results["cfnri.ne"] - z * 100*results["se.cfnri.ne"], dp-2),
                                 " to ", round(100*results["cfnri.ne"] +  z * 100*results["se.cfnri.ne"], 
                                               dp-2), ")", sep = "")) 
  results.matrix[7, ] <- c("cfNRI (dimensionless)", 
                           paste(round(100*results["cfnri"], dp-2), " (", 
                                 round(100*results["cfnri"] - z * 100*results["se.cfnri"], dp-2), 
                                 " to ", round(100*results["cfnri"] + z * 100*results["se.cfnri"], 
                                               dp-2), ")", sep = ""))
  results.matrix[8, ] <- c("NRI and summary statistics","-------------------------")
  results.matrix[9, ] <- c("NRI events (%)", 
                           paste(round(100*results["nri.ev"], dp-2), " (", 
                                 round(100*results["nri.ev"] - z * 100*results["se.nri.ev"], dp-2),
                                 " to ", round(100*results["nri.ev"] + 
                                                 z * 100*results["se.nri.ev"], dp-2), ")", sep = ""))
  results.matrix[10, ] <- c("NRI non-events (%)", 
                            paste(round(100*results["nri.ne"], dp-2), " (",
                                  round(100*results["nri.ne"] - z * 100*results["se.nri.ne"], dp-2),
                                  " to ", round(100*results["nri.ne"] +  z * 100*results["se.nri.ne"],  dp-2), ")", sep = "")) 
  results.matrix[11, ] <- c("NRI (dimensionless)", 
                            paste(round(100*results["nri"], dp-2), " (", 
                                  round(100*results["nri"] - z * 100*results["se.nri"], dp-2), 
                                  " to ", round(100*results["nri"] + z * 100*results["se.nri"], 
                                                dp-2), ")", sep = ""))
  results.matrix[12, ] <- c("Weighted NRI and summary statistics","-------------------------")
  results.matrix[13, ] <- c("wNRI (dimensionless)", 
                            paste(round(100*results["wnri"], dp-2), " (", 
                                  round(100*results["wnri"] - z * 100*results["se.wnri"], dp-2), 
                                  " to ", round(100*results["wnri"] + z * 100*results["se.wnri"], 
                                                dp-2), ")", sep = ""))
  results.matrix[14, ] <- c("IDI and summary statistics","-------------------------")
  results.matrix[15, ] <- c("IDI events", 
                            paste(round(results["idi.ev"], dp), " (", 
                                  round(results["idi.ev"] - z * results["se.idi.ev"], dp), 
                                  " to ", round(results["idi.ev"] + z * results["se.idi.ev"], 
                                                dp), ")", sep = ""))
  results.matrix[16, ] <- c("IDI non-events", 
                            paste(round(results["idi.ne"], dp), " (", 
                                  round(results["idi.ne"] - z * results["se.idi.ne"], dp), 
                                  " to ", round(results["idi.ne"] + z * results["se.idi.ne"], 
                                                dp), ")", sep = ""))
  results.matrix[17, ] <- c("IDI", 
                            paste(round(results["idi"], dp), " (", 
                                  round(results["idi"] - z * results["se.idi"], dp), 
                                  " to ", round(results["idi"] + z * results["se.idi"], 
                                                dp), ")", sep = ""))
  results.matrix[18, ] <- c("Relative IDI (%)", 
                            paste(round(results["relidi"], dp-2), " (", 
                                  round(results["relidi"] - z * results["se.relidi"], dp-2), 
                                  " to ", round(results["relidi"] + z * results["se.relidi"], 
                                                dp-2), ")", sep = ""))
  results.matrix[19, ] <- c("IS (null model)", 
                            paste(round(results["is.x1"], dp), " (", 
                                  round(results["is.x1"] - z * results["se.is.x1"], dp), 
                                  " to ", round(results["is.x1"] + z * results["se.is.x1"], 
                                                dp), ")", sep = ""))
  results.matrix[20, ] <- c("IS (alt model)", 
                            paste(round(results["is.x2"], dp), " (", 
                                  round(results["is.x2"] - z * results["se.is.x2"], dp), 
                                  " to ", round(results["is.x2"] + z * results["se.is.x2"], 
                                                dp), ")", sep = ""))
  results.matrix[21, ] <- c("IP (null model)", 
                            paste(round(results["ip.x1"], dp), " (", 
                                  round(results["ip.x1"] - z * results["se.ip.x1"], dp), 
                                  " to ", round(results["ip.x1"] + z *  results["se.ip.x1"], 
                                                dp), ")", sep = ""))
  results.matrix[22, ] <- c("IP (alt model)", 
                            paste(round(results["ip.x2"], dp), " (", 
                                  round(results["ip.x2"] - z * results["se.ip.x2"], dp), 
                                  " to ", round(results["ip.x2"] + z * results["se.ip.x2"], 
                                                dp), ")", sep = ""))
  results.matrix[23, ] <- c("AUC","-------------------------")
  results.matrix[24, ] <- c("AUC (null model)", 
                            paste(round(results["auc.x1"], dp), " (", 
                                  round(results["auc.x1"] - z * results["se.auc.x1"], dp), 
                                  " to ", round(results["auc.x1"] + z * results["se.auc.x1"], 
                                                dp), ")", sep = ""))
  results.matrix[25, ] <- c("AUC (alt model)", 
                            paste(round(results["auc.x2"], dp), " (", 
                                  round(results["auc.x2"] - z * results["se.auc.x2"], dp), 
                                  " to ", round(results["auc.x2"] +  z * results["se.auc.x2"], 
                                                dp), ")", sep = ""))
  results.matrix[26, ] <- c("difference (P)", round(results["roc.test.x1.x2.pvalue"], dp))
  results.matrix[27, ] <- c("Incidence", round(results["incidence"], dp))
  
  return(results.matrix)
}

PerformNetReclassification <- function(data) {
  #how many people were censored before and after 10 years
  tab <- data %>%
    mutate(tenyear = case_when(
      YearsToI2025 < 10 ~ 0,
      TRUE ~ 1
    )) %>%
    group_by(tenyear) %>%
    filter(completeI2025 == 0) %>%
    dplyr::summarize(count = n())
  
  data10 <- data %>%
    mutate(tenyearevent = case_when(
      YearsToI2025 <= 10 & completeI2025 == 1 ~ 1,
      TRUE ~ 0 #censor events that happen after 10 years
    )) %>%
    select(tenyearevent,
           prs_z, HormoneTherapy, BMI, smokingEver, AlcNow, DiagnosisAge, ThyroidDisease,
           HRTEver, Parity) %>%
    drop_na()
  
  logistic_model_list <- 
    list(basic_model = glm(tenyearevent ~
                             HormoneTherapy +
                             log(BMI) + as.factor(smokingEver) + AlcNow + DiagnosisAge + ThyroidDisease +
                             HRTEver + Parity, 
                           data = data10, family = "binomial"), 
         new_model   = glm(tenyearevent ~
                             prs_z + HormoneTherapy +
                             log(BMI) + smokingEver + AlcNow + DiagnosisAge + ThyroidDisease +
                             HRTEver, 
                           data = data10, family = "binomial"))

  output <- PredictABEL::reclassification(data = as.data.frame(data10), 
                                          cOutcome = 1, 
                   predrisk1 = fitted(logistic_model_list[["basic_model"]]), 
                   predrisk2 = fitted(logistic_model_list[["new_model"]]), 
                   cutoff = c(0, 0.05, 0.1, 1))
  
  output2 <- CI.raplot(fitted(logistic_model_list[["basic_model"]]),
                       fitted(logistic_model_list[["new_model"]]),
                       data10$tenyearevent,
                       s1 = 10,
                       s2 = 1/(0.9),
                       t = c(0, 0.05, 0.1, 1),
                       conf.level = 0.95,
                       n.boot = 4,
                       dp = 4)
  
  return(list("tab" = tab,
              "output" = output,
              "output2" = output2))
}

PerformCompetingRiskAnalysis <- function(data) {

  df <- data.frame(id = data$PersonID,
                   timestart = data$YearsToEntry,
                   timestop = data$time,
                   status = data$censor,
                   event = as.factor(data$competingrisk))
  
  fit <- survfit(Surv(timestart,
                      timestop, 
                      status) ~ 1,
                 data = df, 
                 id = id,
                 etype = event)
  
  #red = incident CAD event
  #green = death from other cause
  compriskplot <- plot(fit, 
                       fun = 'event', 
                       xscale = 1, 
                       xmax = 20,
                       mark.time=FALSE,
                       col = 2:3, 
                       xlab = "Years post-breast cancer diagnosis",
                       ylab = "Proportion")
  
  return(compriskplot)
}

ObtainUniConcordance <- function(data, var) {
  data <- data %>%
    drop_na(c("HormoneTherapy",
              "BMI",
              "smokingEver",
              "AlcNow",
              "ThyroidDisease",
              "Parity",
              "HRTEver",
              "prs_z",
              "DiagnosisAge"))
  
  formula <- paste0("Surv(YearsToEntry, YearsToI2025, completeI2025) ~ ",
                    var)
  fit <- coxph(as.formula(formula),
               data) %>%
    broom::glance()
  
  return(fit)
}

ObtainConcordance <- function(data) {
  #model 1
  fit <- coxph(Surv(YearsToEntry,
                    YearsToI2025,
                    completeI2025) ~ DiagnosisAge,
               data %>% 
                 drop_na(prs_z, DiagnosisAge)) 
  res <- fit %>%
    broom::glance()
  
  fit1a <- coxph(Surv(YearsToEntry,
                      YearsToI2025,
                      completeI2025) ~ prs_z,
                 data %>% 
                   drop_na(prs_z, DiagnosisAge)) 
  res1a <- fit1a %>%
    broom::glance()
  
  fit2 <- coxph(Surv(YearsToEntry,
                     YearsToI2025,
                     completeI2025) ~ prs_z + DiagnosisAge,
                data %>%
                  drop_na(prs_z, DiagnosisAge)) 
  res2 <- fit2 %>%
    broom::glance()
  
  #model 2
  fit3 <- coxph(Surv(YearsToEntry,
                     YearsToI2025,
                     completeI2025) ~ DiagnosisAge + 
                  log(BMI) + as.factor(smokingEver),
                data %>%
                  drop_na(DiagnosisAge, BMI, smokingEver, prs_z)) 
  res3 <- fit3 %>%
    broom::glance()
  
  fit3a <- coxph(Surv(YearsToEntry,
                      YearsToI2025,
                      completeI2025) ~ prs_z,
                 data %>%
                   drop_na(DiagnosisAge, BMI, smokingEver, prs_z)) 
  res3a <- fit3a %>%
    broom::glance()
  
  fit4 <- coxph(Surv(YearsToEntry,
                     YearsToI2025,
                     completeI2025) ~ prs_z + DiagnosisAge + 
                  log(BMI) + as.factor(smokingEver),
                data %>%
                  drop_na(DiagnosisAge, BMI, smokingEver, prs_z)) 
  res4 <- fit4 %>%
    broom::glance()
  
  #model 3
  fit5 <-  coxph(Surv(YearsToEntry,
                      YearsToI2025,
                      completeI2025) ~ HormoneTherapy +
                   log(BMI) + as.factor(smokingEver) + AlcNow + DiagnosisAge +
                   log(IMD),
                 data %>%
                   drop_na(HormoneTherapy,
                           BMI,
                           smokingEver,
                           AlcNow,
                           DiagnosisAge,
                           IMD,
                           prs_z)) 
  res5 <- fit5 %>%
    broom::glance()
  
  fit5a <- coxph(Surv(YearsToEntry,
                      YearsToI2025,
                      completeI2025) ~ prs_z,
                 data %>%
                   drop_na(HormoneTherapy,
                           BMI,
                           smokingEver,
                           AlcNow,
                           DiagnosisAge,
                           IMD,
                           prs_z)) 
  res5a <- fit5a %>%
    broom::glance()
  
  fit6 <-  coxph(Surv(YearsToEntry,
                      YearsToI2025,
                      completeI2025) ~ prs_z + HormoneTherapy +
                   log(BMI) + as.factor(smokingEver) + AlcNow + DiagnosisAge +
                   log(IMD),
                 data %>%
                   drop_na(HormoneTherapy,
                           BMI,
                           smokingEver,
                           AlcNow,
                           DiagnosisAge,
                           IMD,
                           prs_z)) 
  res6 <- fit6 %>%
    broom::glance()
  
  out <- bind_rows(res, res1a, res2, res3, res3a, res4, res5, res5a, res6)
  
  return(out)
}

PlotPRSCumIncCurves <- function(data){
  browser()
  data <- data %>%
    mutate(agescale = YearsToI2025 - YearsToEntry + DiagnosisAge)
  
  ci <- cmprsk::cuminc(
    ftime = data$agescale,
    fstatus = data$completeI2025,
    group = data$prs_z_quintile,
    cencode = 0
  )
  
  ciplotdat <- ci %>%
    list_modify("Tests" = NULL) %>%
    map_dfr(`[`, c("time", "est"), .id = "id") %>%
    mutate(PRS = case_when(
      id == "(-0.253,0.255] 1" ~ 3,
      id == "(-0.828,-0.253] 1" ~ 2,
      id == "(0.255,0.846] 1" ~ 4,
      id == "[-4,-0.828] 1" ~ 1,
      id == "(0.846,3.82] 1" ~ 5
    ))
  
  plot <- ciplotdat %>%
    filter(PRS != 2,
           PRS != 4) %>%
  ggplot(aes(x = time, y = est, col = as.factor(PRS))) + 
    geom_step(lwd = 1.2) + 
    theme_bw() + 
    theme(plot.title = element_text(size = 14),
          legend.position = "bottom") + 
    labs(x = "Age (Years)",
         y = "Cumulative Incidence of CAD") + 
    scale_color_discrete(name = "Quintiles", labels = c("Low (0-20%)",
                                                        "Medium (40-60%)",
                                                        "High (80-100%)")) + 
    scale_x_continuous(limits = c(50, 85)) + 
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                       limits = c(0, 0.5))
  
  ggsave(paste0(analysis_output_path, "/cuminc_curves_prs.png"), 
         plot, device = "png", width = 8, height = 7)
  return(plot)
}

DoAllAnalysis <- function(use_existing_data = T) {
  if (!use_existing_data) {
    data <- MakeFinalDataset()
  } else {
    data <- read_csv(final_data_path)
  }
  
  data_quintiles <- data 
  
  #Plot distributions of covariates
  PlotDistributions(data = data_quintiles)
  
  # Assess Correlation
  AssessCorrelation(data = data_quintiles)
  
  # Make Baseline Tables
  descriptive_data <- PrepareDescriptiveData(data_quintiles)
  baseline_table_vars <- c("groupcause",
                           "completeI20",
                           "completeI2025",
                           "completeI2123",
                           "completeI2425",
                           "completeI2125",
                           "genotyped",
                           "prs_z_quintile")
  map(baseline_table_vars, MakeBaselineTable, data = descriptive_data)
  print("Saved baseline tables!")
  
  missingplot_list <- PerformMissignessAnalysis(data = data_quintiles)
  final_missingplot <- (missingplot_list[[2]] + missingplot_list[[3]]) /
    (missingplot_list[[4]] + missingplot_list[[5]]) /
    (missingplot_list[[6]] + missingplot_list[[7]]) /
    (missingplot_list[[8]] + patchwork::plot_spacer())
  
  # Plot Crude all-cause mortality curves 
  ACM_vars <- c("Chemotherapy",
                "Radiotherapy",
                "HormoneTherapy")
  ACM_var_labels <- c("Chemotherapy",
                      "Radiotherapy",
                      "Hormone Therapy")
  ACM_var_sublabels <- list(chemo = c("Did not receive", "Received"),
                           radio = c("Did not receive", "Received"),
                           horm = c("Did not receive", "Received"))
  ACM_list <- list(ACM_vars, ACM_var_labels, ACM_var_sublabels)
  pmap(ACM_list, PlotACMKMCurves, data = data_quintiles)
  # 
  # Plot Undjusted KM Curves
  KM_vars <- c("Chemotherapy",
               "Radiotherapy",
               "HormoneTherapy",
               "smokingEver",
               "bmi_cat",
               "prs_z_quintile")
  KM_var_labels <- c("Chemotherapy",
                     "Radiotherapy",
                     "Hormone Therapy",
                     "Smoking",
                     "BMI Categories",
                     "Standardized PRS \n Quintiles")
  KM_var_sublabels <- list(chemo = c("Did not receive", "Received"),
                           radio = c("Did not receive", "Received"),
                           horm = c("Did not receive", "Received"),
                           smoking = c("Never-smoker", "Ex-smoker", "Current smoker"),
                           bmi = c("Underweight", "Normal weight", "Overweight", "Obese"),
                           prs = c("Quintile 1", "Quintile 2", "Quintile 3", "Quintile 4", "Quintile 5"))
  KM_list <- list(KM_vars, KM_var_labels, KM_var_sublabels)
  pmap(KM_list, ExploreKM, data = data_quintiles)
  print("Saved Kaplan-Meier curves!")
  
  # Plot PRS Density
  PlotPRSDensity(data = data_quintiles)
  
  relevant_outcomes <- c("I2125", "I2123", "I2425", "I2025")
  model_output <- map_dfr(relevant_outcomes, FitCoxSensitivity, data = data)
  model_output <- model_output %>%
    mutate(outcome = relevant_outcomes)
  write_csv(model_output, paste0(analysis_output_path, "/PRS_models.csv"))
  
  # PlotCoxSensitivity(model_output)
  
  # FitCoxControls(data = data, outcome = "bc")
  # FitCoxControls(data = data, outcome = "VitalStatus")
  
  ltdep_output <- AssessLTDep(data)
  
  # Plot dose response
  quintile_list <- list(
    quintilevar = c("logbmi_quintile", "age_quintile", "prs_z_quintile"),
    var = c("log(BMI)", "DiagnosisAge", "prs_z"),
    label = c("Log(BMI)", "Age at Diagnosis", "Standardized PRS"),
    width = c(0.01, 0.4, 0.05)
  )
  
  dose_list <- quintile_list %>% pmap(AssessQuintileDoseResponse, 
                                      data = data_quintiles)
  dr_plot <- dose_list[[1]] / dose_list[[2]] / dose_list[[3]]
  ggsave(paste0(analysis_output_path_exploratory, "/", "dose-response-cad.png"), 
         dr_plot, device = "png", width = 7, height = 10)
  
  AnalyzeSensitivityExp(data_quintiles)$lrt
  AnalyzeSensitivityExp(data_quintiles)$smoking_ord
  AnalyzeSensitivityExp(data_quintiles)$smoking_cat
  
  print("Finished EDA!")
  
  # Perform univariate modelling
  univariate <- c("DiagnosisAge",
                  "AgeMenarche",
                  "log(IMD)",
                  "ThyroidDisease",
                  "Chemotherapy",
                  "Radiotherapy",
                  "HormoneTherapy",
                  "log(BMI)",
                  "as.factor(smokingEver)",
                  "AlcNow",
                  "Parity",
                  "HRTEver")
  unimodel <- univariate %>% map_dfr(PerformUnivariateModelling, 
                                     data = data_quintiles)
  
  # Plot control dose response
  control_list <- list(
    modelterms = c("prs_z_quintile", "logbmi_quintile", "age_quintile"),
    quintilevar = c("prs_z_quintile", "logbmi_quintile", "age_quintile"),
    var = c("prs_z","log(BMI)", "DiagnosisAge"),
    label = c("PRS","Log(BMI)", "Age at Diagnosis"),
    width = c(0.1, 0.01, 0.4))
  
  allcausedose_list <- control_list %>% pmap(AssessControlDoseResponse, 
                                             outcome = "all",
                                             data = data_quintiles)
  allcause_plot <- allcausedose_list[[1]] / allcausedose_list[[2]] /
    allcausedose_list[[3]]
  ggsave(paste0(analysis_output_path_exploratory, "/dose-response-acm.png"), 
         allcause_plot, device = "png", width = 7, height = 10)
  
  # Perform Sequential Modelling
  PerformSequentialModelling(data = data_quintiles)
  print("Finished sequential modeling!")
  
  # Perform Interaction Modelling
  # PerformInteractionModelling(data = data_quintiles)
  
  #Fit final model
  AnalyzeFinalModel(data = data_quintiles)
  print("Fit final model and explored interactions!")
  
  # bccausedose_list <- control_list %>% pmap(AssessControlDoseResponse, 
  #                                            outcome = "bc",
  #                                            data = data_quintiles)
  # bccause_plot <- bccausedose_list[[1]] / (bccausedose_list[[2]] + 
  #                                            bccausedose_list[[3]])
  PlotAdjKMPRS(data = data_quintiles)
  
  # episodes <- c(1,2,3)
  # timesplit_output <- map_dfr(episodes, AnalyzeTimeSplit, data = data_quintiles)
  # write_csv(timesplit_output, paste0(analysis_output_path, "/split_finalmodel.csv"))
  #
  # therapies <- c("Chemotherapy", "Radiotherapy", "HormoneTherapy")
  # 
  # time_therapy <- expand_grid(therapies, episodes)
  # timesplit_therapy_output <- map2_dfr(time_therapy$therapies,
  #                                      time_therapy$episodes,
  #                                      AnalyzeTimeSplitTherapy,
  #                                      data = data_quintiles)
  # timesplit_therapy_output2 <- map2_dfr(time_therapy$therapies,
  #                                      time_therapy$episodes,
  #                                      AnalyzeTimeSplitTherapy2,
  #                                      data = data_quintiles)
 
  # write_csv(timesplit_therapy_output, 
  #           paste0(analysis_output_path, "/split_therapy_finalmodel.csv"))
  # Plot time split therapy hazard ratios
  # PlotTimeSplitTherapy(data = timesplit_therapy_output2)
  # 
  # AnalyzeTimeSplitInteract(data = data_quintiles)
  
  PerformNetReclassification(data = data_quintiles)$tab
  PerformNetReclassification(data = data_quintiles)$output
  PerformNetReclassification(data = data_quintiles)$output2
  
  PerformCompetingRiskAnalysis(data = data_quintiles)
  
  univariate_C <- c("prs_z",
                    "DiagnosisAge",
                    "ThyroidDisease",
                    "HormoneTherapy",
                    "log(BMI)",
                    "as.factor(smokingEver)",
                    "AlcNow",
                    "Parity",
                    "HRTEver")
  univariate_C_res <- map_dfr(univariate_C, ObtainUniConcordance,
                              data = data_quintiles)
  
  cstat <- ObtainConcordance(data = data_quintiles) %>% 
    mutate(lower = concordance - 1.96 * std.error.concordance,
           upper = concordance + 1.96 * std.error.concordance,
           change = concordance - lag(concordance))
  
  PlotPRSCumIncCurves(data = data_quintiles)
  
  fin_res <- list("data" = data_quintiles,
                  "uni" = unimodel,
                  "cstat" = cstat)
  
  return(fin_res)
}
DoAllAnalysis()
