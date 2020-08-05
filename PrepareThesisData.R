# This script prepares all the data for Lathan Liou's MPhil thesis
library(tidyverse)

#### HIDE PATHS ####
data_path <- "/Users/lathanliou/Desktop/Academic/Cambridge/Thesis/Data"
setwd(data_path)

#onco refers to oncoarray
#cogs refers to cogs genotyping array

lookup_path <- "lathan_cvd/sample_lookup.xlsx"
onco_sample_path <- "lathan_cvd/lathan_cvd_euro_onco_sample_order.txt"
onco_variants_path <- "lathan_cvd/lathan_cvd_euro_onco_imputed_variants_info.txt"
onco_genotypes_path <- "lathan_cvd/lathan_cvd_euro_onco_imputed_dosages_v2a.txt"
cogs_sample_path <- "lathan_cvd/lathan_cvd_euro_icogs_sample_order.txt"
cogs_variants_path <- "lathan_cvd/lathan_cvd_euro_icogs_imputed_variants_info.txt"
cogs_genotypes_path <- "lathan_cvd/lathan_cvd_euro_icogs_imputed_dosages_v2a.txt"
weights_path <- "Abraham2016PGS.txt"
output_path <- "FinalDataPrep"
genotype_output_path <- "genotypes.RData"

#search paths, data was queried in batches of different variables
search1_path <- "SEARCH/SetData_20200321.xlsx"
search2_path <- "SEARCH/SetData_BCAC_Risk factors_20200403.xlsx"
search3_path <- "SEARCH/SetData_BCAC_Risk factors_20200430.xlsx"
search4_path <- "SEARCH/SetData_20200430_1.xlsx"
IMD_path <- "SEARCH/20200727_IMD.xlsx"
hes_path <- "20200610_HES IP operated.xlsx"

#### RUN THIS ####
#dplyr::select SNPs to be used to calculate PRS
PrepareSNPs <- function(LOOKUP_PATH = lookup_path,
                        SHEET,
                        SAMPLE_PATH,
                        VARIANT_PATH,
                        WEIGHTS_PATH = weights_path) {
  weights <- read_tsv(WEIGHTS_PATH, skip=9) %>%
    dplyr::rename(rsid = rsID)
  sample_lookup <- readxl::read_xlsx(LOOKUP_PATH, sheet = SHEET)
  sample_order <- readr::read_delim(SAMPLE_PATH,
                             delim = " ", 
                             col_name = "BCAC_ID") %>% left_join(sample_lookup)
  variants_df <- read_delim(VARIANT_PATH, delim = " ")
  
  variants_cleaned_df <- variants_df %>%
    mutate(rsid = str_split(`1000g_identifier`, ":", simplify=TRUE)[,1],
           palindrome = case_when(a0 == "A" & a1=="T" ~ 1,
                                  a0 == "T" & a1=="A" ~ 1,
                                  a0 == "C" & a1=="G" ~ 1,
                                  a0 == "G" & a1=="C" ~ 1,
                                  TRUE ~ 0)) %>%
    dplyr::select(chr, position, rsid, a0, a1, exp_freq, info, certainty, palindrome) %>%
    # Exclude SNPs with poor imputation
    filter(info > 0.3) %>%
    # For multi-allelic SNPs dplyr::select the most common
    arrange(desc(exp_freq)) %>%
    group_by(rsid) %>%
    dplyr::slice(1) %>%
    # Join to files of SNP weights
    ungroup() %>% left_join(weights, "rsid") %>%
    filter(!is.na(effect_allele)) %>%
    # Generate indicator for palindromic SNPs
    mutate(weight = case_when(a1 == effect_allele & palindrome == 0 ~ effect_weight,
                              a0 == effect_allele & palindrome == 0 ~ -effect_weight,
                              TRUE ~ 0))
  
  return(variants_cleaned_df)
}
PrepareGenotypes <- function(VARIANT_PATH,
                             GENOTYPE_PATH) {
  variants_df <- read_delim(VARIANT_PATH, delim = " ")
  variant_names <- str_split(variants_df$`1000g_identifier`, ":", simplify=TRUE)[,1]
  genotype_df <- read_tsv(GENOTYPE_PATH, col_names = c(variant_names, "X1"))
  return(genotype_df)
}
CalculatePRS <- function(genotype_df, weights, 
                         LOOKUP_PATH = lookup_path, 
                         SHEET, 
                         SAMPLE_PATH) {
  sample_lookup <- readxl::read_xlsx(LOOKUP_PATH, sheet = SHEET)
  sample_order <- readr::read_delim(SAMPLE_PATH,
                                    delim = " ", 
                                    col_name = "BCAC_ID") %>% left_join(sample_lookup)
    #dosage x weight patient by patient
    # genotype_df 8707 x 44188
    # weights : 44188 x 1
    # out: 8707 x 1
    weights_matrix <- matrix(weights, ncol = 1)
    
    prs <- as.matrix(genotype_df) %*% weights_matrix
    
    prs_final <- sample_order %>%
      mutate(prs = prs)
    
  return(prs_final)
}
FinalizePRS <- function(OUTPUT_PATH = output_path,
                        GENOTYPE_OUTPUT_PATH = genotype_output_path,
                        use_existing_data = T,
                        overwrite_data = F) {
  # Read in variant data
  onco_snps <- PrepareSNPs(SHEET = "Oncoarray",
                           SAMPLE_PATH = onco_sample_path,
                           VARIANT_PATH = onco_variants_path)
  cogs_snps <- PrepareSNPs(SHEET = "Cogs",
                           SAMPLE_PATH = cogs_sample_path,
                           VARIANT_PATH = cogs_variants_path)
  
  # Filter so that onco and cogs use same set of SNPs
  snps_in_common <- intersect(onco_snps$rsid, cogs_snps$rsid) 
  snps_in_common_df <- tibble(rsid = snps_in_common)
  
  onco_snps_filtered <- onco_snps %>%
    semi_join(snps_in_common_df, "rsid")
  cogs_snps_filtered <- cogs_snps %>% 
    semi_join(snps_in_common_df, "rsid")
  
  # Read in genotype data
  if (!use_existing_data) {
    onco_patients <- PrepareGenotypes(onco_variants_path,
                                      onco_genotypes_path)
    cogs_patients <- PrepareGenotypes(cogs_variants_path,
                                      cogs_genotypes_path)
  } else {
    load(GENOTYPE_OUTPUT_PATH)
  }
  
  if (overwrite_data) {
    save(onco_patients, cogs_patients, file = GENOTYPE_OUTPUT_PATH)
  }
  
  print("All data loaded!")
  
  #only use genotype dosages that have a matching variant
  onco_patients_filtered <- onco_patients %>%
    dplyr::select(onco_snps_filtered$rsid)
  cogs_patients_filtered <- cogs_patients %>%
    dplyr::select(cogs_snps_filtered$rsid)
  
  # Calculate PRS
  onco_prs <- CalculatePRS(onco_patients_filtered, onco_snps_filtered$weight,
                           SHEET = "Oncoarray",
                           SAMPLE_PATH = onco_sample_path)
  cogs_prs <- CalculatePRS(cogs_patients_filtered, cogs_snps_filtered$weight,
                           SHEET = "Cogs",
                           SAMPLE_PATH = cogs_sample_path)
  
  full_prs <- bind_rows(onco_prs, cogs_prs) %>%
    mutate(prs_z = (prs-mean(prs, na.rm = TRUE))/sd(prs, na.rm = TRUE))
  
  print("PRS has been calculated!")
  
  write_csv(full_prs, paste0(OUTPUT_PATH, "/full_prs.csv", sep = ""))
  
  return(full_prs)
}
MakeSearchData <- function(SEARCH1_PATH = search1_path,
                           SEARCH2_PATH = search2_path,
                           SEARCH3_PATH = search3_path,
                           SEARCH4_PATH = search4_path,
                           IMD_PATH = IMD_path) {
  search_part1 <- readxl::read_excel(SEARCH1_PATH)
  search_part2 <- readxl::read_excel(SEARCH2_PATH)
  search_part3 <- readxl::read_excel(SEARCH3_PATH)
  search_part4 <- readxl::read_excel(SEARCH4_PATH)
  imd <- readxl::read_excel(IMD_PATH)
  
  clean_search_part1 <- search_part1 %>%
    dplyr::select(PersonID, YearsToEntry, YearsToStatus, DiagnosisAge, VitalStatus, Chemotherapy, Radiotherapy, HormoneTherapy,
           TumourGrade, TumourSize, Nodes_excised, Nodes_involved, ER_Status, DetectedByScreeningMammogram, 
           starts_with("CauseDeath")) %>%
    mutate(DiagnosisAge = as.integer(DiagnosisAge),
           VitalStatus = as.factor(VitalStatus),
           Chemotherapy = as.factor(Chemotherapy),
           Radiotherapy = as.factor(Radiotherapy),
           HormoneTherapy = as.factor(HormoneTherapy),
           Nodes_excised = as.integer(Nodes_excised),
           Nodes_involved = as.integer(Nodes_involved))
  
  clean_search_part2 <- search_part2 %>%
    dplyr::select(PersonID, AgeDiagIndex, eduCat, AgeMenarche, Parity, height, weight, BMI, HRTEver, ends_with("Now"),
           smokingEver) %>%
    mutate(AlcNow = case_when(
      AlcBeerNow > 0 & AlcBeerNow != 888 |
        AlcWineNow > 0 & AlcWineNow != 888 |
        AlcFWineNow > 0 & AlcFWineNow != 888 |
        AlcSpiritsNow > 0 & AlcSpiritsNow != 888 ~ 1,
      TRUE ~ 0
    ),
    AgeDiagIndex = as.integer(AgeDiagIndex),
    eduCat = as.integer(eduCat),
    AgeMenarche = as.integer(AgeMenarche),
    Parity = as.integer(Parity),
    BMI = as.numeric(BMI),
    HRTEver = as.factor(HRTEver),
    smokingEver = as.factor(smokingEver),
    AlcNow = as.factor(AlcNow)) %>%
    dplyr::select(-c(AlcBeerNow:AlcSpiritsNow)) %>%
    mutate(Parity = replace_na(Parity, 0),
           AgeDiagIndex = na_if(AgeDiagIndex, 888),
           eduCat = na_if(eduCat, 888),
           AgeMenarche = na_if(AgeMenarche, 888),
           Parity = na_if(Parity, 888),
           HRTEver = na_if(HRTEver, 888),
           smokingEver = na_if(smokingEver, 888))
  
  clean_search_part3 <- search_part3 %>%
    dplyr::select(PersonID, EthnicityClass) %>%
    mutate(EthnicityClass = na_if(EthnicityClass, 888))
  
  full_search_df <- clean_search_part1 %>%
    left_join(clean_search_part2, by = "PersonID") %>%
    left_join(clean_search_part3, by = "PersonID") %>%
    left_join(search_part4, by = "PersonID") %>%
    left_join(imd, by = "PersonID")
  
  #add censor indicator variables
  final_search_df <- full_search_df %>%
    mutate(I2125 = case_when(
      str_detect(CauseDeath1a, "I2[1-5]") |
        str_detect(CauseDeath1b, "I2[1-5]") |
        str_detect(CauseDeath1c, "I2[1-5]") ~ 1,
      TRUE ~ 0
    ),
    I2123 = case_when(
      str_detect(CauseDeath1a, "I2[1-3]") |
        str_detect(CauseDeath1b, "I2[1-3]") |
        str_detect(CauseDeath1c, "I2[1-3]") ~ 1,
      TRUE ~ 0
    ),
    I2425 = case_when(
      str_detect(CauseDeath1a, "I2[4-5]") |
        str_detect(CauseDeath1b, "I2[4-5]") |
        str_detect(CauseDeath1c, "I2[4-5]") ~ 1,
      TRUE ~ 0
    ),
    I20 = case_when(
      str_detect(CauseDeath1a, "I20") |
        str_detect(CauseDeath1b, "I20") |
        str_detect(CauseDeath1c, "I20") ~ 1,
      TRUE ~ 0
    ),
    I2025 = case_when(
      str_detect(CauseDeath1a, "I2[0-5]") |
        str_detect(CauseDeath1b, "I2[0-5]") |
        str_detect(CauseDeath1c, "I2[0-5]") ~ 1,
      TRUE ~ 0
    ),
    bc = case_when(
      str_detect(CauseDeath1a, "C5") |
        str_detect(CauseDeath1b, "C5") |
        str_detect(CauseDeath1c, "C5") ~ 1,
      TRUE ~ 0
    ))
  return(final_search_df)
}
tidy_hes_data <- function(HES_PATH = hes_path, 
                          string, timevariable, variable) {
  hes <- readxl::read_xlsx(HES_PATH)
  tidy_data <- hes %>%
    mutate(variable := case_when(str_detect(DIAG_CONCAT, string) ~ 1,
                                 TRUE ~ 0)) %>%
    filter(DaysFromDiagToOpertn > 0,
           variable == 1) %>%
    group_by(PersonID) %>%
    dplyr::summarize(YearsToEvent = min(DaysFromDiagToOpertn)/365.25) %>%
    ungroup() %>%
    mutate(!!timevariable := as.numeric(str_trim(YearsToEvent)),
           !!variable := 1) %>%
    dplyr::select(-YearsToEvent)
  
  return(tidy_data)
}
MakeHESData <- function(){
  hes_I20 <- tidy_hes_data(string = "I20", 
                           timevariable = "YearsToI20",
                           variable = "hesI20")
  hes_I2125 <- tidy_hes_data(string = "I2[1-5]",
                             timevariable = "YearsToI2125",
                             variable = "hesI2125")
  hes_I2123 <- tidy_hes_data(string = "I2[1-3]",
                             timevariable = "YearsToI2123",
                             variable = "hesI2123")
  hes_I2425 <- tidy_hes_data(string = "I2[4-5]",
                             timevariable = "YearsToI2425",
                             variable = "hesI2425")
  hes_I2025 <- tidy_hes_data(string = "I2[0-5]",
                             timevariable = "YearsToI2025",
                             variable = "hesI2025")
  
  output <- list("I20" = hes_I20, 
                 "I2125" = hes_I2125,
                 "I2123" = hes_I2123, 
                 "I2425" = hes_I2425, 
                 "I2025" = hes_I2025)
  return(output)
}
MakeFinalDataset <- function() {
  PRS <- FinalizePRS()
  print("PRS has been finalized!")
  search <- MakeSearchData()
  print("SEARCH has been finalized!")
  hes_list <- MakeHESData()
  print("HES has been finalized!")
  
  final_thesis_data <- search %>% 
    left_join(PRS %>% dplyr::select(PersonID, prs, prs_z), by = "PersonID") %>%
    left_join(hes_list$`I20`) %>%
    left_join(hes_list$`I2125`) %>%
    left_join(hes_list$`I2123`) %>%
    left_join(hes_list$`I2425`) %>%
    left_join(hes_list$`I2025`) %>%
    replace_na(., list(hes120 = 0,
                       hesI2125 = 0,
                       hesI2123 = 0,
                       hesI2425 = 0,
                       hesI2025 = 0)) %>%
    dplyr::select(PersonID, YearsToEntry, YearsToStatus, YearsToI20,
           YearsToI2125, YearsToI2123, YearsToI2425, YearsToI2025, everything()) %>%
    mutate(YearsToI20 = case_when(
      is.na(YearsToI20) | YearsToStatus < YearsToI20 ~ YearsToStatus,
      TRUE ~ YearsToI20
    ),
    YearsToI2125 = case_when(
      is.na(YearsToI2125) | YearsToStatus < YearsToI2125 ~ YearsToStatus,
      TRUE ~ YearsToI2125
    ),
    YearsToI2123 = case_when(
      is.na(YearsToI2123) | YearsToStatus < YearsToI2123 ~ YearsToStatus,
      TRUE ~ YearsToI2123
    ),
    YearsToI2425 = case_when(
      is.na(YearsToI2425) | YearsToStatus < YearsToI2425 ~ YearsToStatus,
      TRUE ~ YearsToI2425
    ),
    YearsToI2025 = case_when(
      is.na(YearsToI2025) | YearsToStatus < YearsToI2025 ~ YearsToStatus,
      TRUE ~ YearsToI2025
    ),
    completeI20 = case_when(
      I20 == 1 | hesI20 == 1 ~ 1,
      TRUE ~ 0
    ),
    completeI2125 = case_when(
      I2125 == 1 | hesI2125 == 1 ~ 1,
      TRUE ~ 0
    ),
    completeI2123 = case_when(
      I2123 == 1 | hesI2123 == 1 ~ 1,
      TRUE ~ 0
    ),
    completeI2425 = case_when(
      I2425 == 1 | hesI2425 == 1 ~ 1,
      TRUE ~ 0
    ),
    completeI2025 = case_when(
      I2025 == 1 | hesI2025 == 1 ~ 1,
      TRUE ~ 0
    )) %>%
    dplyr::select(-starts_with("I2"), -starts_with("hes"))
  
  #further clean
  final_thesis_data <- final_thesis_data %>% 
    replace_na(., list(Radiotherapy = 0,
                     Chemotherapy = 0,
                     HormoneTherapy = 0)) %>%
    mutate(prs_z_quintile = cut(prs_z, quantile(prs_z, 
                                                probs = 0:5/5,
                                                na.rm = TRUE),
                                include.lowest = TRUE),
           bmi_z = (BMI - mean(BMI, na.rm = TRUE))/sd(BMI, na.rm = TRUE),
           bmi_cat = case_when(
             BMI < 18.5 ~ 1,
             BMI >= 18.5 & BMI < 25 ~ 2,
             BMI >= 25 & BMI < 30 ~ 3,
             TRUE ~ 4
           ),
           competingrisk = case_when(
             completeI2025 == 0 & VitalStatus == 0 ~ 0,
             completeI2025 == 1 ~ 1,
             VitalStatus == 1 & completeI2025 == 0 ~ 2,
             TRUE ~ 0
           ),
           censor = case_when(
             competingrisk == 0 ~ 0,
             TRUE ~ 1
           ),
           time = case_when(
             competingrisk == 0 ~ YearsToStatus,
             competingrisk == 1 ~ YearsToI2025,
             competingrisk == 2 ~ YearsToStatus,
             TRUE ~ YearsToStatus
           ),
           IMD = as.numeric(IMD))
  write_csv(final_thesis_data, paste0(output_path, "/final_thesis_data.csv", sep = ""))
  
  return(final_thesis_data)
}
