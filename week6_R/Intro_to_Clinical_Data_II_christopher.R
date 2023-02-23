setwd("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/analysis_data")

if(!require(survival)){
  install.packages("survival")
}
library(survival)
if(!require(survminer)){
  install.packages("survminer")
}
library(survminer)
if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)

clin_query <- GDCquery(project = "TCGA-BRCA",data.category = "Clinical",file.type = "xml")
GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query,
                            clinical.info = "patient")
clinical <- read.csv("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/analysis_data/brca_clinical_data.csv")
clinical_drug <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query,
                                  clinical.info = "radiation")
sum(is.na(clinical$age_at_initial_pathologic_diagnosis))
# 1. age_at_initial_pathologic_diagnosis
# 2. Discrete
sum(is.na(clinical_drug$days_to_drug_therapy_start))
# 3. days_to_drug_therapy_start    This variable is the days since the diagnosis to the first day of drug use
# 4. Discrete
# 5. The age and days to drug therapy start may have a positive correlation or no correlation
# with each other since the younger patients may want to begin treatment as soon as they can
# but I would not be surprised if there is no correlation. Younger age may be correlated to a higher
# survival rate since younger people are generally more healthy and able to fight the cancer better.
# Less days to drug therapy start may be correlated to higher survival rate since the sooner you start
# treatment, the higher the survival rate since the cancer has less time to develop and increase in size.
nrow(clinical)
nrow(clinical_rad)
nrow(clinical_drug)
clinical_data_merge <- merge(clinical, clinical_drug, by = "bcr_patient_barcode")
plot(x = clinical_data_merge$age_at_initial_pathologic_diagnosis, y = clinical_data_merge$days_to_drug_therapy_start, xlab = "Age", ylab = "Days to Drug Therapy Start", main = "Age vs Days to Drug Therapy Start")
# This plot tells us that there is little association between age and days until the drug therapy starts. I chose
# this plot since both variables are discrete


age_na_mask <- ifelse(!is.na(clinical_data_merge$age_at_initial_pathologic_diagnosis), TRUE, FALSE) # Remove rows without an age
age_cleaned_clinical <- clinical_data_merge[age_na_mask, ]

# Create a new column that turns age at initial diagnosis to a categorical variable
young_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis <= 35, TRUE, FALSE)
middle_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis > 35 & clinical_data_merge$age_at_initial_pathologic_diagnosis <= 50, TRUE, FALSE)
old_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis > 50, TRUE, FALSE)
age_cleaned_clinical$age_at_diagnosis_status <- ifelse(young_mask, "Young", ifelse(middle_mask, "Middle", "Old"))
age_cleaned_clinical$survival_time <- ifelse(is.na(age_cleaned_clinical$days_to_death), age_cleaned_clinical$survival_time <- age_cleaned_clinical$days_to_last_followup, age_cleaned_clinical$survival_time <- age_cleaned_clinical$days_to_death)

# Remove the rows with -Inf and na as the survival time
inf_mask <- ifelse(age_cleaned_clinical$survival_time == "-Inf", FALSE, TRUE)
age_cleaned_clinical <- age_cleaned_clinical[inf_mask, ]
na_mask <- ifelse(is.na(age_cleaned_clinical$survival_time), FALSE, TRUE)
age_cleaned_clinical <- age_cleaned_clinical[na_mask, ]

age_cleaned_clinical$death_event <- ifelse(age_cleaned_clinical$vital_status == "Alive", age_cleaned_clinical$death_event <- FALSE, age_cleaned_clinical$death_event <- TRUE)


# Create survival object and graph it in a KM plot
survival_object <- Surv(time = ifelse(age_cleaned_clinical$vital_status == "Alive", age_cleaned_clinical$days_to_last_followup, age_cleaned_clinical$days_to_death), event = ifelse(age_cleaned_clinical$vital_status == "Alive", FALSE, TRUE))
fit_object <- survfit(survival_object ~ age_cleaned_clinical$age_at_diagnosis_status, data = age_cleaned_clinical)
survplot <- ggsurvplot(fit_object , pval=TRUE, ggtheme =
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend =
                         "right")

jpeg("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/analysis_data/KM_plot1.jpg")
KM_plot <- survplot$plot + theme_bw() + theme(axis.title =
                                                element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=14), legend.text =
                                                element_text(size=12))
KM_plot
dev.off()


# Create a new mask that turns the days to drug therapy start to a categorical variable
early_mask <- ifelse(age_cleaned_clinical$days_to_drug_therapy_start <= 100, TRUE, FALSE)
middle_mask <- ifelse(age_cleaned_clinical$days_to_drug_therapy_start > 100 & age_cleaned_clinical$days_to_drug_therapy_start <= 500, TRUE, FALSE)
late_mask <- ifelse(age_cleaned_clinical$days_to_drug_therapy_start > 500, TRUE, FALSE)
age_cleaned_clinical$days_to_drug_therapy_start <- ifelse(early_mask, "Early", ifelse(middle_mask, "Middle", "Late"))
age_cleaned_clinical$survival_time <- ifelse(is.na(age_cleaned_clinical$days_to_death), age_cleaned_clinical$survival_time <- age_cleaned_clinical$days_to_last_followup, age_cleaned_clinical$survival_time <- age_cleaned_clinical$days_to_death)

# Remove the columns that have -Inf and na as the survival time
inf_mask <- ifelse(age_cleaned_clinical$survival_time == "-Inf", FALSE, TRUE)
age_cleaned_clinical <- age_cleaned_clinical[inf_mask, ]
na_mask <- ifelse(is.na(age_cleaned_clinical$survival_time), FALSE, TRUE)
age_cleaned_clinical <- age_cleaned_clinical[na_mask, ]

# Create a column that turns the vital status to false or true in a new column called "death_event"
age_cleaned_clinical$death_event <- ifelse(age_cleaned_clinical$vital_status == "Alive", age_cleaned_clinical$death_event <- FALSE, age_cleaned_clinical$death_event <- TRUE)

# Create a survival object and graph a KM plot
survival_object <- Surv(time = ifelse(age_cleaned_clinical$vital_status == "Alive", age_cleaned_clinical$days_to_last_followup, age_cleaned_clinical$days_to_death), event = ifelse(age_cleaned_clinical$vital_status == "Alive", FALSE, TRUE))
fit_object <- survfit(survival_object ~ age_cleaned_clinical$days_to_drug_therapy_start, data = age_cleaned_clinical)
survplot <- ggsurvplot(fit_object , pval=TRUE, ggtheme =
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend =
                         "right")

jpeg("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/analysis_data/KM_plot2.jpg")
KM_plot <- survplot$plot + theme_bw() + theme(axis.title =
                                                element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=14), legend.text =
                                                element_text(size=12))
KM_plot
dev.off()

# The age at diagnosis is not statistically significant to make any conclusions or correlations since the p-value is 0.19.
# The days to drug therapy start is very statistically significant since the p-value is less than 0.0001. The KM plot suggests
# that the middle range in drug therapy start correlates to lower survival rate.