setwd("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/midterm_project_Ahn") #Set the working directory to "outputs"
dir.create("outputs")
setwd("outputs")

#Load in required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")
if (!require("maftools", quietly = TRUE))
  BiocManager::install("maftools")
library(BiocManager)
library(TCGAbiolinks)
library(maftools)

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
if(!require(DESeq2)){
  BiocManager::install("DESeq2")
}
library(DESeq2)

if(!require(EnhancedVolcano)){
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

#Load in the clinical data and adjust the column names to the Tumor_Sample_Barcode
clin_query <- GDCquery(project = "TCGA-BRCA",data.category = "Clinical",file.type = "xml")
#GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query,
                            clinical.info = "patient")
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
#Load in maf_object data
maf_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  access = "open",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
#GDCdownload(maf_query)
maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf,
                       clinicalData = clinic,
                       isTCGA = TRUE)

#Load in rna_se object
rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
#GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)


# Number of rows in clinic
nrow(clinic) # 1174

#Check to see if there are any NA in the breast carcinoma estrogen receptor status
sum(is.na(clinic$breast_carcinoma_estrogen_receptor_status))

#Check to see if every patient has an age at initial pathological diagnosis recorded
sum(is.na(clinic$age_at_initial_pathologic_diagnosis))


sum(clinic$breast_carcinoma_estrogen_receptor_status == "Indeterminate") # 50 patients
sum(clinic$breast_carcinoma_estrogen_receptor_status == "") # 2 Patients

#Create a mask that is true for only positive and negative patients (breast carcinoma estrogen receptor status)
positive_negative_mask <- ifelse(clinic$breast_carcinoma_estrogen_receptor_status == "Positive", T, 
                                 ifelse(clinic$breast_carcinoma_estrogen_receptor_status == "Negative", T, F)) 

#Create a new dataframe using the mask
clinic_cleaned <- clinic[positive_negative_mask, ]
nrow(clinic_cleaned) #1122 rows


#Turn the positive and negative values into factors
clinic_cleaned$breast_carcinoma_estrogen_receptor_status <- factor(clinic_cleaned$breast_carcinoma_estrogen_receptor_status)

#Plot the data in a boxplot and save it as a jpeg
par(mar=c(1,1,1,1))
jpeg("Boxplot.jpg")
boxplot(clinic_cleaned$age_at_initial_pathologic_diagnosis ~ clinic_cleaned$breast_carcinoma_estrogen_receptor_status,
        xlab = "Breast Carcinoma Estrogen Receptor Status",
        ylab = "Age at Initial Pathological Diagnosis",
        main = "Age vs. Breast Carcinoma Estrogen Receptor Status")
dev.off()


# Make a KM Plot
sum(is.na(clinic_cleaned$age_at_initial_pathologic_diagnosis)) # Everyone has a reported initial age

# Make a new column that has the survival time of each patient
clinic_cleaned$survival_time <- ifelse(
  is.na(clinic_cleaned$days_to_death), 
  clinic_cleaned$survival_time <- clinic_cleaned$days_to_last_followup, 
  clinic_cleaned$survival_time <- clinic_cleaned$days_to_death)


# Remove any -Inf or NA rows for survival time
inf_mask <- ifelse(clinic_cleaned$survival_time == "-Inf", FALSE, TRUE)
clinic_cleaned <- clinic_cleaned[inf_mask, ]
na_mask <- ifelse(is.na(clinic_cleaned$survival_time), FALSE, TRUE)
clinic_cleaned <- clinic_cleaned[na_mask, ]

# Adds another column that is called the death event
clinic_cleaned$death_event <- ifelse(
  clinic_cleaned$vital_status == "Alive", 
  clinic_cleaned$death_event <- FALSE, 
  clinic_cleaned$death_event <- TRUE)

# Create a survival object
survival_object <- Surv(time = ifelse(clinic_cleaned$vital_status == "Alive", clinic_cleaned$days_to_last_followup, clinic_cleaned$days_to_death), event = ifelse(clinic_cleaned$vital_status == "Alive", FALSE, TRUE))
fit_object <- survfit(survival_object ~ clinic_cleaned$breast_carcinoma_estrogen_receptor_status, data = clinic_cleaned)
survplot <- ggsurvplot(fit_object, pval=TRUE, ggtheme =
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend.title = "Legend", legend = c(0.8,0.1))

png("KM_plot.png")
KM_plot <- survplot$plot + theme_bw() + theme(axis.title =
                                                element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=10), legend.text =
                                                element_text(size=10))
dev.off()


#Create a MAF Plot
sum(is.na(maf_object@clinical.data$breast_carcinoma_estrogen_receptor_status)) # 44 Patients are NA for this 


#Create a subset maf for patients that are positive and negative for breast carcinoma estrogen receptor status
positive_mask <- ifelse(maf_object@clinical.data$breast_carcinoma_estrogen_receptor_status == "Positive", T, F)
positive_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[positive_mask]
positive_maf <- subsetMaf(maf = maf_object, tsb = positive_patient_barcodes)

negative_mask <- ifelse(maf_object@clinical.data$breast_carcinoma_estrogen_receptor_status == "Negative", T, F)
negative_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[negative_mask]
negative_maf <- subsetMaf(maf = maf_object, tsb = negative_patient_barcodes)

# Create a Co-oncoplot for positive and negative patients
coOncoplot(m1 = positive_maf, 
           m2 = negative_maf, 
           m1Name = 'Positive Breast Carcinoma Status', 
           m2Name = 'Negative Breast Carcinoma Status', 
           borderCol = NA)
ggsave("CoOncoplot.png")
dev.off()

# Create a Volcano Plot

# Create new df's for RNA genes and RNA Clinical
rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)

rna_clinical <-rna_se@colData
rna_clinical <- as.data.frame(rna_clinical)


treatments_mask <- !colnames(rna_clinical) %in% c('treatments', 'primary_site', 'disease_type')
rna_clinical <- rna_clinical[, treatments_mask]

# Turn vital status and ajcc pathological stage into factors
rna_clinical$vital_status <- factor(rna_clinical$vital_status)
rna_clinical$ajcc_pathologic_stage <- factor(rna_clinical$ajcc_pathologic_stage)

rna_counts <- rna_se@assays@data$unstranded


# Transfer the Breast Carcinoma Data from the clinic df to the rna_clinical df
# Make a mask that is true for the patients that are in both the rna_clinical df and the clinic df
colnames(rna_clinical) [2] <- "bcr_patient_barcode"
patient_mask <- ifelse(rna_clinical$bcr_patient_barcode %in% clinic$bcr_patient_barcode, T, F)

# Remove the duplicate rows in clinic to change the row names to the bcr_patient_barcode
clinic_nonduplicates <- clinic[!duplicated(clinic), ]
row.names(clinic_nonduplicates) <- clinic_nonduplicates$bcr_patient_barcode

# Apply the patient mask to have only the patients in the clinic df and the rna_clinical df
cleaned_rna_clinical <- rna_clinical[patient_mask, ]

# Create a column in cleaned_rna_clinical that has the breast carcinoma estrogen receptor status data from the clinic df
cleaned_rna_clinical$breast_carcinoma_estrogen_receptor_status <- 
  ifelse(clinic_nonduplicates[cleaned_rna_clinical$bcr_patient_barcode, "breast_carcinoma_estrogen_receptor_status"] == "Positive", T, F)

# Remove any patients with NA
rna_na_mask <- ifelse(!is.na(cleaned_rna_clinical$breast_carcinoma_estrogen_receptor_status), T, F)
cleaned_rna_clinical <- cleaned_rna_clinical[rna_na_mask, ]

# Make the breast carcinoma receptor status variables a factor (T or F)
cleaned_rna_clinical$breast_carcinoma_estrogen_receptor_status <- factor(cleaned_rna_clinical$breast_carcinoma_estrogen_receptor_status)

# Remove any patients that are NA for ajcc_pathological_stage in the cleaned_rna_clinical df and the rna_counts df
na_mask <-  ifelse(is.na(cleaned_rna_clinical$ajcc_pathologic_stage), F, T)
cleaned_rna_clinical <-  cleaned_rna_clinical[na_mask, ]
na_mask <-  ifelse(is.na(rna_clinical$ajcc_pathologic_stage), F, T)
rna_counts <- rna_counts[, na_mask] 

# Remove the genes that have less than 10
row_sums <- rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums < 10, F, T)

# rewrite the rna_counts df, subsetting for only genes with >= 10 total counts
rna_counts <- rna_counts[low_counts_mask,]

#update rna_genes with the low_counts_mas
rna_genes <- rna_genes[low_counts_mask,]

# Run DESeq on the data set
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = cleaned_rna_clinical,
                              design =~ vital_status + ajcc_pathologic_stage + breast_carcinoma_estrogen_receptor_status)


dds_obj <- DESeq(dds) # note: this will likely take a long time (ie 45 minutes to 2 hours)


resultsNames(dds_obj)  # see what comparisons got run

# get the Positive vs. Negative comparison
results <- results(
  dds_obj, format = "DataFrame", contrast = c("breast_carcinoma_estrogen_receptor_status", "TRUE", "FALSE")) 


# Filter out the genes that are not in results
gene_mask <- ifelse(rna_genes$gene_id %in% results@rownames, T, F)
rna_genes_subset <- rna_genes$gene_name[gene_mask]


# Add columns to results
results <- data.frame(rna_genes, results, results$log2FoldChange, results$pvalue, results$padj, -log10(results$padj))
colnames(results) <- c("Gene Name", "Gene ID", "log2foldchange", "p-value", "padj", "-log10(padj)")
keep <- c("gene_name", "gene_id", "log2FoldChange", "pvalue", "padj")
cleaned_results <- results[,names(results) %in% keep]
cleaned_results$"-log10(padj)" <- -log10(cleaned_results$padj)

# Create the threshold for statistical significance
sig_results <- (cleaned_results$"p-value" < 0.05)

# Create a data frame for up-regulated genes
up_reg_results <- cleaned_results[order(cleaned_results$log2FoldChange), ]
up_reg_results <- up_reg_results[(up_reg_results$log2foldchange > 1),]

# Create a data frame for down-regulated genes
neg_log_mask <- !(cleaned_results$log2FoldChange < -1)
down_reg_results <- cleaned_results[neg_log_mask,]
down_reg_results <- down_reg_results[order(down_reg_results$log2FoldChange),]

# Create the volcano plot and save it
rownames(up_reg_results) <- up_reg_results$gene_id
rownames(down_reg_results) <- down_reg_results$gene_id
par(mar=c(1,1,1,1))
jpeg("Enhanced_Volcano_Plot.jpg")
EnhancedVolcano(results, lab = rownames(results), x = "log2FoldChange", y = "pvalue")
dev.off()
