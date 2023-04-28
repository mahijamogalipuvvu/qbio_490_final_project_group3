dir.create("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/final_project_group3/outputs")
setwd("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/final_project_group3/outputs")

if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)

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
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
if(!require(survival)){
  install.packages("survival")
}
library(survival)
if(!require(survminer)){
  install.packages("survminer")
}



# Download TCGA Clinical Data
clin_query <- GDCquery(project = "TCGA-GBM",data.category = "Clinical",file.type = "xml")
GDCdownload(clin_query)
clinical <- GDCprepare_clinic(clin_query,
                            clinical.info = "patient")
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
write.csv(clinical,"clinical.csv")
clinical_drug <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query,
                                  clinical.info = "radiation")


# Download TCGA MAF Data
data_query <- GDCquery(project = "TCGA-GBM",
                       data.category = "Simple Nucleotide Variation",
                       access = "open",
                       workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
                       data.type = "Masked Somatic Mutation")
#GDCdownload(data_query)
maf_data <- GDCprepare(data_query)
maf_object <- read.maf(maf = maf_data, 
                       clinicalData = clinical,
                       isTCGA = TRUE)
mutation_data <- maf_object@data
write.csv(mutation_data,"mutation_data.csv",row.names=F)


# Download TCGA RNA Data
rna_query <- GDCquery(project ="TCGA-GBM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

#GDCdownload(rna_query)

rna_se <- GDCprepare(rna_query)


rna_clinical <-rna_se@colData
treatments_mask <- !colnames(rna_clinical) %in% c('treatments', 'primary_site', 'disease_type')
rna_clinical <- rna_clinical[, treatments_mask]
rna_clinical <- as.data.frame(rna_clinical)
write.csv(rna_clinical, "rna_clinical.csv", row.names = TRUE)

rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
write.csv(rna_genes, "rna_genes.csv", row.names = TRUE)

rna_counts <- rna_se@assays@data$unstranded
rna_counts <- as.data.frame(rna_counts)
write.csv(rna_counts, "rna_counts.csv", row.names = TRUE)


# Create an oncoplot with the top 10 mutated genes in GBM patients
jpeg("oncoplot.jpg")
oncoplot(maf = maf_object,
         top = 10,
         borderCol = NA) 
dev.off()

# Did not use
jpeg("PTEN_lollipop.jpg")
lollipopPlot(maf = maf_object,
             gene = "PTEN")  
dev.off()

# Create a new column in maf_object that tracks survival status
maf_object@clinical.data$Overall_Survival_Status <- ifelse(is.na(maf_object@clinical.data$days_to_death), T, F)  

#mutation_data <- read.csv("mutation_data.csv")

# Create a KM plot for ALDOB
jpeg("ALDOB_KMplot.jpg")
mafSurvival(maf = maf_object,
            genes = "ALDOB", 
            time = "days_to_last_followup",
            Status = "Overall_Survival_Status", 
            isTCGA = TRUE)  
dev.off()
  
# Create a KM plot for SLC2A5
jpeg("SLC2A5_KMplot.jpg")
mafSurvival(maf = maf_object,
            genes = "SLC2A5", 
            time = "days_to_last_followup",
            Status = "Overall_Survival_Status", 
            isTCGA = TRUE)  
dev.off()

# Did not use
jpeg("Somatic_interactions_plot.jpg")
somaticInteractions(maf = maf_object,
                    genes = c('ALDOB', 'SLC2A5','TP53', 'PTEN', 'MUK16','PIK3CA', 'TTN', 'LRP2', 'DST', 'NF1', 'FLG', 'GATA3', 'DMD','CDH1', 'ZFHX4','SPTA1', 'NCOR1','MAP3K1', 'RYR2','SYNE1'),
                    pvalue = c(0.05, 0.1))
dev.off()

#rna_clinical <- read.csv("rna_clinical.csv")
#rna_genes <- read.csv("rna_genes.csv")
#rna_counts <- read.csv("rna_counts.csv")

# Reformat rna_counts and rna_genes to have the gene_id
rna_counts <- rna_counts[,2:176]
colnames(rna_counts) <- rna_clinical$barcode
rownames(rna_counts) <- rna_genes$gene_id
row.names(rna_genes) <- rna_genes$gene_id

# Create the age category with 58 years old as the cutoff
rna_clinical$age_category <- ifelse(rna_clinical$age_at_index > 58, 'Old', 'Young')
rna_clinical$age_category <- factor(rna_clinical$age_category)
rna_clinical$vital_status <- factor(rna_clinical$vital_status)
#rna_clinical$gender <- factor(rna_clinical$gender)

# Remove all the na values
na_mask <-  ifelse(is.na(rna_clinical$gender), F, T)
rna_clinical <-  rna_clinical[na_mask, ]
rna_counts <- rna_counts[, na_mask]

# Check to make sure the na values have all been removed
sum(is.na(rna_clinical$vital_status))
sum(is.na(rna_clinical$gender))
sum(is.na(rna_clinical$age_category))

# Remove the rna rows that have less than 10 counts
row_sums <- rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums < 10, F, T)
rna_counts <- rna_counts[low_counts_mask,]
rna_genes <- rna_genes[low_counts_mask,]


# Run differential analysis on vital status while controlling for age
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design =~ age_category + vital_status)

dds_obj <- DESeq(dds) 

resultsNames(dds_obj)  # see what comparisons got run

# get the Dead vs. Alive comparison
results <- results(dds_obj, format = "DataFrame", contrast = c("vital_status", "Dead", "Alive"))

# Only look at the genes in rna_genes
gene_mask <- ifelse(rna_genes$gene_id %in% results@rownames, T, F)
rna_genes_subset <- rna_genes$gene_name[gene_mask]
results <- data.frame(rna_genes_subset, rownames(results), results$log2FoldChange, results$pvalue, results$padj, -log10(results$padj))

# Create a new dataframe for results with up-regulated and down-regulated genes
colnames(results) <- c("Gene Name", "Gene ID", "log2foldchange", "p-value", "padj", "-log10(padj)")
sig_results <- (results$"p-value" < 0.05)
up_reg_results <- results[order(results$log2foldchange), ]
up_reg_results <- up_reg_results[(up_reg_results$log2foldchange > 1),]

neg_log_mask <- !(results$log2foldchange < -1)
down_reg_results <- results[neg_log_mask,]
down_reg_results <- down_reg_results[order(down_reg_results$log2foldchange),]

rownames(up_reg_results) <- up_reg_results$`Gene ID`
rownames(down_reg_results) <- down_reg_results$`Gene ID`

# Create an enhanced volcano plot
par(mar=c(1,1,1,1))
jpeg("Enhanced_Volcano_Plot.jpg")
EnhancedVolcano(results, lab = results[,"Gene Name"], x = "log2foldchange", y = "p-value", shape = c(1, 4, 23, 25))
dev.off()


# Did not use
colnames(clinical_drug)[ colnames(clinical_drug) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
clinical <- merge(clinical, clinical_drug, by = "Tumor_Sample_Barcode")
sum(is.na(clinical$age_at_initial_pathologic_diagnosis))
young_mask <- ifelse(clinical$age_at_initial_pathologic_diagnosis <= 35, TRUE, FALSE)
middle_mask <- ifelse(clinical$age_at_initial_pathologic_diagnosis > 35 & clinical$age_at_initial_pathologic_diagnosis <= 50, TRUE, FALSE)
old_mask <- ifelse(clinical$age_at_initial_pathologic_diagnosis > 50, TRUE, FALSE)
clinical$age_at_diagnosis_status <- ifelse(young_mask, "Young", ifelse(middle_mask, "Middle", "Old"))
clinical$survival_time <- ifelse(is.na(clinical$days_to_death), clinical$survival_time <- clinical$days_to_last_followup, clinical$survival_time <- clinical$days_to_death)

inf_mask <- ifelse(clinical$survival_time == "-Inf", FALSE, TRUE)
clinical <- clinical[inf_mask, ]
na_mask <- ifelse(is.na(clinical$survival_time), FALSE, TRUE)
clinical <- clinical[na_mask, ]
clinical$vital_status <- ifelse(is.na(clinical$days_to_death), 'Alive', 'Dead')
clinical$death_event <- ifelse(clinical$vital_status == "Alive", clinical$death_event <- FALSE, clinical$death_event <- TRUE)

survival_object <- Surv(time = ifelse(clinical$vital_status == "Alive", clinical$days_to_last_followup, clinical$days_to_death), event = ifelse(clinical$vital_status == "Alive", FALSE, TRUE))
fit_object <- survfit(survival_object ~ clinical$age_at_diagnosis_status, data = clinical)
survplot <- ggsurvplot(fit_object , pval=TRUE, ggtheme =
                         theme(plot.margin = unit(c(1,1,1,1), "cm")), legend =
                         "right")

jpeg('KM_plot_age.jpg')
KM_plot <- survplot$plot + theme_bw() + theme(axis.title =
                                                element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=14), legend.text =
                                                element_text(size=12))
KM_plot
dev.off()





