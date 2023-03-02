setwd("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/analysis_data")

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

clinical <- read.csv("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/analysis_data/brca_clinical_data.csv")
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
write.csv(clinical, "/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/analysis_data/brca_clinical_data.csv", row.names = FALSE)


maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

#GDCdownload(maf_query)

maf <- GDCprepare(maf_query) 

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)

#Makes a new boolean column which groups everyone into white or non-white
maf_object@clinical.data$white_race_category <- ifelse(maf_object@clinical.data$race_list == 'WHITE', T, F)
head(maf_object@clinical.data$white_race_category) #view the top of the column

#Makes a new list of just the white patients' barcode
white_patients_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[maf_object@clinical.data$white_race_category]
white_maf <- subsetMaf(maf = maf_object, tsb = white_patients_barcodes) #Create a new maf object of just the white patients

#Make a new maf of the non-white patients
non_white_patients_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[!(maf_object@clinical.data$white_race_category)]
non_white_maf <- subsetMaf(maf = maf_object, tsb = non_white_patients_barcodes)

#Create a co-oncoplot of the white and non-white patients viewing the top 10 mutated genes
coOncoplot(m1 = white_maf, 
           m2 = non_white_maf,
           genes = maf_object@gene.summary$Hugo_Symbol[1:10],
           m1Name = 'White Patients', 
           m2Name = 'Non-white Patients', 
           borderCol = NA)

ggsave("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/week7_maf/White_vs_Non-white_CoOncoplot.png")
#TP53 is used to encode a protein called p53 which is a tumor suppressor protein. A mutation in this gene could prevent the protein to 
#function properly which can lead to cancer. Maybe TP53 mutations are more common in nnon-white patients since by chance white patients
#tend to have fewer mutations on this gene and pass this on to their white offspring.



#Create a new maf of the patients with a mutated TP53 gene
TP53_maf <- subsetMaf(maf = maf_object, genes = 'TP53')

#Make new lists of the tumor barcode of those white the TP53 mutation and white patients
mut_pats_TP53 <- TP53_maf@clinical.data$Tumor_Sample_Barcode
mut_pats_white <- white_maf@clinical.data$Tumor_Sample_Barcode

#Get the number of TP53 and white patients
num_pats_TP53 <- length(mut_pats_TP53)
num_pats_white <- length(mut_pats_white)

#Get the number of the patients who are both white and have a TP53 mutation
mut_pats_white_TP53 <- intersect(mut_pats_TP53, mut_pats_white)  
num_white_TP53 <- length(mut_pats_white_TP53)  

#Get the number of the patients who have the TP53 mutation and are not white and vice versa
num_TP53_only <- num_pats_TP53 - num_white_TP53
num_white_only <- num_pats_white - num_white_TP53

#Get the number of patients who are non-white and don't have a TP53 mutation
num_neither <- length(maf_object@clinical.data$Tumor_Sample_Barcode) - num_TP53_only - num_white_only - num_white_TP53

#Create a contingency table of TP53 mutations and white patients
contig <- matrix(c(num_white_TP53, 
                   num_white_only,
                   num_TP53_only,
                   num_neither), 
                 nrow=2)
#Display the contingency table
contig
mosaicplot(contig)

#Perform a fisher test on the contingency table
fisher_test <- fisher.test(contig)
fisher_test

#The p-value of the fisher test is 0.03955 which means that the difference between the rate of TP53 mutation in white vs non-white
#patients is statistically significant. The odds ratio is the chances that a random patient will be, in this case, white and negative for TP53
#mutations or non-white and positive for a TP53 mutation.


#Create a lollipop plot of the white and non-white patients and their TP53 mutations
lollipopPlot2(m1 = white_maf, 
              m2 = non_white_maf, 
              m1_name = 'White Patients',
              m2_name = 'Non-white Patients',
              gene = "TP53")
ggsave("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/week7_maf/White_vs_Non-white_TP53_Lollipop_plot.png")

#The mutations seem very similar between white and non-white patients but there does seem to be slightly more frame shift mutations in white patients

#Create a new boolean column of the patients' survival status 
maf_object@clinical.data$Overall_Survival_Status <- ifelse(maf_object@clinical.data$vital_status == 'Alive', T, F)

#Create a survival KM plot based on TP53
mafSurvival(maf = maf_object,
            genes = "TP53", ## pick a gene of your choosing
            time = "days_to_last_followup", ## name of the column in maf_object@clinical.data containing survival time
            Status = "Overall_Survival_Status", ## name of the column that contains a boolean value for death events, you may need to recreate this... 
            isTCGA = TRUE)
ggsave("/Users/christopherahn/Documents/QBIO_490/qbio-490-Christopher-Ahn/week7_maf/TP53_Survival_plot.png")

#There does not seem to be a difference in survival probability between those who have a tp53 mutation and those who do not. This could mean that
#the TP53 mutation itself does not really change the outcome of death or not, but other aspects of the cancer do





