
######The purpose of this script is to merge clinical data to the BRCAMerged-NoClinical######

#increased memory.limit(size = 60000). This is needed for the final transpose
memory.limit(60000)
library(tidyverse)
library(data.table)

#first file should be the already merged data except for the clinical
merged_data <- read.csv("data/BRCAMerged-NoClin.csv")
#merged_data <- read.csv("data/BRCAMergedPract.csv")
#107304 2241

#This file should be the clinical data alone
clinical_data <- read.csv("data/BRCA.clin.merged.csv")
#clinical_data <- read.csv("data/ClinPract.csv")
#3720 1098

#turn data from csv to dataframe for ease of manipulation before transpose
clinical_data <- as.data.frame(clinical_data)
merged_data <- as.data.frame(merged_data)

#transpose merged and clinical for easier merging and editting of clinical's barcode
merged_data <- t(merged_data)
clinical_data <- t(clinical_data)

#setting the column names and removing row names. Row names becomes an actual column
clinical_data <- cbind(Col.Names = row.names(clinical_data), clinical_data)
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1,]
rownames(clinical_data) <- NULL

merged_data <- cbind(Col.Names = row.names(merged_data), merged_data)
colnames(merged_data) <- merged_data[1, ]
merged_data <- merged_data[-1,]
rownames(merged_data) <- NULL

#replacing the "-" in the barcode with . in clinical for easier merge. 
#Then add .01 to end of barcode to match only the .01 from merged data.
#Add "data." as a prefix then capatilize entire string to match merged_data
clinical_data[, "patient.bcr_patient_barcode"] <- gsub("-", ".", clinical_data[, "patient.bcr_patient_barcode"], fixed = TRUE) %>% 
  toupper()
clinical_data[, "patient.bcr_patient_barcode"] <- paste0("Data.", clinical_data[, "patient.bcr_patient_barcode"])
clinical_data[, "patient.bcr_patient_barcode"] <- paste0(clinical_data[, "patient.bcr_patient_barcode"], ".01")

#change clinical's patient barcode row name to match where the barcodes are in the merged data frame
colnames(clinical_data)[colnames(clinical_data) == "patient.bcr_patient_barcode"] <- "Des.GeneSymbol"

#merge the two together, warning will appear but it's okay
all_merged <- merge(merged_data, clinical_data, all = TRUE)
all_merged <- t(all_merged)

#moving row names to an actual column
all_merged <- cbind(all_merged, Row.Names = rownames(all_merged))
colnames(all_merged) <- all_merged[1, ]
all_merged <- all_merged[-1,]
rownames(all_merged) <- NULL

#move Des.Description to front, then move Des.Platform to the front, then move Des.GeneSymbol to the front
all_merged <- subset(all_merged, select = c(Des.Description, Data.TCGA.3C.AAAU.01:Des.GeneSymbol))
all_merged <- subset(all_merged, select = c(Des.Platform, Des.Description:Des.GeneSymbol))
all_merged <- subset(all_merged, select=c(Des.GeneSymbol, Des.Platform:Des.GeneSymbol))
all_merged <- all_merged[, !duplicated(colnames(all_merged))]

#create row numbers
rownames(all_merged) <- c(1:nrow(all_merged))

#eliminate some NA's and properly fill out description
all_merged[, "Des.Platform"][is.na(all_merged[, 'Des.Platform'])] <- "Clinical"
remove_gene_sym <- which(colnames(clinical_data) == "Des.GeneSymbol")
all_merged[, "Des.Description"][is.na(all_merged[, "Des.Description"])] <- colnames(clinical_data)[-remove_gene_sym]

#remove Des.Descrition values from Des.GeneSymbol and replace with "-"
all_merged <- as.data.frame(all_merged)
all_merged <- setDT(all_merged)[all_merged, on = c("Des.GeneSymbol==Des.Description"), Des.GeneSymbol := "-"][]

#sort by Des.Platform
all_merged <- all_merged[order(all_merged[, "Des.Platform"]), ]

write.csv(all_merged, file = "data/BRCAMergeAutomated.csv")




