library(stringr)

### Medicago 
stringtie_dir <-"Medicago_RNA/4stringtie"
stringtie_exp_files <- grep(dir(stringtie_dir,full.names = TRUE),pattern = ".tsv",value = TRUE)

### fpkm merge
fpkm <- c()
tpm <- c()
for (file in stringtie_exp_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE,sep="\t",quote = "")
  temp_file <- temp_file[order(temp_file$Gene.ID,decreasing = T),]
  fpkm <- cbind(fpkm,temp_file[,8])
  tpm <- cbind(tpm,temp_file[,9])
}
fpkm <- as.data.frame(fpkm)
tpm <- as.data.frame(tpm)

rownames(fpkm) = rownames(tpm) <-temp_file[,1]
colnames(fpkm) = colnames(tpm) <- basename(stringtie_exp_files)
colnames(fpkm) = colnames(tpm)<- str_split_fixed(colnames(fpkm),pattern = "\\.",n = 3)[,1]

## tpm  merge
write.csv(fpkm,file = "analysis_result/Mt_FPKM.csv")
write.csv(tpm,file = "analysis_result/Mt_TPM.csv")

#############
## Lotous
#############
library(tidyverse)
stringtie_dir <-"Lotus_RNA/5stringtie_gtf/"
stringtie_exp_files <- grep(dir(stringtie_dir,full.names = TRUE),pattern = ".tsv",value = TRUE)

### fpkm merge
fpkm <- c()
tpm <- c()
for (file in stringtie_exp_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE,sep="\t",quote = "")
  temp_file <- as.data.frame(temp_file%>%group_by(Gene.ID)%>%slice(which.max(TPM)))
  temp_file <- temp_file[order(temp_file$Gene.ID,decreasing = T),]
  fpkm <- cbind(fpkm,temp_file[,8])
  tpm <- cbind(tpm,temp_file[,9])
}
fpkm <- as.data.frame(fpkm)
tpm <- as.data.frame(tpm)

rownames(fpkm) = rownames(tpm) <- temp_file[,1]
colnames(fpkm) = colnames(tpm) <- basename(stringtie_exp_files)
colnames(fpkm) = colnames(tpm)<- str_split_fixed(colnames(fpkm),pattern = "\\.",n = 3)[,1]

## tpm  merge
write.csv(fpkm,file = "Lotus_RNA/analysis_result/Lt_FPKM.csv")
write.csv(tpm,file = "Lotus_RNA/analysis_result/Lt_TPM.csv")

#############
## Gmax
#############
library(tidyverse)
stringtie_dir <-"Gmax_RNA/4stringtie/"
stringtie_exp_files <- grep(dir(stringtie_dir,full.names = TRUE),pattern = ".tsv",value = TRUE)

### fpkm merge
fpkm <- c()
tpm <- c()
for (file in stringtie_exp_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE,sep="\t",quote = "")
  temp_file <- as.data.frame(temp_file%>%group_by(Gene.ID)%>%slice(which.max(TPM)))
  temp_file <- temp_file[order(temp_file$Gene.ID,decreasing = T),]
  fpkm <- cbind(fpkm,temp_file[,8])
  tpm <- cbind(tpm,temp_file[,9])
}
fpkm <- as.data.frame(fpkm)
tpm <- as.data.frame(tpm)

rownames(fpkm) = rownames(tpm) <- temp_file[,1]
colnames(fpkm) = colnames(tpm) <- basename(stringtie_exp_files)
colnames(fpkm) = colnames(tpm)<- str_split_fixed(colnames(fpkm),pattern = "\\.",n = 3)[,1]

## tpm  merge
write.csv(fpkm,file = "Gmax_RNA/analysis_result/GM_FPKM.csv")
write.csv(tpm,file = "Gmax_RNA/analysis_result/Gm_TPM.csv")

#############
## Rice
#############
library(tidyverse)
stringtie_dir <-"Rice_RNA/4stringtie/"
stringtie_exp_files <- grep(dir(stringtie_dir,full.names = TRUE),pattern = ".tsv",value = TRUE)

### fpkm merge
fpkm <- c()
tpm <- c()
for (file in stringtie_exp_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE,sep="\t",quote = "",comment.char = "#")
  temp_file <- as.data.frame(temp_file%>%group_by(Gene.ID)%>%slice(which.max(TPM)))
  temp_file <- temp_file[order(temp_file$Gene.ID,decreasing = T),]
  fpkm <- cbind(fpkm,temp_file[,8])
  tpm <- cbind(tpm,temp_file[,9])
}
fpkm <- as.data.frame(fpkm)
tpm <- as.data.frame(tpm)

rownames(fpkm) = rownames(tpm) <- temp_file[,1]
colnames(fpkm) = colnames(tpm) <- basename(stringtie_exp_files)
colnames(fpkm) = colnames(tpm)<- str_split_fixed(colnames(fpkm),pattern = "\\.",n = 3)[,1]

## tpm  merge
write.csv(fpkm,file = "Rice_RNA/analysis_result/Ri_FPKM.csv")
write.csv(tpm,file = "Rice_RNA/analysis_result/Ri_TPM.csv")

#############
## Maize
#############
library(tidyverse)
stringtie_dir <-"Maize_RNA/4stringtie/"
stringtie_exp_files <- grep(dir(stringtie_dir,full.names = TRUE),pattern = ".tsv",value = TRUE)

### fpkm merge
fpkm <- c()
tpm <- c()
for (file in stringtie_exp_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE,sep="\t",quote = "",comment.char = "#")
  temp_file <- as.data.frame(temp_file%>%group_by(Gene.ID)%>%slice(which.max(TPM)))
  temp_file <- temp_file[order(temp_file$Gene.ID,decreasing = T),]
  fpkm <- cbind(fpkm,temp_file[,8])
  tpm <- cbind(tpm,temp_file[,9])
}
fpkm <- as.data.frame(fpkm)
tpm <- as.data.frame(tpm)

rownames(fpkm) = rownames(tpm) <- temp_file[,1]
colnames(fpkm) = colnames(tpm) <- basename(stringtie_exp_files)
colnames(fpkm) = colnames(tpm)<- str_split_fixed(colnames(fpkm),pattern = "\\.",n = 3)[,1]

## tpm  merge
write.csv(fpkm,file = "Maize_RNA/analysis_result/Zm_FPKM.csv")
write.csv(tpm,file = "Maize_RNA/analysis_result/Zm_TPM.csv")

#############
## Soly 
#############
library(tidyverse)
stringtie_dir <-"Sly_RNA/4stringtie/"
stringtie_exp_files <- grep(dir(stringtie_dir,full.names = TRUE),pattern = ".tsv",value = TRUE)

### fpkm merge
fpkm <- c()
tpm <- c()
for (file in stringtie_exp_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE,sep="\t",quote = "",comment.char = "#")
  temp_file <- as.data.frame(temp_file%>%group_by(Gene.ID)%>%slice(which.max(TPM)))
  temp_file <- temp_file[order(temp_file$Gene.ID,decreasing = T),]
  fpkm <- cbind(fpkm,temp_file[,8])
  tpm <- cbind(tpm,temp_file[,9])
}
fpkm <- as.data.frame(fpkm)
tpm <- as.data.frame(tpm)

rownames(fpkm) = rownames(tpm) <- temp_file[,1]
colnames(fpkm) = colnames(tpm) <- basename(stringtie_exp_files)
colnames(fpkm) = colnames(tpm)<- str_split_fixed(colnames(fpkm),pattern = "\\.",n = 3)[,1]

## tpm  merge
write.csv(fpkm,file = "Sly_RNA/analysis_result/Sl_FPKM.csv")
write.csv(tpm,file = "Sly_RNA/analysis_result/Sl_TPM.csv")
