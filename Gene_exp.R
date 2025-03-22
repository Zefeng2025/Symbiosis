## Medicago

library(stringr)
library(edgeR)
library(clusterProfiler)
library(GOplot)
library(mixOmics)


fc_count_dir <-"~/MyResearch/ATAC_Seq_WET/Medicago_RNA/2feature_count/"
fc_count_files <- grep(dir(fc_count_dir,full.names = TRUE),pattern = ".summary",invert = TRUE,value = TRUE)

fc<-c()
for (file in fc_count_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE)
  fc<-cbind(fc,temp_file[,7])
}

rownames(fc)<-temp_file$Geneid
colnames(fc) <- basename(fc_count_files)
colnames(fc)<- str_split_fixed(colnames(fc),pattern = "\\.fc",n = 2)[,1]
#pheatmap::pheatmap(cor(log(fc+1,2)))

## PCA analysis
group <- str_split_fixed(colnames(fc),pattern = "_",n = 3)[,2]
MyResult.pca <- pca(t(log(fc[rowSums(cpm(fc)>=1)>=3,]+1,2)))
plotIndiv(MyResult.pca,group = group,ellipse = TRUE,title = "Samples",size.axis = 20,size.xlabel = 20,size.ylabel = 20)

## DEG analysis prepare
group_list<-factor(group) # make group list
coldata<-data.frame(row.names = colnames(fc),condition=group_list)  # make coldata

expdat = fc[rowSums(cpm(fc)>=1) >= 2,]  # filter, and 21898/52847 genes were kept
exprSet <- DGEList(counts = expdat, group = group_list) # make dgelist
exprSet <- calcNormFactors(exprSet)  # normlaized factor
exprSet <- estimateCommonDisp(exprSet) # 
exprSet <- estimateTagwiseDisp(exprSet)

## DEG comparative analysis for ck and AM
et <- exactTest(exprSet,pair = c("CK","AM"))
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

logFC_cutoff<-1
tTag$change<-as.factor(ifelse(tTag$FDR<0.01&abs(tTag$logFC)>logFC_cutoff,
                              ifelse(tTag$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))
save(tTag,file = "~/MyResearch/ATAC_Seq_WET/Medicago_RNA/analysis_result/Mt_EdgeR_CKvsAM.Rdata")
this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated gene is ',nrow(tTag[tTag$change =='UP',]) ,
                     '\nThe number of down-regulated gene is ',nrow(tTag[tTag$change =='DOWN',]))
df <- rbind(tTag[tTag$change!="NOT",][1:10,],tTag[tTag$change=="DOWN",][1:10,])

p2<-ggplot(data=tTag,aes(x=logFC,y=-log10(FDR),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(size=15,hjust=0.5),legend.position = c(0.2,0.8))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

up_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(logFC>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(logFC< -1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out <- rbind(up_reg,down_reg) 

## annotation knwon genes
known_genes_mapping <- read_excel("/home/wuzefeng/MyResearch/Medicago_data/1other_data/Medicago_gene2name_mapping.xls")
known_genes_mapping <- known_genes_mapping[,c(1,2)]
known_genes_mapping <- na.omit(known_genes_mapping)
colnames(known_genes_mapping) <- c("Gene_name","Gene_ID")
known_genes_mapping$Gene_ID<- str_replace_all(known_genes_mapping$Gene_ID,pattern = "Medtr",replacement = "MTR_")
degs_out$gene_name <- known_genes_mapping$Gene_name[match(rownames(degs_out),known_genes_mapping$Gene_ID)]
write.csv(degs_out,"~/MyResearch/ATAC_Seq_WET/Medicago_RNA/analysis_result/Mt_CK_AM.deg.csv",na = "")
## functions enrichment
library(gprofiler2)
library(GOplot)
go_result <- gost(query = rownames(degs_out),evcodes = TRUE,organism = "mtruncatula")
save(go_result,file = "~/MyResearch/ATAC_Seq_WET/Medicago_RNA/analysis_result/DEG_AM_GO_gprofile.Rdata")
#load("DEG_GO_gprofile.Rdata")
go_reduce <- go_result$result[,c(10,9,11,3,16)]
colnames(go_reduce) <- c("category","ID","term","adj_pval","genes")
go_deg <- data.frame(ID=rownames(degs_out),logFC=degs_out$logFC)
circ <- circle_dat(go_reduce, go_deg)
GOCircle(circ,nsub = 20,label.size = 3)

##########################################################################
#### DEG comparative analysis for ck and rz
et <- exactTest(exprSet,pair = c("CK","RZ"))
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

logFC_cutoff<-1
tTag$change<-as.factor(ifelse(tTag$FDR<0.01&abs(tTag$logFC)>logFC_cutoff,
                              ifelse(tTag$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))
save(tTag,file = "~/MyResearch/ATAC_Seq_WET/Medicago_RNA/analysis_result/Mt_EdgeR_CKvsRZ.Rdata")
this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated gene is ',nrow(tTag[tTag$change =='UP',]) ,
                     '\nThe number of down-regulated gene is ',nrow(tTag[tTag$change =='DOWN',]))
df <- rbind(tTag[tTag$change!="NOT",][1:10,],tTag[tTag$change=="DOWN",][1:10,])

p2<-ggplot(data=tTag,aes(x=logFC,y=-log10(FDR),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(size=15,hjust=0.5),legend.position = c(0.2,0.8))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

library(tidyverse)
up_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(logFC>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(logFC< -1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out <- rbind(up_reg,down_reg) 

degs_out$gene_name <- known_genes_mapping$Gene_name[match(rownames(degs_out),known_genes_mapping$Gene_ID)]
write.csv(degs_out,"~/MyResearch/ATAC_Seq_WET/Medicago_RNA/analysis_result/Mt_CK_RZ.deg.csv",na = "")

## functions enrichment
library(gprofiler2)
library(GOplot)
go_result <- gost(query = rownames(degs_out),evcodes = TRUE,organism = "mtruncatula")
save(go_result,file = "~/MyResearch/ATAC_Seq_WET/Medicago_RNA/analysis_result/RZ_DEG_GO_gprofile.Rdata")
#load("DEG_GO_gprofile.Rdata")
go_reduce <- go_result$result[,c(10,9,11,3,16)]
colnames(go_reduce) <- c("category","ID","term","adj_pval","genes")
go_deg <- data.frame(ID=rownames(degs_out),logFC=degs_out$logFC)
circ <- circle_dat(go_reduce, go_deg)
GOCircle(circ,nsub = 20,label.size = 3)

#########################################################################################################
## lotus
########################################################################################################
library(stringr)
library(edgeR)

fc_count_dir <-"~/MyResearch/ATAC_Seq_WET/Lotus_RNA/2feature_count/"
fc_count_files <- grep(dir(fc_count_dir,full.names = TRUE),pattern = ".summary",invert = TRUE,value = TRUE)

fc<-c()
for (file in fc_count_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE)
  fc<-cbind(fc,temp_file[,7])
}


rownames(fc)<-temp_file$Geneid
colnames(fc) <- basename(fc_count_files)
colnames(fc)<- str_split_fixed(colnames(fc),pattern = "\\.fc",n = 2)[,1]

## filter by genes
true_genes <- read.table("~/Genomes/Plants/Lotus_japonicus/Lt_genes",stringsAsFactors = F)
fc <- fc[rownames(fc)%in%true_genes$V1,]

#pheatmap::pheatmap(cor(log(fc+1,2)))

## PCA analysis
group <- str_split_fixed(colnames(fc),pattern = "_",n = 3)[,2]
MyResult.pca <- pca(t(log(fc[rowSums(cpm(fc)>=1)>=2,]+1,2)))
plotIndiv(MyResult.pca,group = group,ellipse = TRUE,title = "Samples",size.axis = 20,size.xlabel = 20,size.ylabel = 20)

## DEG analysis prepare
group_list<-factor(group) # make group list
coldata<-data.frame(row.names = colnames(fc),condition=group_list)  # make coldata

expdat = fc[rowSums(cpm(fc)>0) >= 1,]  # filter, and 21898/52847 genes were kept
exprSet <- DGEList(counts = expdat, group = group_list) # make dgelist
exprSet <- calcNormFactors(exprSet)  # normlaized factor
exprSet <- estimateCommonDisp(exprSet) # 
exprSet <- estimateTagwiseDisp(exprSet)

## DEG comparative analysis for ck and AM
et <- exactTest(exprSet,pair = c("CK","AM"))
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

logFC_cutoff<-1
tTag$change<-as.factor(ifelse(tTag$FDR<0.01&abs(tTag$logFC)>logFC_cutoff,
                              ifelse(tTag$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))
save(tTag,file = "~/MyResearch/ATAC_Seq_WET/Lotus_RNA/analysis_result/Lotus_EdgeR_CKvsAM.Rdata")
this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated gene is ',nrow(tTag[tTag$change =='UP',]) ,
                     '\nThe number of down-regulated gene is ',nrow(tTag[tTag$change =='DOWN',]))
df <- rbind(tTag[tTag$change!="NOT",][1:10,],tTag[tTag$change=="DOWN",][1:10,])

p2<-ggplot(data=tTag,aes(x=logFC,y=-log10(FDR),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(size=15,hjust=0.5),legend.position = c(0.2,0.8))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

up_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(logFC>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(logFC< -1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out <- rbind(up_reg,down_reg) 

gene_anno <- read.delim("/home/wuzefeng/Genome.db/gff3/Lt.gene.description.txt",header = F,stringsAsFactors = F,sep="\t",comment.char = "")
degs_out$gene_annotation <- gene_anno$V2[match(rownames(degs_out),gene_anno$V1)]

write.csv(degs_out,file = "~/MyResearch/ATAC_Seq_WET/Lotus_RNA/analysis_result/Lt_CK_AM.deg.csv")
## functions enrichment
## clusterprofiler
degs_out <- read.csv("~/MyResearch/RNA-Seq_WET/Lotus_RNA/analysis_result/Lt_CK_AM.deg.csv")
library(clusterProfiler)
degs_out <- read.csv("~/MyResearch/RNA-Seq_WET/Lotus_RNA/analysis_result/Lt_CK_AM.deg.csv")
gene_term <- read.table("~/Genomes/Plants/Lotus_japonicus/Lotusjaponicus_MG20_v3.0_geneOntology.gaf",stringsAsFactors = FALSE,sep="\t",quote = "")
gene_term <- gene_term[gene_term$V9=="C",]
gene_term$V2 <- str_split_fixed(gene_term$V2,pattern = "\\.",n = 2)[,1]

GO_enrich <- enricher(gene = degs_out$X,pAdjustMethod = "BH",TERM2GENE = gene_term[,c(5,2)],TERM2NAME = gene_term[,c(5,10)])
dotplot(GO_enrich)
write.csv(GO_enrich@result,"MyResearch/RNA-Seq_WET/Lotus_RNA/analysis_result/Lotus.AM_go_CC.csv")
##########################################################################
#### DEG comparative analysis for ck and rz
et <- exactTest(exprSet,pair = c("CK","RZ"))
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

logFC_cutoff<-1
tTag$change<-as.factor(ifelse(tTag$FDR<0.01&abs(tTag$logFC)>logFC_cutoff,
                              ifelse(tTag$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))
save(tTag,file = "~/MyResearch/ATAC_Seq_WET/Lotus_RNA/analysis_result/Lotus_EdgeR_CKvsRZ.Rdata")
this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated gene is ',nrow(tTag[tTag$change =='UP',]) ,
                     '\nThe number of down-regulated gene is ',nrow(tTag[tTag$change =='DOWN',]))
df <- rbind(tTag[tTag$change!="NOT",][1:10,],tTag[tTag$change=="DOWN",][1:10,])

p2<-ggplot(data=tTag,aes(x=logFC,y=-log10(FDR),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(size=15,hjust=0.5),legend.position = c(0.2,0.8))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

library(tidyverse)
up_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(logFC>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(logFC< -1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out <- rbind(up_reg,down_reg) 

gene_anno <- read.delim("/home/wuzefeng/Genome.db/gff3/Lt.gene.description.txt",header = F,stringsAsFactors = F,sep="\t",comment.char = "")
degs_out$gene_annotation <- gene_anno$V2[match(rownames(degs_out),gene_anno$V1)]
write.csv(degs_out,file = "~/MyResearch/ATAC_Seq_WET/Lotus_RNA/analysis_result/Lt_CK_RZ.deg.csv")

## functions enrichment
degs_out <- read.csv("~/MyResearch/RNA-Seq_WET/Lotus_RNA/analysis_result/Lt_CK_RZ.deg.csv")
gene_term <- read.table("~/Genomes/Plants/Lotus_japonicus/Lotusjaponicus_MG20_v3.0_geneOntology.gaf",stringsAsFactors = FALSE,sep="\t",quote = "")
gene_term <- gene_term[gene_term$V9=="C",]
gene_term$V2 <- str_split_fixed(gene_term$V2,pattern = "\\.",n = 2)[,1]

GO_enrich <- enricher(gene = degs_out$X,pAdjustMethod = "BH",TERM2GENE = gene_term[,c(5,2)],TERM2NAME = gene_term[,c(5,10)])
write.csv(GO_enrich@result,"~/MyResearch/RNA-Seq_WET/Lotus_RNA/analysis_result/Lotus.RZ_go_CC.csv")


#load("DEG_GO_gprofile.Rdata")

#8B8378

#########################################################################################################
## G.max
########################################################################################################
library(stringr)
library(edgeR)

fc_count_dir <-"~/MyResearch/ATAC_Seq_WET/Gmax_RNA/2feature_count/"
fc_count_files <- grep(dir(fc_count_dir,full.names = TRUE),pattern = ".summary",invert = TRUE,value = TRUE)

fc<-c()
for (file in fc_count_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE)
  fc<-cbind(fc,temp_file[,7])
}

rownames(fc)<-temp_file$Geneid
colnames(fc) <- basename(fc_count_files)
colnames(fc)<- str_split_fixed(colnames(fc),pattern = "\\.fc",n = 2)[,1]


colnames(fc)<- apply(str_split_fixed(colnames(fc),pattern = "_",n = 4)[,c(1,2,4)],1,function(x)paste(x,collapse =  "_"))
#pheatmap::pheatmap(cor(log(fc+1,2)))

## PCA analysis
group <- str_split_fixed(colnames(fc),pattern = "_",n = 3)[,2]
MyResult.pca <- pca(t(log(fc[rowSums(cpm(fc)>=1)>=3,]+1,2)))
plotIndiv(MyResult.pca,group = group,ellipse = TRUE,title = "Samples",size.axis = 20,size.xlabel = 20,size.ylabel = 20)

## DEG analysis prepare
group_list<-factor(group) # make group list
coldata<-data.frame(row.names = colnames(fc),condition=group_list)  # make coldata

expdat = fc[rowSums(cpm(fc)>=1) >= 2,]  # filter, and 31830/57147 genes were kept
exprSet <- DGEList(counts = expdat, group = group_list) # make dgelist
exprSet <- calcNormFactors(exprSet)  # normlaized factor
exprSet <- estimateCommonDisp(exprSet) # 
exprSet <- estimateTagwiseDisp(exprSet)

## DEG comparative analysis for ck and AM
et <- exactTest(exprSet,pair = c("CK","AM"))
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

logFC_cutoff<-1
tTag$change<-as.factor(ifelse(tTag$FDR<0.01&abs(tTag$logFC)>logFC_cutoff,
                              ifelse(tTag$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))
save(tTag,file = "~/MyResearch/ATAC_Seq_WET/Gmax_RNA/analysis_result/Gm_CK_vs_AM.Rdata")


this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated gene is ',nrow(tTag[tTag$change =='UP',]) ,
                     '\nThe number of down-regulated gene is ',nrow(tTag[tTag$change =='DOWN',]))
df <- rbind(tTag[tTag$change!="NOT",][1:10,],tTag[tTag$change=="DOWN",][1:10,])

p2<-ggplot(data=tTag,aes(x=logFC,y=-log10(FDR),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(size=15,hjust=0.5),legend.position = c(0.2,0.8))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

up_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(logFC>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(logFC< -1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out <- rbind(up_reg,down_reg) 
write.csv(degs_out,"~/MyResearch/ATAC_Seq_WET/Gmax_RNA/analysis_result/DEG_AM_CK.csv",row.names = T,quote = F)

## functions enrichment
## gprofiler
go_result <- gost(query = rownames(degs_out),evcodes = T,organism = "gmax")
save(go_result,file = "~/MyResearch/ATAC_Seq_WET/Gmax_RNA/analysis_result/DEG_CK_AM.GO_gprofile.Rdata")
#load("DEG_GO_gprofile.Rdata")
go_reduce <- go_result$result[,c(10,9,11,3,16)]
colnames(go_reduce) <- c("category","ID","term","adj_pval","genes")
go_deg <- data.frame(ID=rownames(degs_out),logFC=degs_out$logFC)
circ <- circle_dat(go_reduce, go_deg)
pdf("~/MyResearch/ATAC_Seq_WET/Gmax_RNA/analysis_result/Gm_DEG_CK_AM_GO_plot.pdf",width = 10,height = 7)
GOCircle(circ,nsub = 20,label.size = 3)
dev.off()

##########################################################################
#### DEG comparative analysis for ck and rz
et <- exactTest(exprSet,pair = c("CK","RZ"))
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

logFC_cutoff<-1
tTag$change<-as.factor(ifelse(tTag$FDR<0.01&abs(tTag$logFC)>logFC_cutoff,
                              ifelse(tTag$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))
save(tTag,file = "~/MyResearch/ATAC_Seq_WET/Gmax_RNA/analysis_result/Gm_CK_vs_RZ.Rdata")

this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated gene is ',nrow(tTag[tTag$change =='UP',]) ,
                     '\nThe number of down-regulated gene is ',nrow(tTag[tTag$change =='DOWN',]))
df <- rbind(tTag[tTag$change!="NOT",][1:10,],tTag[tTag$change=="DOWN",][1:10,])

p2<-ggplot(data=tTag,aes(x=logFC,y=-log10(FDR),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(size=15,hjust=0.5),legend.position = c(0.2,0.8))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

library(tidyverse)
up_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(logFC>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(logFC< -1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out <- rbind(up_reg,down_reg) 
write.csv(degs_out,"~/MyResearch/ATAC_Seq_WET/Gmax_RNA/analysis_result/DEG_RZ_CK.csv",row.names = T,quote = F)
## functions enrichment
go_result <- gost(query = rownames(degs_out),evcodes = T,organism = "gmax")
save(go_result,file = "~/MyResearch/ATAC_Seq_WET/Gmax_RNA/analysis_result/DEG_CK_RZ.GO_gprofile.Rdata")
#load("DEG_GO_gprofile.Rdata")
go_reduce <- go_result$result[,c(10,9,11,3,16)]
colnames(go_reduce) <- c("category","ID","term","adj_pval","genes")
go_deg <- data.frame(ID=rownames(degs_out),logFC=degs_out$logFC)
circ <- circle_dat(go_reduce, go_deg)
pdf("/home/wuzefeng/MyResearch/ATAC_Seq_WET/Gmax_RNA/analysis_result/Gm_DEG_CK_RZ.pdf",width = 10,height = 7)
GOCircle(circ,nsub = 20,label.size = 3)
dev.off()


#########################################################################################################
## Mazie analysis
########################################################################################################
library(stringr)
library(edgeR)

fc_count_dir <-"~/MyResearch/ATAC_Seq_WET/Maize_RNA//2feature_count/"
fc_count_files <- grep(dir(fc_count_dir,full.names = TRUE),pattern = ".summary",invert = TRUE,value = TRUE)

fc<-c()
for (file in fc_count_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE)
  fc<-cbind(fc,temp_file[,7])
}

rownames(fc)<-temp_file$Geneid
colnames(fc) <- basename(fc_count_files)
colnames(fc)<- str_split_fixed(colnames(fc),pattern = "\\.fc",n = 2)[,1]


colnames(fc)<- apply(str_split_fixed(colnames(fc),pattern = "_",n = 4)[,c(1,2,4)],1,function(x)paste(x,collapse =  "_"))
#pheatmap::pheatmap(cor(log(fc+1,2)))

## PCA analysis
group <- str_split_fixed(colnames(fc),pattern = "_",n = 3)[,2]
MyResult.pca <- pca(t(log(fc[rowSums(cpm(fc)>=1)>=3,]+1,2)))
plotIndiv(MyResult.pca,group = group,ellipse = TRUE,title = "Samples",size.axis = 20,size.xlabel = 20,size.ylabel = 20)

## DEG analysis prepare
group_list<-factor(group) # make group list
coldata<-data.frame(row.names = colnames(fc),condition=group_list)  # make coldata

expdat = fc[rowSums(cpm(fc)>=1) >= 2,]  # filter, and 21898/52847 genes were kept
exprSet <- DGEList(counts = expdat, group = group_list) # make dgelist
exprSet <- calcNormFactors(exprSet)  # normlaized factor
exprSet <- estimateCommonDisp(exprSet) # 
exprSet <- estimateTagwiseDisp(exprSet)

## DEG comparative analysis for ck and AM
et <- exactTest(exprSet,pair = c("CK","AM"))
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

logFC_cutoff<-1
tTag$change<-as.factor(ifelse(tTag$FDR<0.01&abs(tTag$logFC)>logFC_cutoff,
                              ifelse(tTag$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))



this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated gene is ',nrow(tTag[tTag$change =='UP',]) ,
                     '\nThe number of down-regulated gene is ',nrow(tTag[tTag$change =='DOWN',]))
df <- rbind(tTag[tTag$change!="NOT",][1:10,],tTag[tTag$change=="DOWN",][1:10,])

p2<-ggplot(data=tTag,aes(x=logFC,y=-log10(FDR),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(size=15,hjust=0.5),legend.position = c(0.2,0.8))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

up_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(logFC>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(logFC< -1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out <- rbind(up_reg,down_reg) 
write.csv(degs_out,"~/MyResearch/ATAC_Seq_WET/Maize_RNA/analysis_result//DEG_AM_CK.csv",row.names = T,quote = F)
degs_out <- read.csv("~/MyResearch/RNA-Seq_WET/Maize_RNA/analysis_result/DEG_AM_CK.mod.csv",header = T)
## functions enrichment
## gprofiler
go_result <- gost(query = degs_out$X,evcodes = T,organism = "zmays")
save(go_result,file = "~/MyResearch/RNA-Seq_WET/Maize_RNA/analysis_result/DEG_CK_AM.GO_gprofile.Rdata")
write.csv(go_result$result,"~/MyResearch/RNA-Seq_WET/Maize_RNA/analysis_result/DEG_CK_AM.GO.csv")
#load("DEG_GO_gprofile.Rdata")
go_reduce <- go_result$result[,c(10,9,11,3,16)]
colnames(go_reduce) <- c("category","ID","term","adj_pval","genes")
go_deg <- data.frame(ID=degs_out$X,logFC=degs_out$logFC)
circ <- circle_dat(go_reduce, go_deg)
pdf("/home/wuzefeng/MyResearch/ATAC_Seq_WET/Maize_RNA/analysis_result/Zm_DEG_CK_AM.pdf",width = 10,height = 7)
GOCircle(circ,nsub = 20,label.size = 3)
dev.off()

#########################################################################################################
## SLy analysis
########################################################################################################
library(stringr)
library(edgeR)

fc_count_dir <-"~/MyResearch/ATAC_Seq_WET/Sly_RNA/2feature_count/"
fc_count_files <- grep(dir(fc_count_dir,full.names = TRUE),pattern = ".summary",invert = TRUE,value = TRUE)

fc<-c()
for (file in fc_count_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE)
  fc<-cbind(fc,temp_file[,7])
}

rownames(fc)<-temp_file$Geneid
colnames(fc) <- basename(fc_count_files)
colnames(fc)<- str_split_fixed(colnames(fc),pattern = "\\.fc",n = 2)[,1]


colnames(fc)<- apply(str_split_fixed(colnames(fc),pattern = "_",n = 4)[,c(1,2,4)],1,function(x)paste(x,collapse =  "_"))
#pheatmap::pheatmap(cor(log(fc+1,2)))

## PCA analysis
group <- str_split_fixed(colnames(fc),pattern = "_",n = 3)[,2]
MyResult.pca <- pca(t(log(fc[rowSums(cpm(fc)>=1)>=3,]+1,2)))
plotIndiv(MyResult.pca,group = group,ellipse = TRUE,title = "Samples",size.axis = 20,size.xlabel = 20,size.ylabel = 20)

## DEG analysis prepare
group_list<-factor(group) # make group list
coldata<-data.frame(row.names = colnames(fc),condition=group_list)  # make coldata

expdat = fc[rowSums(cpm(fc)>=1) >= 2,]  # filter, and 21898/52847 genes were kept
exprSet <- DGEList(counts = expdat, group = group_list) # make dgelist
exprSet <- calcNormFactors(exprSet)  # normlaized factor
exprSet <- estimateCommonDisp(exprSet) # 
exprSet <- estimateTagwiseDisp(exprSet)

## DEG comparative analysis for ck and AM
et <- exactTest(exprSet,pair = c("CK","AM"))
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

logFC_cutoff<-1
tTag$change<-as.factor(ifelse(tTag$FDR<0.01&abs(tTag$logFC)>logFC_cutoff,
                              ifelse(tTag$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))



this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated gene is ',nrow(tTag[tTag$change =='UP',]) ,
                     '\nThe number of down-regulated gene is ',nrow(tTag[tTag$change =='DOWN',]))
df <- rbind(tTag[tTag$change!="NOT",][1:10,],tTag[tTag$change=="DOWN",][1:10,])

p2<-ggplot(data=tTag,aes(x=logFC,y=-log10(FDR),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(size=15,hjust=0.5),legend.position = c(0.2,0.8))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

up_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(logFC>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(logFC< -1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out <- rbind(up_reg,down_reg) 
write.csv(degs_out,"~/MyResearch/ATAC_Seq_WET/Sly_RNA/analysis_result//DEG_AM_CK.csv",row.names = T,quote = F)

## functions enrichment
## gprofiler
go_result <- gost(query = rownames(degs_out),evcodes = T,organism = "slycopersicum",sources = c("GO:BP", "GO:MF", "GO:CC"))
save(go_result,file = "~/MyResearch/ATAC_Seq_WET/Sly_RNA/analysis_result/DEG_CK_AM.GO_BP_gprofile.Rdata")
#load("DEG_GO_gprofile.Rdata")
go_reduce <- go_result$result[,c(10,9,11,3,16)]
colnames(go_reduce) <- c("category","ID","term","adj_pval","genes")
go_deg <- data.frame(ID=rownames(degs_out),logFC=degs_out$logFC)
circ <- circle_dat(go_reduce, go_deg)
pdf("/home/wuzefeng/MyResearch/ATAC_Seq_WET/Sly_RNA//analysis_result/Sly_DEG_CK_AM_GO.pdf",width = 10,height = 7)
GOCircle(circ,nsub = 19,label.size = 3)
dev.off()

#########################################################################################################
## Rice analysis
########################################################################################################
library(stringr)
library(edgeR)

fc_count_dir <-"~/MyResearch/ATAC_Seq_WET/Rice_RNA/2feature_count/"
fc_count_files <- grep(dir(fc_count_dir,full.names = TRUE),pattern = ".summary",invert = TRUE,value = TRUE)

fc<-c()
for (file in fc_count_files){
  temp_file <- read.table(file,header = TRUE,stringsAsFactors = FALSE)
  fc<-cbind(fc,temp_file[,7])
}

rownames(fc)<-temp_file$Geneid
colnames(fc) <- basename(fc_count_files)
colnames(fc)<- str_split_fixed(colnames(fc),pattern = "\\.fc",n = 2)[,1]


colnames(fc)<- apply(str_split_fixed(colnames(fc),pattern = "_",n = 4)[,c(1,2,4)],1,function(x)paste(x,collapse =  "_"))
#pheatmap::pheatmap(cor(log(fc+1,2)))

## PCA analysis
group <- str_split_fixed(colnames(fc),pattern = "_",n = 3)[,2]
MyResult.pca <- pca(t(log(fc[rowSums(cpm(fc)>=1)>=3,]+1,2)))
plotIndiv(MyResult.pca,group = group,ellipse = TRUE,title = "Samples",size.axis = 20,size.xlabel = 20,size.ylabel = 20)

## DEG analysis prepare
group_list<-factor(group) # make group list
coldata<-data.frame(row.names = colnames(fc),condition=group_list)  # make coldata

expdat = fc[rowSums(cpm(fc)>=1) >= 2,]  # filter, and 21898/52847 genes were kept
exprSet <- DGEList(counts = expdat, group = group_list) # make dgelist
exprSet <- calcNormFactors(exprSet)  # normlaized factor
exprSet <- estimateCommonDisp(exprSet) # 
exprSet <- estimateTagwiseDisp(exprSet)

## DEG comparative analysis for ck and AM
et <- exactTest(exprSet,pair = c("CK","AM"))
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

logFC_cutoff<-1
tTag$change<-as.factor(ifelse(tTag$FDR<0.01&abs(tTag$logFC)>logFC_cutoff,
                              ifelse(tTag$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))



this_title <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up-regulated gene is ',nrow(tTag[tTag$change =='UP',]) ,
                     '\nThe number of down-regulated gene is ',nrow(tTag[tTag$change =='DOWN',]))
df <- rbind(tTag[tTag$change!="NOT",][1:10,],tTag[tTag$change=="DOWN",][1:10,])

p2<-ggplot(data=tTag,aes(x=logFC,y=-log10(FDR),color=change))+
  geom_point(alpha=0.8,size=1.75)+
  geom_hline(yintercept = -log10(0.01),color="#990000",linetype="dashed")+
  geom_vline(xintercept = -1,color="#990000",linetype="dashed")+
  geom_vline(xintercept = 1,color="#990000",linetype="dashed")+
  labs(x="log2 fold change")+ 
  ylab("-log10 pvalue")+
  theme_bw(base_size = 25)+
  theme(plot.title = element_text(size=15,hjust=0.5),legend.position = c(0.2,0.8))+
  scale_color_manual(values=c('orange','gray','steelblue'))+
  
  geom_label_repel(
    data = df,
    aes(label = rownames(df),
        fill=as.factor(change)),
    color="white",
    show.legend = FALSE ,
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = 'grey50',
    max.overlaps = 15
  )+
  scale_fill_manual(values=c("orange","steelblue"))

up_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(desc(logFC),FDR,change)%>%filter(logFC>1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
down_reg <- tTag%>%tibble::rownames_to_column()%>%arrange(logFC,FDR,change)%>%filter(logFC< -1,FDR<0.01)%>%tibble::column_to_rownames(var = "rowname")
degs_out <- rbind(up_reg,down_reg) 
gene_name <- read.table("~/MyResearch/ATAC_Seq_WET/Rice_RNA/analysis_result/rice_name.txt",sep="\t",header = T,quote = "")

degs_out$gene_name <- gene_name$Gene_name[match(rownames(degs_out),gene_name$Gene_ID)]
write.csv(degs_out,"~/MyResearch/ATAC_Seq_WET/Rice_RNA/analysis_result/DEG_AM_CK.csv",row.names = T,quote = F,na = "")

## functions enrichment
## gprofiler
go_result <- gost(query = rownames(degs_out),evcodes = T,organism = "osativa")
save(go_result,file = "~/MyResearch/ATAC_Seq_WET/Rice_RNA/analysis_result/DEG_CK_AM.GO_gprofile.Rdata")
#load("DEG_GO_gprofile.Rdata")
go_reduce <- go_result$result[,c(10,9,11,3,16)]
colnames(go_reduce) <- c("category","ID","term","adj_pval","genes")
go_deg <- data.frame(ID=rownames(degs_out),logFC=degs_out$logFC)
circ <- circle_dat(go_reduce, go_deg)

pdf("/home/wuzefeng/MyResearch/ATAC_Seq_WET/Rice_RNA//analysis_result/Rice_DEG_CK_AM_GO.pdf",width = 10,height = 7)
GOCircle(circ,nsub = 20,label.size = 3)
dev.off()



