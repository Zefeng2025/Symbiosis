library(ggpubr)
library(ggplot2)
library(stringr)
library(readxl)

# Medicago
load("~/MyResearch/RNA-Seq/Medicago_RNA/analysis_result/Mt_EdgeR_CKvsAM.Rdata")
Mt_CK_vs_AM <- tTag
load("~/MyResearch/RNA-Seq/Medicago_RNA/analysis_result/Mt_EdgeR_CKvsRZ.Rdata")
Mt_CK_vs_RZ <- tTag
Mt_CK_vs_RZ$change[Mt_CK_vs_RZ$logFC>1&Mt_CK_vs_RZ$FDR<0.01]<-"DOWN"
Mt_CK_vs_RZ$change[Mt_CK_vs_RZ$logFC< -1&Mt_CK_vs_RZ$FDR<0.01]<-"UP"

Mt_CK_vs_RZ$AM_fc <- Mt_CK_vs_AM$logFC[match(rownames(Mt_CK_vs_RZ),rownames(Mt_CK_vs_AM))]
cor(Mt_CK_vs_RZ$logFC,Mt_CK_vs_RZ$AM_fc) # 0.72

## Glymax
load("~/MyResearch/RNA-Seq/Gmax_RNA/analysis_result/Gm_CK_vs_AM.Rdata")
GM_CK_vs_AM <- tTag
load("~/MyResearch/RNA-Seq/Gmax_RNA/analysis_result/Gm_CK_vs_RZ.Rdata")
GM_CK_vs_RZ <- tTag

GM_CK_vs_RZ$AM_fc <- GM_CK_vs_AM$logFC[match(rownames(GM_CK_vs_RZ),rownames(GM_CK_vs_AM))]
cor(GM_CK_vs_RZ$logFC,GM_CK_vs_RZ$AM_fc) # 0.56

### Lotus
load("~/MyResearch/RNA-Seq/Lotus_RNA/analysis_result/Lotus_EdgeR_CKvsAM.Rdata")
Lt_CK_vs_AM <- tTag
load("~/MyResearch/RNA-Seq/Lotus_RNA/analysis_result/Lotus_EdgeR_CKvsRZ.Rdata")
Lt_CK_vs_RZ <- tTag

Lt_CK_vs_RZ$AM_fc <- Lt_CK_vs_AM$logFC[match(rownames(Lt_CK_vs_RZ),rownames(Lt_CK_vs_AM))]
cor(Lt_CK_vs_RZ$logFC,Lt_CK_vs_RZ$AM_fc) 

## gene name mapping 
known_genes_mapping <- read_excel("~/MyResearch/RNA-Seq/Medicago_RNA/analysis_result/Medicago_gene2name_mapping.xls")
known_genes_mapping <- known_genes_mapping[,c(1,2)]
known_genes_mapping <- na.omit(known_genes_mapping)
colnames(known_genes_mapping) <- c("Gene_name","Gene_ID")
known_genes_mapping$Gene_ID<- str_replace_all(known_genes_mapping$Gene_ID,pattern = "Medtr",replacement = "MTR_")

# gmax
gmax_gene_mapping <- read.csv("~/MyResearch/RNA-Seq/Gmax_RNA/analysis_result/gmax_info/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2022.09.21-03.48.20.02.tsv",sep = "\t")

# lotus 
lotus_gene_mapping <- read_xls("~/MyResearch/RNA-Seq/Lotus_RNA/analysis_result/lotus_info_from_unipot/lotus_gene_name.xls")

## plot
df <- rbind(data.frame(x=Mt_CK_vs_RZ$logFC,y=Mt_CK_vs_RZ$AM_fc,species="M. truncatula"),
            data.frame(x=GM_CK_vs_RZ$logFC,y=GM_CK_vs_RZ$AM_fc,species="G. max"),
            data.frame(x=Lt_CK_vs_RZ$logFC,y=Lt_CK_vs_RZ$AM_fc,species="L. japonicus"))
rownames(df)<- c(rownames(Mt_CK_vs_RZ),rownames(GM_CK_vs_RZ),rownames(Lt_CK_vs_RZ))
df$species<- factor(df$species, levels = c("M. truncatula","G. max","L. japonicus"))
df$name <- known_genes_mapping$Gene_name[match(rownames(df),known_genes_mapping$Gene_ID)]

df$name[df$species=="G. max"] <- gmax_gene_mapping$Gene_Names_primar[match(rownames(df[df$species=="G. max",]),gmax_gene_mapping$Gene_Names_ORF)]
df$name[df$species=="L. japonicus"] <- lotus_gene_mapping$name[match(rownames(df[df$species=="L. japonicus",]),lotus_gene_mapping$id)]


options(ggrepel.max.overlaps = Inf)
ggscatter(df, x = "x", y = "y", 
          add = "reg.line", #loess
          color = "species",
          palette = c("#458B00", "#EE7621", "#5CACEE"),
          facet.by = "species",
          conf.int = TRUE,
          fullrange = TRUE,
          alpha=0.3,
          label =df$name,
          label.select = df$name[df$name%in%c("MtNOD25","MtNodGRP3C","MtLb6","MtNCR306","MtLb5","MtCP1","MtRAM1","MtRAM2","MtRAD1","MtFatM","MtNIN","MtMATE15", " MtMATE14","MtnsLTP96", "MtPNP1","MtPUB1","MtMYB1","MtPt4","MtKunitz8","MtLP18","MtNCR553","MtNodGRP33","MtLb11","MtN18","MtNCR593",
                                              "PT9","PT10","GSTU10","GSTU46","Pht1-13","Chia1","GmNIN","MYB1*","FatM1*","FatM2*","RAD1_1*","RAD1_2*",
                                              "PT8*","LtNIN","Lb6-1*","Lb6-3*","Lb6-2*","PHT14","FatM*", "PLCP","Nodulin-25","ENOD16","COPT6l","ARF2L")], 
          repel = T,
          label.rectangle=T,
          font.label = c(12, "italic")) +
stat_cor(aes(color = species), 
         method = "pearson",
         label.y = 15,
         label.x = -14)+
xlab("Fold change (log) in NFS")+
ylab("Fold change (log) in AMS")+
theme(text = element_text(size = 20),
      strip.text = element_text(face = "italic"),
      legend.position = "none")+
geom_hline(yintercept = c(-1,1), linetype = "dashed", color = "gray") +
geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "gray")

###### AM indueced genes
Ri_CK_vs_AM<-read.csv("~/MyResearch/RNA-Seq/Rice_RNA/analysis_result/DEG_AM_CK.csv",header = T,row.names = 1)
Mz_CK_vs_AM <- read.csv("~/MyResearch/RNA-Seq/Maize_RNA/analysis_result/DEG_AM_CK.csv",header = T,row.names = 1)
Sy_CK_vs_AM <- read.csv("~/MyResearch/RNA-Seq/Sly_RNA/analysis_result/DEG_AM_CK.csv",header = T,row.names = 1)

CK_AM <-  rbind(data.frame(Num=nrow(Mt_CK_vs_AM[Mt_CK_vs_AM$change=="UP",]),species="M. truncatula",change="Up-regulated"),
                data.frame(Num=nrow(Mt_CK_vs_AM[Mt_CK_vs_AM$change=="DOWN",]),species="M. truncatula",change="Down-regulated"),
                data.frame(Num=nrow(GM_CK_vs_AM[GM_CK_vs_AM$change=="UP",]),species="G. max",change="Up-regulated"),
                data.frame(Num=nrow(GM_CK_vs_AM[GM_CK_vs_AM$change=="DOWN",]),species="G. max",change="Down-regulated"),
                data.frame(Num=nrow(Lt_CK_vs_AM[Lt_CK_vs_AM$change=="UP",]),species="L. japonicus",change="Up-regulated"),
                data.frame(Num=nrow(Lt_CK_vs_AM[Lt_CK_vs_AM$change=="DOWN",]),species="L. japonicus",change="Down-regulated"),
                data.frame(Num=nrow(Ri_CK_vs_AM[Ri_CK_vs_AM$change=="UP",]),species="O. sativa",change="Up-regulated"),
                data.frame(Num=nrow(Ri_CK_vs_AM[Ri_CK_vs_AM$change=="DOWN",]),species="O. sativa",change="Down-regulated"),
                data.frame(Num=nrow(Mz_CK_vs_AM[Mz_CK_vs_AM$change=="UP",]),species="Z. mays",change="Up-regulated"),
                data.frame(Num=nrow(Mz_CK_vs_AM[Mz_CK_vs_AM$change=="DOWN",]),species="Z. mays",change="Down-regulated"),
                data.frame(Num=nrow(Sy_CK_vs_AM[Sy_CK_vs_AM$change=="UP",]),species="S. lycopersicum",change="Up-regulated"),
                data.frame(Num=nrow(Sy_CK_vs_AM[Sy_CK_vs_AM$change=="DOWN",]),species="S. lycopersicum",change="Down-regulated"))

CK_AM$Num <-ifelse(CK_AM$change=="Up-regulated",CK_AM$Num,-CK_AM$Num)
ggplot(CK_AM, aes(x = species, y = Num, fill = change)) + 
  geom_bar(stat="identity", position="identity") +
  ylim(- max(CK_AM$Num), max(CK_AM$Num)) +
  scale_y_continuous(breaks = seq(-3000,3000, 1500)) +
  ylab("Gene number")+
  xlab("Species")+
  scale_fill_manual(values = c("#EE7621","#5CACEE"))+
  theme_pubr(base_size = 20,x.text.angle = 45)

#### RZ induced genes
CK_RZ <-  rbind(data.frame(Num=nrow(Mt_CK_vs_RZ[Mt_CK_vs_RZ$change=="UP",]),species="M. truncatula",change="Up-regulated"),
                data.frame(Num=nrow(Mt_CK_vs_RZ[Mt_CK_vs_RZ$change=="DOWN",]),species="M. truncatula",change="Down-regulated"),
                data.frame(Num=nrow(GM_CK_vs_RZ[GM_CK_vs_RZ$change=="UP",]),species="G. max",change="Up-regulated"),
                data.frame(Num=nrow(GM_CK_vs_RZ[GM_CK_vs_RZ$change=="DOWN",]),species="G. max",change="Down-regulated"),
                data.frame(Num=nrow(Lt_CK_vs_RZ[Lt_CK_vs_RZ$change=="UP",]),species="L. japonicus",change="Up-regulated"),
                data.frame(Num=nrow(Lt_CK_vs_RZ[Lt_CK_vs_RZ$change=="DOWN",]),species="L. japonicus",change="Down-regulated"))

CK_RZ$Num <-ifelse(CK_RZ$change=="Up-regulated",CK_RZ$Num,-CK_RZ$Num)
ggplot(CK_RZ, aes(x = species, y = Num, fill = change)) + 
  geom_bar(stat="identity", position="identity") +
  ylim(- max(CK_RZ$Num), max(CK_RZ$Num)) +
  scale_y_continuous(breaks = seq(-3000,3000, 1500)) +
  ylab("Gene number")+
  xlab("Species")+
  scale_fill_manual(values = c("#EE7621","#5CACEE"))+
  theme_pubr(base_size = 20,x.text.angle = 45,margin = 10)

#### veen
library(ggVennDiagram)
DEG_list_Mt<-list(Medicago_AM=rownames(Mt_CK_vs_AM[Mt_CK_vs_AM$change!="NOT",]),
               Medicago_RZ=rownames(Mt_CK_vs_RZ[Mt_CK_vs_RZ$change!="NOT",]))
DEG_list_Gm<-list(Gmax_AM=rownames(GM_CK_vs_AM[GM_CK_vs_AM$change!="NOT",]),
                  Gmax_RZ=rownames(GM_CK_vs_RZ[GM_CK_vs_RZ$change!="NOT",]))
DEG_list_Lt<-list(Lotus_AM=rownames(Lt_CK_vs_AM[Lt_CK_vs_AM$change!="NOT",]),
                  Lotus_RZ=rownames(Lt_CK_vs_RZ[Lt_CK_vs_RZ$change!="NOT",]))


ggv_mt <- ggVennDiagram(DEG_list_Mt,set_size = 5,set_color = "black",category.names = c("M. truncatula (AM)","M. truncatula (RZ)"))+scale_fill_gradient(low = "#66CDAA",high = "#458B00",) 
ggv_gm <- ggVennDiagram(DEG_list_Gm,set_size = 5,set_color = "black",category.names = c("G. max (AM)","G. max (RZ)"))+scale_fill_gradient(low = "#EECFA1",high = "#EE7621",) 
ggv_lt <- ggVennDiagram(DEG_list_Lt,set_size = 5,set_color = "black",category.names = c("L. japonicus (AM)","L. japonicus (RZ)"))+scale_fill_gradient(low = "#87CEFA",high = "#5CACEE",) 
ggarrange(ggv_mt,ggv_gm,ggv_lt,ncol = 1)


################################################
##################### inter-species comparisons between any two species
##################################################

## ANS
one2one_Ortholog <- read.table("~/MyResearch/RNA-Seq/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthologues/Orthologues_Medicago_truncatula/Medicago_truncatula__v__Glycine_max.txt",sep = "\t",header = F,stringsAsFactors = F) 
Mt_CK_vs_AM_with_Ortho <- subset(Mt_CK_vs_AM,rownames(Mt_CK_vs_AM)%in%one2one_Ortholog$V1)
Mt_CK_vs_AM_with_Ortho$Gm <- one2one_Ortholog$V2[match(rownames(Mt_CK_vs_AM_with_Ortho),one2one_Ortholog$V1)]
Mt_CK_vs_AM_with_Ortho$Gm_AM_logFC<-GM_CK_vs_AM$logFC[match(Mt_CK_vs_AM_with_Ortho$Gm,rownames(GM_CK_vs_AM))] # 0.31


one2one_Ortholog <- read.table("~/MyResearch/RNA-Seq/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthologues/Orthologues_Medicago_truncatula/Medicago_truncatula__v__Lotus_japonicus.txt",sep = "\t",header = F,stringsAsFactors = F) 
one2one_Ortholog$V2<-str_split_fixed(one2one_Ortholog$V2,pattern = "\\.",n = 2)[,1]
Mt_CK_vs_AM_with_Ortho <- subset(Mt_CK_vs_AM,rownames(Mt_CK_vs_AM)%in%one2one_Ortholog$V1)
Mt_CK_vs_AM_with_Ortho$Lt <- one2one_Ortholog$V2[match(rownames(Mt_CK_vs_AM_with_Ortho),one2one_Ortholog$V1)]
Mt_CK_vs_AM_with_Ortho$Lt_AM_logFC<-Lt_CK_vs_AM$logFC[match(Mt_CK_vs_AM_with_Ortho$Lt,rownames(Lt_CK_vs_AM))] # 0.05, remove non-degs, 0.26

## NFS
one2one_Ortholog <- read.table("~/MyResearch/RNA-Seq/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthologues/Orthologues_Medicago_truncatula/Medicago_truncatula__v__Glycine_max.txt",sep = "\t",header = F,stringsAsFactors = F) 
Mt_CK_vs_RZ_with_Ortho <- subset(Mt_CK_vs_RZ,rownames(Mt_CK_vs_RZ)%in%one2one_Ortholog$V1)
Mt_CK_vs_RZ_with_Ortho$Gm <- one2one_Ortholog$V2[match(rownames(Mt_CK_vs_RZ_with_Ortho),one2one_Ortholog$V1)]
Mt_CK_vs_RZ_with_Ortho<-subset(Mt_CK_vs_RZ_with_Ortho,abs(Mt_CK_vs_RZ_with_Ortho$logFC)>1)
Mt_CK_vs_RZ_with_Ortho$Gm_RZ_logFC<-GM_CK_vs_RZ$logFC[match(Mt_CK_vs_RZ_with_Ortho$Gm,rownames(GM_CK_vs_RZ))] # -0.07
cor(Mt_CK_vs_RZ_with_Ortho$logFC,Mt_CK_vs_RZ_with_Ortho$Gm_RZ_logFC,use = "complete") # -0.07, remove non-degs, -0.1

one2one_Ortholog <- read.table("~/MyResearch/RNA-Seq/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthologues/Orthologues_Medicago_truncatula/Medicago_truncatula__v__Lotus_japonicus.txt",sep = "\t",header = F,stringsAsFactors = F) 
one2one_Ortholog$V2<-str_split_fixed(one2one_Ortholog$V2,pattern = "\\.",n = 2)[,1]
Mt_CK_vs_RZ_with_Ortho <- subset(Mt_CK_vs_RZ,rownames(Mt_CK_vs_RZ)%in%one2one_Ortholog$V1)
Mt_CK_vs_RZ_with_Ortho$Lt <- one2one_Ortholog$V2[match(rownames(Mt_CK_vs_RZ_with_Ortho),one2one_Ortholog$V1)]
Mt_CK_vs_RZ_with_Ortho<-subset(Mt_CK_vs_RZ_with_Ortho,abs(Mt_CK_vs_RZ_with_Ortho$logFC)>1)
Mt_CK_vs_RZ_with_Ortho$Lt_RZ_logFC<-Lt_CK_vs_RZ$logFC[match(Mt_CK_vs_RZ_with_Ortho$Lt,rownames(Lt_CK_vs_RZ))] # 0.2, remove non-degs, 0.3
cor(Mt_CK_vs_RZ_with_Ortho$logFC,Mt_CK_vs_RZ_with_Ortho$Lt_RZ_logFC,use = "complete")



##########################################
## between species one-to-one orthologous
#########################################

## 1. by tpm

one2ones <- read.table("/home/wuzefeng/MyResearch/RNA-Seq/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthogroups/2sigle_copy.orthogrougp.mod.txt",header = F,stringsAsFactors = F)
one2ones$V2 <- str_split_fixed(one2ones$V2,pattern = "\\.",n = 2)[,1]

Medicago_tpm <- read.csv("Medicago_RNA/analysis_result/Mt_TPM.csv",header = T,row.names = 1)
Lotus_tpm <- read.csv("Lotus_RNA/analysis_result/Lt_TPM.csv",header = T,row.names = 1)
Gm_tpm <- read.csv("Gmax_RNA/analysis_result/Gm_TPM.csv",header = T,row.names = 1)
Ri_tpm <- read.csv("Rice_RNA/analysis_result/Ri_TPM.csv",header = T,row.names = 1)
Zm_tpm <- read.csv("Maize_RNA/analysis_result/Zm_TPM.csv",header = T,row.names = 1)
Sl_tpm <- read.csv("Sly_RNA/analysis_result/Sl_TPM.csv",header = T,row.names = 1)

m1 <- Medicago_tpm[str_split_fixed(one2ones$V4,pattern = "\\|",n = 2)[,2],]
l1 <- Lotus_tpm[str_split_fixed(one2ones$V2,pattern = "\\|",n = 2)[,2],]
g1 <- Gm_tpm[str_split_fixed(one2ones$V1,pattern = "\\|",n = 2)[,2],]
r1 <- Ri_tpm[str_split_fixed(one2ones$V5,pattern = "\\|",n = 2)[,2],]
z1 <- Zm_tpm[str_split_fixed(one2ones$V3,pattern = "\\|",n = 2)[,2],]
s1 <- Sl_tpm[str_split_fixed(one2ones$V6,pattern = "\\|",n = 2)[,2],]

dd<-cbind(m1,l1,g1,r1,z1,s1)
pheatmap::pheatmap(cor(t(t(scale(dd)))))

##2 by fpkm


one2ones <- read.table("/home/wuzefeng/MyResearch/RNA-Seq//1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthogroups/2sigle_copy.orthogrougp.mod.txt",header = F,stringsAsFactors = F)
one2ones$V2 <- str_split_fixed(one2ones$V2,pattern = "\\.",n = 2)[,1]

Medicago_fpkm <- read.csv("Medicago_RNA/analysis_result/Mt_FPKM.csv",header = T,row.names = 1)
Lotus_fpkm <- read.csv("Lotus_RNA/analysis_result/Lt_FPKM.csv",header = T,row.names = 1)
Gm_fpkm <- read.csv("Gmax_RNA/analysis_result/GM_FPKM.csv",header = T,row.names = 1)
Ri_fpkm <- read.csv("Rice_RNA/analysis_result/Ri_FPKM.csv",header = T,row.names = 1)
Zm_fpkm <- read.csv("Maize_RNA/analysis_result/Zm_FPKM.csv",header = T,row.names = 1)
Sl_fpkm <- read.csv("Sly_RNA/analysis_result/Sl_FPKM.csv",header = T,row.names = 1)

m1 <- Medicago_fpkm[str_split_fixed(one2ones$V4,pattern = "\\|",n = 2)[,2],]
l1 <- Lotus_fpkm[str_split_fixed(one2ones$V2,pattern = "\\|",n = 2)[,2],]
g1 <- Gm_fpkm[str_split_fixed(one2ones$V1,pattern = "\\|",n = 2)[,2],]
r1 <- Ri_fpkm[str_split_fixed(one2ones$V5,pattern = "\\|",n = 2)[,2],]
z1 <- Zm_fpkm[str_split_fixed(one2ones$V3,pattern = "\\|",n = 2)[,2],]
s1 <- Sl_fpkm[str_split_fixed(one2ones$V6,pattern = "\\|",n = 2)[,2],]

dd<-cbind(m1,l1,g1,r1,z1,s1)
p1<-pheatmap::pheatmap(cor(t(scale(t(dd)))))

## pca 
group <- str_split_fixed(colnames(dd),pattern = "_",n = 3)[,2]
MyResult.pca <- pca(t(dd))
plotIndiv(MyResult.pca,group = group,ellipse = TRUE,title = "Samples",size.axis = 20,size.xlabel = 20,size.ylabel = 20)

group <- str_split_fixed(colnames(dd),pattern = "_",n = 3)[,1]
MyResult.pca <- pca(t(dd))
p2<-plotIndiv(MyResult.pca,group = group,ellipse = F,title = "Samples",size.axis = 20,size.xlabel = 20,size.ylabel = 20)


## one2ones 
one2ones$ID <- seq(1,nrow(one2ones))
one2ones$Gm_AM <- ifelse(str_split_fixed(one2ones$V1,pattern = "\\|",n = 2)[,2]%in%rownames(GM_CK_vs_AM)[GM_CK_vs_AM$change!="NOT"],1,0) 
one2ones$Lt_AM <- ifelse(str_split_fixed(one2ones$V2,pattern = "\\|",n = 2)[,2]%in%rownames(Lt_CK_vs_AM)[Lt_CK_vs_AM$change!="NOT"],1,0) 
one2ones$Zm_AM <- ifelse(str_split_fixed(one2ones$V3,pattern = "\\|",n = 2)[,2]%in%rownames(Mz_CK_vs_AM)[Mz_CK_vs_AM$change!="NOT"],1,0) 
one2ones$Mt_AM <- ifelse(str_split_fixed(one2ones$V4,pattern = "\\|",n = 2)[,2]%in%rownames(Mt_CK_vs_AM)[Mt_CK_vs_AM$change!="NOT"],1,0) 
one2ones$Ri_AM <- ifelse(str_split_fixed(one2ones$V5,pattern = "\\|",n = 2)[,2]%in%rownames(Ri_CK_vs_AM)[Ri_CK_vs_AM$change!="NOT"],1,0) 
one2ones$Sl_AM <- ifelse(str_split_fixed(one2ones$V6,pattern = "\\|",n = 2)[,2]%in%rownames(Sy_CK_vs_AM)[Sy_CK_vs_AM$change!="NOT"],1,0) 

library(UpSetR)
upset(one2ones[,c(8:13)],nsets = 6,
      main.bar.color = "steelblue",
      nintersects = NA,
      matrix.color = "orange",
      sets.bar.color = "steelblue",
      sets.x.label = "ARGs number \n (One to one orthologs )",
      point.size = 3,
      line.size = 1,
      matrix.dot.alpha = 0.5,
      text.scale = 2)
########################################################
## ALL orthologous groups with degs
## python scripts
aa <- read.table("/home/wuzefeng/MyResearch/RNA-Seq/2cross_species_compare/AMS_orthologous.txt",header = T,row.names = 1)
aa[apply(aa,1,function(x)min(x)>=1),]

dim(aa[apply(aa,1,function(x)sum(x[seq(2,12,2)]>=1)==1),]) # 9286
dim(aa[apply(aa,1,function(x)sum(x[seq(2,12,2)]>=1)==2),]) # 2323
dim(aa[apply(aa,1,function(x)sum(x[seq(2,12,2)]>=1)==3),]) # 908
dim(aa[apply(aa,1,function(x)sum(x[seq(2,12,2)]>=1)==4),]) # 339
dim(aa[apply(aa,1,function(x)sum(x[seq(2,12,2)]>=1)==5),]) # 125
dim(aa[apply(aa,1,function(x)sum(x[seq(2,12,2)]>=1)==6),]) # 28

df <- data.frame(orthogroup=c(9286,2323,908,339,125,28),species_num=c(1,2,3,4,5,6))
library(ggpubr)
ggbarplot(df, x="species_num", y="orthogroup", fill = "steelblue", color = "white", 
                  palette = "jco",#杂志jco的配色 
                  sort.val = "desc",#下降排序 
                  sort.by.groups=FALSE,#不按组排序 
                  x.text.angle=0,
                  label = T)+
  theme(text = element_text(size = 20))



### NFS cross species comparison
#######################################################
## python scripts
aa <- read.table("/home/wuzefeng/MyResearch/RNA-Seq/2cross_species_compare/NFS_orthogous.V2",header = T,row.names = 1)  # 64595
aa <- aa[apply(aa,1,function(x)sum(x[seq(1,5,2)]>=1)==3),] # #13277


on <- nrow(aa[apply(aa,1,function(x)sum(x[seq(2,6,2)]>=1)==1),]) # 3447
tw <- nrow(aa[apply(aa,1,function(x)sum(x[seq(2,6,2)]>=1)==2),]) # 1268
thr<- nrow(aa[apply(aa,1,function(x)sum(x[seq(2,6,2)]>=1)==3),]) # 378


df <- data.frame(orthogroup=c(on,tw,thr),species_num=c(1,2,3))
library(ggpubr)
ggbarplot(df, x="species_num", y="orthogroup", fill = "steelblue", color = "white", 
          palette = "jco",#杂志jco的配色 
          sort.val = "desc",#下降排序 
          sort.by.groups=FALSE,#不按组排序 
          x.text.angle=0,
          label = T)+
  theme(text = element_text(size = 20))


## detailed stats
aa <- read.table("/home/wuzefeng/MyResearch/RNA-Seq/2cross_species_compare/NFS_orthogous",header = T,row.names = 1)  # 64595
for (m in seq(1:3)){
  orthogroups_from_m_legume <- aa[apply(aa,1,function(x)sum(x[seq(1,5,2)]>=1)==m),] # 25217
  orthogroups_from_m_legume_with_NGRs <- orthogroups_from_m_legume[apply(orthogroups_from_m_legume,1,function(x)sum(x[seq(2,6,2)]>=1)==m),] # 1939
  orthogroups_from_m_legume_without_NGRs <- orthogroups_from_m_legume[apply(orthogroups_from_m_legume,1,function(x)sum(x[seq(2,6,2)]>=1)==0),] # 23278
  
  orthogroups_from_m_legume_with_NGRs_has_gene_from_nonlegume <- orthogroups_from_m_legume_with_NGRs[apply(orthogroups_from_m_legume_with_NGRs,1,function(x)sum(x[seq(7,11,2)]>=1)>=1),]
  orthogroups_from_m_legume_with_NGRs_has_no_gene_from_nonlegume <- orthogroups_from_m_legume_with_NGRs[apply(orthogroups_from_m_legume_with_NGRs,1,function(x)sum(x[seq(7,11,2)]>=1)==0),]
  
  orthogroups_from_m_legume_without_NGRs_has_gene_from_nonlegume <- orthogroups_from_m_legume_without_NGRs[apply(orthogroups_from_m_legume_without_NGRs,1,function(x)sum(x[seq(7,11,2)]>=1)>=1),]
  orthogroups_from_m_legume_without_NGRs_has_no_gene_from_nonlegume <-orthogroups_from_m_legume_without_NGRs[apply(orthogroups_from_m_legume_without_NGRs,1,function(x)sum(x[seq(7,11,2)]>=1)==0),]
  
  message(m,"---",
          nrow(orthogroups_from_m_legume),"---",
          nrow(orthogroups_from_m_legume_with_NGRs),"---",
          nrow(orthogroups_from_m_legume_without_NGRs),"---",
          nrow(orthogroups_from_m_legume_with_NGRs_has_gene_from_nonlegume),"---",
          nrow(orthogroups_from_m_legume_with_NGRs_has_no_gene_from_nonlegume),"---",
          nrow(orthogroups_from_m_legume_without_NGRs_has_gene_from_nonlegume),"---",
          nrow( orthogroups_from_m_legume_without_NGRs_has_no_gene_from_nonlegume))
}


##
## combining the two set 
ortho_AMF <- read.table("/home/wuzefeng/MyResearch/RNA-Seq/2cross_species_compare/AMS_orthologous.txt",header = T,row.names = 1)  # 64595
#P1
ortho_AMF_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(8,12,2)]>=1)>=1),] # at least one AMF gene in at least one nonlegume 
ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)>=1),]     # at least one AMF gene in at least one legume


ortho_NFS <- read.table("/home/wuzefeng/MyResearch/RNA-Seq/2cross_species_compare/NFS_orthogous.V2",header = T,row.names = 1)
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=1),] # at least one NRG in at least one legume

a<- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_AMF_in_legume)) # at least one AMF gene in at least one non-legume  and at least one AMF gene in at least one legume
b<-intersect(a,rownames(ortho_NFS_in_legume)) # # at least one AMF gene in at least one non-legume, at least one AMF gene in at least one legume and  at least one NRG gene in at least one legume

### P2  
ortho_AMF_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(8,12,2)]>=1)>=1),]
ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)==0),]
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=1),]

a <- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_AMF_in_legume))
b <- intersect(a,rownames(ortho_NFS_in_legume))

### P3  
ortho_AMF_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(8,12,2)]>=1)==0),]
ortho_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(7,11,2)]>=1)>=1),]

ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)==0),]
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=1),]

a <- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_in_nonlegume))
b<- intersect(a,rownames(ortho_AMF_in_legume))
c <- intersect(b,rownames(ortho_NFS_in_legume))


### P4  
ortho_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(7,11,2)]>=1)==0),]

ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)==0),]
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=1),]

a <- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_AMF_in_legume))
b<- intersect(a,rownames(ortho_NFS_in_legume))



## more strick 
## combining the two set 
ortho_AMF <- read.table("/home/wuzefeng/MyResearch/RNA-Seq/2cross_species_compare/AMS_orthologous.txt",header = T,row.names = 1)  # 64595
#P1
ortho_AMF_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(8,12,2)]>=1)>=2),]
ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)>=2),]


ortho_NFS <- read.table("/home/wuzefeng/MyResearch/RNA-Seq/2cross_species_compare/NFS_orthogous.V2",header = T,row.names = 1)
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=2),]

a<- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_AMF_in_legume))
b<-intersect(a,rownames(ortho_NFS_in_legume))

### P2  
ortho_AMF_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(8,12,2)]>=1)>=2),]
ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)==0),]
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=2),]

a <- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_AMF_in_legume))
b <- intersect(a,rownames(ortho_NFS_in_legume))

### P3  
ortho_AMF_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(8,12,2)]>=1)==0),]
ortho_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(7,11,2)]>=1)>=2),]

ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)==0),]
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=2),]

a <- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_in_nonlegume))
b<- intersect(a,rownames(ortho_AMF_in_legume))
c <- intersect(b,rownames(ortho_NFS_in_legume))


### P4  
ortho_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(7,11,2)]>=1)==0),]

ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)==0),]
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=2),]

a <- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_AMF_in_legume))
b<- intersect(a,rownames(ortho_NFS_in_legume))


## more more strick

## combining the two set 
ortho_AMF <- read.table("/home/wuzefeng/MyResearch/RNA-Seq/2cross_species_compare/AMS_orthologous.txt",header = T,row.names = 1)  # 64595
#P1
ortho_AMF_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(8,12,2)]>=1)>=3),]
ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)>=3),]


ortho_NFS <- read.table("/home/wuzefeng/MyResearch/RNA-Seq/2cross_species_compare/NFS_orthogous.V2",header = T,row.names = 1)
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=3),]

a<- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_AMF_in_legume))
b<-intersect(a,rownames(ortho_NFS_in_legume))

### P2  
ortho_AMF_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(8,12,2)]>=1)>=3),]
ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)==0),]
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=3),]

a <- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_AMF_in_legume))
b <- intersect(a,rownames(ortho_NFS_in_legume))

### P3  
ortho_AMF_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(8,12,2)]>=1)==0),]
ortho_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(7,11,2)]>=1)>=3),]

ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)==0),]
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=3),]

a <- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_in_nonlegume))
b<- intersect(a,rownames(ortho_AMF_in_legume))
c <- intersect(b,rownames(ortho_NFS_in_legume))


### P4  
ortho_in_nonlegume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(7,11,2)]>=1)==0),]

ortho_AMF_in_legume <- ortho_AMF[apply(ortho_AMF,1,function(x)sum(x[seq(2,6,2)]>=1)==0),]
ortho_NFS_in_legume <- ortho_NFS[apply(ortho_NFS,1,function(x)sum(x[seq(2,6,2)]>=1)>=3),]

a <- intersect(rownames(ortho_AMF_in_nonlegume),rownames(ortho_AMF_in_legume))
b<- intersect(a,rownames(ortho_NFS_in_legume))







