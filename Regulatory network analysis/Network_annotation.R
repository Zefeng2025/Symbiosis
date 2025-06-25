# Load network
mt_network <- read.csv("Medicago_regulatory_network.csv")
mt_gene_maps <- read.table("Med.gene_map.txt",header = F)
mt_network$TFs_v3 <- mt_gene_maps$V2[match(mt_network$TFs,mt_gene_maps$V1)]
mt_network$Target_v3 <- mt_gene_maps$V2[match(mt_network$Target.genes,mt_gene_maps$V1)]
mt_network <- subset(mt_network,mt_network$TFs!="MtrunA17_Chr5g0437411")# remove PUB1 from the list of TFs in Medicago 
write.csv(mt_network,"Medicago_regulatory_network2.csv")
