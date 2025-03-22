#! /bin/sh

## orthogroups identification

conda activate orthofinder
orthofinder -f adjust_pep/ -t 10 -M msa  
## extract single-copy gene orthogroups
for m in $(cat /home/wuzefeng/MyResearch/RNA-Seq/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthogroups/Orthogroups_SingleCopyOrthologues.txt); do echo $m; grep $m /home/wuzefeng/MyResearch/RNA-Seq/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthogroups/Orthogroups.txt>>/home/wuzefeng/MyResearch/RNA-Seq/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthogroups/1sigle_copy.orthogrougp.txt;done



## ete3 build species tree.
##  modify orthogroup.txt 
sed -e s'/OG[0-9]*: //g' -e s'/ /\t/g' 1sigle_copy.orthogrougp.txt > 2sigle_copy.orthogrougp.mod.txt
##  build species tree using ete3 (method1)
ete3 build -w clustalo_default-trimal01-none-none -m cog_100-alg_concat_default-fasttree_default -a -a 0proteins/all.pep.fasta --cogs 0proteins/OrthoFinder/Results_Sep27/Orthogroups/2sigle_copy.orthogrougp.mod.txt  -o 2ete3_sptree/ --spname-delimiter "|" -C 10
##  build species tree using ete3 (method2)
ete3 build -w clustalo_default-trimal01-none-none -m cog_100-alg_concat_default-iqtree_bestmodel -a -a 0proteins/all.pep.fasta --cogs 0proteins/OrthoFinder/Results_Sep27/Orthogroups/2sigle_copy.orthogrougp.mod.txt  -o 2ete3_sptree2/ --spname-delimiter "|" -C 14


