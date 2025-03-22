##1. 目蓿
# 1. 比对
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x /home/wuzefeng/Genome.db/genome.db/Medicago_truncatula.MedtrA17_4.0 -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2. 得到reads count进行差异表达分析
mkdir 2feature_count
for m in $(ls 1hisat2_out); do echo $m; featureCounts -a ~/Genome.db/gtf.db/Medicago_truncatula.MedtrA17_4.0.51.gtf -o 2feature_count/${m%.sam}.fc  -p -T 12 1hisat2_out/$m; done

## get Fpkm and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genome.db/gtf.db/Medicago_truncatula.MedtrA17_4.0.51.gtf -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done
#注：stringtie每个样本输出的信息表达信息或基因名字顺序不是一致的


##################################################################################
## 2. 百脉跟Lotus

# 1. 比对
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x ~/Genome.db/genome.db/Lotusjaponicus_MG20_v3.0_genome.fa -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2. 得到reads count进行差异表达分析 (怀疑gtf文件有问题) awk '$3=="gene"' Lotusjaponicus_MG20_v3.0_annotations.gff3 |grep -v "Non-coding"  | grep -v "tRNA" | grep -v "rRNA"| grep -v sequencetype=repeat | awk '$1!="chloro"'| awk '$1!="mito"' | wc -l
cat Lotusjaponicus_MG20_v3.0_annotations.gff3 | awk '$1!="chloro"'| awk '$1!="mito"' | grep -v "Non-coding"   > Lotusjaponicus.gff3 # 过滤非编码和叶绿体、线粒体基因,过滤会丢掉一些gene或exon
用自己脚本转化成gtf文件

mkdir 2feature_count
# for m in $(ls 1hisat2_out); do echo $m; featureCounts -a ~/Genome.db/gff3/Lotusjaponicus.gtf -o 2feature_count/${m%.sam}.fc  -p -T 12 1hisat2_out/$m; done # 所有基因.hisat结果已经删除,只能用bam文件
cd 3sorted_bam
for m in $(ls *.bam); do echo $m; featureCounts -a ~/Genome.db/gff3/Lotus.gtf  -o ../2feature_count/${m%.bam}.fc -p -T 12 $m; done

## get Fpkm and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genomes/Plants/Lotus_japonicus/lotusjaponicus_mg20_v3.0_annotations.gff3  -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done # 重新运行22.10.31,ccamk的问题,通过IGV查看，有些基因的fpkm和tpm统计不正确

for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genomes/Plants/Lotus_japonicus/Lotus.gtf  -e -B -A 5stringtie_gtf/$(basename ${m%.sorted.bam}).genes.tsv -o 5stringtie_gtf/$(basename ${m%.sorted.bam}).gff -p 10; done

#########################################################################################
## 3. 大豆Gmax

# 1. 比对
cd /home/wuzefeng/MyResearch/ATAC_Seq_WET/Gmax_RNA
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x ~/Genome.db/genome.db/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2. 得到reads count进行差异表达分析

mkdir 2feature_count
for m in $(ls 1hisat2_out); do echo $m; featureCounts -a ~/Genome.db/gtf.db/Glycine_max.Glycine_max_v2.1.51.gtf -o 2feature_count/${m%.sam}.fc  -p -T 12 1hisat2_out/$m; done

## get Fpkm and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genome.db/gtf.db/Glycine_max.Glycine_max_v2.1.51.gtf  -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done


## 4. 玉米
# 1. 比对
cd /home/wuzefeng/MyResearch/ATAC_Seq_WET/Maize_RNA
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x ~/Genome.db/genome.db/Zea_mays.B73_RefGen_v4.dna.toplevel.fa -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2. 得到reads count进行差异表达分析

mkdir 2feature_count
for m in $(ls 1hisat2_out); do echo $m; featureCounts -a ~/Genome.db/gtf.db/Zea_mays.B73_RefGen_v4.50.gtf -o 2feature_count/${m%.sam}.fc  -p -T 12 1hisat2_out/$m; done

## get Fpkm and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genome.db/gtf.db/Zea_mays.B73_RefGen_v4.50.gtf   -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done


## 5. 水稻
# 1. 比对
cd /home/wuzefeng/MyResearch/ATAC_Seq_WET/Rice_RNA
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x ~/Genome.db/genome.db/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2. 得到reads count进行差异表达分析

mkdir 2feature_count
for m in $(ls 1hisat2_out); do echo $m; featureCounts -a ~/Genome.db/gtf.db/Oryza_sativa.IRGSP-1.0.51.gtf -o 2feature_count/${m%.sam}.fc  -p -T 12 1hisat2_out/$m; done

## get Fpkm and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genome.db/gtf.db/Oryza_sativa.IRGSP-1.0.51.gtf  -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done
for m in $(ls *.tsv); do echo $m; sed -i 's/\#/_/g' $m; done # some genes have "#" in gene names, and convert them into "_"

## 6. 番茄
# 1. 比对
cd /home/wuzefeng/MyResearch/ATAC_Seq_WET/Sly_RNA
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x ~/Genome.db/genome.db/Solanum_lycopersicum.SL3.0.dna.toplevel.fa -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2. 得到reads count进行差异表达分析

mkdir 2feature_count
for m in $(ls 1hisat2_out); do echo $m; featureCounts -a ~/Genome.db/gtf.db/Solanum_lycopersicum.SL3.0.51.gtf -o 2feature_count/${m%.sam}.fc  -p -T 12 1hisat2_out/$m; done

## get Fpkm and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genome.db/gtf.db/Solanum_lycopersicum.SL3.0.51.gtf  -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done



## 鉴定直系同源基因
conda activate orthofinder
orthofinder  # 构建的物种树有问题
## 提取单拷贝基因group
for m in $(cat /home/wuzefeng/MyResearch/ATAC_Seq_WET/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthogroups/Orthogroups_SingleCopyOrthologues.txt); do echo $m; grep $m /home/wuzefeng/MyResearch/ATAC_Seq_WET/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthogroups/Orthogroups.txt>>/home/wuzefeng/MyResearch/ATAC_Seq_WET/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthogroups/1sigle_copy.orthogrougp.txt;done

### 1. ete3根据单拷贝基因构建物种树（## 线程数太多，t_coffee运行有问题，重启）
## 修改orthogroup.txt格式
sed -e s'/OG[0-9]*: //g' -e s'/ /\t/g' 1sigle_copy.orthogrougp.txt > 2sigle_copy.orthogrougp.mod.txt
## 修改fasta文件header(ete3不需要)
# awk '/^>/{u=split($0,m,"|");$0=">"m[2]}1' all.pep.fasta > all.pep.mod.header.fasta
## 用单拷贝基因建树(需要修改othogroup为tab分割，而且不需要OG的ids，基于clustalomega比对，然后合并文件，用fastatree最大似然树，构树正确！！！)
ete3 build -w clustalo_default-trimal01-none-none -m cog_100-alg_concat_default-fasttree_default -a -a 0proteins/all.pep.fasta --cogs 0proteins/OrthoFinder/Results_Sep27/Orthogroups/2sigle_copy.orthogrougp.mod.txt  -o 2ete3_sptree/ --spname-delimiter "|" -C 10
# -w: 格式为: aligner-trimmer-tester-builder （结果一致）
ete3 build -w clustalo_default-trimal01-none-none -m cog_100-alg_concat_default-iqtree_bestmodel -a -a 0proteins/all.pep.fasta --cogs 0proteins/OrthoFinder/Results_Sep27/Orthogroups/2sigle_copy.orthogrougp.mod.txt  -o 2ete3_sptree2/ --spname-delimiter "|" -C 14

###2. 自己构建物种树 （不管那个流程仍然不对）
2.1 # mattf比对单拷贝基因（也可以选clustal omega; muscle，PRANK等，也可以选近单拷贝基因）
for m in $(ls 0proteins/OrthoFinder/Results_Sep27/Single_Copy_Orthologue_Sequences/*.fa); do mafft --auto --thread 6 $m >2species_tree_self/2.1mafft_align_single_copy_genes/$(basename ${m%.fa}).fasta; done
### 2.2 trim 或者Gblock比对过滤
for m in $(ls 2species_tree_self/2.1mafft_align_single_copy_genes/*.fasta); do echo $m;  trimal -in $m -out 2species_tree_self/2.2trimal_out/$(basename $m) -keepheader -fasta -automated1; done
### 2.3 合并alignment（也可以选seqkit， biopyton脚本）
python concat.align.py

#### 2.3 iqtree 基于super-align构建最大似然树（构建物种树用的少，~30min;树的结果不对，无根树）
iqtree -s 2species_tree_self/2.3merge_trimed_align/superalgn.phy -nt 10 -bb 1000  # 比较慢，先检测546个蛋白模型。

### 2.4 用raxml基于super-align比对结果构建最大似然树（时间长~17h，可能对DNA更合适,结果也不对）
raxmlHPC-PTHREADS -f a -m PROTGAMMAGTR -n raxml_out -# 100 -x 1234 -p 1234 -s 2species_tree_self/2.3merge_trimed_align/superalgn.phy -T 10 -w /home/wuzefeng/MyResearch/ATAC_Seq_WET/1orthofinder/2species_tree_self/2.4raxml_out/ # 注意氨基酸和DNA的-m参数不同

### 2.5 用fasttree基于super-align比对结果构建物种树（速度快）
FastTree 2.3merge_trimed_align/superalgn.phy > species.tree

### 2.6 用基因树推断物种树(先用raxml构建多个基因树+anstral流程等推断物种树)
#2.6.1 修改多序列比对header
for m in $(ls 2species_tree_self/2.2trimal_out/*.fasta); do echo $m;  awk '/^>/{u=split($0,m,"|");$0=m[1]}1' $m >2species_tree_self/2.2trimal_out_mod/$(basename $m); done
#2.6.2 比较费时间
for m in $(ls 2species_tree_self/2.2trimal_out_mod/*.fasta); do echo $m; raxmlHPC-PTHREADS -x 1234 -N 100 -f a -m PROTGAMMAJTTF -n $(basename ${m%.fasta}) -s $m -T 10 -w /home/wuzefeng/MyResearch/ATAC_Seq_WET/1orthofinder/2species_tree_self/2.5gene_tree_raxml/; done #
# Astral 推断物种树 （也不对）
java -jar astra 3.1.1.jar -wq -in [input tree]



########### ATAC-seq 分析
##1. 目蓿 
#bowtie2 建立索引并比对
for m in $(ls *.R1.fastq.gz) ; do echo $m; prefix=${m%.R1.fastq.gz}; bowtie2 -p 12 -x ~/Genome.db/genome.db/Medicago_truncatula.MedtrA17_4.0.dna.toplevel.fa -1 ${prefix}.R1.fastq.gz -2 ${prefix}.R2.fastq.gz  --no-unal --no-mixed --no-discordant -S /media/wuzefeng/软件/1ATAC_seq/2raw_alignement_results/${prefix}.sam >>bt2_mapping_log; done

