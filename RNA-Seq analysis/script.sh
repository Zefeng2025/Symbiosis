##1. Medicago truncatula
# 1. mapping
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x /home/wuzefeng/Genome.db/genome.db/Medicago_truncatula.MedtrA17_4.0 -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2. reads count
mkdir 2feature_count
for m in $(ls 1hisat2_out); do echo $m; featureCounts -a ~/Genome.db/gtf.db/Medicago_truncatula.MedtrA17_4.0.51.gtf -o 2feature_count/${m%.sam}.fc  -p -T 12 1hisat2_out/$m; done

## get FPKM and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genome.db/gtf.db/Medicago_truncatula.MedtrA17_4.0.51.gtf -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done



##################################################################################
## 2. Lotus

# 1. mapping
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x ~/Genome.db/genome.db/Lotusjaponicus_MG20_v3.0_genome.fa -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2. obtain read counts
mkdir 2feature_count
cd 3sorted_bam
for m in $(ls *.bam); do echo $m; featureCounts -a ~/Genome.db/gff3/Lotus.gtf  -o ../2feature_count/${m%.bam}.fc -p -T 12 $m; done

## get Fpkm and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genomes/Plants/Lotus_japonicus/lotusjaponicus_mg20_v3.0_annotations.gff3  -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done 

for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genomes/Plants/Lotus_japonicus/Lotus.gtf  -e -B -A 5stringtie_gtf/$(basename ${m%.sorted.bam}).genes.tsv -o 5stringtie_gtf/$(basename ${m%.sorted.bam}).gff -p 10; done

#########################################################################################
## 3. Glycine. max

# 1. mapping
cd /home/wuzefeng/MyResearch/RNA-Seq/Gmax_RNA
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x ~/Genome.db/genome.db/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2.  obtain read counts

mkdir 2feature_count
for m in $(ls 1hisat2_out); do echo $m; featureCounts -a ~/Genome.db/gtf.db/Glycine_max.Glycine_max_v2.1.51.gtf -o 2feature_count/${m%.sam}.fc  -p -T 12 1hisat2_out/$m; done

## get Fpkm and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genome.db/gtf.db/Glycine_max.Glycine_max_v2.1.51.gtf  -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done


## 4. Zea mays
# 1. mapping
cd /home/wuzefeng/MyResearch/RNA-Seq/Maize_RNA
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x ~/Genome.db/genome.db/Zea_mays.B73_RefGen_v4.dna.toplevel.fa -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2.  obtain read counts

mkdir 2feature_count
for m in $(ls 1hisat2_out); do echo $m; featureCounts -a ~/Genome.db/gtf.db/Zea_mays.B73_RefGen_v4.50.gtf -o 2feature_count/${m%.sam}.fc  -p -T 12 1hisat2_out/$m; done

# get Fpkm and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genome.db/gtf.db/Zea_mays.B73_RefGen_v4.50.gtf   -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done


## 5. Rice
# 1. mapping
cd /home/wuzefeng/MyResearch/RNA-Seq/Rice_RNA
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x ~/Genome.db/genome.db/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2. obtain read counts

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

## 6. tomato
# 1. mapping 
cd /home/wuzefeng/MyResearch/RNA-Seq/Sly_RNA
mkdir 1hisat2_out
for m in $(ls); do echo $m; hisat2 -x ~/Genome.db/genome.db/Solanum_lycopersicum.SL3.0.dna.toplevel.fa -1 $m/*_1.clean.fq.gz -2 $m/*_2.clean.fq.gz --no-mixed --no-discordant --summary-file $m -p 12 --time -S 1hisat2_out/$m.sam ; done

# 2. obtain read counts

mkdir 2feature_count
for m in $(ls 1hisat2_out); do echo $m; featureCounts -a ~/Genome.db/gtf.db/Solanum_lycopersicum.SL3.0.51.gtf -o 2feature_count/${m%.sam}.fc  -p -T 12 1hisat2_out/$m; done

## get Fpkm and TPM using stringtie
# 3. Hisat sam to bam 
mkdir 3sorted_bam
for m in $(ls 1hisat2_out/*.sam); do echo $m; samtools view -u $m | samtools sort -o 3sorted_bam/$(basename ${m%.sam}).sorted.bam -@ 10; done

# 4. Stringtie
mkdir 4stringtie
for m in $(ls 3sorted_bam/*.bam); do echo $m; stringtie $m  -G ~/Genome.db/gtf.db/Solanum_lycopersicum.SL3.0.51.gtf  -e -B -A 4stringtie/$(basename ${m%.sorted.bam}).genes.tsv -o 4stringtie/$(basename ${m%.sorted.bam}).gff -p 10; done




