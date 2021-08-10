#!/bin/bash
BASEDIR="/scratch/kld57880/methods_paper2"

mkdir $BASEDIR/RNA_seq
mkdir $BASEDIR/K9_chip
mkdir $BASEDIR/K4_K27_chip
mkdir $BASEDIR/cutNrun_mods
mkdir $BASEDIR/cutNrun_pol2

##download relevant reference genomes and annotations
##zebrafish:
curl -s http://ftp.ensembl.org/pub/release-103/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz | gunzip -c > $BASEDIR/danio_refseq.fa
curl -s http://ftp.ensembl.org/pub/release-103/gtf/danio_rerio/Danio_rerio.GRCz11.103.gtf.gz | gunzip -c > $BASEDIR/danio_refann.gtf
##spike in genomes for the cutNrun:
curl -s ftp://ftp.ensemblgenomes.org/pub/fungi/release-48/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz | gunzip -c > $BASEDIR/sacc_refseq.fa
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $BASEDIR/ecoli_refseq.fa
#bed file of gene coordinates downloaded from http://genome.ucsc.edu/cgi-bin/hgTables?org=zebrafish
########and saved as $BASEDIR/genes.bed
#gtf file of repeat-masker downloaded from http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1121023327_WLh4AYGK4EabW4ybav66jwsq6NNA&clade=vertebrate&org=Zebrafish&db=danRer11&hgta_group=varRep&hgta_track=refSeqComposite&hgta_table=0&hgta_regionType=genome&position=chr8%3A23%2C404%2C195-23%2C406%2C194&hgta_outputType=gff&hgta_outFileName=
#######repeat-masker gtf was merged with ref_ann.gtf and saved as $BASEDIR/unmasked_ref_ann.gtf
#file with just chr number and the size is saved as $BASEDIR/chrNameLength.txt

########################
####RNA SEQ ANALYSIS####
########################

mkdir $BASEDIR/RNA_seq/raw

wget -O $BASEDIR/RNA_seq/raw/shield_1_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/001/ERR1442601/ERR1442601_1.fastq.gz"
wget -O $BASEDIR/RNA_seq/raw/shield_1_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/001/ERR1442601/ERR1442601_2.fastq.gz"
wget -O $BASEDIR/RNA_seq/raw/shield_2_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/002/ERR1442602/ERR1442602_1.fastq.gz"
wget -O $BASEDIR/RNA_seq/raw/shield_2_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/002/ERR1442602/ERR1442602_2.fastq.gz"
wget -O $BASEDIR/RNA_seq/raw/shield_3_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/003/ERR1442603/ERR1442603_1.fastq.gz"
wget -O $BASEDIR/RNA_seq/raw/shield_3_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/003/ERR1442603/ERR1442603_2.fastq.gz"
wget -O $BASEDIR/RNA_seq/raw/shield_4_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/004/ERR1442604/ERR1442604_1.fastq.gz"
wget -O $BASEDIR/RNA_seq/raw/shield_4_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/004/ERR1442604/ERR1442604_2.fastq.gz"
wget -O $BASEDIR/RNA_seq/raw/shield_5_1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/005/ERR1442605/ERR1442605_1.fastq.gz"
wget -O $BASEDIR/RNA_seq/raw/shield_5_2.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/005/ERR1442605/ERR1442605_2.fastq.gz"

###trimming raw data
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
mkdir $BASEDIR/RNA_seq/trimmed

trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/RNA_seq/trimmed --paired $BASEDIR/RNA_seq/raw/shield_1_1.fastq.gz $BASEDIR/RNA_seq/raw/shield_1_2.fastq.gz
trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/RNA_seq/trimmed --paired $BASEDIR/RNA_seq/raw/shield_2_1.fastq.gz $BASEDIR/RNA_seq/raw/shield_2_2.fastq.gz
trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/RNA_seq/trimmed --paired $BASEDIR/RNA_seq/raw/shield_3_1.fastq.gz $BASEDIR/RNA_seq/raw/shield_3_2.fastq.gz
trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/RNA_seq/trimmed --paired $BASEDIR/RNA_seq/raw/shield_4_1.fastq.gz $BASEDIR/RNA_seq/raw/shield_4_2.fastq.gz
trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/RNA_seq/trimmed --paired $BASEDIR/RNA_seq/raw/shield_5_1.fastq.gz $BASEDIR/RNA_seq/raw/shield_5_2.fastq.gz

##alignment with STAR
module load STAR/2.7.3a-GCC-8.3.0

mkdir $BASEDIR/RNA_seq/genome
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $BASEDIR/RNA_seq/genome --genomeFastaFiles $BASEDIR/danio_refseq.fa --sjdbGTFfile $BASEDIR/danio_refann.gtf --sjdbOverhang 100

mkdir $BASEDIR/RNA_seq/bams

STAR --runThreadN 20 --genomeDir $BASEDIR/RNA_seq/genome --outFileNamePrefix $BASEDIR/RNA_seq/bams/shield_1 --readFilesCommand zcat --readFilesIn $BASEDIR/RNA_seq/trimmed/shield_1_1_val_1.fq.gz $BASEDIR/RNA_seq/trimmed/shield_1_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
STAR --runThreadN 20 --genomeDir $BASEDIR/RNA_seq/genome --outFileNamePrefix $BASEDIR/RNA_seq/bams/shield_2 --readFilesCommand zcat --readFilesIn $BASEDIR/RNA_seq/trimmed/shield_2_1_val_1.fq.gz $BASEDIR/RNA_seq/trimmed/shield_2_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
STAR --runThreadN 20 --genomeDir $BASEDIR/RNA_seq/genome --outFileNamePrefix $BASEDIR/RNA_seq/bams/shield_3 --readFilesCommand zcat --readFilesIn $BASEDIR/RNA_seq/trimmed/shield_3_1_val_1.fq.gz $BASEDIR/RNA_seq/trimmed/shield_3_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
STAR --runThreadN 20 --genomeDir $BASEDIR/RNA_seq/genome --outFileNamePrefix $BASEDIR/RNA_seq/bams/shield_4 --readFilesCommand zcat --readFilesIn $BASEDIR/RNA_seq/trimmed/shield_4_1_val_1.fq.gz $BASEDIR/RNA_seq/trimmed/shield_4_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1
STAR --runThreadN 20 --genomeDir $BASEDIR/RNA_seq/genome --outFileNamePrefix $BASEDIR/RNA_seq/bams/shield_5 --readFilesCommand zcat --readFilesIn $BASEDIR/RNA_seq/trimmed/shield_5_1_val_1.fq.gz $BASEDIR/RNA_seq/trimmed/shield_5_2_val_2.fq.gz --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1

#merging replicates
module load SAMtools/1.10-iccifort-2019.5.281
samtools merge $BASEDIR/RNA_seq/bams/shield_merged.bam $BASEDIR/RNA_seq/bams/shield_1Aligned.sortedByCoord.out.bam $BASEDIR/RNA_seq/bams/shield_2Aligned.sortedByCoord.out.bam $BASEDIR/RNA_seq/bams/shield_3Aligned.sortedByCoord.out.bam $BASEDIR/RNA_seq/bams/shield_4Aligned.sortedByCoord.out.bam $BASEDIR/RNA_seq/bams/shield_5Aligned.sortedByCoord.out.bam

##TPM Calculating to determine which genes are "on" - using TPM cutoff of 0.5
##then making into bed files for bedtools intersect later
##have to cd into desired directory for the TPM calc output to be in the right place
cd /scratch/kld57880/methods_paper/RNA_seq/

module load TPMCalculator/0.0.4
TPMCalculator -b $BASEDIR/RNA_seq/bams/shield_merged.bam -g $BASEDIR/danio_refann.gtf -p -q 20

awk -F'\t' '$7>0.5' shield_merged_genes.out > $BASEDIR/RNA_seq/shield_expressed_genes.out
awk '{print $2 "\t" $3 "\t" $4 }' $BASEDIR/RNA_seq/shield_expressed_genes.out > $BASEDIR/RNA_seq/shield_expressed_genes.bed

grep "exon" shield_merged_genes.ent > shield_merged_exons.ent
awk -F'\t' '$9>0.5' shield_merged_exons.ent > $BASEDIR/RNA_seq/shield_expressed_exons.ent
awk '{print $2 "\t" $5 "\t" $6 }' $BASEDIR/RNA_seq/shield_expressed_exons.ent > $BASEDIR/RNA_seq/shield_expressed_exons.bed

cd ~

########################
####K9 CHIP ANALYSIS####
########################

#raw data in $BASEDIR/K9_chip/raw
#trimming reads
mkdir $BASEDIR/K9_chip/trimmed

for infile in $BASEDIR/K9_chip/raw/*fastq.gz
do
  trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/K9_chip/trimmed $infile
done

#aligning using Bowtie2 and filtering for MAPQ score of at least 20
module load Bowtie2/2.4.1-GCC-8.3.0

bowtie2-build $BASEDIR/danio_refseq.fa $BASEDIR/danio_ref

mkdir $BASEDIR/K9_chip/bams
for infile in $BASEDIR/K9_chip/trimmed/*.gz
do
  base=$(basename ${infile} _trimmed.fq.gz)
  bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/danio_ref -U $infile | samtools view -bq 20 | samtools sort - > $BASEDIR/K9_chip/bams/$base.sort.bam
done

###calling peaks
module load Homer/4.11-foss-2019b
mkdir $BASEDIR/K9_chip/peaks

for infile in $BASEDIR/K9_chip/bams/*.bam
  do base=$(basename ${infile} .bam)
  makeTagDirectory $BASEDIR/K9_chip/peaks/$base.tagdir $infile
done

for infile in $BASEDIR/K9_chip/peaks/K9*.tagdir
do
  base=$(basename ${infile} .tagdir)
  findPeaks $infile -style histone -minDist 1000 -i $BASEDIR/K9_chip/peaks/input_6hpf_1.sort.tagdir -gsize 1.5e9 -o $BASEDIR/K9_chip/peaks/$base.txt
done

for infile in $BASEDIR/K9_chip/peaks/*.txt
do
  base=$(basename ${infile} .txt)
  pos2bed.pl $infile > $BASEDIR/K9_chip/peaks/$base.peaks.bed
done

###intersecting peaks so we only use peaks present in both replicates
module load BEDTools/2.29.2-GCC-8.3.0

bedtools intersect -a $BASEDIR/K9_chip/peaks/K9_6hpf_1.sort.peaks.bed -b $BASEDIR/K9_chip/peaks/K9_6hpf_2.sort.peaks.bed -wa > $BASEDIR/K9_chip/peaks/wt_K9_AB_peaks.bed
bedtools intersect -a $BASEDIR/K9_chip/peaks/K9_6hpf_2.sort.peaks.bed -b $BASEDIR/K9_chip/peaks/K9_6hpf_1.sort.peaks.bed  -wa > $BASEDIR/K9_chip/peaks/wt_K9_BA_peaks.bed
cat $BASEDIR/K9_chip/peaks/wt_K9_AB_peaks.bed $BASEDIR/K9_chip/peaks/wt_K9_BA_peaks.bed | bedtools sort -i stdin |bedtools merge -i stdin > $BASEDIR/K9_chip/peaks/wt_K9_intpeaks.bed

##making into bigwigs
module load deepTools/3.3.1-intel-2019b-Python-3.7.4
mkdir $BASEDIR/K9_chip/bws

for infile in $BASEDIR/K9_chip/bams/*.bam
do
  base=$(basename ${infile} .bam)
  samtools index $infile - > $BASEDIR/K9_chip/bams/$base.bam.bai
done

for infile in $BASEDIR/K9_chip/bams/*.bam
do
  base=$(basename ${infile} .sort.bam)
  bamCoverage -b $infile -bs 10 -p 20 -o $BASEDIR/K9_chip/bws/$base.bw
done

###Make a Black List
###using the K9 chip input to generate a black list of sticky genome regions
for infile in $BASEDIR/K9_chip/peaks/input*.tagdir
do
  base=$(basename ${infile} .tagdir)
  findPeaks $infile -style factor -o $BASEDIR/K9_chip/peaks/$base.txt
done

for infile in $BASEDIR/K9_chip/peaks/input*.txt
do
  base=$(basename ${infile} .txt)
  pos2bed.pl $infile > $BASEDIR/K9_chip/peaks/$base.peaks.bed
done

bedtools intersect -a $BASEDIR/K9_chip/peaks/input_6hpf_1.sort.peaks.bed -b $BASEDIR/K9_chip/peaks/input_6hpf_2.sort.peaks.bed -wa > $BASEDIR/K9_chip/peaks/wt_I_AB_peaks.bed
bedtools intersect -a $BASEDIR/K9_chip/peaks/input_6hpf_2.sort.peaks.bed -b $BASEDIR/K9_chip/peaks/input_6hpf_1.sort.peaks.bed  -wa > $BASEDIR/K9_chip/peaks/wt_I_BA_peaks.bed
cat $BASEDIR/K9_chip/peaks/wt_I_AB_peaks.bed $BASEDIR/K9_chip/peaks/wt_I_BA_peaks.bed | bedtools sort -i stdin |bedtools merge -i stdin > $BASEDIR/K9_chip/peaks/blacklist.bed
###gives us ~500 peaks that are sticky in the genome and we will use this bed file as a BL file in deeptools

########################
##K4/K27 CHIP ANALYSIS##
########################

##downloading the K4/K27 data
module load SRA-Toolkit/2.10.8-centos_linux64

wget -O $BASEDIR/K4_K27_chip/raw/1kcell_I "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-13/SRR6729466/SRR6729466.1"
fastq-dump --split-files --gzip --outdir $BASEDIR/K4_K27_chip/raw $BASEDIR/K4_K27_chip/raw/1kcell_I
wget -O $BASEDIR/K4_K27_chip/raw/1kcell_K27_2 "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-13/SRR6729460/SRR6729460.1"
fastq-dump --split-files --gzip --outdir $BASEDIR/K4_K27_chip/raw $BASEDIR/K4_K27_chip/raw/1kcell_K27_2
wget -O $BASEDIR/K4_K27_chip/raw/1kcell_K27_1 "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-13/SRR6729459/SRR6729459.1"
fastq-dump --split-files --gzip --outdir $BASEDIR/K4_K27_chip/raw $BASEDIR/K4_K27_chip/raw/1kcell_K27_1
wget -O $BASEDIR/K4_K27_chip/raw/1kcell_K4_2 "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-13/SRR6729449/SRR6729449.1"
fastq-dump --split-files --gzip --outdir $BASEDIR/K4_K27_chip/raw $BASEDIR/K4_K27_chip/raw/1kcell_K4_2
wget -O $BASEDIR/K4_K27_chip/raw/1kcell_K4_1 "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-13/SRR6729448/SRR6729448.1"
fastq-dump --split-files --gzip --outdir $BASEDIR/K4_K27_chip/raw $BASEDIR/K4_K27_chip/raw/1kcell_K4_1

##trimming data
for infile in $BASEDIR/K4_K27_chip/raw/*fastq.gz
do
  trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/K4_K27_chip/trimmed $infile
done

##when going to align this data, error message says pairs aren't lined up so I will repair the data before moving forward
mkdir $BASEDIR/K4_K27_chip/repaired

module load BBMap/38.83-GCC-8.3.0
repair.sh in=$BASEDIR/K4_K27_chip/trimmed/1kcell_I_1_trimmed.fq.gz in2=$BASEDIR/K4_K27_chip/trimmed/1kcell_I_2_trimmed.fq.gz out=$BASEDIR/K4_K27_chip/repaired/1kcell_I_1_repaired.fq.gz out2=$BASEDIR/K4_K27_chip/repaired/1kcell_I_2_repaired.fq.gz
repair.sh in=$BASEDIR/K4_K27_chip/trimmed/1kcell_K27_1_1_trimmed.fq.gz in2=$BASEDIR/K4_K27_chip/trimmed/1kcell_K27_1_2_trimmed.fq.gz out=$BASEDIR/K4_K27_chip/repaired/1kcell_K27_1_1_repaired.fq.gz out2=$BASEDIR/K4_K27_chip/repaired/1kcell_K27_1_2_repaired.fq.gz
repair.sh in=$BASEDIR/K4_K27_chip/trimmed/1kcell_K27_2_1_trimmed.fq.gz in2=$BASEDIR/K4_K27_chip/trimmed/1kcell_K27_2_2_trimmed.fq.gz out=$BASEDIR/K4_K27_chip/repaired/1kcell_K27_2_1_repaired.fq.gz out2=$BASEDIR/K4_K27_chip/repaired/1kcell_K27_2_2_repaired.fq.gz
repair.sh in=$BASEDIR/K4_K27_chip/trimmed/1kcell_K4_1_1_trimmed.fq.gz in2=$BASEDIR/K4_K27_chip/trimmed/1kcell_K4_1_2_trimmed.fq.gz out=$BASEDIR/K4_K27_chip/repaired/1kcell_K4_1_1_repaired.fq.gz out2=$BASEDIR/K4_K27_chip/repaired/1kcell_K4_1_2_repaired.fq.gz
repair.sh in=$BASEDIR/K4_K27_chip/trimmed/1kcell_K4_2_1_trimmed.fq.gz in2=$BASEDIR/K4_K27_chip/trimmed/1kcell_K4_2_2_trimmed.fq.gz out=$BASEDIR/K4_K27_chip/repaired/1kcell_K4_2_1_repaired.fq.gz out2=$BASEDIR/K4_K27_chip/repaired/1kcell_K4_2_2_repaired.fq.gz

##now the data is fixed, align using bowtie2 and filer for MAPQ score of 20
mkdir $BASEDIR/K4_K27_chip/bams
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/danio_ref -1 $BASEDIR/K4_K27_chip/repaired/1kcell_I_1_repaired.fq.gz -2 $BASEDIR/K4_K27_chip/repaired/1kcell_I_2_repaired.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/K4_K27_chip/bams/1kcell_I.fs.bam
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/danio_ref -1 $BASEDIR/K4_K27_chip/repaired/1kcell_K27_1_1_repaired.fq.gz -2 $BASEDIR/K4_K27_chip/repaired/1kcell_K27_1_2_repaired.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/K4_K27_chip/bams/1kcell_K27_1.fs.bam
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/danio_ref -1 $BASEDIR/K4_K27_chip/repaired/1kcell_K27_2_1_repaired.fq.gz -2 $BASEDIR/K4_K27_chip/repaired/1kcell_K27_2_2_repaired.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/K4_K27_chip/bams/1kcell_K27_2.fs.bam
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/danio_ref -1 $BASEDIR/K4_K27_chip/repaired/1kcell_K4_1_1_repaired.fq.gz -2 $BASEDIR/K4_K27_chip/repaired/1kcell_K4_1_2_repaired.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/K4_K27_chip/bams/1kcell_K4_1.fs.bam
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/danio_ref -1 $BASEDIR/K4_K27_chip/repaired/1kcell_K4_2_1_repaired.fq.gz -2 $BASEDIR/K4_K27_chip/repaired/1kcell_K4_2_2_repaired.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/K4_K27_chip/bams/1kcell_K4_2.fs.bam

##calling peaks using HOMER
mkdir $BASEDIR/K4_K27_chip/peaks

for infile in $BASEDIR/K4_K27_chip/bams/*.bam
  do base=$(basename ${infile} .bam)
  makeTagDirectory $BASEDIR/K4_K27_chip/peaks/$base.tagdir $infile
done

for infile in $BASEDIR/K4_K27_chip/peaks/*_K*.tagdir
do
  base=$(basename ${infile} .tagdir)
  findPeaks $infile -style histone -minDist 1000 -i $BASEDIR/K4_K27_chip/peaks/1kcell_I.fs.tagdir -gsize 1.5e9 -o $BASEDIR/K4_K27_chip/peaks/$base.txt
done

for infile in $BASEDIR/K4_K27_chip/peaks/*.txt
do
  base=$(basename ${infile} .txt)
  pos2bed.pl $infile > $BASEDIR/K4_K27_chip/peaks/$base.peaks.bed
done

##intersecting peaks so that we only have peaks found in both replicates
bedtools intersect -a $BASEDIR/K4_K27_chip/peaks/1kcell_K27_1.fs.peaks.bed -b $BASEDIR/K4_K27_chip/peaks/1kcell_K27_2.fs.peaks.bed -wa > $BASEDIR/K4_K27_chip/peaks/1kcell_K27_AB_peaks.bed
bedtools intersect -a $BASEDIR/K4_K27_chip/peaks/1kcell_K27_2.fs.peaks.bed -b $BASEDIR/K4_K27_chip/peaks/1kcell_K27_1.fs.peaks.bed  -wa > $BASEDIR/K4_K27_chip/peaks/1kcell_K27_BA_peaks.bed
cat $BASEDIR/K4_K27_chip/peaks/1kcell_K27_AB_peaks.bed $BASEDIR/K4_K27_chip/peaks/1kcell_K27_BA_peaks.bed | bedtools sort -i stdin |bedtools merge -i stdin > $BASEDIR/K4_K27_chip/peaks/1kcell_K27_intpeaks.bed

bedtools intersect -a $BASEDIR/K4_K27_chip/peaks/1kcell_K4_1.fs.peaks.bed -b $BASEDIR/K4_K27_chip/peaks/1kcell_K4_2.fs.peaks.bed -wa > $BASEDIR/K4_K27_chip/peaks/1kcell_K4_AB_peaks.bed
bedtools intersect -a $BASEDIR/K4_K27_chip/peaks/1kcell_K4_2.fs.peaks.bed -b $BASEDIR/K4_K27_chip/peaks/1kcell_K4_1.fs.peaks.bed  -wa > $BASEDIR/K4_K27_chip/peaks/1kcell_K4_BA_peaks.bed
cat $BASEDIR/K4_K27_chip/peaks/1kcell_K4_AB_peaks.bed $BASEDIR/K4_K27_chip/peaks/1kcell_K4_BA_peaks.bed | bedtools sort -i stdin |bedtools merge -i stdin > $BASEDIR/K4_K27_chip/peaks/1kcell_K4_intpeaks.bed

##now intersecting peaks to find differential and shared peaks between K4/K27 to identify bivalent domains
bedtools intersect -a $BASEDIR/K4_K27_chip/peaks/1kcell_K4_intpeaks.bed -b $BASEDIR/K4_K27_chip/peaks/1kcell_K27_intpeaks.bed -v > $BASEDIR/K4_K27_chip/peaks/1kcell_K4_NOT_K27_peaks.bed
bedtools intersect -a $BASEDIR/K4_K27_chip/peaks/1kcell_K27_intpeaks.bed -b $BASEDIR/K4_K27_chip/peaks/1kcell_K4_intpeaks.bed -v > $BASEDIR/K4_K27_chip/peals/1kcell_K27_NOT_K4_peaks.bed

bedtools intersect -a $BASEDIR/K4_K27_chip/peaks/1kcell_K4_intpeaks.bed -b $BASEDIR/K4_K27_chip/peaks/1kcell_K27_intpeaks.bed -wa > $BASEDIR/K4_K27_chip/peaks/1kcell_K4_OVERLAP_K27_peaks.bed
bedtools intersect -a $BASEDIR/K4_K27_chip/peaks/1kcell_K27_intpeaks.bed -b $BASEDIR/K4_K27_chip/peaks/1kcell_K4_intpeaks.bed -wa > $BASEDIR/K4_K27_chip/peaks/1kcell_K27_OVERLAP_K4_peaks.bed
cat $BASEDIR/K4_K27_chip/peaks/1kcell_K4_OVERLAP_K27_peaks.bed $BASEDIR/K4_K27_chip/peaks/1kcell_K27_OVERLAP_K4_peaks.bed | bedtools sort -i stdin |bedtools merge -i stdin > $BASEDIR/K4_K27_chip/peaks/1kcell_bivalent_intpeaks.bed

##making bigwigs for data visualization
mkdir $BASEDIR/K4_K27_chip/bws

for infile in $BASEDIR/K4_K27_chip/bams/*.bam
do
  base=$(basename ${infile} .bam)
  samtools index $infile - > $BASEDIR/K4_K27_chip/bams/$base.bam.bai
done

for infile in $BASEDIR/K4_K27_chip/bams/*.bam
do
  base=$(basename ${infile} .sort.bam)
  bamCoverage -b $infile -bs 10 -p 20 -o $BASEDIR/K4_K27_chip/bws/$base.bw
done

########################
####CUT&RUN ANALYSIS####
#########POL 2##########
########################

##starting with all my raw files in $BASEDIR/raw
mkdir $BASEDIR/cutNrun_pol2/trimmed

#trimming our pol2 reads
trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/cutNrun_pol2/trimmed --paired $BASEDIR/cutNrun_pol2/raw/wt_1_pol2_R1.fastq.gz $BASEDIR/cutNrun_pol2/raw/wt_1_pol2_R2.fastq.gz
trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/cutNrun_pol2/trimmed --paired $BASEDIR/cutNrun_pol2/raw/wt_3_pol2_R1.fastq.gz $BASEDIR/cutNrun_pol2/raw/wt_3_pol2_R2.fastq.gz
trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/cutNrun_pol2/trimmed --paired $BASEDIR/cutNrun_pol2/raw/wt_2_IgG_R1.fastq.gz $BASEDIR/cutNrun_pol2/raw/wt_2_IgG_R2.fastq.gz

#mapping our reads to the zebrafish and the spike in genomes
mkdir $BASEDIR/cutNrun_pol2/bams

bowtie2-build $BASEDIR/ecoli_refseq.fa $BASEDIR/ecoli_ref

bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/danio_ref -1 $BASEDIR/cutNrun_pol2/trimmed/wt_1_pol2_R1_val_1.fq.gz -2 $BASEDIR/cutNrun_pol2/trimmed/wt_1_pol2_R2_val_2.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/cutNrun_pol2/bams/wt_1_pol2.danio_sort.bam
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/danio_ref -1 $BASEDIR/cutNrun_pol2/trimmed/wt_3_pol2_R1_val_1.fq.gz -2 $BASEDIR/cutNrun_pol2/trimmed/wt_3_pol2_R2_val_2.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/cutNrun_pol2/bams/wt_3_pol2.danio_sort.bam
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/danio_ref -1 $BASEDIR/cutNrun_pol2/trimmed/wt_2_IgG_R1_val_1.fq.gz -2 $BASEDIR/cutNrun_pol2/trimmed/wt_2_IgG_R2_val_2.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/cutNrun_pol2/bams/wt_2_IgG.danio_sort.bam

bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/ecoli_ref -1 $BASEDIR/cutNrun_pol2/trimmed/wt_1_pol2_R1_val_1.fq.gz -2 $BASEDIR/cutNrun_pol2/trimmed/wt_1_pol2_R2_val_2.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/cutNrun_pol2/bams/wt_1_pol2.ecoli_sort.bam
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/ecoli_ref -1 $BASEDIR/cutNrun_pol2/trimmed/wt_3_pol2_R1_val_1.fq.gz -2 $BASEDIR/cutNrun_pol2/trimmed/wt_3_pol2_R2_val_2.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/cutNrun_pol2/bams/wt_3_pol2.ecoli_sort.bam
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/ecoli_ref -1 $BASEDIR/cutNrun_pol2/trimmed/wt_2_IgG_R1_val_1.fq.gz -2 $BASEDIR/cutNrun_pol2/trimmed/wt_2_IgG_R2_val_2.fq.gz | samtools view -bq20 | samtools sort - > $BASEDIR/cutNrun_pol2/bams/wt_2_IgG.ecoli_sort.bam

#Now we need to extract all the aligned reads in preperation for spike in normalization
for infile in $BASEDIR/cutNrun_pol2/bams/*.bam
do
  base=$(basename ${infile} .bam)
  bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $BASEDIR/cutNrun_pol2/$base.btb.bed
done

#Now lets do the spike in calibration with RPM normalization
mkdir $BASEDIR/cutNrun_pol2/bedgraphs

sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_pol2/wt_1_pol2.danio_sort.btb.bed $BASEDIR/cutNrun_pol2/wt_1_pol2.ecoli_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_pol2/bedgraphs/wt_1_pol2.norm.bga
sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_pol2/wt_3_pol2.danio_sort.btb.bed $BASEDIR/cutNrun_pol2/wt_3_pol2.ecoli_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_pol2/bedgraphs/wt_3_pol2.norm.bga
sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_pol2/wt_2_IgG.danio_sort.btb.bed $BASEDIR/cutNrun_pol2/wt_2_IgG.ecoli_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_pol2/bedgraphs/wt_2_IgG.norm.bga

#lets make these bedgraphs into bigwigs for data visualization
module load ucsc/359

for infile in $BASEDIR/cutNrun_pol2/bedgraphs/*norm.bga
do
  base=$(basename ${infile} .norm.bga)
  bedSort $infile $BASEDIR/cutNrun_pol2/bedgraphs/$base.norm_sort.bga
done

mkdir $BASEDIR/cutNrun_pol2/bws
for infile in $BASEDIR/cutNrun_pol2/bedgraphs/*_sort.bga
do
  base=$(basename ${infile} .bga)
  bedGraphToBigWig $infile $BASEDIR/chrNameLength.txt $BASEDIR/cutNrun_pol2/bws/$base.bw
done

#calling peaks using HOMER
mkdir $BASEDIR/cutNrun_pol2/peaks

#first turning normalized bedgraphs into bed files for peak calling
for infile in $BASEDIR/cutNrun_pol2/bedgraphs/*.norm.bga
  do base=$(basename ${infile} .norm.bga)
  cat $infile | awk '{print $1 "\t" $2 "\t" $3 "\t" "+" "\t" "+" "\t" "+"}' > $BASEDIR/cutNrun_pol2/peaks/$base.bgato.bed
done

#then making our tag directories and actually calling peaks
for infile in $BASEDIR/cutNrun_pol2/peaks/*bgato.bed
  do base=$(basename ${infile} .bgato.bed)
  makeTagDirectory $BASEDIR/cutNrun_pol2/peaks/$base.BtB.tagdir $infile -format bed
done

for infile in $BASEDIR/cutNrun_pol2/peaks/*.tagdir
do
  base=$(basename ${infile} .tagdir)
  findPeaks $infile -style factor -gsize 1.5e9 -i $BASEDIR/cutNrun_pol2/peaks/wt_2_IgG.BtB.tagdir -o $BASEDIR/cutNrun_pol2/peaks/$base.txt
done

for infile in $BASEDIR/cutNrun_pol2/peaks/*.txt
do
  base=$(basename ${infile} .txt)
  pos2bed.pl $infile > $BASEDIR/cutNrun_pol2/peaks/$base.peaks.bed
done

#intersecting peaks between our two reps
bedtools intersect -a $BASEDIR/cutNrun_pol2/peaks/wt_1_pol2.BtB.peaks.bed -b $BASEDIR/cutNrun_pol2/peaks/wt_3_pol2.BtB.peaks.bed -wa > $BASEDIR/cutNrun_pol2/peaks/wt_pol2_AB_peaks.bed
bedtools intersect -a $BASEDIR/cutNrun_pol2/peaks/wt_3_pol2.BtB.peaks.bed -b $BASEDIR/cutNrun_pol2/peaks/wt_1_pol2.BtB.peaks.bed -wa > $BASEDIR/cutNrun_pol2/peaks/wt_pol2_BA_peaks.bed
cat  $BASEDIR/cutNrun_pol2/peaks/wt_pol2_AB_peaks.bed $BASEDIR/cutNrun_pol2/peaks/wt_pol2_BA_peaks.bed | bedtools sort -i stdin |bedtools merge -i stdin > $BASEDIR/cutNrun_pol2/peaks/wt_pol2_intpeaks.bed

#annotate to find peaks within 1000bp of TSS
annotatePeaks.pl $BASEDIR/cutNrun_pol2/peaks/wt_pol2_intpeaks.bed danRer11 -gtf $BASEDIR/danio_refann.gtf > $BASEDIR/cutNrun_pol2/peaks/wt_pol2_mask_ann.txt

#now filtering for only peaks that are w/i 1000bps of their annotation:
awk -F'\t' 'sqrt($10*$10) <=1000' $BASEDIR/cutNrun_pol2/peaks/wt_pol2_mask_ann.txt > $BASEDIR/cutNrun_pol2/peaks/wt_pol2.1000bp_ann.txt
awk '{print $2 "\t" $3 "\t" $4 }' $BASEDIR/cutNrun_pol2/peaks/wt_pol2.1000bp_ann.txt | tail -n +2 > $BASEDIR/cutNrun_pol2/peaks/wt_pol2.1000bp_ann.bed

#intersecting peaks with expressed genes to determine % association with active transcription
bedtools intersect -a $BASEDIR/cutNrun_pol2/peaks/wt_pol2.1000bp_ann.bed -b $BASEDIR/RNA_seq/shield_expressed_genes.bed -wa > $BASEDIR/cutNrun_pol2/peaks/pol2_genesON.bed

##making pol2 figures
mkdir $BASEDIR/figs
computeMatrix scale-regions -S $BASEDIR/cutNrun_pol2/bws/wt_1_pol2.norm_sort.bw $BASEDIR/cutNrun_pol2/bws/wt_3_pol2.norm_sort.bw  -R $BASEDIR/genes.bed -p 20 -b 2000 -a 2000 -bs=10 --missingDataAsZero -out $BASEDIR/cutNrun_pol2/pol2_genes_reps.gz
plotHeatmap -m $BASEDIR/cutNrun_pol2/pol2_genes_reps.gz -out $BASEDIR/figs/pol2_genes_reps.heatmap.pdf --whatToShow 'heatmap and colorbar' --colorMap=Oranges --regionsLabel genes --samplesLabel pol-II pol-II

computeMatrix scale-regions -S $BASEDIR/cutNrun_pol2/bws/wt_1_pol2.norm_sort.bw  -R $BASEDIR/genes.bed -p 20 -b 2000 -a 2000 -bs=10 --missingDataAsZero -out $BASEDIR/cutNrun_pol2/pol2_genes_rep1.gz
plotProfile -m $BASEDIR/cutNrun_pol2/pol2_genes_rep1.gz -out $BASEDIR/figs/pol2_genes_rep1.profile.pdf --colors orange

########################
####CUT&RUN ANALYSIS####
######HISTONE MODS######
########################

#starting with raw files in $BASEDIR/cutNrun_mods/raw
#trimming the data
for infile in $BASEDIR/cutNrun_mods/raw/*fastq.gz
do
  trim_galore --phred33 --fastqc --illumina --length 20 --output_dir $BASEDIR/cutNrun_mods/trimmed $infile
done

#aligning to danio and yeast genomes
bowtie2-build $BASEDIR/sacc_refseq.fa $BASEDIR/sacc_ref

mkdir $BASEDIR/cutNrun_mods/bams

for infile in $BASEDIR/cutNrun_mods/trimmed/*.gz
do
  base=$(basename ${infile} _trimmed.fq.gz)
  bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/danio_ref -U $infile | samtools view -bq 20 | samtools sort - > $BASEDIR/cutNrun_mods/bams/$base.danio_sort.bam
done

for infile in $BASEDIR/cutNrun_mods/trimmed/*.gz
do
  base=$(basename ${infile} _trimmed.fq.gz)
  bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $BASEDIR/sacc_ref -U $infile | samtools view -bq 20 | samtools sort - > $BASEDIR/cutNrun_mods/bams/$base.sacc_sort.bam
done

#Now we need to extract all the aligned reads
for infile in $BASEDIR/cutNrun_mods/bams/*.bam
do
  base=$(basename ${infile} .bam)
  bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $BASEDIR/cutNrun_mods/$base.btb.bed
done

#Now lets do the spike in calibration with RPM normalization
mkdir $BASEDIR/cutNrun_mods/bedgraphs

sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_mods/wt_K27_1ab.danio_sort.btb.bed $BASEDIR/cutNrun_mods/wt_K27_1ab.sacc_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_mods/bedgraphs/wt_K27_1ab.norm.bga
sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_mods/wt_K27_2ab.danio_sort.btb.bed $BASEDIR/cutNrun_mods/wt_K27_2ab.sacc_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_mods/bedgraphs/wt_K27_2ab.norm.bga
sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_mods/wt_K27_1.danio_sort.btb.bed $BASEDIR/cutNrun_mods/wt_K27_1.sacc_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_mods/bedgraphs/wt_K27_1.norm.bga
sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_mods/wt_K27_2.danio_sort.btb.bed $BASEDIR/cutNrun_mods/wt_K27_2.sacc_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_mods/bedgraphs/wt_K27_2.norm.bga
sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_mods/wt_K4_1.danio_sort.btb.bed $BASEDIR/cutNrun_mods/wt_K4_1.sacc_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_mods/bedgraphs/wt_K4_1.norm.bga
sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_mods/wt_K4_2.danio_sort.btb.bed $BASEDIR/cutNrun_mods/wt_K4_2.sacc_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_mods/bedgraphs/wt_K4_2.norm.bga
sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_mods/wt_K9_1.danio_sort.btb.bed $BASEDIR/cutNrun_mods/wt_K9_1.sacc_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_mods/bedgraphs/wt_K9_1.norm.bga
sh ~/Git2/spike_in_calibration.kd.sh $BASEDIR/cutNrun_mods/wt_K9_2.danio_sort.btb.bed $BASEDIR/cutNrun_mods/wt_K9_2.sacc_sort.btb.bed 100000 bga $BASEDIR/chrNameLength.txt 1 1000 $BASEDIR/cutNrun_mods/bedgraphs/wt_K9_2.norm.bga

#lets make these bedgraphs into bigwigs
for infile in $BASEDIR/cutNrun_mods/bedgraphs/*norm.bga
do
  base=$(basename ${infile} .norm.bga)
  bedSort $infile $BASEDIR/cutNrun_mods/bedgraphs/$base.norm_sort.bga
done

mkdir $BASEDIR/cutNrun_mods/bws
for infile in $BASEDIR/cutNrun_mods/bedgraphs/*_sort.bga
do
  base=$(basename ${infile} .bga)
  bedGraphToBigWig $infile $BASEDIR/chrNameLength.txt $BASEDIR/cutNrun_mods/bws/$base.bw
done

#calling peaks
mkdir $BASEDIR/cutNrun_mods/peaks

for infile in $BASEDIR/cutNrun_mods/bedgraphs/*.norm.bga
  do base=$(basename ${infile} .norm.bga)
  cat $infile | awk '{print $1 "\t" $2 "\t" $3 "\t" "+" "\t" "+" "\t" "+"}' > $BASEDIR/cutNrun_mods/peaks/$base.bgato.bed
done

for infile in $BASEDIR/cutNrun_mods/peaks/*bgato.bed
  do base=$(basename ${infile} .bgato.bed)
  makeTagDirectory $BASEDIR/cutNrun_mods/peaks/$base.BtB.tagdir $infile -format bed
done

##using IgG as input
for infile in $BASEDIR/cutNrun_mods/peaks/*.tagdir
do
  base=$(basename ${infile} .tagdir)
  findPeaks $infile -style histone -minDist 1000 -i $BASEDIR/cutNrun_pol2/peaks/wt_2_IgG.BtB.tagdir -F 6 -gsize 1.5e9 -fdr 0.0001 -o $BASEDIR/cutNrun_mods/peaks/$base.txt
done

for infile in $BASEDIR/cutNrun_mods/peaks/*.txt
do
  base=$(basename ${infile} .txt)
  pos2bed.pl $infile > $BASEDIR/cutNrun_mods/peaks/$base.peaks.bed
done

#intersecting peaks so we only have peaks found in both replicates
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K27_1.BtB.peaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K27_2.BtB.peaks.bed -wa > $BASEDIR/cutNrun_mods/peaks/wt_K27_AB_peaks.bed
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K27_2.BtB.peaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K27_1.BtB.peaks.bed  -wa > $BASEDIR/cutNrun_mods/peaks/wt_K27_BA_peaks.bed
cat $BASEDIR/cutNrun_mods/peaks/wt_K27_AB_peaks.bed $BASEDIR/cutNrun_mods/peaks/wt_K27_BA_peaks.bed | bedtools sort -i stdin |bedtools merge -i stdin > $BASEDIR/cutNrun_mods/peaks/wt_K27_intpeaks.bed

bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K4_1.BtB.peaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K4_2.BtB.peaks.bed -wa > $BASEDIR/cutNrun_mods/peaks/wt_K4_AB_peaks.bed
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K4_2.BtB.peaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K4_1.BtB.peaks.bed  -wa > $BASEDIR/cutNrun_mods/peaks/wt_K4_BA_peaks.bed
cat $BASEDIR/cutNrun_mods/peaks/wt_K4_AB_peaks.bed $BASEDIR/cutNrun_mods/peaks/wt_K4_BA_peaks.bed | bedtools sort -i stdin |bedtools merge -i stdin > $BASEDIR/cutNrun_mods/peaks/wt_K4_intpeaks.bed

bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K9_1.BtB.peaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K9_2.BtB.peaks.bed -wa > $BASEDIR/cutNrun_mods/peaks/wt_K9_AB_peaks.bed
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K9_2.BtB.peaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K9_1.BtB.peaks.bed  -wa > $BASEDIR/cutNrun_mods/peaks/wt_K9_BA_peaks.bed
cat $BASEDIR/cutNrun_mods/peaks/wt_K9_AB_peaks.bed $BASEDIR/cutNrun_mods/peaks/wt_K9_BA_peaks.bed | bedtools sort -i stdin |bedtools merge -i stdin > $BASEDIR/cutNrun_mods/peaks/wt_K9_intpeaks.bed

#annotate peaks with masked gtf
annotatePeaks.pl $BASEDIR/cutNrun_mods/peaks/wt_K27_intpeaks.bed danRer11 -gtf $BASEDIR/danio_refann.gtf > $BASEDIR/cutNrun_mods/peaks/wt_K27_masked_ann.txt
annotatePeaks.pl $BASEDIR/cutNrun_mods/peaks/wt_K4_intpeaks.bed danRer11 -gtf $BASEDIR/danio_refann.gtf > $BASEDIR/cutNrun_mods/peaks/wt_K4_masked_ann.txt
annotatePeaks.pl $BASEDIR/cutNrun_mods/peaks/wt_K9_intpeaks.bed danRer11 -gtf $BASEDIR/danio_refann.gtf > $BASEDIR/cutNrun_mods/peaks/wt_K9_masked_ann.txt

##identifying differential K4 and K27 peaks
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K4_intpeaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K27_intpeaks.bed -v > $BASEDIR/cutNrun_mods/peaks/K4_NOT_K27_peaks.bed
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K27_intpeaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K4_intpeaks.bed -v > $BASEDIR/cutNrun_mods/peaks/K27_NOT_K4_peaks.bed

#identifying bivalent domains
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K4_intpeaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K27_intpeaks.bed -wa -u > $BASEDIR/cutNrun_mods/peaks/K4_OVERLAP_K27_peaks.bed
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K27_intpeaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K4_intpeaks.bed -wa -u > $BASEDIR/cutNrun_mods/peaks/K27_OVERLAP_K4_peaks.bed
cat $BASEDIR/cutNrun_mods/peaks/K4_OVERLAP_K27_peaks.bed $BASEDIR/cutNrun_mods/peaks/K27_OVERLAP_K4_peaks.bed | bedtools sort -i stdin | bedtools merge -i stdin > $BASEDIR/cutNrun_mods/peaks/bivalent_domains.bed

#identifying overlap between CnR bivalent domains and ChIP bivalent domains
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/bivalent_domains.bed -b $BASEDIR/K4_K27_chip/peaks/1kcell_bivalent_intpeaks.bed -wa -u > $BASEDIR/cutNrun_mods/peaks/shield_1kcell_bivalent_overlap.bed
bedtools intersect -a $BASEDIR/K4_K27_chip/peaks/1kcell_bivalent_intpeaks.bed -b $BASEDIR/cutNrun_mods/peaks/bivalent_domains.bed -wa -u > $BASEDIR/cutNrun_mods/peaks/1kcell_shield_bivalent_overlap.bed
cat $BASEDIR/cutNrun_mods/peaks/shield_1kcell_bivalent_overlap.bed $BASEDIR/cutNrun_mods/peaks/1kcell_shield_bivalent_overlap.bed | bedtools sort -i stdin | bedtools merge -i stdin > $BASEDIR/cutNrun_mods/peaks/shield_AND_1kcell_bivalent_domains.bed

bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/bivalent_domains.bed -b $BASEDIR/K4_K27_chip/peaks/1kcell_bivalent_intpeaks.bed -v > $BASEDIR/cutNrun_mods/peaks/shield_NOT_1kcell_bivalent.bed
bedtools intersect -a $BASEDIR/K4_K27_chip/peaks/1kcell_bivalent_intpeaks.bed -b $BASEDIR/cutNrun_mods/peaks/bivalent_domains.bed -v > $BASEDIR/cutNrun_mods/peaks/1kcell_NOT_shield_bivalent.bed

#annotating differential and bivalent domains
annotatePeaks.pl $BASEDIR/cutNrun_mods/peaks/bivalent_domains.bed danRer11 -gtf $BASEDIR/danio_refann.gtf > $BASEDIR/cutNrun_mods/peaks/bivalent_domains_masked_ann.txt
annotatePeaks.pl $BASEDIR/cutNrun_mods/peaks/K27_NOT_K4_peaks.bed danRer11 -gtf $BASEDIR/danio_refann.gtf > $BASEDIR/cutNrun_mods/peaks/K27only_masked_ann.txt
annotatePeaks.pl $BASEDIR/cutNrun_mods/peaks/K4_NOT_K27_peaks.bed danRer11 -gtf $BASEDIR/danio_refann.gtf > $BASEDIR/cutNrun_mods/peaks/K4only_masked_ann.txt

#now filtering for only peaks that are w/i 1000bps of their annotation:
for infile in $BASEDIR/cutNrun_mods/peaks/*_masked_ann.txt
do
  base=$(basename ${infile} _masked_ann.txt)
  awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $BASEDIR/cutNrun_mods/peaks/$base.1000bp_ann.txt
done

#converting into bed file for intersection with TPM Calc output
for infile in $BASEDIR/cutNrun_mods/peaks/*.1000bp_ann.txt
do
  base=$(basename ${infile} .txt)
  awk '{print $2 "\t" $3 "\t" $4 }' $infile | tail -n +2 > $BASEDIR/cutNrun_mods/peaks/$base.bed
done

#intersect differential and bivalent K4/K27 domains that are w/i 1000bps of a TSS with actively expressed exons to determine % assicatied with transcription
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/K27only.1000bp_ann.bed -b $BASEDIR/RNA_seq/shield_expressed_exons.bed -wa > $BASEDIR/cutNrun_mods/peaks/K27only_exonsON.bed
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/K4only.1000bp_ann.bed -b $BASEDIR/RNA_seq/shield_expressed_exons.bed -wa > $BASEDIR/cutNrun_mods/peaks/K4only_exonsON.bed
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/bivalent_domains.1000bp_ann.bed -b $BASEDIR/RNA_seq/shield_expressed_exons.bed -wa > $BASEDIR/cutNrun_mods/peaks/bivalent_exonsON.bed

#pulling K9 peaks that are greater than 1000 bps away to be reannotated with unmasked genome
awk -F'\t' 'sqrt($10*$10) >=1000' $BASEDIR/cutNrun_mods/peaks/wt_K9_masked_ann.txt > $BASEDIR/cutNrun_mods/peaks/wt_K9_moreTHAN1000bp_ann.txt
awk '{print $2 "\t" $3 "\t" $4 }' $BASEDIR/cutNrun_mods/peaks/wt_K9_moreTHAN1000bp_ann.txt > $BASEDIR/cutNrun_mods/peaks/wt_K9_moreTHAN1000bp.bed
annotatePeaks.pl $BASEDIR/cutNrun_mods/peaks/wt_K9_moreTHAN1000bp.bed danRer11 -gtf $BASEDIR/unmasked_ref_ann.gtf > $BASEDIR/cutNrun_mods/peaks/wt_K9_moreTHAN1000bp_unmask_ann.txt

#identifying overlap between K9 CnR and K9 ChIP
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K9_intpeaks.bed -b $BASEDIR/K9_chip/peaks/wt_K9_intpeaks.bed -wa -u > $BASEDIR/cutNrun_mods/peaks/CnR_ChIP_K9_overlap.bed
bedtools intersect -a $BASEDIR/K9_chip/peaks/wt_K9_intpeaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K9_intpeaks.bed -wa -u > $BASEDIR/cutNrun_mods/peaks/ChIP_CnR_K9_overlap.bed
cat $BASEDIR/cutNrun_mods/peaks/CnR_ChIP_K9_overlap.bed $BASEDIR/cutNrun_mods/peaks/ChIP_CnR_K9_overlap.bed | bedtools sort -i stdin | bedtools merge -i stdin > $BASEDIR/cutNrun_mods/peaks/K9_CnR_AND_ChIP_peaks.bed

bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/wt_K9_intpeaks.bed -b $BASEDIR/K9_chip/peaks/wt_K9_intpeaks.bed -v > $BASEDIR/cutNrun_mods/peaks/CnR_NOT_ChIP_K9.bed
bedtools intersect -a $BASEDIR/K9_chip/peaks/wt_K9_intpeaks.bed -b $BASEDIR/cutNrun_mods/peaks/wt_K9_intpeaks.bed -v > $BASEDIR/cutNrun_mods/peaks/ChIP_NOT_CnR_K9.bed

#making cutNrun mods figures
computeMatrix scale-regions -S $BASEDIR/cutNrun_mods/bws/wt_K4_1.norm_sort.bw $BASEDIR/cutNrun_mods/bws/wt_K4_2.norm_sort.bw $BASEDIR/cutNrun_mods/bws/wt_K27_1.norm_sort.bw $BASEDIR/cutNrun_mods/bws/wt_K27_2.norm_sort.bw -R $BASEDIR/genes.bed -p 20 -b 2000 -a 2000 -bs=10 --missingDataAsZero -bl $BASEDIR/K9_chip/peaks/blacklist.bed -out $BASEDIR/cutNrun_mods/K4_K27_genes_reps.gz
plotHeatmap -m $BASEDIR/cutNrun_mods/K4_K27_genes_reps.gz -out $BASEDIR/figs/K4_K27_genes_reps.heatmap.pdf --whatToShow 'heatmap and colorbar' --colorMap Greens Purples --regionsLabel genes --samplesLabel H3K4me3 H3K4me3 H3K27me3 H3K27me3

computeMatrix scale-regions -S $BASEDIR/cutNrun_mods/bws/wt_K4_1.norm_sort.bw $BASEDIR/cutNrun_mods/bws/wt_K27_1.norm_sort.bw -R $BASEDIR/genes.bed -p 20 -b 2000 -a 2000 -bs=10 --missingDataAsZero -out $BASEDIR/cutNrun_mods/K4_K27_genes_rep1.gz
plotProfile -m $BASEDIR/$BASEDIR/cutNrun_mods/K4_K27_genes_rep1.gz -out $BASEDIR/figs/K4_K27_rep1.profile.pdf --perGroup --colors green purple

computeMatrix reference-point -S $BASEDIR/cutNrun_mods/bws/wt_K9_1.norm_sort.bw $BASEDIR/cutNrun_mods/bws/wt_K9_2.norm_sort.bw -R $BASEDIR/cutNrun_mods/peaks/wt_K9_intpeaks.bed --referencePoint center -p 20 -b 5000 -a 5000 -bs=10 --missingDataAsZero -bl $BASEDIR/K9_chip/peaks/blacklist.bed -o $BASEDIR/cutNrun_mods/K9_cNrPeaks_reps.gz
plotHeatmap -m $BASEDIR/cutNrun_mods/K9_cNrPeaks_reps.gz -out $BASEDIR/figs/K9_cNrPeaks_reps_heatmap.pdf --colorMap=Blues --whatToShow 'heatmap and colorbar' --regionsLabel cutNrun_peaks --samplesLabel H3K9me3 H3K9me3

###mappability analysis
#mappability calculated using genmap  - https://github.com/cpockrandt/genmap
#used kmer of 76 and mismatch of 2
#bedgraph output file saved as $BASEDIR/mappability_76bp.bedgraph
#made a bigwig file for data visualization in IGV
bedGraphToBigWig $BASEDIR/mappability_76bp.bedgraph $BASEDIR/chrNameLength.txt $BASEDIR/mappability.bw

#now going to identify regions where mappability score is less than 0.5
awk -F'\t' '$4 <=0.5 {print $1 "\t" $2 "\t" $3 }' $BASEDIR/mappability_76bp.bedgraph > $BASEDIR/mappability_low.bed
bedtools genomecov -bg -i $BASEDIR/mappability_low.bed -g $BASEDIR/chrNameLength.txt > $BASEDIR/mappability_low.bg

#now going to determine peak overlap with these regions of low mappability
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/CnR_NOT_ChIP_K9.bed -b $BASEDIR/mappability_low.bg -wa -u -f 0.3 > $BASEDIR/cutNrun_mods/peaks/CnR_only_LOWmap.bed
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/ChIP_NOT_CnR_K9.bed -b $BASEDIR/mappability_low.bg -wa -u -f 0.3 > $BASEDIR/cutNrun_mods/peaks/ChIP_only_LOWmap.bed
bedtools intersect -a $BASEDIR/cutNrun_mods/peaks/K9_CnR_AND_ChIP_peaks.bed -b $BASEDIR/mappability_low.bg -wa -u -f 0.3 > $BASEDIR/cutNrun_mods/peaks/CnR_ChIP_sharedPeaks_LOWmap.bed


#####these sections from paper revions#####

##making some QC figures to compare replicates
multiBigwigSummary bins -b $BASEDIR/cutNrun_mods/bws/* $BASEDIR/cutNrun_pol2/bws/*pol2.norm* -o $BASEDIR/bw_summ_total.npz -bs 1000 -bl $BASEDIR/K9_chip/peaks/blacklist.bed -p 20
plotCorrelation -in $BASEDIR/bw_summ_total.npz -c spearman -p heatmap --plotNumbers --colorMap Greens -o $BASEDIR/figs/bw_total_summ.heatmap.pdf
plotPCA -in $BASEDIR/bw_summ_total.npz --colors purple orchid green lime navy blue darkorange orange -o $BASEDIR/figs/bw_total_PCA.pdf

###making bws of the mean of the replicates for profile plots
bigwigCompare -b1 $BASEDIR/cutNrun_mods/bws/wt_K27_1.norm_sort.bw -b2 $BASEDIR/cutNrun_mods/bws/wt_K27_2.norm_sort.bw --operation mean -bs 10 -p 20 -o $BASEDIR/cutNrun_mods/bws/wt_K27_meanOFreps.bw
bigwigCompare -b1 $BASEDIR/cutNrun_mods/bws/wt_K4_1.norm_sort.bw -b2 $BASEDIR/cutNrun_mods/bws/wt_K4_2.norm_sort.bw --operation mean -bs 10 -p 20 -o $BASEDIR/cutNrun_mods/bws/wt_K4_meanOFreps.bw
bigwigCompare -b1 $BASEDIR/cutNrun_mods/bws/wt_K9_1.norm_sort.bw -b2 $BASEDIR/cutNrun_mods/bws/wt_K9_2.norm_sort.bw --operation mean -bs 10 -p 20 -o $BASEDIR/cutNrun_mods/bws/wt_K9_meanOFreps.bw
bigwigCompare -b1 $BASEDIR/cutNrun_pol2/bws/wt_1_pol2.norm_sort.bw -b2 $BASEDIR/cutNrun_pol2/bws/wt_3_pol2.norm_sort.bw --operation mean -bs 10 -p 20 -o $BASEDIR/cutNrun_pol2/bws/wt_pol2_meanOFreps.bw

computeMatrix scale-regions -S $BASEDIR/cutNrun_pol2/bws/wt_pol2_meanOFreps.bw  -R $BASEDIR/genes.bed -p 20 -b 2000 -a 2000 -bs=10 --missingDataAsZero -out $BASEDIR/cutNrun_pol2/pol2_genes_meanOFreps.gz
plotProfile -m $BASEDIR/cutNrun_pol2/pol2_genes_meanOFreps.gz -out $BASEDIR/figs/pol2_genes_meanOFreps.profile.pdf --colors orange

computeMatrix scale-regions -S $BASEDIR/cutNrun_mods/bws/wt_K4_meanOFreps.bw $BASEDIR/cutNrun_mods/bws/wt_K27_meanOFreps.bw -R $BASEDIR/genes.bed -p 20 -b 2000 -a 2000 -bs=10 --missingDataAsZero -out $BASEDIR/cutNrun_mods/K4_K27_genes_meanOFreps.gz
plotProfile -m $BASEDIR/cutNrun_mods/K4_K27_genes_meanOFreps.gz -out $BASEDIR/figs/K4_K27_meanOFreps.profile.pdf --perGroup --colors green purple
