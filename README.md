<<<<<<< HEAD
# ewas
Exome-wide association mapping for insecticide resistance associated genes
=======
# ewas_2024
Code for detection of SNP´s and exome modifications in insecticide resistant insects. 
I will put a brief explanation and an example. 
In case of need assistance please write to alejandro.mejiam1@udea.edu.co

###=======================================================================================
#PART 1 (Regular bioinformatic pipeline with some python scripts to streamline the process)
#Quality check and remove any adapaters
Use fastqc for chiecking phred scores

I used cutadapt to remove any potential adapters from illumina sequencing

e.g.
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -a AGATGTGTATAAGAGACAG -A AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT -G CTGTCTCTTATACACATCT -q 30,30 --minimum-length 80 -o 5Feb_Chl_A2_R1_paired.fastq -p 5Feb_Chl_A2_R2_paired.fastq 5Feb_Chlorpyr_A2_S52_L004_R1_001.fastq 5Feb_Chlorpyr_A2_S52_L004_R2_001.fastq

#Map using BWA
Use *bucle_mapp_BWA* its a python script so be aware. The same 

#Samtools sort & index
Use *bucle_samtools.py* for this step. 
Samtools index and then sort

#Convert samtools sorted to pileup to varscan 
Use *bucle_pileup_varscan.py* 

**IMPORTANT** previous to varscan there are some lines that need to be removed. This is already accomplished in *bucle_pileup_varscan.py*. But be aware. 

###=======================================================================================
#PART 2 (Couple of R scripts that... does a lot of things. Read paper in case of need more info Karla et. al. 2021, or 2019)

As S. lozano put it: This pipeline is a series of programs to analyze genome nucleotide variants read counts and to compare the similitude within two replicates. Later, we compare two phenotypically different groups. 

*Required programs.
spliter.r
physmap.r
	physmap.cpp
ven2x2.r
easy_chi2.r
	easy_chi.cpp
	easy_chi_fun.r
annotate.r
replace.r

#The programs are the .r files previously mentioned, however I did some loops for, you guess it, stream line the process. 

The loops are these, in order: 
*splitter_Alejandro_bucle.r
*physmap_Alejandro_bucle.r
*bucle_ven2x2_Alejandro.r
*bucle_easy_chi2_Alejandro.r

The output of all this correspond to .rds files that are easy to handle by R programs. 






