#Pipeline_ewas_linux

#Mpileup
/home/grupobcei/grupobcei/ewas/ewas_Acacias/03._ewas_mapped_subsample/samtools/samtools-1.21/samtools mpileup -f VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta /mnt/disc2/grupobcei/ewas/index/bowtie2_index/HLR1_subsample_bowtie2_sorted.bam > HLR1.pileup
#Using Java version 17 y 24
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -version #Esto funcionó#
/mnt/disc2/grupobcei/java/jdk-24.0.1/bin/java -version
#using gatk
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/gatk/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar
#Using picard
/mnt/disc2/grupobcei/picard/picard.jar
#using samtools que descargué
/mnt/disc2/grupobcei/ewas/ewas_Acacias/03._ewas_mapped_subsample/samtools/samtools-1.21/samtools 
#using bcftools
/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools
#SnpEff
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.jar
#Varscan
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/ewas/ewas_Acacias/Varscan/VarScan.v2.4.6.jar
#PoPoolation2
perl /mnt/disc2/grupobcei/ewas/ewas_Acacias/popoolation2-master/mpileup2sync.pl
#Usando picard para ver integridad de los archivos .sam
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar ValidateSamFile -I LP1_subsample.sam -MODE SUMMARY
/mnt/disc2/grupobcei/ewas/ewas_Acacias/03._ewas_mapped_subsample/samtools/samtools-1.21/samtools view -b -o LP1_subsample.bam LP1_subsample.sam
#tabix, está en el PATH
tabix -p vcf "archivo.vcf.gz"
#Con el reference fasta 
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar ValidateSamFile -I LP1_subsample.sam -MODE SUMMARY -REFERENCE_SEQUENCE /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta
#Plink
/mnt/disc2/grupobcei/ewas/ewas_Acacias/plink_folder/plink

#Tengo que crear estos indices con samtools y picard para poder usar ValidateSamFile y el reference genome
#crear .fai 
/mnt/disc2/grupobcei/ewas/ewas_Acacias/03._ewas_mapped_subsample/samtools/samtools-1.21/samtools faidx /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta

#crear .dict
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar CreateSequenceDictionary -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -O /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.dict


###############################################################################################################11-02-2025######################################################################################################################################
#bbduk.sh de bbmap en el servidor /mnt/disc2/grupobcei/ewas/ewas_Acacias/BBmap_suite/bbmap/
/mnt/disc2/grupobcei/ewas/ewas_Acacias/BBmap_suite/bbmap/bbduk.sh in=/mnt/disc2/grupobcei/ewas/ewas_Acacias/00_fastq/HLR1_R1_001.fastq.gz in2=/mnt/disc2/grupobcei/ewas/ewas_Acacias/00_fastq/HLR1_R2_001.fastq.gz out=HLR1_bbcleaned_1.fastq out2=HLR1_bbcleaned_2.fastq ref=/mnt/disc2/grupobcei/ewas/ewas_Acacias/BBmap_suite/bbmap/resources/truseq.fa.gz tpe tbo

/mnt/disc2/grupobcei/ewas/ewas_Acacias/BBmap_suite/bbmap/bbduk.sh in=/mnt/disc2/grupobcei/ewas/ewas_Acacias/00_fastq/HPR1_R1_001.fastq.gz in2=/mnt/disc2/grupobcei/ewas/ewas_Acacias/00_fastq/HPR1_R2_001.fastq.gz out=HPR1_bbcleaned_1.fastq out2=HPR1_bbcleaned_2.fastq ref=/mnt/disc2/grupobcei/ewas/ewas_Acacias/BBmap_suite/bbmap/resources/truseq.fa.gz tpe tbo
