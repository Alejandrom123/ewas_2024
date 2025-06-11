IIIII Funcionó
IIIII No funcionó
IIIII Problemas

Ubicación
#Samples in /mnt/disc2/grupobcei/ewas/ in 172.16.0.96 (grupobcei) server

Script Python para mapear automáticamente 
"/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/automatizt/bucle_mapp_bowtie.py"

04-06-2025
Use script Python para mapear:
HLR1 y HPR1
 

Al realizar el script por separado, los porcentajes son los mismos. Procedo a hacer el mapeo con el bucle-script. 

	Mapeo con bucle usando bowtie2
	Realicé el index con BWA
	Realicé el bucle del mapeo con BWA, hay que correrlo mañana

10-06-2025
IA y pipeline Johana Tejada
Muestras HLR1 (HLR2) y HPR1 (HPR2)

Primero exclusivamente con HLR1 bowtie2:
Pasar a bam:
samtools view -S -b HLR1_subsample_bowtie2.sam > HLR1_subsample_bowtie2.bam
samtools sort -o HLR1_sbs_bwt2_sortd.bam HLR1_subsample_bowtie2.bam
samtools index HLR1_sbs_bwt2_sortd.bam

Añadir Read Groups
Primero exclusivamente con HLR1 bwa:
Pasar a bam:
samtools view -S -b HLR1_subsample_BWA.sam> HLR1_subsample_BWA.bam
samtools sort -o HLR1_sbs_BWA_sortd.bam HLR1_subsample_BWA.bam
samtools index HLR1_sbs_BWA_sortd.bam


#Read group para HLR1 con bowtie2
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar AddOrReplaceReadGroups I=HLR1_sbs_bwt2_sortd.bam O=HLR1_bwt2_rg.bam SO=coordinate CREATE_INDEX=true RGID=HLR1 RGLB=lib1 RGPL=illumina RGPU=HLR1 RGSM=sample1

#Read group para HLR1 con BWA
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar AddOrReplaceReadGroups I=HLR1_sbs_BWA_sortd.bam O=HLR1_BWA_rg.bam SO=coordinate CREATE_INDEX=true RGID=HLR1 RGLB=lib1 RGPL=illumina RGPU=HLR1 RGSM=sample1

Marcar duplicados usando picard 
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar MarkDuplicates I=HLR1_BWA_rg.bam O=HLR1_BWA_-nodups.bam M=HLR1_BWA.metrics REMOVE_DUPLICATES=TRUE

/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar MarkDuplicates I=HLR1_bwt2_rg.bam O=HLR1_bowtie2_-nodups.bam M=HLR1_bowtie2.metrics REMOVE_DUPLICATES=TRUE
#Variant calling using GATK 
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java "-Xmx30G" -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar HaplotypeCaller -I HLR1_BWA_-nodups.bam -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -O HLR1.variants.vcf --native-pair-hmm-threads 30

#Gatk filtering
nohup /mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/gatk/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar VariantFiltration -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -V HLR1.variants.vcf --filter-name FAIL --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -O HLR1_filt_10x.vcf

/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools view HLR1_filt_10x.vcf --threads 30 -Oz -O HLR1_filt_10x.vcf.gz

## Generate before/after plots [Only necessary if you want to look at the data; requires R libraries ggplot2, reshape, gplots, gsalib]
nohup /mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/gatk/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar AnalyzeCovariates -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -V HLR1.variants.vcf -before HLR1before_recal_data.table -after HLR1_postRecal_data.table -plots HLR1_recalibration_plots.pdf

# SnpEff for SNP annotation
Genoma de referencia en snpeff
Aedes_aegypti_lvpagwg

java -jar /data1/softwares/snpEff/snpEff.jar Aedes_aegypti_lvpagwg -c /data1/softwares/snpEff/snpEff.config tu_archivo.vcf > tu_archivo_annotated.vcf

/mnt/disc2/grupobcei/java/jdk-24.0.1/bin/java -jar /mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.jar Aedes_aegypti_lvpagwg -c /mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.config HLR1_filt_10x.vcf > HLR1_GATK_annotated.vcf



########################################################################
## Index the BAM file using Picard Tools
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar BuildBamIndex INPUT= HLR1_BWA_-nodups.bam

/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar BuildBamIndex INPUT= HLR1_bowtie2_-nodups.bam 

Realineación local alrededor de indels (opcional pero recomendado)
gatk --java-options "-Xmx4g" RealignerTargetCreator -R reference.fasta -I HLR1_BWA_-nodups.bam -O realignment_targets.list

/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -Xmx30g -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar RealignerTargetCreator -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -I HLR1_BWA_-nodups.bam -O HLR1_BWA_realignment_targets.list

/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -Xmx30g -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar BaseRecalibrator -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -I HLR1_BWA_-nodups.bam -O HLR1_BWA_recal_data.table

