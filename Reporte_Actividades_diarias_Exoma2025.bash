#########################################################################
######################### Indicaciones iniciales##############################
#########################################################################

IIIII Funcionó
IIIII No funcionó
IIIII Problemas

https://bio.tools/VCF-Server
http://diseasegps.sjtu.edu.cn/VCF-Server?lan=eng

Elementos para abrir para trabajar en exoma: 
Reporte_Actividades_diarias_Exoma2025 (word) y el .sh con este mismo nombre lo completo al final del día para subir los cambios hechos en GIT

El pipeline de exoma en el que estoy trabajando 
C:\Users\User\Desktop\exoma\Script_exoma.bash

A fecha de 17-06-2025 lo que estoy trabajando está en el servidor de biología
/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/

Ubicación de lo de EWAS en general
#Samples in /mnt/disc2/grupobcei/ewas/ in 172.16.0.96 (grupobcei) server


#########################################################################
################################################################
#########################################################################



Script Python para mapear automáticamente 
"/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/automatizt/bucle_mapp_bowtie.py"

04-06-2025
Use script Python para mapear:
HLR1 y HPR1
 

Al realizar el script por separado, los porcentajes son los mismos. Procedo a hacer el mapeo con el bucle-script. 

	Mapeo con bucle usando bowtie2
	Realicé el index con BWA
	Realicé el bucle del mapeo con BWA, hay que correrlo mañana

#########################################################################
###############################10-06-2025#################################
#########################################################################

IA y pipeline Johana Tejada
Muestras HLR1 (HLR2) y HPR1 (HPR2)

Primero exclusivamente con HLR1 bowtie2:
Pasar a bam:
samtools view -S -b HLR1_subsample_bowtie2.sam > HLR1_subsample_bowtie2.bam
samtools sort -o HLR1_sbs_bwt2_sortd.bam HLR1_subsample_bowtie2.bam
samtools index HLR1_sbs_bwt2_sortd.bam

Añadir Read Groups
Mapeo con BWA
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

# SnpEff for SNP annotation

Genoma de referencia en snpeff
Aedes_aegypti_lvpagwg

java -jar /data1/softwares/snpEff/snpEff.jar Aedes_aegypti_lvpagwg -c /data1/softwares/snpEff/snpEff.config tu_archivo.vcf > tu_archivo_annotated.vcf

/mnt/disc2/grupobcei/java/jdk-24.0.1/bin/java -jar /mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.jar Aedes_aegypti_lvpagwg -c /mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.config HLR1_filt_10x.vcf > HLR1_GATK_annotated.vcf 


#########################################################################
###############################12-06-2025#################################
#########################################################################

Realicé ordenamiento de script en Script_exoma
Y realicé samtools bucle para view, sort y luego index. 
Se deja corriendo

#########################################################################
###############################17-06-2025#################################
#########################################################################


Comando luego de samtools es read groups: 
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar AddOrReplaceReadGroups I=HLR1_sbs_BWA_sortd.bam O=HLR1_BWA_rg.bam SO=coordinate CREATE_INDEX=true RGID=HLR1 RGLB=lib1 RGPL=illumina RGPU=HLR1 RGSM=sample1


#########################################################################
###############################18-06-2025#################################
#########################################################################


Creé el archivo de métricas usando samtools stat:
nohup bash -c 'for file in *.sam; do samtools stats "$file" > "${file%.sam}.stats.txt"; done' > nohup.out 2>&1 &  

y luego para unir:
cat *stats.txt > merged_BWA_stats.txt

samtools stats HLR1_subsample_BWA.sam > HLR1_bwa_stats.txt

En archivo métricas C:\Users\User\Desktop\exoma\outputs_mobaxterm\Metricas_Mapeo_Exoma.xlsx
Están las métricas del mapeo

La muestra HLR1 parece que está mala.  

#########################################################################
###############################19-06-2025#################################
#########################################################################

Hay un problema con HLR1. Cuando mapeo con BWA, dice que el header está malo. 
Estoy mirando si se comparte con HLS1. 

Pasar a BAM a ver sí samtools lo reconoce más fácilmente. 

--- 
Estos son los comandos: 
samtools view -bS HLR1_BWA_19-06-25.sam > HLR1_BWA_19-06-25.bam 
samtools flagstat HLR1_BWA_19-06-25.bam 
samtools stats HLR1_BWA_19-06-25.bam > HLR1_BWA_19-06-25.stats.txt

Volví a cargar la muestra HLR1 desde mi tera y funcionó. HLR1 percentage of properly paired reads (%):	82.1

/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/





#########################################################################
###############################Restantes#################################
#########################################################################

## Index the BAM file using Picard Tools
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar BuildBamIndex INPUT= HLR1_BWA_-nodups.bam

/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar BuildBamIndex INPUT= HLR1_bowtie2_-nodups.bam 

Realineación local alrededor de indels (opcional pero recomendado)
gatk --java-options "-Xmx4g" RealignerTargetCreator -R reference.fasta -I HLR1_BWA_-nodups.bam -O realignment_targets.list

/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -Xmx30g -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar RealignerTargetCreator -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -I HLR1_BWA_-nodups.bam -O HLR1_BWA_realignment_targets.list

/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -Xmx30g -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar BaseRecalibrator -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -I HLR1_BWA_-nodups.bam -O HLR1_BWA_recal_data.table

