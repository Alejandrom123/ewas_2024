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

## Index the BAM file using Picard Tools
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar BuildBamIndex INPUT= HLR1_BWA_-nodups.bam

/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar BuildBamIndex INPUT= HLR1_bowtie2_-nodups.bam

Realineación local alrededor de indels (opcional pero recomendado)
gatk --java-options "-Xmx4g" RealignerTargetCreator -R reference.fasta -I HLR1_BWA_-nodups.bam -O realignment_targets.list

/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -Xmx30g -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar RealignerTargetCreator -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -I HLR1_BWA_-nodups.bam -O HLR1_BWA_realignment_targets.list

/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -Xmx30g -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar BaseRecalibrator -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -I HLR1_BWA_-nodups.bam -O HLR1_BWA_recal_data.table

Genome resequencing and genome-wide polymorphisms in mosquito vectors Aedes aegypti and Aedes albopictus from south India 







gatk --java-options "-Xmx30g" IndelRealigner -R reference.fasta -I HLR1_BWA_-nodups.bam -targetIntervals realignment_targets.list -O HLR1_BWA_realigned.bam

