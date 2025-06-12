#Script exoma usando muestras de Aedes aegypti por Alejandro Mejía Muñoz
#Mapeador BWA
#Llamado de variantes GATK 
#Anotación de variantes SnpEff

#Index usando BWA
bwa index VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta

#Mapeo con BWA
bwa mem -t 30 -M /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta _R1_subsample.fastq _R2_subsample.fastq > _subsample_BWA.sam
    #Mapeo Automático BWA
    #EN PYTHON -> bucle_mapp_BWA.py 

#Pasar de SAM a BAM
samtools view -S -b HLR1_subsample_BWA.sam > HLR1_subsample_BWA.bam
samtools sort -o HLR1_sbs_BWA_sortd.bam HLR1_subsample_BWA.bam
samtools index HLR1_sortd_index.bam
    #samtools automático 
    #

#Read group #Adiciona una serie de columnas necesarias para los pasos subsecuentes. "Assigns all the reads in a file to a single new read-group"
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar AddOrReplaceReadGroups I=HLR1_sbs_BWA_sortd.bam O=HLR1_BWA_rg.bam SO=coordinate CREATE_INDEX=true RGID=HLR1 RGLB=lib1 RGPL=illumina RGPU=HLR1 RGSM=sample1

#Marcar duplicados usando picard
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar MarkDuplicates I=HLR1_BWA_rg.bam O=HLR1_BWA_-nodups.bam M=HLR1_BWA.metrics REMOVE_DUPLICATES=TRUE

#Llamada de Variantes usando GATK
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java "-Xmx30G" -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar HaplotypeCaller -I HLR1_BWA_-nodups.bam -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -O HLR1.variants.vcf --native-pair-hmm-threads 30

#Filtrado usando GATK
nohup /mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/gatk/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar VariantFiltration -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -V HLR1.variants.vcf --filter-name FAIL --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -O HLR1_filt_10x.vcf

#Comprimir las variantes en el caso de que sea necesario para aplicaciones upstream
/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools view HLR1_filt_10x.vcf --threads 30 -Oz -O HLR1_filt_10x.vcf.gz

#SnpEff para la anotación de SNPs
/mnt/disc2/grupobcei/java/jdk-24.0.1/bin/java -jar /mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.jar Aedes_aegypti_lvpagwg -c /mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.config HLR1_filt_10x.vcf > HLR1_GATK_annotated.vcf 

