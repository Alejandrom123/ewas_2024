#!/bin/bash

#----------------------#
# INDEXACIÓN DEL GENOMA
#----------------------#
bwa index VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta

#----------------------#
# MAPEAMIENTO POR MUESTRA (ejemplo para HLR1, repetir para cada muestra)
#----------------------#
bwa mem -t 30 -M VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta HLS1_R1_subsample.fastq HLS1_R2_subsample.fastq > HLS1_subsample_BWA.sam

samtools view -S -b HLS1_subsample_BWA.sam > HLS1_subsample_BWA.bam
samtools sort -o HLS1_sbs_BWA_sortd.bam HLS1_subsample_BWA.bam
samtools index HLS1_sbs_BWA_sortd.bam

#----------------------#
# READ GROUPS Y DUPLICADOS (ejemplo para HLR1, repetir para cada muestra)
#----------------------#
java -jar /mnt/disc2/grupobcei/picard/picard.jar AddOrReplaceReadGroups I=HLS1_sbs_BWA_sortd.bam O=HLS1_BWA_rg.bam SO=coordinate CREATE_INDEX=true RGID=HLR1 RGLB=lib1 RGPL=illumina RGPU=HLR1 RGSM=HLR1

java -jar /mnt/disc2/grupobcei/picard/picard.jar MarkDuplicates I=HLS1_BWA_rg.bam O=HLS1_BWA_-nodups.bam M=HLS1_BWA.metrics REMOVE_DUPLICATES=TRUE

#----------------------#
# LLAMADO DE VARIANTES MODO GVCF (para cada muestra)
#----------------------#
for sample in HLS1 HLS2 SP1 SP2; do
  /mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java "-Xmx30G" -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar HaplotypeCaller \
    --native-pair-hmm-threads 30 \
    -I ${sample}_BWA_nodups.bam \
    -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta \
    -O ${sample}.g.vcf.gz \
    -ERC GVCF
done

#----------------------#
# COMBINAR GVCFs (si son pocas muestras)
#----------------------#
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java "-Xmx30G" -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar CombineGVCFs \
  -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta \
  --variant HLS1.g.vcf.gz \
  --variant HLS2.g.vcf.gz \
  --variant SP1.g.vcf.gz \
  --variant SP2.g.vcf.gz \
  -O cohort.g.vcf.gz

nohup /mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java "-Xmx30G" -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar CombineGVCFs -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta --variant HLS1.g.vcf.gz --variant HLS2.g.vcf.gz --variant SP1.g.vcf.gz --variant SP2.g.vcf.gz -O cohort.g.vcf.gz

#----------------------#
# JOINT GENOTYPING
#----------------------#
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java "-Xmx30G" -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar GenotypeGVCFs \
  --native-pair-hmm-threads 30 \  
  -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta \
  -V cohort.g.vcf.gz \
  -O cohort_joint.vcf.gz
  
nohup /mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java "-Xmx30G" -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar GenotypeGVCFs --native-pair-hmm-threads 30 -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -V cohort.g.vcf.gz -O cohort_joint.vcf.gz

#----------------------#
# INDEXAR VCF FINAL
#----------------------#
bgzip -f cohort_joint.vcf #Si no está comprimido
tabix -p vcf cohort_joint.vcf.gz

#----------------------#
# CONVERTIR A PLINK Y GWAS
#----------------------#
/mnt/disc2/grupobcei/ewas/ewas_Acacias/plink_folder/plink --vcf cohort_joint.vcf.gz --make-bed --out /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/plink --aec

awk '{
if ($2 == ".") {
  $2 = $1 ":" $4 ":" NR
}
print $0
}' /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/plink.bim > /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/plink_fixed.bim

# Identificar SNPs duplicados
/mnt/disc2/grupobcei/ewas/ewas_Acacias/plink_folder/plink --bfile /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/plink_fixed --list-duplicate-vars ids-only --out /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/duplicados --aec

# Excluir SNPs duplicados
/mnt/disc2/grupobcei/ewas/ewas_Acacias/plink_folder/plink --bfile /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/plink_fixed --exclude /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/duplicados.dupvar --make-bed --out /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/plink_final --aec

# Crear archivo de fenotipos (ejemplo)
cat > /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/phenotype.txt << EOF
FID IID Phenotype
HLS1 HLS1 2
HLS2 HLS2 2
SP1 SP1 1
SP2 SP2 1
EOF

# Ejecutar GWAS
/mnt/disc2/grupobcei/ewas/ewas_Acacias/plink_folder/plink --bfile /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/plink_final --pheno /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/phenotype.txt --assoc --out /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/resultados_gwas --allow-no-sex --aec

# Filtrar SNPs significativos
awk '$9 < 5e-5' /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/resultados_gwas.assoc > /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/snps_significativos_0.00005.txt
awk '$9 < 5e-2' /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/resultados_gwas.assoc > /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/snps_significativos_0.05.txt

#----------------------#
# VISUALIZACIONES (ajusta según tu flujo)
#----------------------#
