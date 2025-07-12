#Script exoma usando muestras de Aedes aegypti por Alejandro Mejía Muñoz
#Mapeador BWA
#Modificaciones con samtools
#Llamado de variantes GATK 
#Merge de VCFs con bcftools 
#Index con tabix
#GWAS con PLINK

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
    #Bucle está en "/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/bucle_samtools.py"

#Read group #Adiciona una serie de columnas necesarias para los pasos subsecuentes. "Assigns all the reads in a file to a single new read-group"
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar AddOrReplaceReadGroups I=HLR1_sbs_BWA_sortd.bam O=HLR1_BWA_rg.bam SO=coordinate CREATE_INDEX=true RGID=HLR1 RGLB=lib1 RGPL=illumina RGPU=HLR1 RGSM=sample1

#Marcar duplicados usando picard
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar MarkDuplicates I=HLR1_BWA_rg.bam O=HLR1_BWA_-nodups.bam M=HLR1_BWA.metrics REMOVE_DUPLICATES=TRUE

#Llamada de Variantes usando GATK
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java "-Xmx30G" -jar /home/administrador/programas/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar HaplotypeCaller -I HLR1_BWA_-nodups.bam -R /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta -O HLR1.variants.vcf --native-pair-hmm-threads 30

#################Unir VCFs de muestras diferentes diferentes################

# Comprimir los archivos VCF
bgzip *.vcf

#Indexar HLS1, HLS2, SP1, SP2 
bcftools index *vcf.gz

#Bucle para index 
for f in *.vcf.gz; do
        /mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools index "$f"
done

#Merge
/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools merge --force-samples HLS1.variants.vcf.gz HLS2.variants.vcf.gz SP1.variants.vcf.gz SP2.variants.vcf.gz -Oz -o HLS1_HLS2_SP1_SP2.variants.vcf.gz

##Modificar el header
   # Write out the header to be modified
/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools view -h HLS1_HLS2_SP1_SP2.variants.vcf.gz > header.txt

   # Edit the header using your favorite text editor
   nano header.txt

   # Reheader the file
/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools reheader -h header.txt -o reHLS1_HLS2_SP1_SP2.variants.vcf.gz HLS1_HLS2_SP1_SP2.variants.vcf.gz

###El archivo VCF ya está listo###

#########################################################################################################################
##################Plink - (association_gwas_MOD.sh) Suministrado por Cristial Velarde mod Alejandro Mejía################
#########################################################################################################################


#Indexar archivos vcf con tabix
tabix -c vcf "Archivo.vcf.gz" 

#usar plink
/mnt/disc2/grupobcei/ewas/ewas_Acacias/plink_folder/plink --vcf HLS1_HLS2_SP1_SP2_re.variants.vcf.gz --make-bed --out /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/plink --aec
#El --aec omite el formato std de nombres de cromosomas y deja otros formatos disponibles para que corra sin errores

#####################
###script en BASH####
#####################
#Este script cambia el archivo plink para que no haya ambiguedades. 
awk '{
  if ($2 == ".") {
    $2 = $1 ":" $4 ":" NR
  }
  print $0
}' plink.bim > plink_fixed.bim

#Identificación y SNPs duplicados con PLINK
/mnt/disc2/grupobcei/ewas/ewas_Acacias/plink_folder/plink --bfile plink_fixed --list-duplicate-vars ids-only --out duplicados --aec

#A partir del duplicados.dupvar crear un archivo excluyendo los SNPs duplicados
/mnt/disc2/grupobcei/ewas/ewas_Acacias/plink_folder/plink --bfile plink_fixed --exclude duplicados.dupvar --make-bed --out plink_final --aec

#Crear archivo de fenotipos
nano phenotype.txt
"""
FID IID Phenotype
HLS1 HLS1 2
HLS2 HLS2 2
SP1 SP1 1
SP2 SP2 1
"""

#Ejecutar el GWAS
/mnt/disc2/grupobcei/ewas/ewas_Acacias/plink_folder/plink --bfile plink_final --pheno phenotype.txt --assoc --out resultados_gwas --allow-no-sex --aec

#### Filtrar los snps significativos
awk '$9 < 5e-5' resultados_gwas.assoc > snps_significativos_0.00005.txt
awk '$9 < 5e-2' resultados_gwas.assoc > snps_significativos_0.05.txt

####BASH para visualizaciones
outdir="/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/visualizaciones/"
cd ${outdir}
export OUTDIR_ABS=$(pwd)
