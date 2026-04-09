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
    #Bucle está en "/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/bucle_samtools.py"

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

    #El bucle para lo anterior se llama bucle.variantes.py -> Funcionó excepcionalmente desde que lo corrí el 19-06-2025
    
#mpileup (Analiza los archivos para mirar cobertura y las bases encontradas comparado con una referencia) ## AL parecer ya no se hace smatools mpileup sino bcftools mpileup
samtools mpileup -f /mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta _BWA_sortd.bam > .pileup
samtools mpileup -B XX_f2.bam YY_f2.bam ZZ_f2.bam -o all_pools_ordered.mpileup # el -o era el output en versiones anteriores

#PoPoolation2 a partir de mpileup
perl /mnt/disc2/grupobcei/ewas/ewas_Acacias/popoolation2-master/mpileup2sync.pl --input muestra.mpileup --output muestra.sync

#Convert mpileup file to readcounts using Varscan
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/ewas/ewas_Acacias/Varscan/VarScan.v2.4.6.jar readcounts .pileup --min-coverage 25 --min-base-qual 30 --output-file .readcounts

#########################################################################################################################
##################Plink - (association_gwas_MOD.sh) Suministrado por Cristial Velarde mod Alejandro Mejía################
#########################################################################################################################

#Modificado para uso manual en servidor biología UDEA

#Indexar archivos vcf con tabix
tabix -p vcf "Archivo.vcf.gz" 

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

####R para visualizaciones

echo -e "\e[32m################################################################################################\e[0m"
echo -e "\e[32m#### Generando gráficos de Manhattan y QQ ===> \e[31m$(date)\e[0m \e[32m                ####\e[0m"
echo -e "\e[32m################################################################################################\e[0m"
mkdir -p ${outdir}/visualizaciones
cd ${outdir}

export OUTDIR_ABS=$(pwd)

R --vanilla <<EOF
if (!requireNamespace("qqman", quietly = TRUE)) {
    install.packages("qqman", repos = "https://cloud.r-project.org")
}

library(qqman)

assoc_path <- file.path(Sys.getenv("OUTDIR_ABS"), "resultados_gwas", "resultados_gwas.qassoc")
out_dir <- file.path(Sys.getenv("OUTDIR_ABS"), "visualizaciones")

assoc <- read.table(assoc_path, header = TRUE)
assoc <- assoc[!is.na(assoc\$P) & assoc\$P > 0, ]

png(filename = file.path(out_dir, "manhattan_plot.png"), width = 1000, height = 500)
manhattan(assoc,
          chr = "CHR",
          bp = "BP",
          snp = "SNP",
          p = "P",
          main = "Manhattan Plot",
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-5),
          col = c("blue4", "red3"))
dev.off()

png(filename = file.path(out_dir, "qq_plot.png"), width = 500, height = 500)
qq(assoc\$P, main = "QQ Plot")
dev.off()
EOF

echo ""
echo -e "\e[32m################################################################################################\e[0m"
echo -e "\e[32m#### Exportando resultados a MySQL ===> \e[31m$(date)\e[0m \e[32m                       ####\e[0m"
echo -e "\e[32m################################################################################################\e[0m"

cut -f1-2,5,6 ${outdir}/plink/chr1_plink_final.fam > ${outdir}/plink/individuos.tsv
cut -f2,1,4,5,6 ${outdir}/plink/chr1_plink_final.bim > ${outdir}/plink/variantes.tsv
awk 'NR>1 {OFS="\t"; print $3, $9, $6, $7, $8, $10}' ${outdir}/resultados_gwas/resultados_gwas.qassoc > ${outdir}/plink/asociaciones.tsv

if ! command -v mysql &> /dev/null; then
    echo -e "\e[33mMySQL no está instalado. Instalando...\e[0m"
    sudo apt update && sudo apt install -y mysql-server
else
    echo -e "\e[32mMySQL ya está instalado.\e[0m"
fi

# Iniciar el servicio MySQL
sudo systemctl start mysql

# Ejecutar comandos SQL como root desde script (NO modo interactivo)
sudo mysql --local-infile=1 -u root <<EOF
DROP DATABASE IF EXISTS gwasdb;
CREATE DATABASE gwasdb;
USE gwasdb;

CREATE TABLE individuos (
    fid VARCHAR(50),
    iid VARCHAR(50),
    sexo INT,
    fenotipo FLOAT
);

CREATE TABLE variantes (
    snp_id VARCHAR(100),
    chr VARCHAR(10),
    pos INT,
    alelo1 VARCHAR(10),
    alelo2 VARCHAR(10)
);

CREATE TABLE asociaciones (
    snp_id VARCHAR(100),
    pval FLOAT,
    beta FLOAT,
    se FLOAT,
    tstat FLOAT,
    r2 FLOAT
);

LOAD DATA LOCAL INFILE '${outdir}/plink/individuos.tsv'
INTO TABLE individuos
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';

LOAD DATA LOCAL INFILE '${outdir}/plink/variantes.tsv'
INTO TABLE variantes
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';

LOAD DATA LOCAL INFILE '${outdir}/plink/asociaciones.tsv'
INTO TABLE asociaciones
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';

SELECT v.snp_id, v.chr, v.pos, v.alelo1, v.alelo2,
       a.pval, a.beta, a.se, a.tstat, a.r2
INTO OUTFILE '/var/lib/mysql-files/resultados_concatenados.tsv'
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
FROM variantes v
JOIN asociaciones a ON v.snp_id = a.snp_id;
EOF

# Mover archivo desde carpeta segura al destino final
sudo cp /var/lib/mysql-files/resultados_concatenados.tsv ${outdir}/plink/
sudo chown $(whoami) ${outdir}/plink/resultados_concatenados.tsv

echo -e "\e[32mArchivo final creado en:\e[0m ${outdir}/plink/resultados_concatenados.tsv"
