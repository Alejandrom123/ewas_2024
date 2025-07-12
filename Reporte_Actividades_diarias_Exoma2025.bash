Indicaciones Iniciales
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

Más cosas iniciales
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

10-06-2025
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

12-06-2025
#########################################################################
###############################12-06-2025#################################
#########################################################################

Realicé ordenamiento de script en Script_exoma
Y realicé samtools bucle para view, sort y luego index. 
Se deja corriendo

17-06-2025
#########################################################################
###############################17-06-2025#################################
#########################################################################


Comando luego de samtools es read groups: 
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar AddOrReplaceReadGroups I=HLR1_sbs_BWA_sortd.bam O=HLR1_BWA_rg.bam SO=coordinate CREATE_INDEX=true RGID=HLR1 RGLB=lib1 RGPL=illumina RGPU=HLR1 RGSM=sample1

18-06-2025
#########################################################################
###############################18-06-2025#################################
#########################################################################


Creé el archivo de métricas usando samtools stat:
nohup bash -c 'for file in *.sam; do samtools stats "$file" > "${file%.sam}.stats.txt"; done' > nohup.out 2>&1 &  

nohup bash -c 'for file in *.sam; do samtools flagstats "$file" > "${file%.sam}.stats.txt"; done' > nohup.out 2>&1 & 


y luego para unir:
cat *stats.txt > merged_BWA_stats.txt

samtools stats HLR1_subsample_BWA.sam > HLR1_bwa_stats.txt

En archivo métricas C:\Users\User\Desktop\exoma\outputs_mobaxterm\Metricas_Mapeo_Exoma.xlsx
Están las métricas del mapeo

La muestra HLR1 parece que está mala.  

19-06-2025 
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


Aquí va el bucle variantes. bucle_variantes.py

24-06-2025 Aquí está el script largo con el bucle_variantes.py
#########################################################################
###############################24-06-2025#################################
#########################################################################
Input a perplexity
1.	“Imagen de el output vcf de las columnas que está en HLR1.variants_excel.xlsx” 
Ubicación: "C:\Users\User\Desktop\exoma\outputs_mobaxterm\HLR1_y varios para entrenarme\HLR1.variants_excel.xlsx"
2.	El bucle de bucle_variantes.py
"C:\Users\User\Desktop\exoma\bucles_exoma\bucle_variantes.py"
Pregunta a perplexity
En un artículo encontré justo lo que quiero hacer: Our analysis compared the frequency of the alternate allele at each polymorphic site between the two different phenotypes. The alternate allele consisted of an allele not present in the reference AaegL5 genome assembly. In this study, we performed three pairwise comparisons: 1) kdr vs recovered, 2) recovered vs dead, and 3) kdr vs dead. "A series of R scripts were used to split the genomic sites in the three chromosomes. A list of common sites within biological replicates and common sites between the phenotypes was then generated. The following script selected polymorphic sites (SNPs) and generated tables for allele counts in each of the four libraries (e.g., kdr1, kdr2, dead1, dead2). Allele frequencies for the alternant allele were calculated for each phenotype, and a goodness of fit test identified SNPs with consistent proportions within replicates (p > 0.05). Then, we built contingency tables and calculated the heterogeneity χ2 with n—1 degrees of freedom to compare the proportion of the alternate allele between the phenotypes (the probability derived from this analysis was -log10 transformed to provide a “LOD” value). Additionally, we calculated the expected heterozygosity (Hexp) of each site where and n is the number of alternate nucleotides at a site. We applied a Benjamini-Hochberg correction for false discovery rate [18] for each chromosome separately (α = 0.01)." Basado en lo que hice en "bucle_variantes.py", cómo puedo hacer para obtener algo similar a lo planteado en el párrafo anterior. Puedo seguir desde el archivo VCF o debo de devolverme?



Perfecto, ese artículo describe exactamente el análisis que necesitas. Basándote en tu script bucle_variantes.py, puedes continuar directamente desde los archivos VCF que ya generaste en la carpeta 05_SnpEffAnnotation. No necesitas devolverté.
Aquí te proporciono un pipeline completo en Python y R para replicar el análisis del artículo:
________________________________________
Pipeline para análisis de frecuencias alélicas HPR vs SP
Paso 1: Script Python para procesar VCFs y extraer información alélica
python
import pandas as pd
import numpy as np
from scipy import stats
import os
import re

def parse_vcf_allele_info(vcf_file, sample_name):
    """Extrae información alélica de archivos VCF individuales"""
    variants = []
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
                
            chrom, pos, id_var, ref, alt, qual, filter_field, info = parts[:8]
            
            # Extraer información del campo INFO
            ac = re.search(r'AC=(\d+)', info)
            af = re.search(r'AF=([\d.]+)', info)
            an = re.search(r'AN=(\d+)', info)
            dp = re.search(r'DP=(\d+)', info)
            
            if ac and af and an and dp:
                variants.append({
                    'CHROM': chrom,
                    'POS': int(pos),
                    'REF': ref,
                    'ALT': alt,
                    'QUAL': float(qual),
                    'AC': int(ac.group(1)),  # Alternate allele count
                    'AF': float(af.group(1)),  # Alternate allele frequency
                    'AN': int(an.group(1)),  # Total allele number
                    'DP': int(dp.group(1)),  # Depth
                    'SAMPLE': sample_name
                })
    
    return pd.DataFrame(variants)

def process_group_vcfs(group_samples, group_name):
    """Procesa todos los VCFs de un grupo (ej: HPR1, HPR2)"""
    all_variants = []
    
    for sample in group_samples:
        vcf_file = f"05_SnpEffAnnotation/{sample}_GATK_annotated.vcf"
        if os.path.exists(vcf_file):
            print(f"Procesando {vcf_file}")
            sample_variants = parse_vcf_allele_info(vcf_file, sample)
            all_variants.append(sample_variants)
        else:
            print(f"Archivo no encontrado: {vcf_file}")
    
    if all_variants:
        combined = pd.concat(all_variants, ignore_index=True)
        combined['GROUP'] = group_name
        return combined
    return pd.DataFrame()

def main():
    # Definir grupos basándose en tu tabla RG_info
    groups = {
        'HPR': ['HPR1', 'HPR2'],  # Resistentes a permetrina
        'SP': ['SP1', 'SP2'],     # Sin presión
        'HPS': ['HPS1', 'HPS2'],  # Susceptibles a permetrina
        'HLR': ['HLR1', 'HLR2'],  # Resistentes a lambdacialotrina
        'HLS': ['HLS1', 'HLS2']   # Susceptibles a lambdacialotrina
    }
    
    # Procesar cada grupo
    group_data = {}
    for group_name, samples in groups.items():
        print(f"Procesando grupo {group_name}")
        group_data[group_name] = process_group_vcfs(samples, group_name)
    
    # Combinar datos para análisis
    all_data = pd.concat(group_data.values(), ignore_index=True)
    all_data.to_csv('all_variants_by_group.csv', index=False)
    print("Datos guardados en all_variants_by_group.csv")
    
    return group_data

if __name__ == "__main__":
    group_data = main()
________________________________________
Paso 2: Script R para análisis estadístico siguiendo el artículo
r
library(dplyr)
library(ggplot2)

# Cargar datos procesados
all_variants <- read.csv("all_variants_by_group.csv")

# Función para comparación por pares siguiendo el artículo
pairwise_allele_comparison <- function(data, group1, group2) {
  
  # Filtrar datos para los dos grupos
  group1_data <- data[data$GROUP == group1, ]
  group2_data <- data[data$GROUP == group2, ]
  
  # Encontrar sitios comunes entre grupos y réplicas
  common_sites <- intersect(
    paste(group1_data$CHROM, group1_data$POS, sep="_"),
    paste(group2_data$CHROM, group2_data$POS, sep="_")
  )
  
  print(paste("Sitios comunes entre", group1, "y", group2, ":", length(common_sites)))
  
  results <- data.frame()
  
  for(site in common_sites) {
    site_parts <- strsplit(site, "_")[[1]]
    chrom <- site_parts[1]
    pos <- as.numeric(site_parts[2])
    
    # Datos del sitio para cada grupo
    g1_site <- group1_data[group1_data$CHROM == chrom & group1_data$POS == pos, ]
    g2_site <- group2_data[group2_data$CHROM == chrom & group2_data$POS == pos, ]
    
    if(nrow(g1_site) > 0 & nrow(g2_site) > 0) {
      
      # Calcular frecuencias alélicas por grupo
      g1_ac_total <- sum(g1_site$AC)
      g1_an_total <- sum(g1_site$AN)
      g1_af <- ifelse(g1_an_total > 0, g1_ac_total / g1_an_total, 0)
      
      g2_ac_total <- sum(g2_site$AC)
      g2_an_total <- sum(g2_site$AN)
      g2_af <- ifelse(g2_an_total > 0, g2_ac_total / g2_an_total, 0)
      
      # Test de bondad de ajuste dentro de réplicas (consistencia)
      if(nrow(g1_site) > 1 & nrow(g2_site) > 1) {
        # Chi-cuadrado para consistencia dentro de grupos
        g1_consistency <- tryCatch({
          chisq.test(cbind(g1_site$AC, g1_site$AN - g1_site$AC))$p.value
        }, error = function(e) NA)
        
        g2_consistency <- tryCatch({
          chisq.test(cbind(g2_site$AC, g2_site$AN - g2_site$AC))$p.value
        }, error = function(e) NA)
      } else {
        g1_consistency <- NA
        g2_consistency <- NA
      }
      
      # Tabla de contingencia para comparar entre grupos
      contingency_table <- matrix(c(
        g1_ac_total, g1_an_total - g1_ac_total,
        g2_ac_total, g2_an_total - g2_ac_total
      ), nrow = 2, byrow = TRUE)
      
      # Chi-cuadrado de heterogeneidad
      chi_test <- tryCatch({
        chisq.test(contingency_table)
      }, error = function(e) list(statistic = NA, p.value = NA))
      
      # Calcular LOD score (-log10(p-value))
      lod_score <- ifelse(is.na(chi_test$p.value), NA, -log10(chi_test$p.value))
      
      # Calcular heterocigosidad esperada (Hexp)
      # Hexp = 1 - sum(pi^2) donde pi es la frecuencia de cada alelo
      p_ref <- 1 - ((g1_af + g2_af) / 2)  # Frecuencia del alelo referencia
      p_alt <- (g1_af + g2_af) / 2        # Frecuencia del alelo alternativo
      hexp <- 1 - (p_ref^2 + p_alt^2)
      
      # Guardar resultados
      results <- rbind(results, data.frame(
        CHROM = chrom,
        POS = pos,
        REF = g1_site$REF[1],
        ALT = g1_site$ALT[1],
        G1_AF = g1_af,
        G2_AF = g2_af,
        G1_consistency_p = g1_consistency,
        G2_consistency_p = g2_consistency,
        Chi_squared = chi_test$statistic,
        P_value = chi_test$p.value,
        LOD_score = lod_score,
        Hexp = hexp,
        Comparison = paste(group1, "vs", group2)
      ))
    }
  }
  
  return(results)
}

# Realizar comparaciones por pares
comparisons <- list(
  c("HPR", "SP"),   # Resistentes permetrina vs Sin presión
  c("HPR", "HPS"),  # Resistentes vs Susceptibles permetrina
  c("HLR", "HLS"),  # Resistentes vs Susceptibles lambdacialotrina
  c("HLR", "SP")    # Resistentes lambdacialotrina vs Sin presión
)

all_results <- data.frame()

for(comp in comparisons) {
  print(paste("Analizando comparación:", comp[1], "vs", comp[2]))
  comp_results <- pairwise_allele_comparison(all_variants, comp[1], comp[2])
  all_results <- rbind(all_results, comp_results)
}

# Filtrar SNPs con consistencia dentro de réplicas (p > 0.05)
consistent_snps <- all_results[
  (is.na(all_results$G1_consistency_p) | all_results$G1_consistency_p > 0.05) &
  (is.na(all_results$G2_consistency_p) | all_results$G2_consistency_p > 0.05), 
]

print(paste("SNPs consistentes dentro de réplicas:", nrow(consistent_snps)))

# Aplicar corrección Benjamini-Hochberg por cromosoma (α = 0.01)
consistent_snps$P_adjusted <- NA

for(chrom in unique(consistent_snps$CHROM)) {
  chrom_data <- consistent_snps$CHROM == chrom & !is.na(consistent_snps$P_value)
  if(sum(chrom_data) > 0) {
    consistent_snps$P_adjusted[chrom_data] <- p.adjust(
      consistent_snps$P_value[chrom_data], 
      method = "BH"
    )
  }
}

# Identificar SNPs significativos (α = 0.01)
significant_snps <- consistent_snps[
  !is.na(consistent_snps$P_adjusted) & consistent_snps$P_adjusted < 0.01, 
]

print(paste("SNPs significativos después de corrección FDR:", nrow(significant_snps)))

# Guardar resultados
write.csv(all_results, "allele_frequency_analysis_all.csv", row.names = FALSE)
write.csv(consistent_snps, "allele_frequency_analysis_consistent.csv", row.names = FALSE)
write.csv(significant_snps, "allele_frequency_analysis_significant.csv", row.names = FALSE)

print("Análisis completado. Archivos guardados:")
print("- allele_frequency_analysis_all.csv")
print("- allele_frequency_analysis_consistent.csv") 
print("- allele_frequency_analysis_significant.csv")
________________________________________
Paso 3: Manhattan Plot
r
# Manhattan plot para visualizar resultados
library(ggplot2)

create_manhattan_plot <- function(data, title) {
  # Preparar datos para Manhattan plot
  data$CHROM_NUM <- as.numeric(gsub("AaegL5_", "", data$CHROM))
  data <- data[!is.na(data$LOD_score) & !is.na(data$CHROM_NUM), ]
  
  # Calcular posiciones cumulativas
  data <- data[order(data$CHROM_NUM, data$POS), ]
  data$pos_cumulative <- 0
  
  cumulative_pos <- 0
  for(chrom in unique(data$CHROM_NUM)) {
    chrom_indices <- data$CHROM_NUM == chrom
    data$pos_cumulative[chrom_indices] <- data$POS[chrom_indices] + cumulative_pos
    cumulative_pos <- max(data$pos_cumulative[chrom_indices])
  }
  
  # Crear plot
  p <- ggplot(data, aes(x = pos_cumulative, y = LOD_score)) +
    geom_point(aes(color = factor(CHROM_NUM)), alpha = 0.7, size = 1) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.01), color = "red", linetype = "solid") +
    labs(
      title = title,
      x = "Chromosome",
      y = "-log10(p-value) [LOD score]",
      color = "Chromosome"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  return(p)
}

# Crear Manhattan plots para cada comparación
for(comp_name in unique(significant_snps$Comparison)) {
  comp_data <- significant_snps[significant_snps$Comparison == comp_name, ]
  if(nrow(comp_data) > 0) {
    p <- create_manhattan_plot(comp_data, paste("Manhattan Plot:", comp_name))
    ggsave(paste0("manhattan_plot_", gsub(" ", "_", comp_name), ".png"), 
           p, width = 12, height = 6, dpi = 300)
  }
}
________________________________________
Resumen del Pipeline
1.	Usa tus VCFs existentes - No necesitas devolverté
2.	Extrae información alélica de cada grupo (HPR, SP, etc.)
3.	Aplica exactamente la metodología del artículo:
o	Compara frecuencias de alelos alternativos entre fenotipos
o	Test de bondad de ajuste para consistencia dentro de réplicas
o	Tablas de contingencia y chi-cuadrado de heterogeneidad
o	Cálculo de LOD scores (-log10 p-values)
o	Corrección Benjamini-Hochberg por cromosoma
4.	Genera Manhattan plots para visualizar resultados
¿Quieres que ajuste algún parámetro específico o que te ayude con la interpretación de los resultados?
1.	https://pplx-res.cloudinary.com/image/private/user_uploads/12711040/d34ee22c-b784-4aec-925b-7c538626e6e2/image.jpg
2.	https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/12711040/0fb1e7a2-89bf-4f30-bbed-6a2bb22518f3/bucle_variantes.py
3.	https://www.reddit.com/r/bioinformatics/comments/1kijnyu/comparing_variant_call_data_in_a_vcf_file_with/
4.	https://cran.r-project.org/web/packages/contingencytables/contingencytables.pdf
5.	https://tuos-bio-data-skills.github.io/intro-stats-book/goodness-of-fit-tests.html
6.	https://popgen.nescent.org/StartSNP.html
7.	https://search.r-project.org/CRAN/refmans/simpleFDR/html/simFDR.html
8.	https://bookdown.org/mike/data_analysis/sec-false-discovery-rate.html
9.	https://bioinformatics.stackexchange.com/questions/15292/allele-count-and-allele-frequency-in-vcf-files
10.	https://pubmed.ncbi.nlm.nih.gov/35891781/
11.	https://speciationgenomics.github.io/filtering_vcfs/
12.	https://pmc.ncbi.nlm.nih.gov/articles/PMC10394849/
13.	https://www.datacamp.com/tutorial/chi-square-test-r
14.	https://cran.r-project.org/web/packages/fsthet/fsthet.pdf
15.	https://www.goldenhelix.com/blog/alternate-allele-frequency-vcf-file-format/
16.	https://biostatistics.letgen.org/mikes-biostatistics-book/inferences-categorical-data/chi-square-test-goodness-of-fit/
17.	https://manpages.ubuntu.com/manpages/xenial/man1/vcftools.1.html
18.	https://cran.r-project.org/web/packages/vcfR/vcfR.pdf
19.	https://www.bioconductor.org/packages/release/bioc/vignettes/TVTB/inst/doc/Introduction.html
20.	https://tuos-bio-data-skills.github.io/intro-stats-book/contingency-tables.html
21.	https://bookdown.org/hhwagner1/LandGenCourse_book/WE_3.html
22.	https://github.com/ANGSD/angsd/issues/248
23.	https://adegenet.r-forge.r-project.org/files/montpellier/practical-MVAintro.1.0.pdf
24.	https://vcftools.sourceforge.net/man_latest.html
25.	https://github.com/vcflib/vcflib/blob/master/README.md
26.	https://www.biorxiv.org/content/10.1101/2022.02.11.479970v1.full.pdf
27.	https://www.pnas.org/doi/pdf/10.1073/pnas.1114759108?download=true
28.	https://search.r-project.org/CRAN/refmans/GARCOM/html/vcf_counts_annot.html
29.	https://iastate.pressbooks.pub/quantitativegenetics/chapter/gene-frequencies/
30.	https://faculty.cnr.ncsu.edu/fikretisik/wp-content/uploads/sites/3/2015/06/SNP-marker-Data-Analysis-CTGN-Presentation.pdf
31.	https://pmc.ncbi.nlm.nih.gov/articles/PMC3154648/
32.	https://www.scielo.br/j/aabc/a/9dQtPV96bBcVHptpPdSpvtL/?format=pdf&lang=en
33.	https://www.biologysimulations.com/post/how-to-use-chi-squared-to-test-for-hardy-weinberg-equilibrium
34.	https://www.mun.ca/biology/scarr/Midterm_HWP.html
35.	https://researchguides.library.vanderbilt.edu/c.php?g=156859&p=3057842
36.	https://github.com/green-striped-gecko/dartR/blob/master/R/gl.report.heterozygosity.r
37.	https://search.r-project.org/CRAN/refmans/rCNV/html/h.zygosity.html
38.	https://www.youtube.com/watch?v=921IgQ6T6_4
39.	https://vcftools.github.io/man_0112a.html
40.	https://github.com/disulfidebond/VCF_Parsing_Analysis
41.	https://pmc.ncbi.nlm.nih.gov/articles/PMC9870988/
42.	https://toolshed.g2.bx.psu.edu/repository/display_tool?changeset_revision=cf2af5c3118c&render_repository_actions_for=tool_shed&repository_id=5f7d83aa4f577607&tool_config=%2Fsrv%2Ftoolshed%2Fmain%2Fvar%2Fdata%2Frepos%2F000%2Frepo_417%2Fallele-counts.xml

26-06-2025
#########################################################################
###############################26-06-2025#################################
#########################################################################

Voy a hacer mpileup y las modificaciones con varscan que me Saul Lozano específico
Hay que instalar Varscan

samtools mpileup -f /data/black_lab/ReferenceSeqs/test_Aaegl5_map/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta Col_perm_A1_replica1_sorted.bam > Col_perm_A1_replica1.pileup
Bucle mpileup

Ya hice el bucle en Python, se llama bucle_variantes.py

28-06-2025
#########################################################################
###############################28-06-2025#################################
#########################################################################

Hice un archivo de metadata. 

Según IA: 
Modelos estadísticos:
•	Para EWAS, usa regresión logística en PLINK ajustando por componentes principales (PCAs) para controlar ancestría6:

plink --vcf annotated.vcf --pheno pheno.txt --logistic --covar pca_covariates.txt  --out ewas_results

 Corrección y visualización
•	Corrección múltiple:
Aplica corrección de Bonferroni o FDR a los p-valores usando R:

results <- read.table("ewas_results.assoc.logistic", header=TRUE) results$FDR <- p.adjust(results$P, method="fdr")

30-06-2025
#########################################################################
###############################30-06-2025#################################
#########################################################################

Ya tengo los readcoutns con varscan?
Verificar. 

Tengo varias opciones. 

Pipeline Saul	Recomendado NSTC	Protocolo.io 
 Sortd.bam	BWA	Pileup.bam
Samtools mpileup 	Variant Calling: FreeBayes/Varscan (Correct for ploidy and pooled settings)	Pairwise FST:
Poolfstat (R)
Nozeros.pileup	Allele Frequency Extraction: Bcftools/Variant Caller (Allele Frequency Extraction)	Heatmap in R
Convertir mpileup to Readcounts: Varscan	Statistical Comparison: PoPoolation2  (CMH test of Fisher’s Exact Test in R)	Nucleotide diversity:
PoPoolation (Repetir para cada pool)
t test 
Linear regression 
Modifying Varscan readcount files. Readcounts_mod (Chrom-y t/)	Filtering Significant SNPs: Based on p-values, FDR correction 	
FORTRAN	Annotation: SnpEff or ANNOVAR for functional impact	
Visualization: 
Manhattan plots, volcano plots	Visualization: 
Manhattan plots, volcano plots	

TAREAS
MÉTRICAS DE BWA
MERGE DE LOS ARCHIVOS DE GATK 
Continuar con el pipeline
-	Hacerlo con el nozero.pileup 
•	perl /mnt/disc2/grupobcei/ewas/ewas_Acacias/popoolation2-master/mpileup2sync.pl --input HLR1_no_zeros.pileup --ouput HLR1.sync
•	perl /mnt/disc2/grupobcei/ewas/ewas_Acacias/popoolation2-master/mpileup2sync.pl --input HLR2_no_zeros.pileup --ouput HLR2.sync
•	perl /mnt/disc2/grupobcei/ewas/ewas_Acacias/popoolation2-master/mpileup2sync.pl --input SP1_no_zeros.pileup --ouput SP1.sync
•	perl /mnt/disc2/grupobcei/ewas/ewas_Acacias/popoolation2-master/mpileup2sync.pl --input SP2_no_zeros.pileup --ouput SP2.sync 

-	Si no da con ese hacer las correcciones que dice en protocols.io y seguir con ese protocolo

Need an mpileup for all pools and one for each pool individually.  

Trabajar con los .mpileup -> .sync -> Ejecutar analisis 
Trabajar con los archivos .readcounts

01-07-2025 Recapitulación SAM->BAM->samtools view -q 20	-> samtools sort -> samtools mpileup -B -Q 0 XX.mpileup -> mpileup2sync.jar -> HLR1_HLR2_SP1_SP2.sync -> poolfstat_R.r
#########################################################################
###############################01-07-2025#################################
#########################################################################

Hay que agregar un -B -Q 0  al código de samtools mpileup. 

Recapitulación: 
1.	Realicé el paso de samtools view y sort con los archivos .BAM y un calidad de 20. Esto está en script_original_popoolation2.bash. 
2.	Luego con el sort creé el mpileup 
3.	Luego con el mpileup creé el .sync
4.	Luego ese lo utilicé como insumo para poolfstat_R.r (Script en R) 
LA EVIDENCIA MUESTRAS QUE DEBIDO A LA BAJA COBERTURA, 



Thanks — this is great progress! Your relaxed filters worked:
✅ 745,071 SNPs retained after filtering — that's a strong dataset.
________________________________________
❗ Error Explanation
You got this error:
Error in is.data.frame(x) : 
  no slot of name "snp.read.counts" for this object of class "pooldata"
That means this line is failing:
write.csv(pool_file@snp.read.counts, file = "SNP_read_counts.csv", row.names = FALSE)
Because @snp.read.counts is not a valid slot in the new pooldata object — it likely doesn’t exist anymore in the version of the poolfstat package you're using.
________________________________________
✅ Fix: Check valid slots
You can see the available slots with:
slotNames(pool_file)
You’ll probably see something like:
[1] "refallele.readcount" "readcoverage" "snp.names" ...
To get SNP data or read counts, use the correct slots:
🔁 Instead of @snp.read.counts, try:
write.csv(pool_file@refallele.readcount, file = "refallele_readcounts.csv", row.names = FALSE)
write.csv(pool_file@readcoverage, file = "total_coverage.csv", row.names = FALSE)
These two matrices give you:
•	refallele.readcount: Number of reads supporting the reference allele for each SNP in each pool
•	readcoverage: Total read depth per SNP per pool
________________________________________
✅ Updated Code Snippet for Output
# Guardar matrices de lectura si existen
write.csv(pool_file@refallele.readcount, file = "refallele_readcounts.csv", row.names = FALSE)
write.csv(pool_file@readcoverage, file = "total_coverage.csv", row.names = FALSE)
If you want to join them into one table (optional):
readcount_df <- as.data.frame(pool_file@refallele.readcount)
coverage_df <- as.data.frame(pool_file@readcoverage)
combined_df <- cbind(readcount_df, coverage_df)
write.csv(combined_df, file = "SNP_readcounts_and_coverage.csv", row.names = FALSE)
________________________________________
✅ Summary
Fix	Action
❌ @snp.read.counts	🔄 Replace with @refallele.readcount and @readcoverage
✅	Use slotNames(pool_file) to explore valid components
💡	You can save per-SNP metrics or coverage matrices now
Would you like to also:
•	Export the list of SNP positions and alleles?
•	Get a summary table per chromosome?
Let me know!

JULIO 02-07-2025 Resumen en Excel de pipelines
o	El resumen de los pipelines está en inputs_outputs
o	Lo que estoy haciendo es una mezcla de samtools con el bam mapeado con BWA unos pocos filtros de ese samtools. 
-	Luego se crea un mpileup
•	Luego con mpileup2sync. Jar se crea el sync
•	Luego usando R con el archivo poolfstat_R.r
Esto está bueno: 
https://www.cyted.org/assets/img/redes/511/evento/1_TALLER_FENO_AXEL.pdf
diapo 65 para crear los archivos .ped

JULIO 07-07-2025 Pipeline GWAS Cristian Velarde modificado por Alejandro Mejía
Empieza desde los VCFs merged. -> el script original se llama (association_gwas_MOD.sh) el scrpt mío está en script_exoma.bash. 
Index – tabix
Plin.bim
AWK -> plink_fixed.bim
	Duplicados.dupbar
		Plink_final + Archivo phenotype.txt
			PLINK GWAS: resultados.gwas.assoc
			
JULIO 08-07-2025 Joint Genotyping
Se creó el archivo “Script_GATK_jointgenotyping_PLINK” para llamado de variantes desde 02_MarkDuplicates. 
En este se usa el llamado de variantes con el modo G y luego se combinan los gVCFs y luego se hace el joint genotyping. 

JULIO 09-07-2025 Repetición llamado de variantes con GATK jointgenotyping
Se crea combine_joint_GWAS.bash en servidor /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/02_MarkDuplicates/ para realizar: 
-	Combinado de gvcfs
-	Joint de gvcf
-	Comprimido
-	Index
-	Pling
-	GWAS 
Los resultados deberían de estar en /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/resultados_gwas.assoc > /mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/resultados_plink/
Lo que hay que hacer mañana: 
-	Ver si sí hay variantes y con esta nueva aproximación sí funciona.
•	Al parecer necesito mínimo 10 muestras para calcular el coeficiente de inbreeding. REVISAR
 
-	Mirar sí el genotipo sigue estando raro. 
-	Mirar lo de karla en el drive 
-	Mirar lo de llegar hasta especies en lo de microbiota. 
JULIO 10-07-2025 Reunión Omar Triana. Repetición llamado de variantes con GATK jointgenotyping
-	Los readgroups están malos. (RGSM). 
-	Hay que volver a crear los archivos .bam 
o	Se habló con Omar. Conclusiones: 
-	Mirar las bacterias burkholderia y Pseudomonas viridiflava. 
-	Hacer un merge entre los datos de RNAseq upregulated y las variantes encontradas. Mirar sí hay solapamiento.
JULIO 11-07-2025 Tareas pendientes y recapitulación. 
Se debe de hacer el joint genotyping nuevamente, basado en el RGSM (Read Groups Nuevo). Porque está tomando las dos réplicas como una sola. Hay que hacer prácticamente todo el “Pipeline GATK” y modificar el bucles_variantes.py. NO ES PRIORIDAD. 
Hay que hacer un merge. Ya sea con Python o con R de los archivos VCF y los genes Upregulated o downregualted. Pero primero hacer el filtro para solo los
Cómo leer un archivo VCF en Python o en R? R/ como tengo los archivos en .txt y en xlsx puedo hacerlo con R o con Python. 
Cómo quitar dejar solamente el Header del archivo VCF?R/ como ya es un archivo .txt es más fácil.
El código GATK_snpeff.py tiene la forma de poner los archivos snpeff.txt. 
Cómo hacer para filtrar únicamente las variantes que pasaron el filtro?

/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools view -f PASS HLR1_GATK_annotated.vcf -o HLR1_GATK_annotated_PASS.vcf







Restantes 
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


#########################################################################
#########################Espacio para cosas temporales#######################
#########################################################################

