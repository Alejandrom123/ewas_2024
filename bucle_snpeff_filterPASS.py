import os
import glob
import subprocess
import shutil

# Rutas a herramientas
BCFTOOLS = "/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools"
JAVA = "/mnt/disc2/grupobcei/java/jdk-24.0.1/bin/java"
SNPEFF_JAR = "/mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.jar"
SNPEFF_DB = "Aedes_aegypti_lvpagwg"
SNPEFF_CONFIG = "/mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.config"

# Carpeta de salida para los resultados
output_dir = "05_SnpEffAnnotation"
os.makedirs(output_dir, exist_ok=True)

for vcf_file in glob.glob("*_GATK_annotated.vcf"):
    sample = vcf_file.replace("_GATK_annotated.vcf", "")
    pass_vcf = os.path.join(output_dir, f"{sample}_PASS.vcf")
    html_report = os.path.join(output_dir, f"{sample}_snpEff.html")
    genes_report = os.path.join(output_dir, f"{sample}_snpEff_genes.txt")

    print(f"Filtrando variantes PASS en {vcf_file} ...")
    # 1. Filtrar variantes con FILTER=PASS
    cmd_bcftools = [
        BCFTOOLS, "view", "-f", "PASS", vcf_file, "-o", pass_vcf
    ]
    subprocess.run(cmd_bcftools, check=True)

    print(f"Anotando variantes con SnpEff para {pass_vcf} ...")
    # 2. Anotar con SnpEff el archivo filtrado
    cmd_snpeff = [
        JAVA, "-jar", SNPEFF_JAR,
        SNPEFF_DB,
        "-c", SNPEFF_CONFIG,
        "-stats", html_report,
        pass_vcf
    ]
    annotated_vcf = os.path.join(output_dir, f"{sample}_PASS_annotated.vcf")
    with open(annotated_vcf, "w") as out_vcf:
        subprocess.run(cmd_snpeff, stdout=out_vcf, check=True)

    # 3. Renombrar el archivo snpEff_genes.txt generado al nombre específico
    if os.path.exists("snpEff_genes.txt"):
        shutil.move("snpEff_genes.txt", genes_report)

    print(f"Listos: {pass_vcf}, {annotated_vcf}, {html_report}, {genes_report}")

print("¡Filtrado y anotación completados para todos los archivos!")
