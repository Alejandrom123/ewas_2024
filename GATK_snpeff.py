import os
import glob
import subprocess
import shutil

# Rutas
SNPEFF_JAR = "/mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.jar"
JAVA = "/mnt/disc2/grupobcei/java/jdk-24.0.1/bin/java"
SNPEFF_DB = "Aedes_aegypti_lvpagwg"
SNPEFF_CONFIG = "/mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.config"

# Carpeta de salida para los resultados
output_dir = "05_SnpEffAnnotation"
os.makedirs(output_dir, exist_ok=True)

# Procesar todos los archivos *_filt_10x.vcf
for vcf_file in glob.glob("*_filt_10x.vcf"):
    sample = vcf_file.replace("_filt_10x.vcf", "")
    annotated_vcf = os.path.join(output_dir, f"{sample}_GATK_annotated.vcf")
    html_report = os.path.join(output_dir, f"{sample}_snpEff.html")
    genes_report = os.path.join(output_dir, f"{sample}_snpEff_genes.txt")

    print(f"Procesando {vcf_file}...")

    # Comando SnpEff
    cmd = [
        JAVA, "-jar", SNPEFF_JAR,
        SNPEFF_DB,
        "-c", SNPEFF_CONFIG,
        "-stats", html_report,
        vcf_file
    ]

    with open(annotated_vcf, "w") as out_vcf:
        subprocess.run(cmd, stdout=out_vcf, check=True)

    # Mover el archivo snpEff_genes.txt generado al nombre específico
    if os.path.exists("snpEff_genes.txt"):
        shutil.move("snpEff_genes.txt", genes_report)

    print(f"Listos: {annotated_vcf}, {html_report}, {genes_report}")

print("¡Anotación SnpEff completada para todos los archivos!")
