import os
import subprocess

# Directorios de salida
output_dirs = [
    "01_AddOrReplaceReadGroups",
    "02_MarkDuplicates",
    "03_HaplotypeCaller",
    "04_VariantFiltration",
    "05_SnpEffAnnotation"
]

for d in output_dirs:
    os.makedirs(d, exist_ok=True)

# Rutas a herramientas (ajustar según tu entorno)
PICARD = "/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/picard/picard.jar"
GATK = "/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java -jar /mnt/disc2/grupobcei/gatk/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar"
SNPEFF = "/mnt/disc2/grupobcei/java/jdk-24.0.1/bin/java -jar /mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.jar"
REFERENCE = "/mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta"
SNPEFF_DB = "Aedes_aegypti_lvpagwg"

# Tabla RG info
RG_info = {
    "HLR1": {"RGID": "HLR1", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "HLR1", "RGSM": "HLR1"},
    "HLR2": {"RGID": "HLR2", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "HLR2", "RGSM": "HLR2"},
    "HLS1": {"RGID": "HLS1", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "HLS1", "RGSM": "HLS1"},
    "HLS2": {"RGID": "HLS2", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "HLS2", "RGSM": "HLS2"},
    "LL1": {"RGID": "LL1", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "LL1", "RGSM": "LL1"},
    "LL2": {"RGID": "LL2", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "LL2", "RGSM": "LL2"},
    "HPR1": {"RGID": "HPR1", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "HPR1", "RGSM": "HPR1"},
    "HPR2": {"RGID": "HPR2", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "HPR2", "RGSM": "HPR2"},
    "HPS1": {"RGID": "HPS1", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "HPS1", "RGSM": "HPS1"},
    "HPS2": {"RGID": "HPS2", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "HPS2", "RGSM": "HPS2"},
    "LP1": {"RGID": "LP1", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "LP1", "RGSM": "LP1"},
    "LP2": {"RGID": "LP2", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "LP2", "RGSM": "LP2"},
    "SP1": {"RGID": "SP1", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "SP1", "RGSM": "SP1"},
    "SP2": {"RGID": "SP2", "RGLB": "lib1", "RGPL": "illumina", "RGPU": "SP2", "RGSM": "SP2"},
}

samples = list(RG_info.keys())

def run_command(cmd):
    print(f"Ejecutando: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def main():
    for sample in samples:
        print(f"Procesando muestra: {sample}")

        input_bam = f"{sample}_BWA_sortd.bam"
        if not os.path.isfile(input_bam):
            print(f"Archivo no encontrado: {input_bam}, saltando...")
            continue

        rg_bam = f"01_AddOrReplaceReadGroups/{sample}_BWA_rg.bam"
        nodups_bam = f"02_MarkDuplicates/{sample}_BWA_nodups.bam"
        metrics = f"02_MarkDuplicates/{sample}_BWA.metrics"
        variants = f"03_HaplotypeCaller/{sample}.variants.vcf"
        filt_variants = f"04_VariantFiltration/{sample}_filt_10x.vcf"
        annotated_variants = f"05_SnpEffAnnotation/{sample}_GATK_annotated.vcf"

        # 1. AddOrReplaceReadGroups
        cmd_rg = [
            *PICARD.split(),
            "AddOrReplaceReadGroups",
            f"I={input_bam}",
            f"O={rg_bam}",
            "SO=coordinate",
            "CREATE_INDEX=true",
            f"RGID={RG_info[sample]['RGID']}",
            f"RGLB={RG_info[sample]['RGLB']}",
            f"RGPL={RG_info[sample]['RGPL']}",
            f"RGPU={RG_info[sample]['RGPU']}",
            f"RGSM={RG_info[sample]['RGSM']}"
        ]
        run_command(cmd_rg)

        # 2. MarkDuplicates
        cmd_md = [
            *PICARD.split(),
            "MarkDuplicates",
            f"I={rg_bam}",
            f"O={nodups_bam}",
            f"M={metrics}",
            "REMOVE_DUPLICATES=true",
            "CREATE_INDEX=true"
        ]
        run_command(cmd_md)

        # 3. HaplotypeCaller
        cmd_hc = [
            *GATK.split(),
            "HaplotypeCaller",
            "-I", nodups_bam,
            "-R", REFERENCE,
            "-O", variants,
            "--native-pair-hmm-threads", "10"
        ]
        run_command(cmd_hc)

        # 4. VariantFiltration
        cmd_vf = [
            *GATK.split(),
            "VariantFiltration",
            "-R", REFERENCE,
            "-V", variants,
            "--filter-name", "FAIL",
            "--filter-expression", "QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10 || MQRankSum < -12.5 || ReadPosRankSum < -8.0",
            "-O", filt_variants
        ]
        run_command(cmd_vf)

        # 5. SNPeff annotation
        cmd_snp = [
            *SNPEFF.split(),
            SNPEFF_DB,
            "-c", "/mnt/disc2/grupobcei/ewas/ewas_Acacias/snpEff/snpEff.config",
            filt_variants
        ]
        with open(annotated_variants, "w") as out_vcf:
            print(f"Ejecutando: {' '.join(cmd_snp)} > {annotated_variants}")
            subprocess.run(cmd_snp, stdout=out_vcf, check=True)

        print(f"Procesamiento completado para {sample}\n")

if __name__ == "__main__":
    main()
