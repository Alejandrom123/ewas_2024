import os
import subprocess

# Lista de muestras según tu pipeline
samples = [
    "HLR1", "HLR2", "HLS1", "HLS2", "LL1", "LL2",
    "HPR1", "HPR2", "HPS1", "HPS2", "LP1", "LP2", "SP1", "SP2"
]

# Rutas
REFERENCE = "/mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta"
JAVA_EXE = "/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java"
VARSCAN_JAR = "/mnt/disc2/grupobcei/ewas/ewas_Acacias/Varscan/VarScan.v2.4.6.jar"

# Carpetas de salida por paso
output_dirs = {
    "06_mpileup": "06_mpileup",
    "07_nozeros": "07_nozeros",
    "08_readcounts": "08_readcounts",
    "09_readcounts_mod": "09_readcounts_mod"
}

# Crear carpetas si no existen
for d in output_dirs.values():
    os.makedirs(d, exist_ok=True)

def run_command(cmd):
    print(f"Ejecutando: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def process_samples():
    for sample in samples:
        input_bam = f"{sample}_BWA_sortd.bam"
        pileup_file = os.path.join(output_dirs["06_mpileup"], f"{sample}.pileup")
        no_zeros_pileup = os.path.join(output_dirs["07_nozeros"], f"{sample}_no_zeros.pileup")
        readcounts_file = os.path.join(output_dirs["08_readcounts"], f"{sample}.readcounts")
        readcounts_mod_file = os.path.join(output_dirs["09_readcounts_mod"], f"{sample}.readcounts")

        if not os.path.isfile(input_bam):
            print(f"Archivo no encontrado: {input_bam}, saltando...")
            continue

        # 1. Ejecutar samtools mpileup
        print(f"Generando mpileup para {sample}")
        cmd_mpileup = ["samtools", "mpileup", "-f", REFERENCE, input_bam]
        with open(pileup_file, "w") as out:
            subprocess.run(cmd_mpileup, stdout=out, check=True)

        # 2. Eliminar líneas con 0 reads (columna 4 igual a 0)
        print(f"Eliminando líneas con 0 reads en {pileup_file}")
        with open(pileup_file, 'r') as infile, open(no_zeros_pileup, 'w') as outfile:
            for line in infile:
                cols = line.strip().split()
                if len(cols) > 3 and cols[3] != '0':
                    outfile.write(line)

        # 3. Convertir mpileup a readcounts con Varscan
        print(f"Generando archivo readcounts para {sample}")
        cmd_varscan = [
            JAVA_EXE, "-jar", VARSCAN_JAR, "readcounts", no_zeros_pileup,
            "--min-coverage", "25", "--min-base-qual", "30", "--output-file", readcounts_file
        ]
        run_command(cmd_varscan)

        # 4. Modificar archivo readcounts (reemplazar cromosomas y tabs)
        print(f"Modificando archivo readcounts para {sample}")
        replacements = {"AaegL5_1": "1", "AaegL5_2": "2", "AaegL5_3": "3"}
        with open(readcounts_file, 'r') as infile, open(readcounts_mod_file, 'w') as outfile:
            for line in infile:
                for old, new in replacements.items():
                    line = line.replace(old, new)
                line = line.replace("\t", " ")
                outfile.write(line)

if __name__ == "__main__":
    process_samples()
