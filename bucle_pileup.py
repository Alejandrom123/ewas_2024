import os
import subprocess

# Ruta al genoma de referencia
REFERENCE = "/mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta"

# Lista de muestras según tu pipeline
samples = [
    "HLR1", "HLR2", "HLS1", "HLS2", "LL1", "LL2",
    "HPR1", "HPR2", "HPS1", "HPS2", "LP1", "LP2", "SP1", "SP2"
]

def run_mpileup_for_samples():
    for sample in samples:
        input_bam = f"{sample}_BWA_sortd.bam"
        output_pileup = f"{sample}.pileup"
        if not os.path.isfile(input_bam):
            print(f"Archivo no encontrado: {input_bam}, saltando...")
            continue
        cmd = [
            "samtools", "mpileup",
            "-f", REFERENCE,
            input_bam
        ]
        print(f"Ejecutando: {' '.join(cmd)} > {output_pileup}")
        with open(output_pileup, "w") as out:
            subprocess.run(cmd, stdout=out, check=True)

if __name__ == "__main__":
    run_mpileup_for_samples()
