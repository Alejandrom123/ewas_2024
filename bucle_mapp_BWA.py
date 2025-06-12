import os
from glob import glob
import subprocess

INPUT_DIR = "/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample"
REFERENCE = "/mnt/disc2/grupobcei/ewas/index/VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta"
THREADS = 30
OUTPUT_SUFFIX = "_subsample_BWA.sam"

def encontrar_pares(directorio):
    archivos_r1 = glob(os.path.join(directorio, "*_R1_subsample.fastq"))
    archivos_r2 = glob(os.path.join(directorio, "*_R2_subsample.fastq"))
    
    dict_r1 = {}
    dict_r2 = {}
    
    for f in archivos_r1:
        base = os.path.basename(f).replace("_R1_subsample.fastq", "")
        dict_r1[base] = f
    
    for f in archivos_r2:
        base = os.path.basename(f).replace("_R2_subsample.fastq", "")
        dict_r2[base] = f
    
    pares = {}
    for base in dict_r1:
        if base in dict_r2:
            pares[base] = (dict_r1[base], dict_r2[base])
        else:
            print(f"⚠️ No se encontró archivo R2 para la muestra {base}")
    
    return pares

if __name__ == "__main__":
    muestras = encontrar_pares(INPUT_DIR)
    
    for base, (r1, r2) in muestras.items():
        archivo_salida = os.path.join(INPUT_DIR, f"{base}{OUTPUT_SUFFIX}")
        
        if os.path.exists(archivo_salida):
            print(f"⏩ Saltando {base} (archivo existente)")
            continue
        
        comando = [
            "bwa", "mem",
            "-t", str(THREADS),
            "-M",
            REFERENCE,
            r1,
            r2
        ]
        
        print(f"🔍 Procesando: {base}")
        print(f"   R1: {os.path.basename(r1)}")
        print(f"   R2: {os.path.basename(r2)}")
        
        try:
            with open(archivo_salida, "w") as out_sam:
                subprocess.run(comando, check=True, stdout=out_sam)
            print(f"✅ {base} completado\n")
        except subprocess.CalledProcessError as e:
            print(f"❌ Error en {base}: {e}\n")
            if os.path.exists(archivo_salida):
                os.remove(archivo_salida)
