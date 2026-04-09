import os #Manipula rutas- directorios e incluso ejecuta comandos
import subprocess #ejecuta otras cosas, otros "Procesos"
from glob import glob 

INPUT_DIR = "/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/"
INDEX_PATH = "/mnt/disc2/grupobcei/ewas/index/bowtie2_index/aedes_genome_index"
THREADS = 30
OUTPUT_SUFFIX = "_subsample_bowtie2.sam"

def encontrar_pares(directorio):
    # Listar archivos R1 y R2
    archivos_r1 = glob(os.path.join(directorio, "*_R1_subsample.fastq"))
    archivos_r2 = glob(os.path.join(directorio, "*_R2_subsample.fastq"))
    
    # Crear diccionarios para emparejar por nombre base
    dict_r1 = {}
    dict_r2 = {}
    
    for f in archivos_r1:
        base = os.path.basename(f).replace("_R1_subsample.fastq", "")
        dict_r1[base] = f
    
    for f in archivos_r2:
        base = os.path.basename(f).replace("_R2_subsample.fastq", "")
        dict_r2[base] = f
    
    # Emparejar solo las bases que tengan ambos archivos
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
            "bowtie2",
            "--threads", str(THREADS),
            "-x", INDEX_PATH,
            "-1", r1,
            "-2", r2,
            "-S", archivo_salida
        ]
        
        print(f"🔍 Procesando: {base}")
        print(f"   R1: {os.path.basename(r1)}")
        print(f"   R2: {os.path.basename(r2)}")
        
        try:
            subprocess.run(comando, check=True)
            print(f"✅ {base} completado\n")
        except subprocess.CalledProcessError as e:
            print(f"❌ Error en {base}: {e}\n")
            if os.path.exists(archivo_salida):
                os.remove(archivo_salida)
