import os
import subprocess

def run_samtools_pipeline():
    # Listar todos los archivos que terminan en _subsample_BWA.sam
    sam_files = [f for f in os.listdir('.') if f.endswith('_subsample_BWA.sam')]

    for sam_file in sam_files:
        # Extraer prefijo antes de _subsample_BWA.sam
        prefix = sam_file.replace('_subsample_BWA.sam', '')

        bam_file = f"{prefix}_subsample_BWA.bam"
        sorted_bam = f"{prefix}_BWA_sortd.bam"
        index_file = f"{prefix}_sortd_index.bam"

        print(f"Procesando muestra: {sam_file}")
        
        # 1. Convertir SAM a BAM
        cmd_view = ['samtools', 'view', '-S', '-b', sam_file, '-o', bam_file]
        print(f"Ejecutando: {' '.join(cmd_view)}")
        subprocess.run(cmd_view, check=True)

        # 2. Ordenar BAM
        cmd_sort = ['samtools', 'sort', '-o', sorted_bam, bam_file]
        print(f"Ejecutando: {' '.join(cmd_sort)}")
        subprocess.run(cmd_sort, check=True)

        # 3. Indexar BAM ordenado
        cmd_index = ['samtools', 'index', '-o', index_file, sorted_bam]
        print(f"Ejecutando: {' '.join(cmd_index)}")
        subprocess.run(cmd_index, check=True)

        print(f"Procesamiento completado para {sam_file}\n")

if __name__ == "__main__":
    run_samtools_pipeline()