#Script extraído de https://sourceforge.net/p/popoolation2/wiki/Tutorial/

#Remover reads mapeados ambiguamente 
#ESto se hace a partir del sam pero lo borré para hacer espacio
#Con el bam: 

samtools view -q 20 -b HLR1_BWA_sortd.bam | samtools sort -o HLR1_filtered_sorted.bam
samtools view -q 20 -b HLR2_BWA_sortd.bam | samtools sort -o HLR2_filtered_sorted.bam
samtools view -q 20 -b SP1_BWA_sortd.bam | samtools sort -o SP1_filtered_sorted.bam
samtools view -q 20 -b SP2_BWA_sortd.bam | samtools sort -o SP2_filtered_sorted.bam

#Crear los mpileups
/mnt/disc2/grupobcei/ewas/ewas_Acacias/03._ewas_mapped_subsample/samtools/samtools-1.21/samtools mpileup -B -Q 0 HLR1_filtered_sorted.bam HLR2_filtered_sorted.bam > HLR1_HLR2.mpileup
/mnt/disc2/grupobcei/ewas/ewas_Acacias/03._ewas_mapped_subsample/samtools/samtools-1.21/samtools mpileup -B HLR1_filtered_sorted.bam > HLR1.mpileup
/mnt/disc2/grupobcei/ewas/ewas_Acacias/03._ewas_mapped_subsample/samtools/samtools-1.21/samtools mpileup -B HLR2_filtered_sorted.bam > HLR2.mpileup
##Este lo hice y funcionó
/mnt/disc2/grupobcei/ewas/ewas_Acacias/03._ewas_mapped_subsample/samtools/samtools-1.21/samtools mpileup -B -Q 0 HLR1_filtered_sorted.bam HLR2_filtered_sorted.bam SP1_filtered_sorted.bam SP2_filtered_sorted.bam > HLR1_HLR2_SP1_SP2.mpileup

#Crear los sync
# El original: perl /mnt/disc2/grupobcei/ewas/ewas_Acacias/popoolation2-master/mpileup2sync.pl --input muestra.mpileup --ouput muestra.sync

#El mejorado
/mnt/disc2/grupobcei/java/jdk-17.0.12/bin/java "-Xmx30G" -jar /mnt/disc2/grupobcei/ewas/ewas_Acacias/popoolation2-master/mpileup2sync.jar --input HLR1_HLR2_SP1_SP2.mpileup --output HLR1_HLR2_SP1_SP2.sync --fastq-type sanger --min-qual 20 --threads 30

#Luego de esto, el resultado HLR1_HLR2_SP1_SP2.sync lo puse a correr en R con el script "poolfstat_R.r"