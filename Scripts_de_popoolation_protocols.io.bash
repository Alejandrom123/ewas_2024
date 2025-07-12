#Pasos para popoolation
#En popoolation https://www.protocols.io/view/pool-seq-pipeline-stanford-et-al-e6nvwk372vmk/v1?step=10.1
#Los pasos son: Filter mapped reads

java -Xmx50g -XX:ParallelGCThreads=5 -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=$filename \
O=/path/to/outfile/directory/${filename2}_nodup.bam \
M=/path/to/outfile/directory/${filename2}_nodup.txt \
REMOVE_DUPLICATES=TRUE \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 MAX_RECORDS_IN_RAM=50000 ASSUME_SORT_ORDER=coordinate SORTING_COLLECTION_SIZE_RATIO=0.1

#Filtrar por calidad 

samtools view -@ 8 -q20 $filename \
-o /path/to/q20/file/directory/${filename2}_q20.bam

#Quitar los reads que no están pareados 

samtools view -@ 8 -f2 $filename \
-o /path/to/f2/file/directory/${filename2}_q20.bam/${filename2}_f2.bam

#Validar los Reads - Picard

java -Xmx50g -XX:ParallelGCThreads=5 -jar $EBROOTPICARD/picard.jar ValidateSamFile \
I=$filename\
MODE=SUMMARY \
MAX_OPEN_TEMP_FILES=1000 \
MAX_RECORDS_IN_RAM=50000

#Ahí sí crear el mpileup - Samtools 
###Warning!! Hay que hacerlo aparentemente con versiones actuales de Samtools###
    #Para todos los pools?
samtools mpileup -B XX_f2.bam YY_f2.bam ZZ_f2.bam \
-o all_pools_ordered.mpileup

    #Para las muestras individuales
samtools mpileup -B $filename -o ${filename2}.mpileup


