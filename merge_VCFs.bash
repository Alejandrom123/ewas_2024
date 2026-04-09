#Unir VCFs de muestras diferentes diferentes

# Comprimir los archivos VCF
bgzip *.vcf

#Indexar HLS1, HLS2, SP1, SP2 
bcftools index *vcf.gz

#Bucle para index 
for f in *.vcf.gz; do
        /mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools index "$f"
done

#Merge
/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools merge --force-samples HLS1.variants.vcf.gz HLS2.variants.vcf.gz SP1.variants.vcf.gz SP2.variants.vcf.gz -Oz -o HLS1_HLS2_SP1_SP2.variants.vcf.gz

##Modificar el header
   # Write out the header to be modified
/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools view -h HLS1_HLS2_SP1_SP2.variants.vcf.gz > header.txt

   # Edit the header using your favorite text editor
   nano header.txt

   # Reheader the file
/mnt/disc2/grupobcei/ewas/ewas_Acacias/bcftools/bcftools-1.22/bcftools reheader -h header.txt -o reHLS1_HLS2_SP1_SP2.variants.vcf.gz HLS1_HLS2_SP1_SP2.variants.vcf.gz

###El archivo VCF ya está listo###




