rm(list = ls())

library(data.table)
library(doParallel)
library(foreach)

# Fuente de funciones externas (ajusta la ruta si es necesario)
source("/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/pileup_varscan/09_readcounts_mod/easy_chi2_fun.r")

## Parámetros
numberOfProcesors <- 20
markThreshold <- 0.05

## Sufijos cromosomas
chromosomes <- c("c1", "c2", "c3")

## Define aquí el prefijo (XXX) de la muestra que quieres procesar
samplePrefix <- "HLR1"  # <--- Cambia aquí por el prefijo deseado, ejemplo "HLS"

# Rutas base del directorio output
baseOutputPath <- "/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/pileup_varscan/09_readcounts_mod/data/output/"

for (chr in chromosomes) {
  cat(sprintf("Procesando muestra: %s cromosoma: %s\n", samplePrefix, chr))
  
  polySites_file <- paste0(baseOutputPath, samplePrefix, "vSP_avd_", chr, ".rds")
  txtOutputFile <- paste0(baseOutputPath, samplePrefix, "vSP_ezchi_", chr, ".chi")
  rawEzchiResults_file <- paste0(baseOutputPath, samplePrefix, "vSP_raw_ezchi_", chr, ".rds")
  rdsEzChiResults_file <- paste0(baseOutputPath, samplePrefix, "vSP_ezchi_", chr, ".rds")
  
  # Verificar que el archivo exista, si no existe salta al siguiente
  if (!file.exists(polySites_file)) {
    warning(sprintf("Archivo no encontrado: %s. Se omite esta iteración.\n", polySites_file))
    next
  }
  
  polySites <- readRDS(polySites_file)
  
  # Ajustes preliminares (eliminar columnas y modificar ref)
  polySites$ref <- polySites$ref1
  polySites$ref1 <- NULL
  polySites$ref2 <- NULL
  polySites$ref3 <- NULL
  polySites$ref4 <- NULL
  
  polySites$chrom1 <- NULL
  polySites$chrom2 <- NULL
  polySites$chrom3 <- NULL
  polySites$chrom4 <- NULL
  
  polySites$sumDepth1 <- NULL
  polySites$sumDepth2 <- NULL
  polySites$sumDepth3 <- NULL
  polySites$sumDepth4 <- NULL
  
  nLines <- nrow(polySites)
  
  cl <- parallel::makeCluster(numberOfProcesors)
  doParallel::registerDoParallel(cl)
  
  rawEzChiResults <- foreach(i = 1:nLines, .combine = 'rbind') %dopar% {
    GetEasyChiEstimates(polySite = polySites[i, ])
  }
  
  stopCluster(cl)
  
  saveRDS(rawEzChiResults, file = rawEzchiResults_file, compress = FALSE)
  
  rm(polySites)
  
  ezChiResults <- as.data.table(rawEzChiResults, key = "nucPosition")
  rm(rawEzChiResults)
  
  ezChiResults[, totalProb := 1 - pchisq(q = totalChiSqr, df = totalDegFreedom), by = nucPosition]
  ezChiResults[, group1Prob := 1 - pchisq(q = group1ChiSqr, df = group1DegFreedom), by = nucPosition]
  ezChiResults[, group2Prob := 1 - pchisq(q = group2ChiSqr, df = group2DegFreedom), by = nucPosition]
  
  bhThreshold <- -log10(GetBenjaminiHochberThreshold(ezChiResults$totalProb))
  
  ezChiResults <- ezChiResults[ezChiResults$lod > bhThreshold]
  
  ezChiResults[, inconsistency := MarkInconsistency(
    chi1 = group1ChiSqr,
    degFreedom1 = group1DegFreedom,
    inconsistencyMark1 = "1*",
    chi2 = group2ChiSqr,
    degFreedom2 = group2DegFreedom,
    inconsistencyMark2 = "2*",
    markThreshold = markThreshold),
    by = nucPosition]
  
  ezChiResults[, alleles := GetAllelesLabel(
    nucPosition = nucPosition,
    refNucleotide = refNuc, 
    As = As,
    Cs = Cs,
    Gs = Gs,
    Ts = Ts,
    Is = Is,
    Ds = Ds),
    by = nucPosition]
  
  saveRDS(ezChiResults, file = rdsEzChiResults_file, compress = FALSE)
  
  nRemainingSites <- nrow(ezChiResults)
  
  ezChiResults$alleles <- sprintf("%-6s", ezChiResults$alleles)
  ezChiResults$alleles <- gsub(" ", "_", ezChiResults$alleles)
  
  fileConn <- file(txtOutputFile, open = "wt")
  
  header <- paste(
    "SNPID,MUTATION,FREQ(ALIVE),FREQ(DEAD),LOD,HET(ALL),HET(ALIVE),HET(DEAD),",
    "FREQ(A),FREQ(C),FREQ(G),FREQ(T),FREQ(I),FREQ(D),CHISQ(ALL),CHISQ(ALIVE),",
    "CHISQ(DEAD),DF(ALL),DF(ALIVE),DF(DEAD)", sep = "")
  writeLines(header, fileConn)
  
  printFormatTemp <- data.frame(
    nucPosition  = "%10.0f",
    alleles = "%6s",
    group1AltAllFreq = "%8.5f",
    group2AltAllFreq = "%8.5f",
    lod = "%7.2f",
    group1Heteroz = "%8.5f",
    group2Heteroz = "%8.5f",
    totalHeteroz = "%8.5f",
    As = "%8.0f",
    Cs = "%7.0f",
    Gs = "%8.0f",
    Ts = "%8.0f",
    Is = "%8.0f",
    Ds = "%8.0f",
    group1ChiSqr = "%11.5f",
    group2ChiSqr = "%10.5f",
    totalChiSqr = "%10.5f",
    group1DegFree = "%5.0f",
    group2DegFree = "%5.0f",
    totalDegFree = "%5.0f",
    inconsistency = "%3s"
  )
  
  printFormat <- paste(printFormatTemp[, ], collapse = ",")
  
  for (i in 1:nRemainingSites) {
    line <- sprintf(printFormat,
                    ezChiResults$nucPosition[i],    #1
                    ezChiResults$alleles[i],        #2
                    ezChiResults$group1AltAllFreq[i],#3
                    ezChiResults$group2AltAllFreq[i],#4
                    ezChiResults$lod[i],            #5
                    ezChiResults$totalHeteroz[i],   #6
                    ezChiResults$group1Heteroz[i],  #7
                    ezChiResults$group2Heteroz[i],  #8
                    ezChiResults$As[i],             #9
                    ezChiResults$Cs[i],             #10
                    ezChiResults$Gs[i],             #11
                    ezChiResults$Ts[i],             #12
                    ezChiResults$Is[i],             #13
                    ezChiResults$Ds[i],             #14
                    ezChiResults$group1ChiSqr[i],  #15
                    ezChiResults$group2ChiSqr[i],  #16
                    ezChiResults$totalChiSqr[i],   #17
                    ezChiResults$group1DegFree[i], #18
                    ezChiResults$group2DegFree[i], #19
                    ezChiResults$totalDegFree[i],  #20
                    ezChiResults$inconsistency[i]  #21
    )
    writeLines(line, fileConn)
  }
  
  close(fileConn)
  
  cat(sprintf("Procesamiento terminado para muestra: %s, cromosoma: %s\n\n", samplePrefix, chr))
}

cat("Procesamiento completo para la muestra y cromosomas indicados.\n")
