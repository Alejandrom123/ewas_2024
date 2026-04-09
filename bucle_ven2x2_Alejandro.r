rm(list = ls())
library(data.table)

# Prefijos para archivos de las dos muestras que vas a comparar
phen1_1_baseName <- "LP1"  # Ejemplo: "HLS1"
phen1_2_baseName <- "LP2"  # Ejemplo: "HLS2"

# Etiquetas personalizables que aparecerán en la salida (pueden o no ser iguales a los prefijos)
phen1_1Name <- "LP1"
phen1_2Name <- "LP2"

# Etiquetas para tus muestras de referencia constantes
phen2_1Name <- "SP1"
phen2_2Name <- "SP2"

# Sufijos para cromosomas (archivos)
suffixes <- c("c1", "c2", "c3")

# Rutas base donde están tus archivos
base_path_phen1 <- "/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/pileup_varscan/09_readcounts_mod/data/midput/"
base_path_phen2 <- "/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/pileup_varscan/09_readcounts_mod/data/input/"

for (suf in suffixes) {
  cat("Procesando cromosoma:", suf, "\n")

  # Construir rutas para los archivos .rds según sufijos y prefijos
  phen1_1_file <- paste0(base_path_phen1, phen1_1_baseName, "_", suf, ".rds")
  phen1_2_file <- paste0(base_path_phen1, phen1_2_baseName, "_", suf, ".rds")
  phen2_1_file <- paste0(base_path_phen2, phen2_1Name, "_", suf, ".rds")
  phen2_2_file <- paste0(base_path_phen2, phen2_2Name, "_", suf, ".rds")
  
  # Leer los archivos
  phen1_1 <- readRDS(phen1_1_file)
  phen1_2 <- readRDS(phen1_2_file)
  phen2_1 <- readRDS(phen2_1_file)
  phen2_2 <- readRDS(phen2_2_file)
  
  # Definir nombres para los archivos de salida, usando el sufijo del cromosoma
  rdsFile <- paste0("/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/pileup_varscan/09_readcounts_mod/data/output/", phen1_1_baseName, "vSP_avd_", suf, ".rds")
  textFile <- paste0("/mnt/disc2/grupobcei/ewas/ewas_Acacias/02.1_ewas_trimmed_subsample/pileup_varscan/09_readcounts_mod/data/output/", phen1_1_baseName, "vSP_avd_", suf, ".tbl")

  # Mergear muestras de cada grupo por posición
  p1_1xp1_2 <- merge(x = phen1_1, y = phen1_2, by = "pos", suffixes = c("1", "2"))
  p2_1xp2_2 <- merge(x = phen2_1, y = phen2_2, by = "pos", suffixes = c("3", "4"))
  
  # Merge final entre los dos grupos
  p1Vp2 <- merge(x = p1_1xp1_2, y = p2_1xp2_2, by = "pos")
  
  rm(p1_1xp1_2)
  rm(p2_1xp2_2)
  
  # Filtrar sitios polimórficos: filas en las que la suma de valores 0 es <= 19
  pol <- p1Vp2[rowSums(p1Vp2 == 0) <= 19, ]
  
  # Guardar el objeto pol en formato RDS
  saveRDS(pol, file = rdsFile)
  
  # Escribir archivo texto con formato personalizado
  nRows <- nrow(pol)
  fileConn <- file(textFile, open = "wt")
  
  for (i in 1:nRows) {
    line1 <- sprintf("%10i", pol$pos[i])
    line2 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
                     phen1_1Name, pol$a1[i], pol$c1[i], pol$g1[i], pol$t1[i], pol$i1[i], pol$d1[i])
    line3 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
                     phen1_2Name, pol$a2[i], pol$c2[i], pol$g2[i], pol$t2[i], pol$i2[i], pol$d2[i])
    line4 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
                     phen2_1Name, pol$a3[i], pol$c3[i], pol$g3[i], pol$t3[i], pol$i3[i], pol$d3[i])
    line5 <- sprintf("%s%20i%10i%10i%10i%10i%10i",
                     phen2_2Name, pol$a4[i], pol$c4[i], pol$g4[i], pol$t4[i], pol$i4[i], pol$d4[i])
    
    writeLines(c(line1, line2, line3, line4, line5), fileConn)
  }
  
  close(fileConn)
  
  cat("Terminado procesamiento para cromosoma:", suf, "\n\n")
}

cat("Procesamiento completo para todos los cromosomas.\n")
