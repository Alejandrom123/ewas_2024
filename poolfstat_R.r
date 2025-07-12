# Carga la librería poolfstat
library(poolfstat)

# Intentar cargar y procesar el archivo .sync
pool_file <- tryCatch({
  popsync2pooldata(
    sync.file = "./HLR1_HLR2_SP1_SP2.sync",
    poolsizes = c(15, 15, 15, 15),
    poolnames = c("HLR1", "HLR2", "SP1", "SP2"),
    min.rc = 1,
    min.cov.per.pool = 2,       # más permisivo
    max.cov.per.pool = 1000,    # permite más cobertura
    min.maf = 0.0001,           # variantes aún más raras
    noindel = TRUE
  )
}, error = function(e) {
  cat("❌ Error al crear el objeto PoolData:\n", e$message, "\n")
  NULL
})

# Verifica si pool_file fue creado con éxito
if (!is.null(pool_file)) {

  # Guardar el objeto PoolData
  save(pool_file, file = "pool_file_relaxed.RData")

  # Mostrar cuántos SNPs fueron retenidos
  cat("✅ SNPs retenidos después del filtrado:", pool_file@nsnp, "\n")

  # (Opcional) Guardar tabla de recuento por SNP
  write.csv(pool_file@snp.read.counts, file = "SNP_read_counts.csv", row.names = FALSE)

  # Calcular FST pareado (Anova)
  PW_fst <- compute.pairwiseFST(
    pool_file,
    method = "Anova",
    min.cov.per.pool = 2,
    max.cov.per.pool = 1000,
    min.maf = 0.0001,
    output.snp.values = TRUE
  )

  # Guardar resultados
  save(PW_fst, file = "PW_fst_relaxed.RData")
  write.csv(PW_fst@PairwiseFSTmatrix, file = "PW_fst_relaxed.csv")

  # Mensaje final
  cat("🎉 Análisis de FST completado exitosamente con filtros relajados.\n")

} else {
  cat("⚠️ Se abortó el análisis debido a errores en el archivo .sync o en los filtros.\n")
}

### Lo anterior generó un error pero se obtuvo pool_file_relaxed.RData con lo que se puede continuar

# Cargar objeto pool_file desde archivo .RData
load("pool_file_relaxed.RData")
slotNames(pool_file)
#Probablemente te devolverá algo como: [1] "refallele.readcount" "readcoverage" "snp.names" "chr.info" ...

# Guardar los recuentos del alelo de referencia
write.csv(pool_file@refallele.readcount, file = "refallele_readcounts.csv", row.names = FALSE)

# Guardar la cobertura total por SNP y por pool
write.csv(pool_file@readcoverage, file = "total_coverage.csv", row.names = FALSE)

# (Opcional) Guardar ambas matrices combinadas
readcount_df <- as.data.frame(pool_file@refallele.readcount)
coverage_df <- as.data.frame(pool_file@readcoverage)
combined_df <- cbind(readcount_df, coverage_df)
write.csv(combined_df, file = "SNP_readcounts_and_coverage.csv", row.names = FALSE)


#####################Aquí todo el resto###############################
# Cargar el objeto ya procesado
load("pool_file_relaxed.RData")

# Verifica qué slots tiene
print(slotNames(pool_file))

# Guardar los datos válidos
write.csv(pool_file@refallele.readcount, file = "refallele_readcounts.csv", row.names = FALSE)
write.csv(pool_file@readcoverage, file = "total_coverage.csv", row.names = FALSE)

# Guardar matriz combinada (opcional)
combined_df <- cbind(as.data.frame(pool_file@refallele.readcount),
                     as.data.frame(pool_file@readcoverage))
write.csv(combined_df, file = "SNP_readcounts_and_coverage.csv", row.names = FALSE)

cat("✅ Exportación corregida realizada sin volver a correr todo el análisis.\n")

#######################################################################
##########Código completo corregido####################################
#######################################################################

