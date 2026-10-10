# ------- Cargar librerias ----------
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scuttle)
  library(scater)
  library(scDblFinder)
  library(ggplot2)
})

# Cargar datos filtrados por células vacías y ambiente (EmptyDrops y SoupX)
sce <- readRDS("../datos/sce_limpio_50_auto.rds")
figures <- "../figuras"
dir_salida <- "../datos"

# Usado en proceso DblFinder
set.seed(6)

cat("Clase del objeto:", class(sce), "\n")
cat("Dimensiones:", dim(sce), "\n")

if (is(sce, "dgCMatrix")) {
  sce <- SingleCellExperiment(assays = list(counts = sce))
}

head(rownames(sce))

anotacion <- readRDS("../datos/anotacion_genes.rds")

rowData(sce)$ID     <- rownames(sce)
rowData(sce)$Symbol <- anotacion$Symbol[match(rownames(sce), anotacion$ID)]

cat("\nGenes sin símbolo:", sum(is.na(rowData(sce)$Symbol)), "de", nrow(sce), "\n")

# Guardar simbolos del dataset para usarlos dentro del análisis
simbolos <- rowData(sce)$Symbol

# ============= MITOCONDRIALES ==========
# Obtener lista de genes mitocondriales con prefijo MT
mitocondriales <- grep("^MT-", simbolos)
cat("\nGenes mitocondriales encontrados:", length(mitocondriales), "\n")
cat(simbolos[mitocondriales], "\n")

# Obtener lista de genes ribosomales con prefijo RP
ribosomales <- grep("^RP[SL][0-9]+[A-Z]?[0-9]*$|^RPLP[0-9]$|^RPSA$", simbolos)

sce <- addPerCellQCMetrics(sce, subsets = list(MT = mitocondriales, RIBO = ribosomales))

cat("\nCuentas totales por célula: \n")
print(summary(sce$sum))

cat("\nMétricas de porcentaje de genes mitocondriales por célula: \n")
print(summary(sce$subsets_MT_percent))

# Guardar imagenes de QC mitocondrial
# Cuentas totales
png(file.path(figures, "MT_cuentas_totales.png"), width = 1600, height = 1200, res = 300, type = "cairo")
print(
  plotColData(sce, y = "sum") + 
  scale_y_log10() + labs(y = "cuentas totales (log)")
)
dev.off()

# Porcentaje mitocondrial
png(file.path(figures, "MT_porcentajes.png"), width = 1600, height = 1200, res = 300, type = "cairo")
print(
  plotColData(sce, y = "subsets_MT_percent") +
  labs(y = "% mitocondrial")
)
dev.off()

### Comparación de umbrales mitocondriales ###
# a) Umbral fijo (mayoría tutoriales), descarta los TRUE
fijo <- !(sce$subsets_MT_percent < 5 & sce$detected > 200)

# b) MAD adaptativo (solo tenemos un lote), descarta TRUE
qc_mt <- perCellQCFilters(sce, sub.fields = "subsets_MT_percent")
mad_mt <- qc_mt$discard

cat("\nDescarte de células por umbral: \n")
comp_mit = data.frame(
  criterio = c("fijo de tutorial", "MAD adaptativo"),
  descarta = c(sum(fijo), sum(mad_mt)),
  pct_descartado = round(100 * c(sum(fijo), sum(mad_mt)) / ncol(sce), 1)
)
print(comp_mit)

# Desarrollo de cada métrica
cat("\nDesglose por métrica (MAD): \n")
print(colSums(as.matrix(qc_mt[, c("low_lib_size", "low_n_features", "high_subsets_MT_percent")])))

# Concordancia entre ambos criterios
cat(" \nComparación fijo VS MAD: \n")
print(table(fijo = fijo, MAD = mad_mt))

# ============= RIBOSOMALES ==========
cat("\nGenes ribosomales encontrados:", length(ribosomales), "\n")
cat(rownames(sce)[ribosomales], "\n")

cat("\nMétricas de porcentaje de genes ribosomales por célula: \n")
print(summary(sce$subsets_RIBO_percent))

# Porcentaje ribosomal
png(file.path(figures, "RIBO_porcentajes.png"), width = 1600, height = 1200, res = 300, type = "cairo")
print(
  plotColData(sce, y = "subsets_RIBO_percent") +
    labs(y = "% ribosomal")
)
dev.off()

# Relación MT vs RIBO
png(file.path(figures, "MT_vs_RIBO.png"), width = 1600, height = 1200, res = 300, type = "cairo")
print(
  plotColData(sce, x = "subsets_RIBO_percent", y = "subsets_MT_percent") +
    labs(x = "% ribosomal", y = "% mitocondrial")
)
dev.off()

# Valores atípicos altos por MAD
mad_ribo <- isOutlier(sce$subsets_RIBO_percent, type = "higher", nmads = 3)

cat("\nCélulas con % ribosomal atípicamente alto (MAD, 3 nmads): ",
    sum(mad_ribo), " (", round(100 * sum(mad_ribo) / ncol(sce), 1), "%)\n", sep = "")

# Comparación con descarte mitocondrial
cat("\nSolapamiento con descarte mitocondrial (MAD): \n")
print(table(RIBO_alto = mad_ribo, MAD_mt = mad_mt))

# ======= DOBLETES ========
# Clasifica cada célula como singlet/doublet
sce <- scDblFinder(sce)

# Singlets y doublets encontrados en db
cat("\nDobletes encontrados: \n")
print(table(sce$scDblFinder.class))

# Guardar imagen de células singlets y doublets
png(file.path(figures, "dobletes.png"), width = 1600, height = 1200, res = 300, type = "cairo")
print(plotColData(sce, x = "sum", y = "scDblFinder.score",
            colour_by = "scDblFinder.class") +
            scale_x_log10() + labs(x = "cuentas totales (log)", y = "puntuación")
)
dev.off()

# Obtener objeto de dobletes
dobletes = sce$scDblFinder.class == "doublet"

# ===== FILTRADO FINAL =========
# Descartar por MAD mitocondriales y dobletes
sce$descartar <- mad_mt | dobletes

cat("\nCélulas descartadas por criterio:\n")
print(table(MAD = mad_mt, doblete = dobletes))

# Guardar dataset filtrado
sce_filtrado <- sce[, !sce$descartar]

# Imprimir resultados finales de células tras filtrado
cat("\nAntes:", ncol(sce), "celulas | Después:", ncol(sce_filtrado), "\n")

# Salvar dataset final después de QC
saveRDS(sce_filtrado, file.path(dir_salida, "sce_limpio_MT_dbl.rds"))

