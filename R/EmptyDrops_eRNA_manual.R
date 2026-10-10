getwd()

# -------- Instalar librerias si hace falta --------

lib_extra <- "~/R/x86_64-pc-linux-gnu-library/4.6"

## BiocManager se requiere para varios paquetes
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}

#paquetes_cran <- c("ggplot2", "Seurat", "Matrix", "SoupX")
#paquetes_bioc <- c("SingleCellExperiment", "scuttle", "scater", "scDblFinder", "DropletUtils", "biomaRt")

#instalar_si_falta <- function(paquetes, origen = c("CRAN", "Bioconductor")) {
#  origen <- match.arg(origen)
#  faltantes <- paquetes[!sapply(paquetes, requireNamespace, quietly = TRUE)]
#  if (length(faltantes) > 0) {
#    message("Instalando desde ", origen, ": ", paste(faltantes, collapse = ", "))
#    if (origen == "CRAN") {
#      install.packages(faltantes)
#    } else {
#      BiocManager::install(faltantes, update = FALSE, ask = FALSE)
#    }
#  } else {
#    message("Todos los paquetes de ", origen, " ya estaban instalados.")
#  }
#}

# instalar_si_falta(paquetes_cran, "CRAN")
# instalar_si_falta(paquetes_bioc, "Bioconductor")


# -------- Cargar librerias -------

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ggplot2)
  library(Seurat)
  library(Matrix)
  library(DropletUtils)
  library(SoupX)
  library(scran)
})

.libPaths(c(lib_extra, .libPaths()))
suppressPackageStartupMessages({
  library(biomaRt)
})

# -------- Cargar datos -------

sce <- readRDS("/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-A.rds")
figures <- "/mnt/data/bioinfo3/ysosa/proyecto-equipo-A/figuras"

set.seed(6)

# -------- Crear la curva de rango --------
totales <- colSums(counts(sce))
br <- barcodeRanks(counts(sce))

knee       <- metadata(br)$knee
inflection <- metadata(br)$inflection
c(codo = knee, inflexion = inflection)

orden_br <- order(br$rank)

# png(file.path(figures, "Curva_de_rango.png"), width = 1600, height = 1200, res = 300, type = "cairo")
plot(br$rank[orden_br], br$total[orden_br], log = "xy",
     xlab = "Rango del código de barras", ylab = "Total de cuentas",
     main = "Curva de rango", type = "l", lwd = 2)
abline(h = knee,       col = "aquamarine2",  lty = 2, lwd = 2)
abline(h = inflection, col = "coral1", lty = 2, lwd = 2)
legend("bottomleft", bty = "n", lty = 2, lwd = 2,
       col = c("aquamarine2", "coral1"),
       legend = c(paste("codo:", knee),
                  paste("inflexión:", inflection)))
dev.off()

# -------- Utilizar EmptyDrops --------
emptyD_100 <- emptyDrops(counts(sce))
llamadas_emptyD_100 <- !is.na(emptyD_100$FDR) & emptyD_100$FDR <= 0.001

emptyD_50 <- emptyDrops(counts(sce), lower = 50, niters = 10000)
llamadas_emptyD_50 <- !is.na(emptyD_50$FDR) & emptyD_50$FDR <= 0.001

sce_100 <- sce[, llamadas_emptyD_100]
cat("\nLas dimensiones después del filtrado por emptyDrops utilizando lower = 100 son las siguientes: \n")
dim(sce_100) 

sce_50 <- sce[, llamadas_emptyD_50]
cat("\n\nLas dimensiones después del filtrado por emptyDrops utilizando lower = 50 son las siguientes: \n")
dim(sce_50)

# ========== Filtrado Ambiental ============

# Ver el perfil del ambiente
soup_channel_100 <- SoupChannel(tod = counts(sce), toc = counts(sce_100))
# head(soup_channel_100$soupProfile[order(-soup_channel_100$soupProfile$est), ], 20)

perfil <- soup_channel_100$soupProfile
perfil$symbol <- rowData(sce)$Symbol[match(rownames(perfil), rownames(sce))]

# Top 20 con symbol
cat("\nEl perfil del ambiente con límite 100\n")
head(perfil[order(-perfil$est), ], 20)

# limit 50
soup_channel_50 <- SoupChannel(tod = counts(sce), toc = counts(sce_50))

perfil2 <- soup_channel_50$soupProfile
perfil2$symbol <- rowData(sce)$Symbol[match(rownames(perfil2), rownames(sce))]

# Top 20 con symbol
cat("\n\nEl perfil2 del ambiente con límite 50\n")
head(perfil2[order(-perfil2$est), ], 20)

# ==== Genes de Hemoglobina e inmunoglobulina para eliminar el eRNA de fondo manualmente ====

# Le pasamos la lista de genes de hemoglobina y hemoglobulina que pueden ser un marcador de fondo
genes_manual <- list(HB = c("HBB", "HBA1", "HBA2"),
              IG = c("IGKC", "IGHG1", "IGHA1", "JCHAIN"))

# Creamos una función para checar que los IDs de ENSEMBL estén dentro de nuestro objeto SCE
mapear_genes <- function(simbolos, sce){
  ensembl_ids <- rownames(sce)[rowData(sce)$Symbol %in% simbolos]
  if (length(ensembl_ids) == 0){
    stop("No se encontraron IDs de Ensembl para: ", paste(simbolos, collapse = ", "))
  }
  ensembl_ids
}

genes_ensembl <- lapply(genes_manual, mapear_genes, sce = sce)

cat("\nA continuación se muestran los genes marcadores (Ensembl ID) a usar para estimar la contaminación: \n")
genes_ensembl

# -------- Verificación cruzada: símbolo original vs. símbolo en el dataset --------
verificacion_genes <- do.call(rbind, lapply(names(genes_ensembl), function(grupo) {
  ids <- genes_ensembl[[grupo]]
  data.frame(
    grupo          = grupo,
    ensembl_id     = ids,
    symbol_dataset = rowData(sce)[ids, "Symbol"],
    stringsAsFactors = FALSE
  )
}))

cat("\nVerificación cruzada de los genes marcadores: \n")
verificacion_genes

set.seed(6)
clusters_rapidos <- quickCluster(sce_50)
names(clusters_rapidos) <- colnames(sce_50)
soup_channel_50 <- setClusters(soup_channel_50, clusters_rapidos)

usar <- estimateNonExpressingCells(soup_channel_50, nonExpressedGeneList = genes_ensembl, clusters = clusters_rapidos)

# -------- Calcular fracción de contaminación usando ambos grupos --------
soup_channel_50 <- calculateContaminationFraction(soup_channel_50, genes_ensembl, useToEst = usar, forceAccept = TRUE)

cat("\nLa fracción de contaminación de la primer célula es: ", round(soup_channel_50$metaData$rho[1], 4), "\n")
cat("\nLa fracción de contaminación en promedio dentro del dataset es: ", round(mean(soup_channel_50$metaData$rho), 4), "\n")

# Descontaminar dataset original
limpio_50 <- adjustCounts(soup_channel_50, roundToInt = TRUE, verbose = 0)
head(limpio_50)

# =========== Guardar el dataset final ==========

dir_salida <- "../datos"

saveRDS(limpio_50, file.path(dir_salida, "sce_limpio_50_manual.rds"))
