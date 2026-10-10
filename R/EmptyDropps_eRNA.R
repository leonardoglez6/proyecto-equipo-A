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

# ==== Genes de Hemoglobina para eliminar el eRNA de fondo ====

# La anotación fue realizada con ENSEMBL110 por lo que se conecto a Ensembl de esta versión
#ensembl_110_genes <- useEnsembl(biomart = "genes",
#                          dataset = "hsapiens_gene_ensembl",
#                          host = "https://jul2023.archive.ensembl.org")

# Se consulta a ENSEMBL para extraer información de los genes indicados para hemoglobina
#hb_biomart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
#                                   "chromosome_name", "gene_biotype"),
#                    filters = "external_gene_name",
#                    values = c("HBA1", "HBA2", "HBB", "HBD", "HBE1",
#                               "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"),
#                    mart = ensembl_110_genes)

# Si no funciona la base de datos de ENSEMBL están los IDs manuales
hb_biomart <- data.frame(
  symbol      = c("HBA1", "HBA2", "HBB", "HBD", "HBE1",
                  "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"),
  ensembl_id  = c("ENSG00000206172", # HBA1
                  "ENSG00000188536", # HBA2
                  "ENSG00000244734", # HBB
                  "ENSG00000223609", # HBD
                  "ENSG00000213931", # HBE1
                  "ENSG00000213934", # HBG1
                  "ENSG00000196565", # HBG2
                  "ENSG00000206177", # HBM
                  "ENSG00000086506", # HBQ1
                  "ENSG00000130656"),# HBZ
  stringsAsFactors = FALSE
)

cat("\nA continuación se muestran los genes de la hemoglobina que pueden causar señal de fondo en nuestro dataset: \n")
hb_biomart

# Cruzar contra los genes presentes en el dataset para ver que coincidan los IDs
hemoglobina <- intersect(hb_biomart$ensembl_id, rownames(sce))

if (length(hemoglobina) == 0) {
  stop("No se encontraron genes de hemoglobina en rownames(sce). Revisar IDs.")
}

# Verificación cruzada: símbolo reportado por biomaRt vs. símbolo en rowData del dataset
verificacion_hb <- data.frame(
  ensembl_id      = hemoglobina,
  symbol_biomart  = hb_biomart$symbol[match(hemoglobina, hb_biomart$ensembl_id)],
  symbol_dataset  = rowData(sce)[hemoglobina, "Symbol"]
)
cat("\nA continuación vemos los IDs de ENSEMBL de los genes de hemoglobina junto a los IDs de nuestro dataset para asegurarnos que sean los mismos: \n")
verificacion_hb

# Pasar lista de genes de hemoglobina
usar <- estimateNonExpressingCells(soup_channel_50, nonExpressedGeneList = list(HB = hemoglobina), clusters = FALSE, )

# Calcular fracción de contaminación
soup_channel_50 <- calculateContaminationFraction(soup_channel_50, list(HB = hemoglobina), useToEst = usar, forceAccept = TRUE)
# Calcular fracción de primera célula
cat("\nLa fracción de contaminación de la primer célula es: ", round(soup_channel_50$metaData$rho[1], 4), "\n")
# Calcular media para todas las células
cat("\nLa fracción de contaminación en promedio dentro del dataset es: ", round(mean(soup_channel_50$metaData$rho), 4), "\n")

# Descontaminar dataset original
limpio_50 <- adjustCounts(soup_channel_50, roundToInt = TRUE, verbose = 0)
head(limpio_50)


