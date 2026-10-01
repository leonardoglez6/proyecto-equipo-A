getwd()

# -------- Instalar librerias si hace falta --------

## BiocManager se requiere para varios paquetes
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}

#paquetes_cran <- c("ggplot2", "Seurat", "Matrix", "SoupX")
#paquetes_bioc <- c("SingleCellExperiment", "scuttle", "scater", "scDblFinder", "DropletUtils")

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

#instalar_si_falta(paquetes_cran, "CRAN")
#instalar_si_falta(paquetes_bioc, "Bioconductor")


# -------- Cargar librerias -------

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ggplot2)
  library(Seurat)
  library(Matrix)
  library(DropletUtils)
  library(SoupX)
})


# -------- Cargar datos -------

sce <- readRDS("/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-A.rds")
figures <- "/mnt/data/bioinfo3/ysosa/proyecto-equipo-A/figuras"

# Set seed
set.seed(6)

# -------- Crear la curva de rango --------
totales <- colSums(counts(sce))
br <- barcodeRanks(counts(sce))

knee       <- metadata(br)$knee
inflection <- metadata(br)$inflection
c(codo = knee, inflexion = inflection)

orden_br <- order(br$rank)

png(file.path(figures, "Curva_de_rango.png"), width = 1600, height = 1200, res = 300, type = "cairo")
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
