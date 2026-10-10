suppressPackageStartupMessages(library(SingleCellExperiment))

sce_orig <- readRDS("/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-A.rds")

cat("Clase:", class(sce_orig), "\n")
cat("Columnas de rowData:", colnames(rowData(sce_orig)), "\n")
stopifnot("No existe la columna Symbol en rowData" = "Symbol" %in% colnames(rowData(sce_orig)))

anotacion <- data.frame(ID = rownames(sce_orig),
                        Symbol = rowData(sce_orig)$Symbol)

cat("Genes:", nrow(anotacion), "| sin símbolo:", sum(is.na(anotacion$Symbol)), "\n")
cat("Ejemplo de IDs:", head(anotacion$ID, 3), "\n")

saveRDS(anotacion, "../datos/anotacion_genes.rds")
cat("Terminado")
