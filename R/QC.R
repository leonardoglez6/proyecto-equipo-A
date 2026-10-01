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
  library(scuttle)
  library(scater)
  library(scDblFinder)
  library(ggplot2)
  library(Seurat)
  library(Matrix)
  library(DropletUtils)
  library(SoupX)
})


# -------- Cargar datos -------

sc_data <- readRDS("/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-A.rds")
cat("\n\nA continuación se muestra la clase de nuestro objeto (dataset):\n")
class(sc_data) #Checamos si es SCE o Seurat


# -------- PREGUNTA 1 | Son cuentas crudas? Si hay decimales, no lo son. -------
cat(" -------- PREGUNTA 1 | Son cuentas crudas? Si hay decimales, no lo son. -------")

if(all(counts(sc_data)@x == round(counts(sc_data)@x, 0))){ #Esta linea de comando compara la matriz de cuentas donde no hay 0s para que no muera por el poder de cómputo y checa si el valor redondeado da lo mismo
  cat("\n Las cuentas son crudas.\n")
} else{
  cat("\n Las cuentas son procesadas.\n")
}


# -------- PREGUNTA 2 | ¿Las células ya están filtradas? --------
cat("# -------- PREGUNTA 2 | ¿Las células ya están filtradas? --------\n")

# Número de genes y barcodes
n_genes = nrow(sc_data)
n_barcodes = ncol(sc_data)

# Revisar si existen métricas de QC por célula (si tiene columnas tendría información nueva)
tiene_QC <- length(colData(sc_data)) > 0

# Revisar el origen del objeto
origen <- metadata(sc_data)$derivado_de
# Este nombre es bastante común en datos de 10xGenomics y Cell Ranger
es_raw <- grepl("raw_feature_bc_matrix", origen, ignore.case = TRUE)

# Evaluación manual
cat("Número de genes:", n_genes, "\n")
cat("Número de barcodes:", n_barcodes, "\n")
cat("¿Tiene métricas de QC por célula?:", tiene_QC, "\n")
cat("Origen:", origen, "\n\n")

# Evaluar con código
# Se asume que es raro que un dataset filtrado tenga más de unos pocos miles de barcodes
if (es_raw && n_barcodes > 100000) {
	cat("Conclusión: Las células NO parecen estar filtradas.\n")
	cat("- El objeto contiene", n_barcodes, "barcodes.\n")
	cat("- El objeto viene de 'raw_feature_bc_matrix'.\n")
  
  # Existen métricas de calidad generadas?
  if (!tiene_QC) {
	cat("- No contiene métricas de QC..\n")
  }
  
  # Si no se puede confirmar, se debería reportar para preguntar a personas correspondientes
} else if (n_barcodes > 100000) {
	cat("Conclusión: Las células probablemente NO están filtradas.\n")
	cat("- El número de barcodes es muy elevado:", n_barcodes, "\n")
} else {
	cat("Conclusión: No se puede determinar solamente con estas lecturas.\n")
}
cat("\n\n")

# -------- PREGUNTA 3 | ¿Dónde están los metadatos de lote, donante y condición? --------
cat("# -------- PREGUNTA 3 | ¿Dónde están los metadatos de lote, donante y condición? --------")

verificar <- function(nombre, genes, celulas, suma, tiene_lote) {
  data.frame(objeto = nombre, genes = genes, celulas = celulas,
             suma = suma, lote = tiene_lote)
}

resultados <- verificar(
  "SingleCellExperiment",
  nrow(sc_data), ncol(sc_data),
  sum(counts(sc_data)),
  "lote" %in% colnames(colData(sc_data))
)

cat("\n\n")
resultados

# Metadatos de las columnas/células
cat("\n\nA continuación se muestran los metadatos de las células (columnas):\n")
colData(sc_data)
int_colData(sc_data)

# Metadatos de las filas/genes
cat("\n\nA continuación se muestran los metadatos de los genes (filas):\n")
rowData(sc_data)

# Metadatos
cat("\n\nA continuación se muestran los metadatos del dataset:\n")
metadata(sc_data)
int_metadata(sc_data)


# -------- PREGUNTA 4 | ¿Qué versión de la anotación de genes se usó?
cat("\n -------- PREGUNTA 4 | ¿Qué versión de la anotación de genes se usó?----------\n")
head(rowData(sc_data))

cat("\nNotamos que la anotación proviene de Ensembl\n")
cat("Buscando en la ruta del dataset, observamos la procedencia de los datos\n")

texto <- readLines("/mnt/data/bioinfo3/compartido/10x/PROCEDENCIA.md")

inicio <- grep("^## refdata", texto)
fin <- grep("^## Verificación", texto)

cat(texto[inicio:(fin - 1)], sep = "\n")
