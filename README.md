# Proyecto final — equipo A

> Plantilla del módulo de análisis de células individuales, LCG 2027-1.

## Qué es este repositorio

El análisis del conjunto de datos que se le asignó a nuestro equipo, de principio a
fin, de forma que **cualquier persona pueda reproducirlo desde cero**. Esa es la
condición que se evalúa (criterio A2 de la rúbrica, 30 %): no que el resultado sea
bonito, sino que corra en otra máquina sin intervención manual.

## Integrantes

| Nombre | Qué hizo |
|---|---|
|Leonardo Daniel González López | Copiar plantilla de proyecto, Reproducibilidad, Pregunta 1 y 3, y Control de Calidad (EmptyDrops, SoupX) y reporte.qmd |
|Yuliana Denisse Sosa Gómez | Leer los datos, Pregunta 2 y 4, y Control de calidad (EmptyDrops, mitocondrial, ribosomal y dobletes) y README |

Esta tabla no es burocracia: el criterio B de la rúbrica es una respuesta
individual dirigida, y saber quién trabajó en qué parte es lo que la hace justa.

## Cómo reproducir este análisis

```bash
git clone https://github.com/leonardoglez6/proyecto-equipo-A.git
cd proyecto-equipo-A
quarto render reporte.qmd
```

**Los datos no están en el repositorio** y no deben estarlo: `.gitignore` los
excluye. Declara aquí de dónde salen.

| Dato | Dónde vive | Cómo se obtuvo |
|---|---|---|
| Dataset Equipo A | `/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-A.rds` | [Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Manual (channel 5) from 10X Genomics](https://www.10xgenomics.com/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-manual-channel-5-3-1-standard-3-1-0) |

### ¿Cómo se corrieron los códigos en orden?

1. El código `preguntas.R` ayuda a contestar las 4 preguntas iniciales para poder entender el conjunto de datos.
2. El código `EmptyDrops_eRNA_auto.R` ayuda a realizar el primer filtro de calidad. Incluye: filtrado de gotas vacías (EmptyDrops) y filtrado del perfil ambiental con opción automática (SoupX)
   * Input: `proyecto-equipo-A.rds`. Output: `sce_limpio_50_auto.rds`
   * Este código se puede correr con Slurm, con el código `eD_SX_auto.sh`.
   * También, se tiene la versión `EmptyDrops_eRNA_manual.R`, que se puede correr en Slurm con `eD_SX_manual.sh`. En este se deben escoger los genes de forma manual.
4. El código `anotacion_genes.rds` ayudará a obtener las anotaciones que se perdieron al generar el nuevo RDS.
5. El código `MT_RP_Dbl.R` ayuda a realizar el segundo filtro de calidad. Incluye: filtrado por proporción mitocondrial, por proporción ribosomal y por dobletes.
   * Input: `sce_limpio_50_auto.rds` y `anotacion_genes.rds` (anotaciones). Output: sce_limpio_MT_dbl.rds
   * Este código se puede correr con Slurm, con el código `MT_job.sh`.

## Estructura

| Carpeta | Qué va aquí |
|---|---|
| `reporte.qmd` | El reporte ejecutable. Es el entregable |
| `R/` | Funciones que uses en más de un lugar |
| `datos/` | Vacía a propósito: los datos no se versionan |
| `figuras/` | Salida del reporte. Tampoco se versiona |

```
proyecto-equipo-A 			# Carpeta principal del proyecto
├── datos								# Datos generados durante los scripts
│   ├── anotacion_genes.rds			# Anotación de datos (symbols)
│   ├── proyecto-equipo-A.rds -> /mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-A.rds	# Dataset original
│   ├── sce_limpio_50_auto.rds			# Dataset obtenido de EmptyDrops_eRNA_auto.R	
│   ├── sce_limpio_50_manual.rds		# Dataset obtenido de EmptyDrops_eRNA_manual.R
│   └── sce_limpio_MT_dbl.rds			# Dataset obtenido de MT_RP_Dbl.R
├── figuras
│   ├── autoEstCont_rho.png 			# Estimaciones de filtrado ambiental automático
│   ├── Curva_de_rango.png				# Umbral fijo para células vacías
│   ├── dobletes.png					# Dobletes y singlets
│   ├── MT_cuentas_totales.png			# Cuentas totales por célula
│   ├── MT_porcentajes.png				# Proporción mitocondrial
│   ├── MT_vs_RIBO.png					# Porcentaje mitocondrial VS ribosomal
│   └── RIBO_porcentajes.png			# Proporción ribosomal
├── outs								# Outlogs de códigos (exitosos y fallidos)
│   └── ...
├── R									# Códigos utilizados
│   ├── anotacion.R					# Recuperar anotación de RDS
│   ├── eD_SX_auto.sh					# Bash que corre EmptyDrops_eRNA_auto.R
│   ├── eD_SX_manual.sh				# Bash que corre EmptyDrops_eRNA_manual.R
│   ├── EmptyDrops_eRNA_auto.R			# QC: EmptyDrops y SoupX (genes escogidos manualmente)
│   ├── EmptyDrops_eRNA_manual.R		# QC: EmptyDrops y SoupX (genes escogidos por código)
│   ├── MT_job.sh						# Bash que corre MT_RP_Dbl.R
│   ├── MT_RP_Dbl.R					# QC: Mitocondriales, Ribosomales y doublets				
│   └── preguntas.R					# Código para resolver las 4 preguntas
├── README.md							# Este archivo
└── reporte.qmd						# Reporte generado con interpretaciones
```

## Bitácora de decisiones

Cada vez que tomes una decisión de análisis —un umbral, un método, un filtro—
anótala aquí con su razón. **La tabla de decisiones es el criterio A3 de la
rúbrica** (15 %), y llenarla al final, de memoria, se nota.

| Fecha | Decisión | Por qué | Qué se probó antes |
|---|---|---|---|
| 1-octubre-2026 | Utilizar filtrado de cuentas vacías con `lower = 50` (límite inferior de UMIs) de EmptyDrops | Decidimos utilizar este filtro debido a que el umbral fijo suele perder todas las células pequeñas por el poco RNA. Usar EmptyDrops es más confiable debido a que es una prueba estadística que compara directamente con el perfil ambiental. Al probar con lower = 100 y lower = 50 se obtuvieron números de células muy parecidos (3999 y 4085, respectivamente), lo cual nos llevó a confiar en el resultado que conservó un poco más de células para no perder información (`lower = 50`). | Se probó el umbral fijo (se generó la imagen de Curva de Rango), EmptyDrops con `lower = 100` y EmptyDrops con `lower = 50`. |
| 6-octubre-2026 | Utilizar filtrado automático de genes de ambiente con SoupX | Decidimos utilizar esta opción debido a que al filtrar por genes de Hemoglobina salió un porcentaje muy alto de contaminación (66.85%), lo cual no es común. Por lo tanto, nos dimos cuenta que este dataset ya estaba filtrado de eritrocitos, lo que ocasionaba este error. Por lo que teníamos 2 opciones: escoger los genes manualmente o realizar una estimación automática con el paquete de SoupX. Para el filtrado manual, se escogieron algunos genes recomendados de Hemoglobina e Immunoglobina, que dieron una contaminación de 2.74%. Por otro lado, para la estimación automática, se realizó primero un clustering rápido y después la selección mediante la función autoEstCont() y la estimación batesiana de la fracción de contaminación (rho), lo cual resultó en una contaminación de 1%. Decidimos escoger este último resultado debido a que se escogen los genes marcadores respecto a los clusters y posee menos sesgos de los que podríamos tener nosotros al escoger genes | Filtrado con muchos genes de Hemoglobina, filtrado con pocos genes de hemoglobina e immunoglobulina, filtrado automático con clustering rápido. |
| 9-octubre-2026 | Filtrado mitocondrial con MAD | Se encontraron 13 genes mitocondriales, de los cuales los porcentajes de genes mitocondriales se encontraron con una mediana de 4.76% y una media de 8.7%. Se probó un umbral fijo de tutorial y un MAD adaptativo global. El umbral fijo descartaba el 47.8% de las células, mientras el MAD adaptativo solo 27.2%. Decidimos realizar el filtrado con el MAD adaptativo, debido a que se adapta a nuestros datos al marcar células atípicas que se apartan más de 3 desviaciones estándar y no elimina tantas células como el umbral estándar. | Umbral estándar y MAD mitocondrial. | 
| 9-octubre-2026 | No realizar filtrado ribosomal | Se encontraron 86 genes ribosomales, de los cuales los porcentajes de genes ribosomales se encontraron con una mediana de 21.41% y una media de 19.51%. Se probó un MAD adaptativo global, pero no se encontró ninguna célula con un porcentaje ribosomal atípicamente alto (+3 desviaciones estándar), por lo tanto no se pudo filtrar en base a ello. | MAD ribosomal. Se investigó y se vio que rangos entre 15-40% son normales. |
| 9-octubre-2026 | Filtrado por dobletes con scDblFinder() | Es un paquete robusto, preciso y eficiente para eliminar dobletes. Se eligió un umbral con score de 0.37 de forma automática, lo cual es un umbral razonable que detectó 252 dobletes y minimiza errores de clasificación. | Nada |

Al final, de las 4085 células que se habían conservado de EmptyDrops. Quedamos con 2784 células tras el filtrado mitocondrial y de dobletes. 
