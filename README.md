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
|Leonardo Daniel González López | Copiar plantilla de proyecto, Pregunta 1 y 3, y Control de Calidad |
|Yuliana Denisse Sosa Gómez | Leer los datos, Pregunta 2 y 4, y Control de calidad |

Esta tabla no es burocracia: el criterio B de la rúbrica es una respuesta
individual dirigida, y saber quién trabajó en qué parte es lo que la hace justa.

## Cómo reproducir este análisis

```bash
git clone https://github.com/leonardoglez6/proyecto-equipo-A.git
cd <el-repositorio>
quarto render reporte.qmd
```

**Los datos no están en el repositorio** y no deben estarlo: `.gitignore` los
excluye. Declara aquí de dónde salen.

| Dato | Dónde vive | Cómo se obtuvo |
|---|---|---|
| Dataset Equipo A | `/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-A.rds` | [Peripheral blood mononuclear cells (PBMCs) from a healthy donor - Manual (channel 5) from 10X Genomics](https://www.10xgenomics.com/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-manual-channel-5-3-1-standard-3-1-0) |

## Estructura

| Carpeta | Qué va aquí |
|---|---|
| `reporte.qmd` | El reporte ejecutable. Es el entregable |
| `R/` | Funciones que uses en más de un lugar |
| `datos/` | Vacía a propósito: los datos no se versionan |
| `figuras/` | Salida del reporte. Tampoco se versiona |

## Bitácora de decisiones

Cada vez que tomes una decisión de análisis —un umbral, un método, un filtro—
anótala aquí con su razón. **La tabla de decisiones es el criterio A3 de la
rúbrica** (15 %), y llenarla al final, de memoria, se nota.

| Fecha | Decisión | Por qué | Qué se probó antes |
|---|---|---|---|
| | | | |
