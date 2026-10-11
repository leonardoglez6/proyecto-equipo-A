# Software y referencia de 10x Genomics

Instalado el 25 de agosto de 2026 para el módulo de single-cell RNA-seq
(Bioinformática y Estadística 3, LCG-UNAM, semestre 2027-1).

## cellranger-10.1.0/

Cell Ranger 10.1.0, binario para Linux x86_64. 636 MB comprimido.

Descargado de https://www.10xgenomics.com/support/software/cell-ranger/downloads
El enlace de descarga es firmado y requiere cuenta en 10x Genomics más aceptar
la licencia; caduca en horas, así que no se puede reproducir con una URL fija.

Licencia: 10x Genomics End User License Agreement. Uso académico permitido.

## refdata-gex-GRCh38-2024-A/

Referencia transcriptómica humana GRCh38, versión 2024-A. 11.5 GB comprimida.

https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz

URL pública, sin registro. Construida por 10x Genomics a partir de GENCODE.

## Verificación

`cellranger testrun` y una corrida real de `cellranger count` sobre el conjunto
de la sesión 03. Los .tar originales se borraron tras extraer.

PENDIENTE: contrastar las sumas de comprobación que publica 10x. No se hizo al
descargar.

## Sumas de comprobación

cellranger-10.1.0.tar.xz
  md5 local: bb88407acb40cd9e4cd2749b6d774743
  calculado el 25 de agosto de 2026
  md5 publicado por 10x: bb88407acb40cd9e4cd2749b6d774743
  VERIFICADO: coinciden. Contrastado el 25 de agosto de 2026.

El .tar.xz se CONSERVA a proposito: el enlace de descarga es firmado y caduca,
asi que este archivo es la unica via de reinstalar sin volver a pedirlo a 10x.

La referencia GRCh38 y los FASTQ ya no se pueden contrastar contra las sumas
publicadas: sus .tar se borraron tras extraer. En su lugar hay un CHECKSUMS.md5
por carpeta, linea base local para detectar corrupcion futura.
  Verificar con:  md5sum -c CHECKSUMS.md5
