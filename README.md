# Receta para el análisis de Metagenómica

Tutorial de Metagenómica Shotgun (WGS): del read crudo a los genomas ensamblados a partir de metagenomas (MAGs)

## Objetivos:

El siguiente tutorial tiene como objetivo introducir al alumno en el análisis bioinformático de
experimentos de secuenciación metagenómica (shotgun / WGS) provenientes de distintos ambientes.

A diferencia del práctico de RNA-Seq, aquí **no** tenemos un organismo de referencia único: cada muestra es
una comunidad microbiana completa, por lo que el desafío es **reconstruir quién está** (composición
taxonómica) y **qué puede hacer** (potencial funcional) esa comunidad.

Se realizarán las siguientes tareas:

- Control de Calidad de los reads con [__FastQC__](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
- Filtrado y recorte de reads paired-end con [__fastp__](https://github.com/OpenGene/fastp).
- (Opcional, muestras de hospedero humano) Remoción de reads del hospedero con [__Bowtie2__](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) / [__Hostile__](https://github.com/bede/hostile) usando la referencia humana completa __T2T-CHM13__.
- Perfil taxonómico basado en reads con [__Kraken2__](https://github.com/DerrickWood/kraken2) y re-estimación de abundancias con [__Bracken__](https://github.com/jenniferlu717/Bracken).
- Visualización interactiva de la composición con [__Krona__](https://github.com/marbl/Krona/wiki).
- Ensamblado _de novo_ del metagenoma con [__MEGAHIT__](https://github.com/voutcn/megahit).
- Evaluación del ensamblado con [__metaQUAST__](https://quast.sourceforge.net/metaquast) y [__seqkit__](https://bioinf.shenwei.me/seqkit/).
- Predicción de genes con [__Prodigal__](https://github.com/hyattpd/Prodigal) (modo `meta`).
- Reconstrucción de genomas (_binning_) con [__MetaBAT2__](https://bitbucket.org/berkeleylab/metabat/src/master/) y evaluación de la calidad de los MAGs con [__CheckM__](https://github.com/Ecogenomics/CheckM).
- (Opcional) Asignación taxonómica de los MAGs con [__GTDB-Tk__](https://github.com/Ecogenomics/GTDBTk).


## Materiales:

Una excelente fuente de datos metagenómicos públicos es el [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra)
de NCBI y su espejo europeo, el [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home).
Un buen catálogo de estudios y muestras curadas es [MGnify](https://www.ebi.ac.uk/metagenomics/) de EBI.

En este práctico cada grupo analizará **una muestra metagenómica de un ambiente distinto**. Todas son
Illumina, paired-end (`R1` + `R2`) y de tipo shotgun (WGS / `METAGENOMIC`). Cada grupo trabajará de forma
independiente con su propia muestra.

Ustedes utilizarán los siguientes datos:

|Grupo| Ambiente | Accesión (Run) | Link |
| ----- | ----- | ----- | ----- |
|1| Intestino humano | `SRR39501602` | [ENA](https://www.ebi.ac.uk/ena/browser/view/SRR39501602) |
|2| Cavidad oral humana | `ERR17497937` | [ENA](https://www.ebi.ac.uk/ena/browser/view/ERR17497937) |
|3| Piel humana | `SRR37938837` | [ENA](https://www.ebi.ac.uk/ena/browser/view/SRR37938837) |
|4| Marino (agua de mar) | `ERR16838384` | [ENA](https://www.ebi.ac.uk/ena/browser/view/ERR16838384) |
|5| Agua dulce | `ERR17336719` | [ENA](https://www.ebi.ac.uk/ena/browser/view/ERR17336719) |
|6| Aguas residuales | `SRR39106523` | [ENA](https://www.ebi.ac.uk/ena/browser/view/SRR39106523) |
|7| Fuente termal | `ERR16138956` | [ENA](https://www.ebi.ac.uk/ena/browser/view/ERR16138956) |
|8| Sedimento | `DRR1069364` | [ENA](https://www.ebi.ac.uk/ena/browser/view/DRR1069364) |
|9| Rizosfera | `SRR38834441` | [ENA](https://www.ebi.ac.uk/ena/browser/view/SRR38834441) |
|10| Compost / suelo | `ERR17020363` | [ENA](https://www.ebi.ac.uk/ena/browser/view/ERR17020363) |

Los datos ya se encuentran en el servidor y no es necesario descargarlos nuevamente. Se ponen los links a
disposición para que puedan identificar de dónde se obtienen los datos en las bases de datos internacionales
y para que puedan descargar metadata y el artículo asociado a cada investigación.

Si quisieran descargarlos por su cuenta, se usa la suite [sra-tools](https://github.com/ncbi/sra-tools) de NCBI,
específicamente [prefetch](https://github.com/ncbi/sra-tools/wiki/HowTo:-Download-and-verify-data) y
[fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump):

```bash
prefetch <ACCESION> --output-directory .
fasterq-dump <ACCESION> --split-files --threads 4 --outdir .
gzip <ACCESION>_1.fastq <ACCESION>_2.fastq
```

### Ubicación de los reads en el servidor:

Los reads de cada grupo están alojados en el servidor, en la carpeta de su ambiente:

```
/mnt/biostore/dipBG/metagenomics/
├── 01_intestino_humano/     SRR39501602_1.fastq.gz  SRR39501602_2.fastq.gz
├── 02_oral_humano/          ERR17497937_1.fastq.gz  ERR17497937_2.fastq.gz
├── 03_piel_humana/          SRR37938837_1.fastq.gz  SRR37938837_2.fastq.gz
├── 04_marino/               ERR16838384_1.fastq.gz  ERR16838384_2.fastq.gz
├── 05_agua_dulce/           ERR17336719_1.fastq.gz  ERR17336719_2.fastq.gz
├── 06_aguas_residuales/     SRR39106523_1.fastq.gz  SRR39106523_2.fastq.gz
├── 07_fuente_termal/        ERR16138956_1.fastq.gz  ERR16138956_2.fastq.gz
├── 08_sedimento/            DRR1069364_1.fastq.gz   DRR1069364_2.fastq.gz
├── 09_rizosfera/            SRR38834441_1.fastq.gz  SRR38834441_2.fastq.gz
└── 10_compost_suelo/        ERR17020363_1.fastq.gz  ERR17020363_2.fastq.gz
```

Convención de nombres: `<ACCESION>_1.fastq.gz` es el read forward (`R1`) y `<ACCESION>_2.fastq.gz` el reverse (`R2`).

Cada grupo trabajará en su carpeta de home creando una carpeta de trabajo, por ejemplo:

```bash
mkdir -p ~/METAGENOMICS/reads
cp /mnt/biostore/dipBG/metagenomics/<SU_CARPETA>/*.fastq.gz ~/METAGENOMICS/reads/
cd ~/METAGENOMICS/reads
```

### Bases de datos de referencia:

A diferencia de RNA-Seq, en metagenómica no usamos un genoma único, sino **bases de datos** que agrupan miles
de genomas. Ya están descargadas en el servidor (¡ocupan mucho espacio y demoran en construirse!):

| Base de datos | Uso | Ubicación en el servidor |
| ---- | ---- | ---- |
| Kraken2 DB (Standard/PlusPF) | Clasificación taxonómica de reads | `/mnt/biostore/dipBG/kraken2_db/` |
| Genoma humano **T2T-CHM13** (`chm13v2.0`, con chrY; índice Bowtie2) | Remoción de hospedero (muestras humanas) | `/mnt/biostore/dipBG/HostRef/` |
| GTDB (release para GTDB-Tk) | Taxonomía de MAGs (opcional) | `/mnt/biostore/dipBG/gtdbtk_db/` |

> Nota para el/la docente: ajuste estas rutas a las de su servidor. Si no dispone de la Kraken2 DB Standard
> (~50-100 GB), puede usar la [MiniKraken / k2_standard_08gb](https://benlangmead.github.io/aws-indexes/k2).

### Software:

Todo el software está disponible vía [conda](https://docs.conda.io/) / [bioconda](https://bioconda.github.io/).

- Control de calidad: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
- Filtrado/recorte: [fastp](https://github.com/OpenGene/fastp).
- Remoción de hospedero: [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) o [Hostile](https://github.com/bede/hostile) (envoltorio que usa Bowtie2/minimap2 con referencia curada) + [SAMtools](https://www.htslib.org/).
- Perfil taxonómico: [Kraken2](https://github.com/DerrickWood/kraken2) + [Bracken](https://github.com/jenniferlu717/Bracken) + [Krona](https://github.com/marbl/Krona/wiki).
- Ensamblado: [MEGAHIT](https://github.com/voutcn/megahit) (alternativa: [metaSPAdes](https://github.com/ablab/spades)).
- Evaluación del ensamblado: [metaQUAST](https://quast.sourceforge.net/metaquast) + [seqkit](https://bioinf.shenwei.me/seqkit/).
- Predicción de genes: [Prodigal](https://github.com/hyattpd/Prodigal).
- Binning y MAGs: [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) + [CheckM](https://github.com/Ecogenomics/CheckM) + [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) (opcional).

### Activar el ambiente

El ambiente del curso ya está creado en el servidor con **todo el software** de la lista anterior. Actívelo al
inicio de cada sesión (y déjelo activo durante todo el práctico):

```bash
conda activate metagenomica
```

Sabrá que está activo porque el prompt cambia a `(metagenomica) usuario@servidor:$`. Puede comprobar que las
herramientas responden, por ejemplo:

```bash
fastqc --version
kraken2 --version
```

> **GTDB-Tk (paso 8, opcional)** vive en un **ambiente aparte** llamado `gtdbtk` (sus dependencias son pesadas y
> se aíslan para no interferir con el resto). Solo al momento de correr GTDB-Tk cambie de ambiente:
>
> ```bash
> conda deactivate          # sale de metagenomica
> conda activate gtdbtk     # entra al ambiente de GTDB-Tk
> ```
>
> Al terminar ese paso, vuelva con `conda deactivate && conda activate metagenomica`.

Para salir del ambiente al final de la sesión:

```bash
conda deactivate
```


# Inicio del Práctico

## 1. Quality Check de los Reads

Partiremos con FastQC para revisar nuestros reads: los tamaños, las calidades y el contenido de adaptadores.
Como son datos paired-end, revisamos `R1` y `R2`.

```bash
fastqc --threads 6 *.fastq.gz
```

Nos traemos los resultados a nuestro computador para evaluarlos:

```bash
scp <usuario>@<servidor>:METAGENOMICS/reads/*_fastqc.html .
```

> Reemplace `<usuario>` por su usuario/grupo y `<servidor>` por la IP o nombre del servidor del curso.

Revise, para cada muestra: el gráfico **Per base sequence quality** (dónde cae la calidad), el
**Adapter Content** y el **Per base sequence content**. En metagenómica **no** esperamos duplicados bajos ni
GC uniforme (¡es una mezcla de muchos organismos!), así que algunos módulos aparecerán en rojo y eso es normal.

Recordemos la tabla de calidad [Phred](https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf):

| Valor Phred | Probabilidad de base errónea | Precisión |
| ----- | ---- | ---- |
| 10 | 1 en 10 | 90% |
| 20 | 1 en 100 | 99% |
| 30 | 1 en 1000 | 99.9% |
| 40 | 1 en 10000 | 99.99% |

## 2. Filtrado y recorte de los reads (fastp)

Usaremos **fastp**, que en un solo comando recorta adaptadores, filtra por calidad y largo, y genera un
reporte `HTML` + `JSON`. Al ser paired-end, entrega los dos archivos (`R1` y `R2`) a la vez para mantener el
pareo de las lecturas.

Para una muestra:

```bash
fastp \
  -i <ACCESION>_1.fastq.gz  -I <ACCESION>_2.fastq.gz \
  -o <ACCESION>_1.clean.fastq.gz -O <ACCESION>_2.clean.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 6 \
  --html <ACCESION>_fastp.html --json <ACCESION>_fastp.json
```

| Parámetro | Descripción |
| ---- | ---- |
| `-i` / `-I` | archivos de entrada `R1` / `R2` |
| `-o` / `-O` | archivos de salida filtrados `R1` / `R2` |
| `--detect_adapter_for_pe` | detecta automáticamente los adaptadores en modo paired-end |
| `--qualified_quality_phred` | calidad Phred mínima aceptada por base (aquí `20` = 99% de precisión) |
| `--length_required` | descarta reads más cortos que este largo (bp) después del recorte |
| `--thread` | número de hebras |

Como este trabajo es repetitivo (10 muestras si lo hiciéramos para todos los grupos), en lugar de correr una
por una podemos usar un pequeño bucle sobre todos los pares de archivos de la carpeta:

```bash
for r1 in *_1.fastq.gz; do
    acc=${r1%_1.fastq.gz}
    fastp -i ${acc}_1.fastq.gz -I ${acc}_2.fastq.gz \
          -o ${acc}_1.clean.fastq.gz -O ${acc}_2.clean.fastq.gz \
          --detect_adapter_for_pe --qualified_quality_phred 20 \
          --length_required 50 --thread 6 \
          --html ${acc}_fastp.html --json ${acc}_fastp.json
done
```

> Este bucle es el equivalente a los _wrappers_ (`run_trim.pl`, etc.) del práctico de RNA-Seq: automatiza la
> tarea repetitiva. Tráigase los reportes con `scp` y revise cuántos reads sobrevivieron el filtrado.

## 3. (Opcional) Remoción de reads del hospedero — solo muestras humanas

Las muestras de **intestino, oral y piel** (grupos 1, 2 y 3) contienen ADN humano. Para un análisis
microbiano limpio conviene eliminar esos reads mapeándolos contra el genoma humano y **quedándonos con los que
NO mapean**. Los grupos de ambientes (marino, suelo, etc.) pueden **saltarse este paso**.

### ¿Qué referencia humana usar? T2T-CHM13, no GRCh38

Aunque en el práctico de RNA-Seq usamos GRCh38 (porque ahí importa la **anotación de genes** para contar), para
**descontaminar** conviene la referencia **más completa posible**. La razón es simple: _todo read humano que no
esté en la referencia no alinea y se "cuela" como falso read microbiano_. GRCh38 todavía tiene gaps en
centrómeros, ADN satélite y duplicaciones segmentales; el [__T2T-CHM13__](https://www.science.org/doi/10.1126/science.abj6987)
(consorcio Telomere-to-Telomere) es el primer genoma humano **completo, sin gaps**, y por eso remueve más ADN
humano y deja un set microbiano más limpio.

> :warning: **Cromosoma Y:** la línea CHM13 original es 46,XX y **no incluye chrY**. Como las muestras humanas
> pueden provenir de hombres, use la versión que **incorpora el chrY de HG002** (`chm13v2.0`), o agregue ese
> chrY a la referencia; de lo contrario, los reads del cromosoma Y se escaparían sin removerse.

En el servidor la referencia `chm13v2.0` y su índice de Bowtie2 ya están construidos.

# Importante: el índice de Bowtie2 de T2T-CHM13 ya está construido en el servidor (demora mucho). Esta sección es para que sepan cómo se hace.

Construcción del índice (NO ejecutar, solo referencia):

```bash
bowtie2-build /mnt/biostore/dipBG/HostRef/chm13v2.0.fa /mnt/biostore/dipBG/HostRef/chm13v2
```

### Opción A — Bowtie2 (didáctica: se ve el concepto)

Remoción del hospedero (los grupos 1-3 sí ejecutan esto):

```bash
bowtie2 -x /mnt/biostore/dipBG/HostRef/chm13v2 \
        -1 <ACCESION>_1.clean.fastq.gz -2 <ACCESION>_2.clean.fastq.gz \
        --un-conc-gz <ACCESION>_nohost.%.fastq.gz \
        -p 8 -S /dev/null
```

`--un-conc-gz` guarda los pares que **no** alinearon al humano en `<ACCESION>_nohost.1.fastq.gz` y
`<ACCESION>_nohost.2.fastq.gz`; esos son los reads microbianos que usaremos de aquí en adelante. (El
`%` en el nombre es reemplazado por bowtie2 con `1` y `2`; si se omite, los archivos salen como
`<ACCESION>_nohost.fastq.1.gz`, que no calza con los pasos siguientes.)

| Parámetro | Descripción |
| ---- | ---- |
| `-x` | prefijo del índice del genoma del hospedero (T2T-CHM13) |
| `-1` / `-2` | reads limpios `R1` / `R2` |
| `--un-conc-gz` | escribe los pares que NO alinearon (lo que queremos conservar) |
| `-S /dev/null` | descartamos el SAM de los reads que sí alinearon (no nos interesan) |

### Opción B — Hostile (recomendada)

[Hostile](https://github.com/bede/hostile) es una herramienta dedicada a la descontaminación de hospedero. Por
debajo usa **Bowtie2** (reads cortos) o **minimap2** (largos), pero con una **referencia humana ya curada y
enmascarada** en las regiones homólogas a microbios —para no borrar por error reads microbianos reales— y que
ya incluye el chrY. Resuelve en una sola línea lo que en la Opción A hacemos a mano:

```bash
hostile clean \
        --fastq1 <ACCESION>_1.clean.fastq.gz \
        --fastq2 <ACCESION>_2.clean.fastq.gz \
        --index human-t2t-hla \
        --threads 8 \
        --out-dir <ACCESION>_nohost
```

La primera vez descarga el índice `human-t2t-hla` (basado en T2T-CHM13). Entrega los FASTQ limpios y un
reporte con cuántos reads se removieron.

> **Completitud vs. especificidad:** T2T-CHM13 aporta *completitud* (atrapa más ADN humano); el *enmascarado*
> de Hostile aporta *especificidad* (no elimina microbios que se parecen a humano). Son dos mejoras
> complementarias, y por eso Hostile combina ambas.

> De aquí en adelante, si su grupo hizo remoción de hospedero use los archivos limpios (`*_nohost.1/2.fastq.gz`
> de Bowtie2, o los de la carpeta de salida de Hostile); si no, use los `*_1.clean.fastq.gz` /
> `*_2.clean.fastq.gz`.

## 4. Perfil taxonómico basado en reads (Kraken2 + Bracken)

Aquí respondemos **¿quiénes están en mi muestra?** sin necesidad de ensamblar. Kraken2 asigna cada read a un
taxón comparando k-mers contra la base de datos.

Es recomendable usar `screen` porque Kraken2 carga la base de datos completa en memoria:

```bash
screen -S kraken
```

> :information_source: **Grupos 1-3 (muestras humanas):** de aquí en adelante reemplacen `<ACCESION>_1.clean.fastq.gz` /
> `<ACCESION>_2.clean.fastq.gz` por sus archivos sin hospedero (`*_nohost.1/2.fastq.gz` de Bowtie2, o los de la
> carpeta de salida de Hostile). Los demás grupos siguen usando los `*.clean.fastq.gz`.

```bash
kraken2 --db /mnt/biostore/dipBG/kraken2_db \
        --paired <ACCESION>_1.clean.fastq.gz <ACCESION>_2.clean.fastq.gz \
        --threads 8 \
        --report <ACCESION>.kreport \
        --output <ACCESION>.kraken
```

| Parámetro | Descripción |
| ---- | ---- |
| `--db` | ruta a la base de datos de Kraken2 |
| `--paired` | indica que la entrada son dos archivos paired-end |
| `--report` | archivo resumen con el % de reads por taxón (el que más usaremos) |
| `--output` | asignación read-por-read (archivo grande, se puede omitir con `--output -`) |

Kraken2 tiende a **sobre-estimar** la cantidad de especies. **Bracken** re-estima las abundancias reales a un
nivel taxonómico dado (ej. especie `S` o género `G`):

```bash
bracken -d /mnt/biostore/dipBG/kraken2_db \
        -i <ACCESION>.kreport \
        -o <ACCESION>.bracken \
        -r 150 -l S -t 10
```

| Parámetro | Descripción |
| ---- | ---- |
| `-r` | largo de los reads (150 para 2×150) |
| `-l` | nivel taxonómico: `S` especie, `G` género, `P` phylum |
| `-t` | umbral mínimo de reads para considerar un taxón |

### Visualización con Krona

```bash
ktImportTaxonomy -q 2 -t 3 <ACCESION>.kraken -o <ACCESION>.krona.html
```

Ábralo en su navegador: es un gráfico circular interactivo de la composición taxonómica.

### En su informe deberá mostrar el gráfico de Krona y describir los phyla/géneros dominantes de su ambiente, e indicar si tienen sentido biológico (por ejemplo, cianobacterias en muestras marinas, o _Bacteroides_/_Firmicutes_ en intestino).

## 5. Ensamblado de novo del metagenoma (MEGAHIT)

Ahora reconstruimos fragmentos largos de genoma (_contigs_) uniendo los reads. Usaremos **MEGAHIT**, que es
rápido y eficiente en memoria, ideal para metagenomas.

```bash
screen -S assembly
```

```bash
megahit -1 <ACCESION>_1.clean.fastq.gz -2 <ACCESION>_2.clean.fastq.gz \
        -t 12 \
        --min-contig-len 1000 \
        -o <ACCESION>_megahit
```

| Parámetro | Descripción |
| ---- | ---- |
| `-1` / `-2` | reads limpios `R1` / `R2` |
| `-t` | número de hebras |
| `--min-contig-len` | descarta contigs más cortos que este largo (1000 bp es un buen mínimo) |
| `-o` | carpeta de salida (¡debe NO existir previamente!) |

El ensamblado final queda en `<ACCESION>_megahit/final.contigs.fa`.

> Alternativa (más lenta pero a veces más contigua): **metaSPAdes**
> `spades.py --meta -1 R1 -2 R2 -t 12 -o <ACCESION>_spades`

## 6. Evaluación del ensamblado (metaQUAST / seqkit)

Estadísticas rápidas de los contigs:

```bash
seqkit stats -a <ACCESION>_megahit/final.contigs.fa
```

Reporte más completo (N50, largo total, contig más largo, etc.):

```bash
metaquast.py <ACCESION>_megahit/final.contigs.fa -o <ACCESION>_quast --max-ref-number 0 -t 8
```

### En su informe deberá reportar, para su ensamblado: número de contigs, largo total (bp), N50 y el contig más largo. Discuta qué tan bien ensambló su muestra y por qué (pista: complejidad/diversidad de la comunidad y profundidad de secuenciación).

## 7. Predicción de genes (Prodigal)

Sobre los contigs predecimos los genes codificantes de proteínas en modo metagenómico (`-p meta`):

```bash
prodigal -i <ACCESION>_megahit/final.contigs.fa \
         -a <ACCESION>_proteins.faa \
         -d <ACCESION>_genes.fna \
         -o <ACCESION>_genes.gff -f gff \
         -p meta
```

| Parámetro | Descripción |
| ---- | ---- |
| `-a` | proteínas predichas (aminoácidos) |
| `-d` | genes predichos (nucleótidos) |
| `-o` / `-f gff` | coordenadas de los genes en formato GFF |
| `-p meta` | modo metagenómico (muestras con múltiples organismos) |

Cuente cuántos genes predijo:

```bash
grep -c ">" <ACCESION>_proteins.faa
```

## 8. Reconstrucción de genomas (binning) → MAGs

El _binning_ agrupa los contigs que probablemente pertenecen al **mismo organismo** usando su composición y su
cobertura. El producto son los **MAGs** (Metagenome-Assembled Genomes). Este es el análogo metagenómico de
"recuperar genomas completos" a partir de la mezcla.

Primero necesitamos la **cobertura** de cada contig, mapeando los reads de vuelta al ensamblado:

```bash
# índice del ensamblado
bowtie2-build <ACCESION>_megahit/final.contigs.fa <ACCESION>_contigs
# mapeo de los reads limpios a los contigs
bowtie2 -x <ACCESION>_contigs -1 <ACCESION>_1.clean.fastq.gz -2 <ACCESION>_2.clean.fastq.gz -p 8 | \
    samtools sort -@ 8 -o <ACCESION>.bam -
samtools index <ACCESION>.bam
```

Calculamos la profundidad por contig y corremos MetaBAT2:

```bash
jgi_summarize_bam_contig_depths --outputDepth <ACCESION>_depth.txt <ACCESION>.bam
metabat2 -i <ACCESION>_megahit/final.contigs.fa \
         -a <ACCESION>_depth.txt \
         -o <ACCESION>_bins/bin \
         -m 1500 -t 8
```

Cada archivo `<ACCESION>_bins/bin.*.fa` es un MAG candidato.

### Evaluación de la calidad de los MAGs (CheckM)

CheckM estima la **completitud** y **contaminación** de cada bin usando genes marcadores de copia única:

```bash
checkm lineage_wf -x fa -t 8 <ACCESION>_bins <ACCESION>_checkm
```

Un MAG se considera de **alta calidad** si tiene completitud > 90% y contaminación < 5%
(criterios [MIMAG](https://www.nature.com/articles/nbt.3893)).

### En su informe deberá reportar cuántos MAGs recuperó y una tabla con completitud y contaminación de cada uno, indicando cuáles cumplen criterios de calidad media/alta.

### (Opcional) Taxonomía de los MAGs con GTDB-Tk

```bash
gtdbtk classify_wf --genome_dir <ACCESION>_bins -x fa --out_dir <ACCESION>_gtdbtk --cpus 8
```

Esto le dará el nombre (dominio → especie) de cada genoma reconstruido según la taxonomía GTDB.

---

# Resumen del flujo

```
reads crudos (R1/R2)
   │  FastQC
   ▼
fastp  ──►  (Bowtie2: remoción de hospedero, solo muestras humanas)
   │
   ├──► Kraken2 + Bracken + Krona  ──►  ¿QUIÉNES están?  (perfil taxonómico)
   │
   └──► MEGAHIT (ensamblado)
             │  metaQUAST / seqkit (calidad del ensamblado)
             ├──► Prodigal  ──►  ¿QUÉ pueden hacer?  (genes)
             └──► MetaBAT2 + CheckM (+ GTDB-Tk)  ──►  MAGs (genomas reconstruidos)
```

# Entrega del informe

Su informe debe incluir, para su ambiente asignado:

1. **Calidad**: resumen de FastQC/fastp (reads antes/después del filtrado).
2. **Composición taxonómica**: gráfico de Krona + descripción de los taxones dominantes.
3. **Ensamblado**: tabla con nº de contigs, largo total, N50 y contig más largo.
4. **Genes**: número de genes predichos por Prodigal.
5. **MAGs**: número de genomas recuperados y su completitud/contaminación (CheckM); taxonomía si usó GTDB-Tk.
6. **Discusión**: ¿los resultados tienen sentido biológico para su ecosistema? Contraste con literatura.
