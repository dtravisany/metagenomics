#!/usr/bin/env bash
# =============================================================================
# Práctico de metagenómica — 10 muestras de AMBIENTES DISTINTOS
# Illumina, paired-end, shotgun (WGS / METAGENOMIC). Una accesión por grupo.
# Descarga en el symlink ./metagenomics -> /mnt/biostore/dipBG/metagenomics (TrueNAS)
# Requiere: sra-tools (prefetch, fasterq-dump)  ->  ya instalados en /usr/bin
# =============================================================================
set -euo pipefail

# Destino = symlink que está junto a este script (apunta al TrueNAS)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$SCRIPT_DIR/metagenomics"        # -> /mnt/biostore/dipBG/metagenomics/
THREADS="${THREADS:-4}"
TMPDIR_SRA="${TMPDIR_SRA:-$DEST/.tmp}" # temporales de fasterq-dump

# grupo|ambiente|accesión   (formato: etiqueta:accesión)
SAMPLES=(
  "01_intestino_humano:SRR39501602"   # Grupo 1  - human gut     - NovaSeq 6000 - ~191 MB
  "02_oral_humano:ERR17497937"        # Grupo 2  - human oral    - NovaSeq 6000 - ~148 MB
  "03_piel_humana:SRR37938837"        # Grupo 3  - human skin    - HiSeq 2500   - ~787 MB
  "04_marino:ERR16838384"             # Grupo 4  - marine        - NextSeq 2000 - ~376 MB
  "05_agua_dulce:ERR17336719"         # Grupo 5  - freshwater    - NovaSeq 6000 - ~335 MB
  "06_aguas_residuales:SRR39106523"   # Grupo 6  - wastewater    - NextSeq 500  - ~552 MB
  "07_fuente_termal:ERR16138956"      # Grupo 7  - hot spring    - NovaSeq 6000 - ~184 MB
  "08_sedimento:DRR1069364"           # Grupo 8  - sediment      - MiSeq        - ~583 MB
  "09_rizosfera:SRR38834441"          # Grupo 9  - rhizosphere   - MiSeq        - ~528 MB
  "10_compost_suelo:ERR17020363"      # Grupo 10 - compost/soil  - MiSeq        - ~121 MB
)

mkdir -p "$DEST" "$TMPDIR_SRA"
echo ">> Destino: $DEST"
echo ">> Hilos:   $THREADS"

for entry in "${SAMPLES[@]}"; do
  label="${entry%%:*}"
  acc="${entry##*:}"
  out="$DEST/$label"
  echo ""
  echo "=================================================================="
  echo ">> $label   ($acc)"
  echo "=================================================================="

  # Saltar si ya está descargada (ambos R1 y R2 comprimidos)
  if [[ -s "$out/${acc}_1.fastq.gz" && -s "$out/${acc}_2.fastq.gz" ]]; then
    echo "   ya existe, se omite."
    continue
  fi

  mkdir -p "$out"
  # 1) Descargar el .sra
  prefetch "$acc" --output-directory "$out" --max-size 5g
  # 2) Convertir a FASTQ paired (R1/R2)
  fasterq-dump "$acc" \
      --split-files \
      --threads "$THREADS" \
      --temp "$TMPDIR_SRA" \
      --outdir "$out"
  # 3) Comprimir
  gzip -f "$out/${acc}_1.fastq" "$out/${acc}_2.fastq"
  # 4) Limpiar el .sra descargado (opcional, ahorra espacio)
  rm -rf "$out/$acc"
done

rm -rf "$TMPDIR_SRA"
echo ""
echo ">> LISTO. Resumen:"
du -sh "$DEST"/*/ 2>/dev/null || true
