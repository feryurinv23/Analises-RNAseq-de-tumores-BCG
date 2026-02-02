#!/usr/bin/env bash
set -euo pipefail

############################
# CONFIG (edite se quiser) #
############################
THREADS=12
READ_LENGTH=150                      # usado apenas para sjdbOverhang (= READ_LENGTH-1)
BASE_DIR="$HOME/Documentos/Nina"

FASTQ_DIR="$BASE_DIR/FastQ/WTNIBCG"
REF_DIR="$BASE_DIR/GCF_000001635.27"
FASTA="$REF_DIR/GCF_000001635.27_GRCm39_genomic.fna"
GTF="$REF_DIR/genomic.gtf"

GENOME_DIR="$BASE_DIR/STAR_index_GRCm39"
OUT_DIR="$BASE_DIR/STAR_out"

SJDB_OVERHANG=$((READ_LENGTH-1))

#########################
# CHECAGENS INICIAIS    #
#########################
echo "==> Checando arquivos e pastas..."

command -v STAR >/dev/null 2>&1 || { echo "ERRO: STAR não encontrado no PATH."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERRO: samtools não encontrado no PATH."; exit 1; }

[[ -d "$FASTQ_DIR" ]] || { echo "ERRO: FASTQ_DIR não existe: $FASTQ_DIR"; exit 1; }
[[ -f "$FASTA" ]] || { echo "ERRO: FASTA não encontrado: $FASTA"; exit 1; }
[[ -f "$GTF" ]] || { echo "ERRO: GTF não encontrado: $GTF"; exit 1; }

mkdir -p "$GENOME_DIR" "$OUT_DIR"

#########################
# 1) GERAR ÍNDICE STAR   #
#########################
# STAR cria arquivos como Genome, SA, SAindex etc.
if [[ -f "$GENOME_DIR/Genome" && -f "$GENOME_DIR/SA" ]]; then
  echo "==> Índice STAR já existe em: $GENOME_DIR (pulando genomeGenerate)"
else
  echo "==> Gerando índice STAR em: $GENOME_DIR"
  STAR --runThreadN "$THREADS" \
       --runMode genomeGenerate \
       --genomeDir "$GENOME_DIR" \
       --genomeFastaFiles "$FASTA" \
       --sjdbGTFfile "$GTF" \
       --sjdbOverhang "$SJDB_OVERHANG"
fi

#########################
# 2) ALIGNMENT          #
#########################
echo "==> Iniciando alinhamento dos FASTQs em: $FASTQ_DIR"
cd "$FASTQ_DIR"

shopt -s nullglob
R1_FILES=(*_R1_001.fastq.gz)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  echo "ERRO: Não encontrei arquivos *_R1_001.fastq.gz em $FASTQ_DIR"
  exit 1
fi

for r1 in "${R1_FILES[@]}"; do
  base="${r1%_R1_001.fastq.gz}"
  r2="${base}_R2_001.fastq.gz"

  if [[ ! -f "$r2" ]]; then
    echo "AVISO: Não encontrei o par R2 para $r1 (esperado: $r2). Pulando..."
    continue
  fi

  prefix="$OUT_DIR/${base}_"
  bam_out="${prefix}Aligned.sortedByCoord.out.bam"

  # Se já existe BAM final, pula para evitar reprocessar
  if [[ -f "$bam_out" ]]; then
    echo "==> Já existe BAM para $base (pulando): $bam_out"
    continue
  fi

  echo "==> Alinhando: $base"
  STAR --runThreadN "$THREADS" \
       --genomeDir "$GENOME_DIR" \
       --readFilesIn "$r1" "$r2" \
       --readFilesCommand zcat \
       --outFileNamePrefix "$prefix" \
       --outSAMtype BAM SortedByCoordinate \
       --quantMode GeneCounts
done

#########################
# 3) INDEXAR BAMs       #
#########################
echo "==> Indexando BAMs com samtools..."
cd "$OUT_DIR"
for bam in *_Aligned.sortedByCoord.out.bam; do
  [[ -f "$bam" ]] || continue
  if [[ -f "${bam}.bai" ]]; then
    echo "   OK (já indexado): $bam"
  else
    samtools index "$bam"
    echo "   Indexado: $bam"
  fi
done

#########################
# 4) RESUMO DE MAPEAMENTO
#########################
echo "==> Gerando resumo de mapeamento (STAR_Log.final.summary.txt)"
SUMMARY="$OUT_DIR/STAR_Log.final.summary.txt"
{
  echo -e "Sample\tUniquelyMapped(%)\tMultiMapped(%)\tTooManyLoci(%)\tUnmapped(%)"
  for log in *_Log.final.out; do
    [[ -f "$log" ]] || continue
    sample="${log%_Log.final.out}"

    uniq=$(grep -m1 "Uniquely mapped reads %" "$log" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
    multi=$(grep -m1 "% of reads mapped to multiple loci" "$log" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
    tooMany=$(grep -m1 "% of reads mapped to too many loci" "$log" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
    unmapped=$(grep -m1 "% of reads unmapped: too many mismatches" "$log" | awk -F'|' '{gsub(/ /,"",$2); a=$2}') || true
    # unmapped total: soma 3 linhas "unmapped" (mismatches + too short + other)
    u1=$(grep -m1 "% of reads unmapped: too many mismatches" "$log" | awk -F'|' '{gsub(/ /,"",$2); print $2}' | tr -d '%')
    u2=$(grep -m1 "% of reads unmapped: too short" "$log" | awk -F'|' '{gsub(/ /,"",$2); print $2}' | tr -d '%')
    u3=$(grep -m1 "% of reads unmapped: other" "$log" | awk -F'|' '{gsub(/ /,"",$2); print $2}' | tr -d '%')
    # se algum estiver vazio, vira 0
    u1=${u1:-0}; u2=${u2:-0}; u3=${u3:-0}
    unmapped_total=$(awk -v a="$u1" -v b="$u2" -v c="$u3" 'BEGIN{printf "%.2f", a+b+c}')

    echo -e "${sample}\t${uniq}\t${multi}\t${tooMany}\t${unmapped_total}%"
  done
} > "$SUMMARY"

echo "==> Finalizado!"
echo "Saídas em: $OUT_DIR"
echo "Resumo: $SUMMARY"

