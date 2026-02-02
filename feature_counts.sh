#!/usr/bin/env bash
set -euo pipefail

############################
# CONFIG
############################
THREADS=12

BAM_DIR="/home/linoma/Documentos/Nina/STAR_out"
GTF="/home/linoma/Documentos/Nina/GCF_000001635.27/genomic.gtf"

OUT_DIR="/home/linoma/Documentos/Nina/counts_featureCounts_GeneID"
GTF_GENEID="/home/linoma/Documentos/Nina/GCF_000001635.27/genomic.geneid_as_gene_id.gtf"

############################
# PRE-CHECKS
############################
command -v featureCounts >/dev/null 2>&1 || { echo "ERRO: featureCounts não encontrado. Instale: conda install -c bioconda subread"; exit 1; }
[[ -d "$BAM_DIR" ]] || { echo "ERRO: BAM_DIR não existe: $BAM_DIR"; exit 1; }
[[ -f "$GTF" ]] || { echo "ERRO: GTF não encontrado: $GTF"; exit 1; }

mkdir -p "$OUT_DIR"

############################
# 1) Criar GTF com gene_id = GeneID
############################
echo "==> Criando GTF com gene_id = GeneID..."
awk -F'\t' 'BEGIN{OFS="\t"}
{
  attr=$9
  if (match(attr, /GeneID:[0-9]+/)) {
    gid=substr(attr, RSTART, RLENGTH)
    sub(/^GeneID:/,"",gid)

    if (attr ~ /gene_id "[^"]+"/) {
      sub(/gene_id "[^"]+"/, "gene_id \"" gid "\"", attr)
    } else {
      attr=attr "; gene_id \"" gid "\""
    }
    $9=attr
  }
  print
}' "$GTF" > "$GTF_GENEID"

echo "GTF_GENEID: $GTF_GENEID"

############################
# 2) Rodar featureCounts (todos os BAMs)
############################
echo "==> Rodando featureCounts em todos os BAMs..."
cd "$BAM_DIR"

# lista BAMs (garante ordem fixa)
ls -1 *_Aligned.sortedByCoord.out.bam > "$OUT_DIR/bams.list"
echo "Nº de BAMs: $(wc -l < "$OUT_DIR/bams.list")"

featureCounts \
  -T "$THREADS" \
  -a "$GTF_GENEID" \
  -o "$OUT_DIR/featureCounts_counts.txt" \
  -g gene_id \
  -t exon \
  -p -B -C \
  $(cat "$OUT_DIR/bams.list")

echo "OK: $OUT_DIR/featureCounts_counts.txt"
echo "Summary: $OUT_DIR/featureCounts_counts.txt.summary"

############################
# 3) Matriz final limpa (GeneID + counts)
############################
echo "==> Gerando counts_matrix_GeneID.tsv..."
cd "$OUT_DIR"

awk -F'\t' 'BEGIN{OFS="\t"}
NR==2{
  printf $1
  for(i=7;i<=NF;i++){
    gsub(/^.*\//,"",$i)
    sub(/_Aligned\.sortedByCoord\.out\.bam$/,"",$i)
    printf OFS $i
  }
  printf "\n"
  next
}
NR>2{
  printf $1
  for(i=7;i<=NF;i++) printf OFS $i
  printf "\n"
}' featureCounts_counts.txt > counts_matrix_GeneID.tsv

echo "OK: $OUT_DIR/counts_matrix_GeneID.tsv"
echo "Colunas (NF) únicas:"
awk -F'\t' '{print NF}' counts_matrix_GeneID.tsv | sort -nu | head

############################
# 4) Mapa GeneID -> Symbol/Biotype (limpo)
############################
echo "==> Gerando geneid_to_symbol.clean.tsv..."
OUTMAP="$OUT_DIR/geneid_to_symbol.clean.tsv"

awk -F'\t' 'BEGIN{OFS="\t"; print "GeneID","Symbol","Biotype"}
$3=="gene"{
  gid=""; sym=""; bio=""
  if (match($9, /gene_id "[^"]+"/)) { gid=substr($9,RSTART+9,RLENGTH-10) }
  if (match($9, /gene "[^"]+"/))    { sym=substr($9,RSTART+6,RLENGTH-7) }
  if (match($9, /gene_biotype "[^"]+"/)) { bio=substr($9,RSTART+13,RLENGTH-14) }

  gsub(/"/,"",sym)
  gsub(/"/,"",bio)

  if (gid!="") print gid, sym, bio
}' "$GTF_GENEID" > "$OUTMAP"

echo "OK: $OUTMAP"

echo "==> Finalizado!"
echo "Counts (GeneID): $OUT_DIR/counts_matrix_GeneID.tsv"
echo "Mapa GeneID->Symbol: $OUTMAP"

