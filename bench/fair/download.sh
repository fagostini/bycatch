#!/bin/bash
# alternative assemblies for simulation (never the indexed ones) + NCBI taxdump
set -euo pipefail
cd /var/tmp/bench_fair/genomes_alt
while read -r name acc; do
  [ -s $name.fa ] && continue
  curl -sS -f --retry 5 -m 3600 -o $name.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA"
  unzip -p $name.zip 'ncbi_dataset/data/*/*.fna' > $name.fa
  rm $name.zip
  echo "$(date +%T) $name $acc $(grep -c '>' $name.fa) seqs $(du -h $name.fa | cut -f1)"
done <<'L'
cod_alt GCA_010882105.1
pacific_cod GCF_031168955.1
human_alt GCF_009914755.1
mouse_alt GCA_921999005.2
lsinapis_alt GCA_949711785.1
drosophila_alt GCA_002300595.1
arabidopsis_alt GCA_965638225.1
yeast_alt GCA_000181435.1
mtb_alt GCF_000008585.1
ecoli_alt GCF_000008865.2
L
cd ..; mkdir -p taxdump; curl -sS -f --retry 5 -o taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz -C taxdump nodes.dmp names.dmp; rm taxdump.tar.gz
echo "$(date +%T) download done"
