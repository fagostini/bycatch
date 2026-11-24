#!/bin/bash
# ART read pools per organism x read length x seed, from assemblies that are NOT the indexed ones
# (phiX excepted: only one sequence exists). Pools are oversampled; mix.py draws exact counts.
set -euo pipefail
W=/var/tmp/bench_fair; cd $W
export PATH=$W/env/bin:$PATH
A=$W/genomes_alt; P=$W/pools; mkdir -p $P
D=/home/agostif/Projects/bycatch/data

# off-panel genomes kept compressed in data/: ART needs plain FASTA
[ -s $A/dog.fa ] || zcat $D/dog.fa.gz > $A/dog.fa
[ -s $A/frog.fa ] || zcat $D/genomes/common_frog.fa.gz > $A/frog.fa
[ -s $A/phix.fa ] || cp $D/phix.fa $A/phix.fa

# organism  fasta  pairs needed (max over mixtures; see mix.py)
cat > $P/orgs.tsv <<L
cod	cod_alt.fa	1998000
human	human_alt.fa	3889
mouse	mouse_alt.fa	3889
lsinapis	lsinapis_alt.fa	3889
drosophila	drosophila_alt.fa	3889
arabidopsis	arabidopsis_alt.fa	3889
yeast	yeast_alt.fa	3889
mtb	mtb_alt.fa	3889
ecoli	ecoli_alt.fa	3889
phix	phix.fa	3889
pacific_cod	pacific_cod.fa	5000
dog	dog.fa	5000
frog	frog.fa	5000
L

job() {  # org fasta need L seed
  local org=$1 fa=$A/$2 need=$3 L=$4 seed=$5 out=$P/${1}_L${4}_s${5}
  [ -s ${out}.done ] && return 0
  local G=$(grep -v '>' $fa | tr -d '\n' | wc -c)
  # oversample: x1.2 for the host (rounding only), x3 for the rest (small contigs round to 0 reads)
  local over=3; [ $org = cod ] && over=1.2
  local f=$(python3 -c "print(f'{max($need*$over*2*$L/$G, 1e-6):.8f}')")
  art_illumina -ss HS25 -i $fa -p -l $L -f $f -m 250 -s 80 -rs $((seed*100000 + L*100 + ${#org})) -na -q -o ${out}_ > ${out}.log 2>&1
  local got=$(( $(wc -l < ${out}_1.fq) / 4 ))
  echo "$org L$L s$seed pairs=$got need=$need" > ${out}.done
}
export -f job; export A P

for s in 1 2 3; do for L in 50 100 150; do
  while IFS=$'\t' read -r org fa need; do echo "$org $fa $need $L $s"; done < $P/orgs.tsv
done; done | xargs -P 8 -L 1 bash -c 'job $0 $1 $2 $3 $4'
# throughput pool: 10 M host pairs, 2x150, seed 9 (mid mixture)
job cod cod_alt.fa 9900000 150 9
for o in human mouse lsinapis drosophila arabidopsis yeast mtb ecoli; do job $o ${o}_alt.fa 77778 150 9; done
job phix phix.fa 77778 150 9
for o in pacific_cod dog frog; do job $o $o.fa 50000 150 9; done
cat $P/*.done | sort > $P/pools_summary.txt
echo "$(date +%T) pools done"
