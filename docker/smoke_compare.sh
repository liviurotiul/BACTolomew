#!/usr/bin/env bash
# Run the pipeline's heavier tools on a tiny input set in two tool images and compare the
# outputs. Use it before publishing a new image version (docker/PUBLISHING.md step 4): the
# lock check only reads conda-meta, so it cannot see files removed after install (static
# libs, docs); real runs can.
#
# Usage: bash docker/smoke_compare.sh <image-a> <image-b> <inputs-dir> <work-dir>
#   inputs-dir holds:
#     gff/*.gff3  three or four Bakta GFF3 files (with their ##FASTA section), for Panaroo
#     fna/*.fna   the same genomes' contigs, for chewBBACA and Snippy
#     ref.fna     a Snippy reference
#     aln.fa      a small core alignment (about 10 sequences, a few hundred kb), for
#                 RAxML-NG, IQ-TREE and Gubbins
#   work-dir    an empty folder for the outputs (<work-dir>/<a|b>/<tool>)
# Every container runs with --network none, as the calling user, with the two folders
# mounted at their own paths. Nothing is pushed or changed in the images.
set -uo pipefail
IMAGE_A="${1:?usage: smoke_compare.sh <image-a> <image-b> <inputs-dir> <work-dir>}"
IMAGE_B="${2:?usage}"
IN="$(cd "${3:?usage}" && pwd)"
WORK="${4:?usage}"
mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
THREADS="${THREADS:-4}"
FIRST_FNA="$(ls "$IN"/fna/*.fna | head -1)"

run_in() {  # run_in <image> <out-dir> <env> <command...>
    local image="$1" out="$2" env="$3"; shift 3
    docker run --rm --network none --user "$(id -u):$(id -g)" \
        -e HOME="$out/.home" -v "$IN:$IN:ro" -v "$WORK:$WORK" -w "$out" \
        "$image" conda run --no-capture-output -n "$env" "$@"
}

run_all() {  # run_all <image> <tag>
    local image="$1" w="$WORK/$2"
    mkdir -p "$w/.home" "$w/raxml" "$w/iqtree" "$w/gubbins" "$w/chew"
    echo "[smoke] $image -> $w"
    run_in "$image" "$w" phylogeny snippy --outdir "$w/snippy" --ref "$IN/ref.fna" \
        --ctgs "$FIRST_FNA" --cpus "$THREADS" --force > "$w/snippy.log" 2>&1
    echo "[smoke]   snippy rc=$?"
    run_in "$image" "$w/raxml" phylogeny raxml-ng --search --msa "$IN/aln.fa" --model GTR+G \
        --seed 1 --threads 2 --tree pars{2} --prefix "$w/raxml/r" > "$w/raxml.log" 2>&1
    echo "[smoke]   raxml-ng rc=$?"
    run_in "$image" "$w/iqtree" gubbins iqtree2 -s "$IN/aln.fa" -m GTR+G -seed 1 -T 2 \
        --prefix "$w/iqtree/i" > "$w/iqtree.log" 2>&1
    echo "[smoke]   iqtree2 rc=$?"
    run_in "$image" "$w/gubbins" gubbins run_gubbins.py --prefix g --tree-builder raxml \
        --seed 1 --threads 2 --iterations 3 "$IN/aln.fa" > "$w/gubbins.log" 2>&1
    echo "[smoke]   gubbins rc=$?"
    run_in "$image" "$w" panaroo panaroo -i "$IN"/gff/*.gff3 -o "$w/panaroo" \
        --clean-mode strict -t "$THREADS" > "$w/panaroo.log" 2>&1
    echo "[smoke]   panaroo rc=$?"
    run_in "$image" "$w/chew" chewbbaca chewBBACA.py CreateSchema -i "$IN/fna" -o "$w/chew/schema" \
        --cpu "$THREADS" > "$w/chew_schema.log" 2>&1
    echo "[smoke]   chewBBACA CreateSchema rc=$?"
    run_in "$image" "$w/chew" chewbbaca chewBBACA.py AlleleCall -i "$IN/fna" \
        -g "$w/chew/schema/schema_seed" -o "$w/chew/allelecall" --cpu "$THREADS" > "$w/chew_call.log" 2>&1
    echo "[smoke]   chewBBACA AlleleCall rc=$?"
}

same() {  # same <label> <relative path>; sorted compare for row-order-free tables
    local a="$WORK/a/$2" b="$WORK/b/$2"
    if [ ! -s "$a" ] || [ ! -s "$b" ]; then
        echo "MISSING  $1 ($2)"
    elif cmp -s <(sort "$a") <(sort "$b"); then
        echo "SAME     $1"
    else
        echo "DIFFERS  $1 ($2)"
    fi
}

run_all "$IMAGE_A" a
run_all "$IMAGE_B" b

echo "[smoke] $IMAGE_A (a) vs $IMAGE_B (b)"
same "snippy snps.tab" snippy/snps.tab
same "raxml-ng best tree" raxml/r.raxml.bestTree
same "iqtree2 tree" iqtree/i.treefile
same "gubbins final tree" gubbins/g.final_tree.tre
same "gubbins recombination" gubbins/g.recombination_predictions.gff
same "panaroo presence/absence" panaroo/gene_presence_absence.Rtab
same "chewBBACA alleles" chew/allelecall/results_alleles.tsv
