#!/usr/bin/env bash
# Run the pipeline's heavier tools on a tiny input set in two tool images and compare the
# outputs. Use it before publishing a new image version (docker/PUBLISHING.md step 4): the
# lock check only reads conda-meta, so it cannot see files removed after install (static
# libs, docs); real runs can.
#
# Usage: bash docker/smoke_compare.sh <image-a> <image-b> <inputs-dir> <work-dir>
#   runs both images here and compares (<work-dir>/a and <work-dir>/b);
#        bash docker/smoke_compare.sh --run <image> <inputs-dir> <out-dir>
#   runs one image into <out-dir> (for example the arm64 image on an Apple Silicon Mac);
#        bash docker/smoke_compare.sh --compare <out-dir-a> <out-dir-b>
#   compares two such folders made on different machines from the same inputs.
#   inputs-dir holds:
#     gff/*.gff3  three or four Bakta GFF3 files (with their ##FASTA section), for Panaroo
#     fna/*.fna   the same genomes' contigs, for chewBBACA and Snippy
#     ref.fna     a Snippy reference
#     aln.fa      a small core alignment (about 10 sequences, a few hundred kb), for
#                 RAxML-NG, IQ-TREE, Gubbins and FastBAPS (an image without the fastbaps
#                 env, 0.2.0 and older, reports it MISSING)
#   mlst and QUAST also run on the fna files (QUAST in the quast env, or in assembly_qc in
#   images before 0.3.0).
#   work-dir    an empty folder for the outputs (<work-dir>/<a|b>/<tool>)
# Every container runs with --network none, as the calling user, with the two folders
# mounted at their own paths. Nothing is pushed or changed in the images.
set -uo pipefail
USAGE="usage: smoke_compare.sh <image-a> <image-b> <inputs-dir> <work-dir> | --run <image> <inputs-dir> <out-dir> | --compare <dir-a> <dir-b>"
MODE=both
case "${1:-}" in
    --run) MODE=run; shift ;;
    --compare) MODE=compare; shift ;;
esac
THREADS="${THREADS:-4}"

run_in() {  # run_in <image> <out-dir> <env|quast> <command...>
    local image="$1" out="$2" env="$3"; shift 3
    docker run --rm --network none --user "$(id -u):$(id -g)" \
        -e HOME="$out/.home" -v "$IN:$IN:ro" -v "$WORK:$WORK" -w "$out" \
        "$image" sh -c 'e="$1"; shift; [ -d "/opt/conda/envs/$e" ] || [ "$e" != quast ] || e=assembly_qc; exec conda run --no-capture-output -n "$e" "$@"' \
        sh "$env" "$@"
}

run_all() {  # run_all <image> <out-dir>
    local image="$1" w="$2"
    mkdir -p "$w/.home" "$w/raxml" "$w/iqtree" "$w/gubbins" "$w/chew" "$w/fastbaps"
    docker run --rm --network none "$image" uname -m > "$w/arch.txt"
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
    run_in "$image" "$w" mlst mlst --threads "$THREADS" "$IN"/fna/*.fna 2> "$w/mlst.log" \
        | sed "s#$IN/fna/##" > "$w/mlst.tsv"
    echo "[smoke]   mlst rc=${PIPESTATUS[0]}"
    run_in "$image" "$w" quast quast.py -o "$w/quast" -t "$THREADS" "$IN"/fna/*.fna > "$w/quast.log" 2>&1
    echo "[smoke]   quast rc=$?"
    run_in "$image" "$w/fastbaps" fastbaps Rscript -e 'suppressPackageStartupMessages(library(fastbaps)); set.seed(1);
        s <- optimise_prior(import_fasta_sparse_nt(commandArgs(TRUE)[1]), type = "optimise.symmetric");
        write.table(multi_res_baps(s, levels = 2), commandArgs(TRUE)[2], sep = "\t", quote = FALSE, row.names = FALSE)' \
        "$IN/aln.fa" "$w/fastbaps/clusters.tsv" > "$w/fastbaps.log" 2>&1
    echo "[smoke]   fastbaps rc=$?"
}

same() {  # same <label> <relative path>; sorted compare for row-order-free tables
    local a="$DIR_A/$2" b="$DIR_B/$2"
    if [ ! -s "$a" ] || [ ! -s "$b" ]; then
        echo "MISSING  $1 ($2)"
    elif cmp -s <(sort "$a") <(sort "$b"); then
        echo "SAME     $1"
    else
        echo "DIFFERS  $1 ($2)"
    fi
}

if [ "$MODE" = run ]; then
    IMAGE="${1:?$USAGE}"
    IN="$(cd "${2:?$USAGE}" && pwd)"
    mkdir -p "${3:?$USAGE}"
    FIRST_FNA="$(ls "$IN"/fna/*.fna | head -1)"
    WORK="$(cd "$3" && pwd)"
    run_all "$IMAGE" "$WORK"
    exit 0
fi
if [ "$MODE" = both ]; then
    IMAGE_A="${1:?$USAGE}"
    IMAGE_B="${2:?$USAGE}"
    IN="$(cd "${3:?$USAGE}" && pwd)"
    mkdir -p "${4:?$USAGE}"
    FIRST_FNA="$(ls "$IN"/fna/*.fna | head -1)"
    WORK="$(cd "$4" && pwd)"
    run_all "$IMAGE_A" "$WORK/a"
    run_all "$IMAGE_B" "$WORK/b"
    DIR_A="$WORK/a" DIR_B="$WORK/b"
    echo "[smoke] $IMAGE_A (a) vs $IMAGE_B (b)"
else
    DIR_A="$(cd "${1:?$USAGE}" && pwd)"
    DIR_B="$(cd "${2:?$USAGE}" && pwd)"
    echo "[smoke] $DIR_A ($(cat "$DIR_A/arch.txt" 2>/dev/null)) vs $DIR_B ($(cat "$DIR_B/arch.txt" 2>/dev/null))"
fi
same "snippy snps.tab" snippy/snps.tab
same "raxml-ng best tree" raxml/r.raxml.bestTree
same "iqtree2 tree" iqtree/i.treefile
same "gubbins final tree" gubbins/g.final_tree.tre
same "gubbins recombination" gubbins/g.recombination_predictions.gff
same "panaroo presence/absence" panaroo/gene_presence_absence.Rtab
same "chewBBACA alleles" chew/allelecall/results_alleles.tsv
same "mlst" mlst.tsv
same "quast report" quast/report.tsv
same "fastbaps clusters" fastbaps/clusters.tsv
