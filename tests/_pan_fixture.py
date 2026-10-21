"""Pangenome outputs for the synthetic run of the golden tests (window plan, Phase 9).

``write_pan_outputs(output_dir)`` writes a small tree, a Panaroo gene
presence table, fastGEAR clades, chewBBACA allele calls and MOB-suite plasmid
clusters for the samples S1 to S6, so every standalone pan panel has data.
With ``snippy_reference`` the tree is a snippy tree: that sample's leaf is
named "Reference" and run_config.json names it as the snippy reference.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

SAMPLES = ("S1", "S2", "S3", "S4", "S5", "S6")

TREE = "((S1:0.10,S2:0.20):0.05,((S3:0.10,S4:0.12):0.04,S5:0.30):0.03,S6:0.40);\n"

# fastGEAR: two clades (lineages), three clusters.
CLADES = {"S1": ("1", "1"), "S2": ("1", "1"), "S3": ("2", "2"), "S4": ("2", "2"), "S5": ("2", "3"), "S6": ("1", "3")}

# Panaroo genes: (gene, annotation, samples with the gene).
GENES = (
    ("geneCore", "core gene", SAMPLES),
    ("geneA", "accessory A", ("S1", "S2", "S6")),
    ("geneB", "accessory B", ("S3", "S4")),
    ("geneC", "accessory C", ("S4", "S5")),
    ("geneD", "accessory D", ("S2",)),
)

# chewBBACA allele calls (one row per sample, five loci).
ALLELES = {
    "S1": (1, 1, 1, 1, 1),
    "S2": (1, 1, 1, 2, 1),
    "S3": (2, 2, 1, 1, 3),
    "S4": (2, 2, 1, 1, 4),
    "S5": (2, 3, 2, 1, 4),
    "S6": (1, 1, 3, 3, 1),
}

# MOB-suite plasmids: (sample, plasmid id, primary cluster, replicon, mobility).
PLASMIDS = (
    ("S1", "S1:AA001", "AA001", "rep_cluster_1", "conjugative"),
    ("S4", "S4:AA001", "AA001", "rep_cluster_1", "conjugative"),
    ("S2", "S2:AA002", "AA002", "rep_cluster_2", "mobilizable"),
)


def write_pan_outputs(output_dir: str | Path, snippy_reference: Optional[str] = None) -> Path:
    root = Path(output_dir)
    (root / "trees").mkdir(parents=True, exist_ok=True)
    tree = TREE
    if snippy_reference:
        tree = tree.replace(f"{snippy_reference}:", "Reference:")
        config = {"inputs": {"snippy_reference": f"/references/{snippy_reference}.fasta"}}
        (root / "run_config.json").write_text(json.dumps(config), encoding="utf-8")
    (root / "trees" / "core_gene_tree.nwk").write_text(tree, encoding="utf-8")

    panaroo = root / "panaroo" / "out"
    panaroo.mkdir(parents=True, exist_ok=True)
    header = ["Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", *SAMPLES]
    lines = [",".join(header)]
    for gene, annotation, present in GENES:
        cells = [f"{gene}_{sample}" if sample in present else "" for sample in SAMPLES]
        lines.append(",".join([gene, "", annotation, str(len(present)), str(len(present)), *cells]))
    (panaroo / "gene_presence_absence.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fastgear = root / "fastgear"
    fastgear.mkdir(parents=True, exist_ok=True)
    rows = ["Sample\tFastGearLineage\tFastGearCluster"]
    rows += [f"{sample}\t{lineage}\t{cluster}" for sample, (lineage, cluster) in CLADES.items()]
    (fastgear / "sample_clades.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    chewbbaca = root / "chewbbaca" / "out"
    chewbbaca.mkdir(parents=True, exist_ok=True)
    rows = ["FILE\t" + "\t".join(f"locus{i}.fasta" for i in range(1, 6))]
    rows += [f"{sample}\t" + "\t".join(str(v) for v in values) for sample, values in ALLELES.items()]
    (chewbbaca / "results_alleles.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    mobsuite = root / "mobsuite"
    mobsuite.mkdir(parents=True, exist_ok=True)
    rows = ["Sample\tPlasmidID\tRepliconTypes\tPredictedMobility\tPrimaryClusterID\tSecondaryClusterID"]
    rows += [f"{s}\t{pid}\t{rep}\t{mob}\t{primary}\t-" for s, pid, primary, rep, mob in PLASMIDS]
    (mobsuite / "plasmids.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    return root
