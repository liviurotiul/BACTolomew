# Validation plan: BACTolomew against published results

Goal: show that BACTolomew gives the same results as well-known, published analyses, on several species. "Same" means the same biological conclusion within a stated tolerance, not identical files (tool and database versions always differ).

Status (2026-10-19): plan only, nothing downloaded or run. Full research notes (about 25 other candidates, why they ranked lower, all sources): `docs/VALIDATION_DATASETS_RESEARCH_2026-10-19.md`.

## 1. Datasets

| # | Dataset | Species | Genomes | Published results we compare against | Priority |
|---|---|---|---|---|---|
| A | Holt et al. 2016, Global Clone 1 | *Acinetobacter baumannii* | 44 assemblies | RAxML tree (Gubbins-masked), Gubbins recombination blocks, Pasteur + Oxford STs, KL/OCL, resistance genes and phenotypes | 1 |
| B | KASPAH, Gorrie et al. 2022 | *Klebsiella pneumoniae* complex | 328 assemblies (subset of ~60 is fine) | Panaroo output on the same assemblies, STs and genotypes (Table S1), Vitek2 phenotypes, Microreact tree | 1 |
| C | Mixao et al. 2023 (ReporTree data) | *Listeria monocytogenes* | ~1,874 assemblies (subset ~50) | chewBBACA cgMLST profiles (Pasteur scheme, 1,748 loci), STs | 2 |
| D | van der Putten et al. 2022, simulated | *E. coli* (in silico evolved) | 96 per simulation | True tree, plus RF/KC scores of 17 published workflows | 2 |
| E | Lees et al. 2018 "which_tree", simulated | pneumococcus-like | 96 assemblies | True tree, true gene sets, published IQ-TREE / RAxML trees with their distance to truth | 3 |
| F | GenomeTrakr Salmonella Bareilly 1203NYJAP-1 | *Salmonella enterica* | 23 assemblies (18 outbreak, 5 outgroup) | Reference tree, outbreak membership | 3 |

Plus a check that needs no new run (section 4, step 0): the reference run against PubMLST and NCBI Pathogen Detection.

### Where to get them

- **A. Holt 2016.** Microb Genom 2016;2(2):e000052, doi:10.1099/mgen.0.000052. Results: figshare doi:10.4225/49/5690B341A0FDB (`GC1.tree`, `GubbinsOutput.gff` on A1 coordinates, `K_OC.csv`); Supp Tables 1 and 2 (STs, KL, OCL, AbaR, phenotypes, N50). Genomes: 22 ENA WGS sets FBWP01 to FBXJ01, closed A1 (CP010781) and D36 (CP012952), plus 20 public GenBank genomes listed in Supp Table 1. Two leaf names contain spaces: rename.
- **B. KASPAH.** Nat Commun 2022;13:3017, doi:10.1038/s41467-022-30717-6. figshare doi:10.26180/16811344 (assemblies 577 MB, Panaroo 3.6 GB, Table S1), CC BY 4.0.
- **C. Mixao.** Zenodo 7230003 (v1.0.1): `Lm_public.zip` (1.7 GB, can be read by range requests to extract a subset), `Lm_public/profiles/Lm_profile.tsv`, `Lm_metadata.xlsx`. Assemblies named by run accession.
- **D. van der Putten.** Zenodo 10.5281/zenodo.5036179 (`assemblies.tar` 3.2 GB, `full_data.csv` with RF/KC per workflow); true tree in https://github.com/niekh-13/phylogenetic_workflows `input/tree_div100.nwk`.
- **E. Lees 2018.** Wellcome Open Res 2018;3:33, doi:10.12688/wellcomeopenres.14265.2. figshare 10.6084/m9.figshare.5483461 (`simulations.tar.bz2`, 562 MB), trees 10.6084/m9.figshare.5483464, code https://github.com/johnlees/which_tree.
- **F. GenomeTrakr.** Timme et al., PeerJ 2017;5:e3893. https://github.com/WGS-standards-and-analysis/datasets `datasets/Salmonella_enterica_1203NYJAP-1.tsv` (GCA, sha256 per isolate; pin the superseded `.1` versions, they are still on the NCBI FTP). Tree: https://api.opentreeoflife.org/v3/study/ot_301/tree/tree1.tre.
- **Step 0 sources.** PubMLST: https://rest.pubmlst.org/db/pubmlst_abaumannii_isolates. NCBI Pathogen Detection AMRFinderPlus results: https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Acinetobacter/latest_kmer/AMR/ (`PDG000000010.*.amr.metadata.tsv`; store the snapshot with its PDG version, it changes daily). Shelenkov et al. 2024 *A. baumannii* cgMLST/ST/AMR tables: Zenodo 19186935 (newest version of 10143378, CC BY-SA 4.0, fetch, do not bundle).

## 2. What each dataset checks

| BACTolomew step | A Holt | B KASPAH | C Mixao | D van der Putten | E Lees | F Salmonella |
|---|---|---|---|---|---|---|
| MLST | yes | yes | yes | | | |
| Panaroo | consistency | yes (same input) | | | yes (true genes) | |
| Core alignment + IQ-TREE | yes | yes | | yes (true tree) | yes (true tree) | yes |
| Snippy tree | yes | yes | | yes | yes | yes |
| fastGEAR | lineages | | | | | |
| Gubbins | yes (blocks) | | | | | |
| chewBBACA cgMLST | | | yes | | | |
| AMRFinderPlus | genes, phenotype | phenotype | | | | |

No good external benchmark: Bakta annotation, CheckM2, the DIAMOND searches of user-added databases, MOB-suite (the gene categories come from Bakta, AMRFinderPlus and MOB-suite through `bakta_evidence`). Smoke tests only (CDS counts on closed genomes such as ACICU GCF_005519135.2 and AB5075-UW GCF_000963815.1).

## 3. Metrics and what counts as reproduced

- **MLST:** exact ST for every genome with all 7 alleles called. A mismatch is acceptable only if it comes from a newer PubMLST allele missing from our frozen DB.
- **Trees:** prune both to the same leaves, collapse weak branches (UFBoot < 95, bootstrap < 70), normalized Robinson-Foulds (ete3 or dendropy) and quartet distance (tqDist). Real data: named lineages and clades identical, nRF ≤ 0.2. Simulated: our distance to the true tree within the spread of the published IQ-TREE/RAxML trees (E) or among the better half of the 17 workflows (D). Branch lengths are not compared. Holt and KASPAH trees come from read mapping, so compare both our core tree and our Snippy tree.
- **Recombination (A):** map Holt's blocks (A1 coordinates) to our core alignment via `gene_map.tsv` and the Bakta coordinates of A1 (or pick A1 as the ordering genome); base-level recall and precision on core genes. Large blocks (≥ 10 kb) recovered on the same clade.
- **Pangenome (B, E):** same thresholds for core/soft-core/shell/cloud; map clusters (DIAMOND or blastn, ≥ 95% identity, ≥ 80% coverage); core count within 5%, median per-genome Jaccard ≥ 0.95.
- **cgMLST (C, step 0):** pairwise allele distances on loci called by both sides: Spearman ≥ 0.98, median absolute difference ≤ 2 for pairs ≤ 20 alleles, adjusted Rand index ≥ 0.95 for single-linkage clusters at ≤ 9 alleles. Allele IDs compared only for base-scheme alleles (INF alleles are local).
- **AMR:** against NCBI PD, ≥ 98% per-genome agreement of element symbols including point mutations, each difference explained by DB version. Against phenotypes: sensitivity and specificity reported, not pass/fail.

## 4. Steps

0. **No new run (about 1 h):** compare the reference run `runs/sample_78eaa9126850/run_78eaa9126850` (GCA names) with PubMLST STs, NCBI PD AMRFinderPlus results and the Shelenkov cgMLST table for the shared genomes.
1. **Run A (Holt, 44 genomes, about 3 h).** Snippy reference A1 (CP010781).
2. **Run B (KASPAH subset of ~60, about 3.5 h; full 328 overnight).**
3. **Run C (Listeria subset of ~50).** Needs the Listeria Pasteur scheme (chewie-NS) to match their profiles; the installed cgMLST.org Listeria scheme gives distances only.
4. **Runs D, E, F** (phylogeny truth). See the cgMLST caveat below.
5. Write `compare.py` per dataset producing one table: published result, our result, within tolerance or not. Record tool and DB versions from each run's `run_report.json` and `docker/locks/`.

Each run in its own run folder; never touch the reference run.

## 5. Caveats and open decisions

- **cgMLST is never skippable**, so every run needs a scheme for its species. Installed: *A. baumannii*, *L. monocytogenes* (and others in `DB/cgmlst/`). *K. pneumoniae* and *E. coli* are auto-installed from cgMLST.org during the run if listed there (check). *Salmonella* is not on cgMLST.org (EnteroBase): F needs an imported scheme or runs outside the app. For the simulated sets (D, E) chewBBACA results are meaningless; check whether a scheme exists for their species, otherwise compare those trees with a harness outside the app.
- **Assemblies only.** Strong *A. baumannii* sets that publish reads only (Higgins 2017, same cgMLST scheme, ST2 isolates; Schultz 2016, GC2 with Gubbins blocks) need an assembly step we do not have.
- **No Gram-positive** dataset beyond Listeria; a *S. aureus* or *E. faecium* set with published results is still to be found.
- **Where results live:** ship accession lists, fetch scripts, `compare.py` and our expected summaries (e.g. `benchmarks/<set>/`), never genomes or CC BY-SA tables. INSDC data and CC BY data may be redistributed with attribution; check the Holt and Lees figshare licenses before bundling anything.
- **AMR DB version:** pin to the NCBI PD version for exact agreement, or keep ours and report drift.
