# Published benchmark datasets for validating BACTolomew

Date: 2026-10-19. Each candidate's existence, its genome access and its published results were checked by a second pass. Where the two passes disagreed, the corrections from the second pass are used here.

## 1. Short answer

**Partly.** No published dataset covers the whole pipeline with accepted results, and none comes close for *A. baumannii* ST2/ST492. Four datasets together cover most steps, and all four offer **assemblies**, so nothing has to be reassembled:

1. **Holt et al. 2016, *A. baumannii* GC1 (44 genomes).** Same species. It publishes a tree, Gubbins blocks, STs, KL/OCL types and phenotypes as machine-readable files.
2. **Shelenkov et al. 2024, *A. baumannii* typing tables (a ST2/ST492 subset).** It has cgMLST profiles on the same 2,390 cgMLST.org loci (`ACICU_RS...`) that our chewBBACA uses, plus Pasteur STs and AMR genes. Combined with the **NCBI Pathogen Detection** AMRFinderPlus table for the same GCA accessions, it gives a same-input check of MLST, cgMLST and AMR on your own STs.
3. **Lees et al. 2018 "which_tree" simulation (96 genomes).** It is the only option here with a known **true tree** and true gene content. Use it for core alignment, IQ-TREE and Panaroo.
4. **KASPAH *K. pneumoniae* (328 genomes).** The published Panaroo output was built from exactly the assemblies that are deposited, so it is the cleanest real-data Panaroo comparison.

## 2. Recommended datasets

| | Holt 2016 GC1 | Shelenkov 2024 + NCBI PD | Lees 2018 which_tree | KASPAH (Gorrie 2022) |
|---|---|---|---|---|
| **Citation** | Holt KE et al., Microb Genom 2016;2(2):e000052, doi:10.1099/mgen.0.000052 | Shelenkov A, Mikhaylova Y, Akimkin V, Data 2024;9(2):22, doi:10.3390/data9020022; Feldgarden M et al., Microb Genom 2022;8(6), doi:10.1099/mgen.0.000832 | Lees JA et al., Wellcome Open Res 2018;3:33, doi:10.12688/wellcomeopenres.14265.2 | Gorrie CL et al., Nat Commun 2022;13:3017, doi:10.1038/s41467-022-30717-6 |
| **Organism** | *A. baumannii* GC1 (ST1, ST19, ST20, ST81) | *A. baumannii*, subset to ST2 and ST492 | Simulated *S. pneumoniae*-like genomes on a 96-tip *Listeria* tree topology | *K. pneumoniae* species complex |
| **Genomes / assemblies** | 44 (the tree leaves), all as assemblies: 22 ENA WGS sets FBWP01 to FBXJ01, closed A1 (CP010781) and D36 (CP012952), plus 20 public GenBank genomes | 17,546 GCA assemblies in total: 35 are ST492 and about 11,100 are ST2 | 96 assemblies of simulated reads (Velvet), plus the 96 true genomes | 328 assemblies (Unicycler; hybrid where long reads exist) |
| **Published results** | `GC1.tree` (RAxML, Newick, bootstraps); `GubbinsOutput.gff` (about 455 blocks on A1 coordinates); `K_OC.csv`; Supp Tables 1 and 2 with Pasteur and Oxford ST, KL, OCL, AbaR type, phenotype string and N50 | `table_summary` (IC, Pasteur ST, OXA-51-like, KL, OCL), `table_cgmlst` (about 2,390 `ACICU_RS` loci), `table_amr`; NCBI PD `*.amr.metadata.tsv` (AMRFinderPlus genotypes per GCA) and SNP cluster trees | `RealTree.nwk` (true tree), true gene MSAs (1,992), true per-genome gene sets; inferred trees from RAxML, IQ-TREE, FastTree and parsnp with RF and KC distances to the truth | Panaroo v1.1.2 output on the same assemblies (`panaroo.tar.gz`), Table S1 (ST, genotypes), Vitek2 phenotypes with EUCAST 2020 breakpoints, a Microreact tree (FastTree on read-mapping SNVs) |
| **Steps validated** | MLST, core tree and Snippy tree topology, Gubbins, fastGEAR (lineages 1 and 2), AMR (AbaR content, phenotype) | MLST, chewBBACA (distances, and allele IDs where they come from the base scheme), AMRFinderPlus genes and point mutations | core_alignment, core_tree, Snippy tree, Panaroo (against true gene content) | Panaroo, MLST, AMRFinderPlus (against phenotype), core tree (topology only) |
| **Download** | figshare doi:10.4225/49/5690B341A0FDB (article 2061936); assemblies from ENA (e.g. `ftp.ebi.ac.uk/pub/databases/ena/wgs/public/fbw/FBWR01.fasta.gz`) and NCBI | Zenodo 10143378 (use the newest version: 19186935, which adds VFDB and replicon tables); `https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Acinetobacter/latest_kmer/AMR/`; genomes by GCA through NCBI Datasets | figshare 10.6084/m9.figshare.5483461 (`simulations.tar.bz2`, 562 MB); trees at 10.6084/m9.figshare.5483464; code at https://github.com/johnlees/which_tree | figshare 10.26180/16811344 (assemblies 577 MB, Panaroo 3.6 GB, Table S1), CC BY 4.0 |
| **Size / time** (Bakta full DB at about 2.5 min per genome) | about 2 h of Bakta, about 3 h end to end | about 60 genomes (35 ST492 plus about 25 ST2): about 2.5 h of Bakta, 3.5 h end to end | 96 genomes: about 4 h of Bakta, 5 to 6 h end to end; the MLST, chewBBACA and AMR steps are meaningless and chewBBACA may stop the run | 328 genomes: about 14 h of Bakta, about 18 h end to end (overnight; setup installs only the full Bakta DB) |

## 3. Coverage of each step

✔ = a direct comparison with published or true results; ~ = qualitative or tool-concordance only; blank = nothing to compare.

| BACTolomew step | Holt GC1 | Shelenkov + NCBI PD | Lees which_tree | KASPAH | Best available if none |
|---|---|---|---|---|---|
| MLST | ✔ Pasteur + Oxford | ✔ Pasteur | (meaningless) | ✔ | PubMLST for the reference run's own GCAs |
| Bakta / annotation | | | ~ true CDS sets | ~ Prokka gene counts | **No good benchmark.** Smoke test: CDS counts against RefSeq PGAP on ACICU GCF_005519135.2 and AB5075-UW GCF_000963815.1 |
| QUAST / CheckM2 | ~ N50 and contigs in Supp Tables | | | | **Weak.** Closed genomes should be about 100% complete; for contamination, the Panaroo paper's contaminated simulation (Zenodo 3599800) |
| Panaroo | ~ consistency only | | ✔ true gene content (gene sets must be mapped) | ✔ same-input Panaroo | Panaroo paper simulations: true presence/absence for 21 sets |
| core alignment + IQ-TREE | ✔ topology against GC1.tree | ~ NCBI SNP trees | ✔ **true tree** | ~ Microreact topology | van der Putten 2022 leaderboard (RF/KC, IQ-TREE) |
| fastGEAR | ~ GC1 lineages 1 and 2 | | (no recombination simulated) | | SPARC: 15 sequence clusters (qualitative) |
| Gubbins | ✔ blocks in a GFF (coordinates need mapping) | | (no recombination simulated) | | **No set with true recombination.** Optional: Schultz 2016 GC2 blocks (reads only) |
| Snippy | ✔ topology (they mapped reads to A1) | | ✔ true tree | | Gen-FS simulated *Salmonella* truth VCF (reads only) |
| chewBBACA cgMLST | ~ consistency only | ✔ same 2,390 loci | (would fail or be meaningless) | | Higgins 2017 S4 Table (same scheme, ST2 isolates, reads only) |
| AMRFinderPlus genes + mutations | ~ AbaR genes, phenotype string | ✔ NCBI runs the same tool | | ✔ against phenotype | Javkar 2021 *A. baumannii* carbapenem phenotypes (349 assemblies) |
| DIAMOND marker DBs | | | | | **No external equivalent** (in-house databases) |
| MOB-suite | | | | ~ via the hybrid closed genomes | **No good benchmark for *A. baumannii*.** Optional: PHA4GE 2022 (closed genome against shovill drafts of the same isolate, Zenodo 5604579) |

## 4. How each comparison would work, and what counts as "reproduced"

General rule: tool and database versions differ from every published set, so "reproduced" means **the same biological conclusion within a stated tolerance**, not identical files. Record the DB versions of both sides for each comparison.

**Trees (Holt GC1, Lees, KASPAH)**
- Prune both trees to the same leaves, then compute the **normalized Robinson-Foulds** distance (ete3 or `dendropy`) and the **quartet distance** (tqDist).
- For real data, first collapse branches with UFBoot below 95 or bootstrap below 70 on both sides.
- Holt GC1:
  - Membership of lineage 1 and lineage 2 must match exactly.
  - Normalized RF on the resolved splits should be 0.2 or less, which is a reasonable target because their tree was built by read mapping to A1 with Gubbins on a whole-genome alignment, not on a core-gene alignment.
  - Branch lengths are not compared.
- Lees (true tree):
  - Score our tree by RF and KC distance to `RealTree.nwk`.
  - It counts as reproduced if it lands within the spread of the published IQ-TREE and RAxML core-gene trees in `trees.tar.gz`. A tree clearly worse than theirs is a failure.

**Recombination (Holt GC1 Gubbins)**
- Their blocks are on A1 genome coordinates; ours are on the ordered core alignment.
- To map them, make A1 (closed) the ordering genome, or use `gene_map.tsv` together with the Bakta coordinates of A1.
- Then compute **base-level recall and precision** of the recombined positions, restricted to core genes, plus the overlap of blocks per branch.
- Expected outcome: large blocks (10 kb or more) are recovered and placed on the matching clade. Our r/m will differ by design, because we run on core genes only and per lineage.

**Pangenome (KASPAH, Lees, Panaroo simulations)**
- Compare core, soft-core, shell and cloud counts at the same thresholds. For KASPAH, compare per species, *K. pneumoniae* sensu stricto first.
- Map our `pan_genome_reference` clusters to theirs, or to the true genes, with DIAMOND or blastn at 95% identity and 80% coverage or more.
- Then compute the per-genome **Jaccard** of presence vectors on the mapped clusters.
- Target: core count within 5%, and median per-genome Jaccard of 0.95 or more. Bakta against Prokka and a different Panaroo version account for most of the remainder.

**cgMLST (Shelenkov, Higgins)**
- Compute pairwise allele distances on the loci both sides call, and compare them by **Spearman correlation** (target 0.98 or more) and by the median absolute difference among pairs at 20 alleles or fewer (target 2 or fewer).
- **Cluster concordance**: adjusted Rand index of single-linkage clusters at 9 or fewer alleles (target 0.95 or more). Higgins reports 0 to 8 alleles within a pulsotype.
- Allele IDs: compare them only for alleles in the base cgMLST.org scheme. INF alleles are local.

**MLST**
- Exact ST match is required for every genome whose 7 alleles are all called. A mismatch is acceptable only when it traces to a newer PubMLST allele that is missing from our frozen 2025-12-29 DB.

**AMR**
- Against NCBI PD: per-genome agreement of element symbols, point mutations included. Target 98% or more. Each discordance must be explained by the AMRFinderPlus DB version, which NCBI records per isolate. For an exact check, rerun with the same DB version NCBI used.
- Against phenotype (KASPAH, Holt, optionally Javkar): sensitivity and specificity of predicted carbapenem resistance (OXA-23/24/58-like, NDM, ISAba1 upstream of OXA-51-like where called) against R/S calls. Report these as a description, not as a pass/fail, because AST methods differ.

**QC**
- Our N50 and contig count should match the published per-assembly values exactly (it is the same file), apart from QUAST's 500 bp contig cutoff.

## 5. Other candidates and why they ranked lower

- **Higgins 2017 *A. baumannii* cgMLST (53 isolates, PRJEB12082):** best for chewBBACA with the same scheme and ST2 isolates, but reads only (64 runs to map to 53 isolates). It becomes first tier if we add an assembly step.
- **Schultz 2016 Vietnam GC2 (147 isolates):** ST2, MICs, a Gubbins GFF and a tree, but reads only, and the tree and blocks cover only the 66 GC2 isolates (lane IDs with `#`).
- **van der Putten 2022 (in silico *E. coli*):** an RF/KC leaderboard for 17 workflows with a true tree, but it has no Panaroo entry and the 3 GB assembly archives were not unpacked. It is the alternative to Lees.
- **SPARC *S. pneumoniae* (615 assemblies, Dryad core tree and COGs):** strong for fastGEAR (15 clusters) and the tree, but large and off-species; recombination blocks exist only as figures.
- **Panaroo paper simulations (Zenodo 3599800):** true presence/absence, but the contigs must be extracted from Prokka GFFs and a gene-mapping script written. A good add-on for Panaroo.
- **Panaroo paper *Klebsiella* 319 and TB outbreak:** consistency checks only. AMRFinderPlus has no TB organism option.
- **Gen-FS *Salmonella* Bareilly (23 assemblies):** tiny and clonal, with a polytomy reference tree; it tests cluster recovery only. The other Gen-FS sets are mostly reads.
- **Gen-FS simulated *Salmonella*:** a true SNP VCF for Snippy, but reads only.
- **Ahrenfeldt 2017 evolved *E. coli*:** a true tree, but the genomes are reference-mapped pseudo-assemblies, which most steps cannot use.
- **Mixao *Listeria* chewBBACA profiles (Zenodo 7230003):** an exact allele-by-allele chewBBACA regression check, but off-species and chewBBACA 2.8.5 against itself.
- **Palma, Lüth, pyMLST, Mellmann ring trial:** reads only, and the results are mostly PDF tables or ad hoc schemes.
- **chewBBACA tutorial (*S. agalactiae*):** a tool-integration regression test with chewBBACA 2.x expected outputs, not biology.
- **PubMLST *A. baumannii*:** a live curated database, not a frozen benchmark. Its cgMLST v1 is marked experimental. Use it for an MLST check on the reference run.
- **Javkar 2021 *A. baumannii* carbapenem cohort (349 assemblies):** the best *A. baumannii* phenotype set; a good AMR add-on, but it has no tree file.
- **CDC/FDA AR Bank ACI panel:** MICs are in the BioSample antibiograms, but only about 8 of 51 isolates have assemblies.
- **BV-BRC AMR table, NARMS AMRFinder set, ResFinder 4.0:** phenotype only and heterogeneous; NARMS and ResFinder have no *Acinetobacter*.
- **EuSCAPE *K. pneumoniae*:** 1,771 assemblies, but its supplements (MICs, clusters) were not reachable.
- **Neil 2026 IC2 chromosomes:** the closest published method stack (Bakta, Panaroo, IQ-TREE, Gubbins) and IC2, but there is no machine-readable tree or group table, chromosomes only.
- **Koong 2021 GC1 188:** tables are in a PDF and the tree is a figure; Holt 2016 is the machine-readable core of it.
- **MOB-suite 2018 benchmark, Paganini 2021:** the truth would have to be re-derived from closed genomes, and the drafts need assembly.
- **Bakta benchmark genomes:** Bakta v1.1 outputs, one genome per species.
- **Kaptive KL/OCL sets:** BACTolomew has no KL/OCL step.
- **Gubbins PMEN1/ST239, ClonalFrameML data, CG258 (Wyres/Bowers), fastGEAR paper, Horesh *E. coli*, PIRATE *S. aureus*:** these are alignments only, figures only, or off-method, and cannot be fed to an assembly-only pipeline.

## 6. Practical plan

**Tier 0 (no new run, about 1 h of scripting).** The reference run `run_78eaa9126850` already uses GCA accessions.
1. Compare `mlst/results.csv` with PubMLST Pasteur STs (`rest.pubmlst.org/db/pubmlst_abaumannii_isolates`).
2. Compare `amrfinder/<s>.tsv` with NCBI PD `PDG000000010.*.amr.metadata.tsv`.
3. Compare `chewbbaca/out/results_alleles.tsv` with Shelenkov `table_cgmlst` on the overlapping genomes. Many of the 35 GenBank ST492 genomes are probably in the run already; check the overlap.

**Tier 1 (the shippable "validation set", about 100 genomes, 6 to 7 h on the reference machine).** Two runs, kept separate so that each tree comparison has an exact leaf set:
- **Run A: Holt GC1, 44 genomes.** Set `phylo.outgroup` or compare unrooted; set the Snippy reference to A1 (CP010781). Validates MLST, both trees, Gubbins, fastGEAR and AMR.
- **Run B: *A. baumannii* ST2/ST492, about 60 genomes.** All 35 ST492, plus about 25 ST2 picked across IC2 diversity from the Shelenkov table, plus ACICU GCF_005519135.2 as a GC2 anchor. Validates MLST, chewBBACA, AMRFinderPlus and fastGEAR/Gubbins on your own lineages.

**Tier 2 (method truth, run once for the paper).**
- **Run C:** Lees which_tree, 96 genomes. Needs a way to skip MLST, chewBBACA and AMR (see question 1).
- **Run D:** KASPAH, 328 genomes, overnight.
- **Optional:** the Panaroo simulation with true presence/absence.

**What to ship.** Do not ship genomes; the sizes are too large.

What to ship:
- `benchmarks/<set>/accessions.tsv` (accession, exact version, sha256, original name, BACTolomew-safe name)
- a fetch script (NCBI Datasets, ENA, figshare, Zenodo by DOI)
- a frozen comparator snapshot, or the DOI plus file checksum to fetch it
- a `compare.py` that writes one pass/fail table using the tolerances in section 4
- the expected-output summary of our own run, with tool and DB versions from `run_report.json` and `docker/locks/tool_versions.json`

Track renames: Holt has 2 leaf names with spaces, Schultz and the Panaroo sets use `#` lane IDs, and several GCA `.1` versions are superseded, so pin the exact version.

**Licensing.**
- **INSDC data** (ENA, GenBank, SRA) is free to redistribute, so the genomes and accession lists are fine.
- **KASPAH** figshare is CC BY 4.0: attribution only.
- **Shelenkov** Zenodo is **CC BY-SA 4.0**: tables derived from it must carry the same license, so fetch it at test time rather than bundling it.
- **Holt** figshare and **Lees** figshare: the license was not checked (figshare data usually defaults to CC BY 4.0). Confirm before bundling; otherwise fetch by DOI.
- **NCBI PD** tables are US government public domain, but they change daily, so store the snapshot file with its PDG version (e.g. `PDG000000010.1980`).
- **Journal supplements** (Higgins S4, Javkar Table S1): PLOS and Frontiers are CC BY. Cite them and fetch them rather than redistributing modified copies.

## 7. Open questions

1. Should the pipeline get a way to skip steps per run (MLST, chewBBACA, AMR) so simulated benchmarks like Lees can run end to end, or should Run C be driven by a small harness outside the app?
2. Is adding a read-assembly step (Shovill) in scope? It would unlock Higgins 2017 (the same cgMLST scheme as ours) and Schultz 2016 GC2 as first-tier *A. baumannii* benchmarks.
3. Where should the validation set live: in the repo (`benchmarks/`, with fetch and compare scripts) and run on demand, or only in the paper's supplement?
4. For AMR, should the comparison pin our AMRFinderPlus DB to the version NCBI PD used (exact concordance), or keep our 2025-03-25.1 DB and report the differences as DB drift?

## Sources (main)
- Holt 2016: https://pmc.ncbi.nlm.nih.gov/articles/PMC5320584/ ; https://api.figshare.com/v2/articles/2061936
- Shelenkov 2024: https://zenodo.org/records/10143378
- NCBI PD: https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Acinetobacter/latest_kmer/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9455714/
- Lees 2018: https://pmc.ncbi.nlm.nih.gov/articles/PMC5930550/ ; https://api.figshare.com/v2/articles/5483461 ; https://github.com/johnlees/which_tree
- KASPAH: https://api.figshare.com/v2/articles/16811344 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9156735
- van der Putten 2022: https://zenodo.org/records/5036179 ; https://github.com/niekh-13/phylogenetic_workflows
- Higgins 2017: https://pmc.ncbi.nlm.nih.gov/articles/PMC5464626/
- Schultz 2016: https://pmc.ncbi.nlm.nih.gov/articles/PMC5320574/ ; https://api.figshare.com/v2/articles/2064357
- Panaroo simulations: https://zenodo.org/records/3599800 ; https://github.com/gtonkinhill/panaroo_manuscript
- SPARC: https://pmc.ncbi.nlm.nih.gov/articles/PMC4622223/ ; https://datadryad.org/dataset/doi:10.5061/dryad.t55gq
- Javkar 2021: https://www.frontiersin.org/articles/10.3389/fmicb.2021.714284
- PHA4GE: https://zenodo.org/records/5604579
- PubMLST: https://rest.pubmlst.org/db/pubmlst_abaumannii_seqdef/schemes
