# Phylogeny rework (2026-10-10): implementation plan

> Status (2026-10-21): implemented, including Rework 2 (Gubbins per fastGEAR lineage) and the
> branch length sets; the notes and reruns below are the log. Section 2 describes the code
> before the rework, not the current code. Current pipeline steps: `core_alignment`,
> `fastgear_clades`, `gubbins`, `core_tree` (IQ-TREE) and `phylogeny` (Snippy tree only).

Recombination-aware core tree from a genome-ordered Panaroo core alignment, the same
automatic rules on any dataset, every choice recorded in the run, and a GUI that shows
the recombination-filtered tree by default.

Two lanes work in parallel: the **pipeline lane** (`src/pipeline`) and the **GUI lane**
(`src/annotation_explorer_polars`). Section 4 is the contract between them.

## 1. Goal and scope

Goal: one core alignment, ordered by a real genome, feeds Gubbins, the core tree and
fastGEAR. The default tree is IQ-TREE (ModelFinder, 1000 UFBoot, constant-site counts) on the Gubbins
recombination-filtered sites, rooted by a recorded rule. The Snippy tree stays as a
secondary tree labeled "not recombination-corrected".

Review items covered (from `pipeline_review/assessment.json`):

| Item | Problem | Covered by |
|---|---|---|
| P09 | No tree is recombination-corrected; dead self-recursive `_resolve_core_tree_alignment`; misleading log line | 3.3 core tree on Gubbins output; resolver deleted; log line rewritten |
| P10 | GUI shows the Snippy tree, cannot select the Gubbins tree | 4 contract, GUI default = recombination-filtered tree, source selector |
| P30 | No branch support, fixed model, logs discarded | 3.3 IQ-TREE MFP + UFBoot, raw outputs kept |
| P39 | Gubbins runs on a gene concatenation in arbitrary Panaroo order | 3.1 genome ordering plus gene map |
| P62 | Gubbins hits the 5-iteration cap, no seed, nothing configurable | 3.2 seed, iterations, convergence parsed and surfaced |
| P63 | Tree, Gubbins and fastGEAR use different Panaroo alignments | 3.1 one shared resolver and one ordered file for all three |
| P64 | No rooting decision; GUI roots by a fastGEAR display heuristic | 3.5 midpoint or outgroup, recorded, shown in the GUI |
| P77 | Mixed tree names, `Reference` tip, corrected tree outside `trees/` | 3.6 `trees/*.nwk`, `trees/trees.json`, tip renamed |

Out of scope (user decisions): rerunning Snippy, Snippy reference choice and masking
(P29, P88), a Panaroo-graph consensus gene order (possible later upgrade).

## 2. Current behavior (verified in code and on the reference run)

Reference run `runs/sample_78eaa9126850/run_78eaa9126850`, 153 genomes.

Gubbins (`src/pipeline/GubbinsWorkflow.py`)
- Input: first existing of `PANAROO_ALIGNMENT_CANDIDATES` (13-20, 51-64), i.e. the
  unfiltered `panaroo/out/core_gene_alignment.aln` (2,627,569 bp, 2676 gene blocks in
  Panaroo order). It is copied and sanitized to `gubbins/input/core_gene_alignment.aln`
  (86-169; 193 non-ACGTN- characters set to N on this run).
- Command (215-222): only `--prefix core_gene --threads N`. No `--seed`, so RAxML seeds are
  random per iteration (`command.stdout.log` lines 191, 245: `-p 4239`, `-p 2088`); default
  `--iterations 5`; default converge method `weighted_robinson_foulds`; no outgroup.
- Reference run: Gubbins 3.4.3, `command.stdout.log` line 1258 "Maximum number of
  iterations (5) reached.", 1105.9 s total at 64 threads. The step is still "success".
- `summary.json` (241-257) records paths and counts only: no convergence, no parameters.
- Outputs on the run: 28,672 filtered polymorphic sites, 4129 recombination blocks, SNPs
  inside/outside recombination 88,918 / 32,894 (r/m as SNP ratio 2.70), final tree 153 tips.
- A stale `gubbins/tmp1tj7krj1/` temp folder from June is still in the run.

Phylogeny (`src/pipeline/PhylogenyWorkflow.py`)
- `_resolve_core_tree_alignment` (97-105) has no caller and calls itself.
- The core tree (468-484) uses `_TREE_ALIGNMENT_CANDIDATES` (21-24), whose first entry is
  Panaroo's **filtered** `core_gene_alignment_filtered.aln` (2,324,240 bp, 2471 genes). So
  the core tree, Gubbins and fastGEAR each read a different alignment (P63). The task
  brief called it "unfiltered"; the code and file listing say filtered.
- `_run_raxml_tree` (108-149): `raxml-ng --search --model GTR+G --seed 1` in a temp dir,
  only `.raxml.bestTree` kept; no support, model, log or likelihood survive.
- `trees/core_gene_tree.nwk` and `trees/snippy_raxml.tree` (475-476): both 153 tips, no
  internal labels. The Snippy tree has a `Reference` tip for GCA_028622775.1 (the reference
  is never run as a sample, 388-389; `snippy/reference.txt` holds its path).
- Snippy and the core tree run in one function and one step (`run_phylogeny`, 449-531).

Runner (`src/pipeline/Runner.py`)
- `PHYLOGENY_OUT_FILES` 180-183; thread plan 200-227 (`gubbins` 64, `phylogeny_raxml` 64).
- Gubbins step 1617-1673: reused only when Panaroo was reused and the outputs validate.
  No input fingerprint, no parameter check.
- Phylogeny step 1686: reused only when `not gubbins_regenerated`, so **any Gubbins rerun
  reruns Snippy** (hours), although Snippy does not read Gubbins output. Line 1699 logs
  "rebuilding phylogeny from filtered alignment", which is false (P09).
- fastGEAR step 1740-1762: reused when Panaroo was reused and three files exist. No input
  fingerprint.
- There is no run-level warnings list in `run_state.json` / `run_report.json`.

fastGEAR (`src/pipeline/FastGearWorkflow.py`)
- Own candidate list (16-23) and resolver (57-68): unfiltered alignment, copied to
  `fastgear/input/core_gene_alignment.aln` (71-85, a second 409 MB copy). On the reference
  run: 2 lineages (6 and 147 genomes), 4 clusters; about 2 h wall time (input staged
  17:37, result 19:34).

Validation (`src/pipeline/RunValidation.py`)
- `expected_phylogeny_paths` 447-452 requires both old tree names; `validate_gubbins_outputs`
  596-622; the legacy Gubbins input check reads `gubbins/input/core_gene_alignment.aln`
  (266-274). `RESULT_OPTIONAL_STEPS` (49-51) is the mechanism for steps old runs lack.
- The GUI deep validation (`data.py` via `validate_completed_run`) uses the same
  validators, so old runs must keep validating after the change.

Tools (`src/pipeline/ToolEnv.py`, container `bactolomew-tools:0.1.0`)
- IQ-TREE 2.4.0 (`iqtree2`) exists only in the `gubbins` env; `_tool_env_key` (126-148)
  has no entry for it. `minimap2` is in `phylogeny` and `assembly_qc`.
- Checked in the container: `-m MFP+ASC -B 1000 --seed 1` works (picked SYM+ASC on a
  12-taxon test); `--fast` refuses `-B` ("Ultrafast bootstrap does not work with -fast"),
  `--fast --alrt 1000` works.

GUI
- `model/pan_data.py` 146-197 `resolve_raxml_tree_path`: Snippy tree first, then
  `core_gene_tree.nwk`, then any `*raxml*` file. 200-248: sources RAxML SNP / RAxML core
  gene / Panaroo; the Gubbins tree is never offered (P10).
- `panaroo_core_accessory.py` 948-972 `tree_sources` / `_canonical_tree`: first source wins;
  there is no source selector in the UI. 53-60: branch-length unit guessed from the
  source name ("raxml"). 1021: the hidden leaf-order tree maps the `Reference` leaf.
- `panels_controller.py` 85-102: every pan tab is built with `tree_root_mode=ROOT_GROUPS`
  (root on the edge that keeps fastGEAR clades together, `phylo.py` 243-248) (P64).
- `run_expectations.py` 17-24, 382: old file names hard-coded for in-progress views.
- `tests/test_retired_shims.py` 116 asserts that `set_tree_source`, `tree_source_key`,
  `tree_source_label`, `tree_source_path` do not exist: the new selector must use other
  names.

Ordering feasibility (scratch probe on the reference run)
- Most contiguous genome by `assembly_qc/summary.tsv`: GCA_003931755.1 (1 contig,
  N50 3,943,275; ties at 1 contig broken by N50).
- Its Bakta GFF3 `ID` equals the Panaroo gene ids. Full alignment: 2610 of 2676 genes
  placed, 65 absent from it, 1 refound gene without a GFF entry, 6 multi-copy clusters.
  Filtered alignment: 2411 of 2471 placed.
- Cohort diversity on the full alignment: 67,806 variable sites, mean pairwise
  p-distance 0.0030, max 0.0153 (median 6034 SNPs).

## 3. Design

### 3.1 Ordered core alignment (new step `core_alignment`)

New module `src/pipeline/CoreAlignment.py`, output folder `core_alignment/`.

- **Source alignment (one resolver for all consumers, P63):** `phylo.core_alignment`
  = `"full"` (default, `core_gene_alignment.aln` + `core_alignment_header.embl`) or
  `"filtered"` (`core_gene_alignment_filtered.aln` + `core_alignment_filtered_header.embl`).
  Default full: Gubbins is the step that masks recombinant or misaligned dense SNP
  clusters, and it keeps the before/after comparison about ordering only. Genes that
  Panaroo's entropy filter removes are flagged (`in_panaroo_filtered`, from membership in
  the filtered header) so hotspots on them can be read with care. `.gz` variants are
  accepted. The three local candidate tuples and resolvers (Gubbins, Phylogeny, fastGEAR)
  are deleted.
- **Gene blocks:** parsed from the Panaroo EMBL header (`FT feature a..b` + `/label=`
  cluster`.aln`). A missing header fails the step with a clear message.
- **Placement of a cluster in a cohort genome:** ids from that genome's column of
  `gene_presence_absence.csv`, split on `;`; positions from `bakta/<s>/<s>.gff3` (`ID`,
  seqid, start, end, strand). Ids on one contig and strand within 5 kb of each other are
  fragments of one gene (one placement spanning them); otherwise the cluster is
  multi-copy and the copy with the lowest (contig rank, start) is used. `_refound` ids and
  ids without a GFF row do not place the cluster. Counts of fragmented, multi-copy and
  refound cells are recorded.
- **Ordering genome (auto):** candidates are cohort samples present in the alignment,
  ranked by contiguity from `AssemblyQualityWorkflow._fasta_metrics` on every genome
  (all contigs; QUAST's `QuastContigs` ignores contigs under 500 bp, so it is not used):
  fewest contigs, then highest N50, then name. Candidates are checked in that order and
  the first with placed core genes >= `phylo.ordering_min_placed` (default 0.95) wins; if
  none of the first 20 pass, the best placed of them wins with warning
  `ordering_many_unplaced`. `ordering.json` records the evaluated candidates (at most 20:
  name, contigs, N50, placed, total).
- **Override `phylo.order_reference`:**
  - a cohort sample name: use that sample (error at run start if unknown);
  - a path to an external genome FASTA: map each core cluster's representative
    (`pan_genome_reference.fa`, else the first ungapped sequence of
    `aligned_gene_sequences/<cluster>.aln.fas`) with `minimap2 -c -x asm20
    --secondary=no` (phylogeny env); keep the best hit with query coverage >= 0.8; strand
    from the hit.
- **Order and orientation:** contigs longest first (`##sequence-region`, then name),
  genes by start. A gene on the ordering genome's minus strand is reverse-complemented
  (A/T, C/G swapped, `N` and `-` kept; neutral for the GTR family), so adjacent genes are
  contiguous in genome direction and a tract crossing two genes on opposite strands stays
  one column range. Unplaced core genes are appended at the end in Panaroo order, forward,
  and counted. Intergenic regions are absent: Gubbins window sizes (100 to 10,000) are in
  concatenation coordinates, not genome coordinates (documented in `ordering.json` notes).
- **Writing:** stream one sequence at a time, slice blocks in the new order, sanitize
  (uppercase, `ACGTN-` kept, anything else to N, as Gubbins staging does now). Memory
  O(alignment length).
- **Diversity pass (same stream):** a numpy `4 x L` int32 count array of A/C/G/T per
  column (about 42 MB on the reference run). Mean pairwise p-distance = mean over columns
  with >= 2 called bases of 1 - sum c_a (c_a - 1) / (m (m - 1)). Variable column indices
  are kept; max pairwise p-distance = max over pairs of differences on variable columns
  (both called) divided by the number of columns with >= 2 called bases, computed with
  numpy on the `n x V` variable-site matrix (O(n^2 V), no genome cap). Per-taxon gap
  fraction is recorded, so the taxa Gubbins will drop (above `--filter-percentage`) are
  predicted and logged here.
- **Outputs:**
  - `core_alignment/core_gene_alignment.ordered.aln`
  - `core_alignment/gene_map.tsv`, one row per gene in the new order: `order`, `cluster`,
    `aln_start`, `aln_end` (1-based inclusive, ordered alignment), `length`, `placed`,
    `reversed` (true when reverse-complemented), `ordering_contig`, `ordering_start`,
    `ordering_end`, `ordering_strand`, `ordering_locus_tag`, `copies_in_ordering_genome`,
    `fragmented`, `in_panaroo_filtered`, `panaroo_aln_start`, `panaroo_aln_end`.
  - `core_alignment/ordering.json` (`schema_version: 1`, `algorithm_version`): method
    (`auto_most_contiguous` | `config_sample` | `external_genome`), ordering genome (name or
    path + sha256), its contigs and N50, candidate table, genes total / placed / appended /
    multi-copy / fragmented / refound / reversed, inputs (see 3.7 stamps), sanitized
    character count, sequence count, alignment length, output sha256, diversity (mean,
    max, variable sites, per-taxon gap fraction above the Gubbins filter), notes.
- `map_alignment_interval(gene_map, start, end)` returns the genes a block overlaps with
  gene-relative and ordering-genome coordinates (inverted for reversed genes); used for
  Gubbins blocks (3.2) and fastGEAR recombination coordinates.

### 3.2 Gubbins on the ordered alignment

- Input: `core_alignment/core_gene_alignment.ordered.aln` directly (no `gubbins/input`
  copy). Prefix stays `core_gene` (manuscript scripts read
  `gubbins/core_gene.final_tree.tre`). Gubbins 3.4.3 in the container.
- Command adds `--seed <phylo.seed>` (default 1), `--iterations <phylo.gubbins_iterations>`
  (default 20), `--converge-method <phylo.gubbins_converge_method>` (default
  `weighted_robinson_foulds`; `robinson_foulds` is the documented, less strict
  alternative), `--filter-percentage <phylo.gubbins_filter_percentage>` (default 25),
  `--tree-builder <phylo.gubbins_tree_builder>` (default `raxml`), `--outgroup` when
  `phylo.outgroup` is set, `--threads` from the thread plan. `--invariant-site-correction`
  stays off (Gubbins' default) and is recorded. Convergence needs an exact repeat of a
  tree, which a fixed seed makes reachable only at a fixed thread count (RAxML pthreads
  results can change with `-T`), so threads are part of the params.
- Run into a temp folder inside `gubbins/`, then replace the folder content (no partial
  outputs on failure); the stale `tmp*` folder goes with it.
- **Failure** (tool exit non-zero, missing outputs): the step fails and the run fails,
  as for every other tool. There is no silent fallback to an unfiltered tree.
- **Convergence parse** from stdout: `Convergence after N iterations: Tree observed
  before.` or `...: Recombinations observed before.` gives `converged: true,
  iterations_run: N, convergence_reason`; `Maximum number of iterations (N) reached.`
  gives `converged: false`; neither gives `converged: null` (treated as not converged).
- **Non-convergence:** step status success; `logger.warning`; warning
  `gubbins_not_converged` in step info, run-level `warnings` (3.7) and `trees.json`; the
  GUI shows a chip.
- **Diversity warnings (before the call):** `gubbins_high_diversity` when mean >
  `phylo.gubbins_diversity_warn` (0.01) or max > `phylo.gubbins_max_diversity_warn`
  (0.02). Reference run: 0.0030 / 0.0153, no warning.
- **Lineage structure warning (after the call):** from `core_gene.per_branch_statistics.csv`,
  share of `Number of SNPs Outside Recombinations` on the single branch with the most of
  them; `gubbins_lineage_structure` when above `phylo.gubbins_longest_branch_warn` (0.5).
  Reference run: about 0.81 (one branch separating a 6-genome lineage), so it warns. The
  message says Gubbins is designed for closely related isolates and advises running it per
  lineage (fastGEAR or PopPUNK clusters) for such cohorts. Never a refusal.
- **Left-out samples:** taxa in the input but not in the final tree are
  `samples_left_out`; warning `gubbins_samples_left_out` (count and names).
- **Skip below `phylo.gubbins_min_genomes` (default 4):** no tool call; `summary.json` with
  `ran: false`, `skip_reason`; step success with `skipped: true`; warning
  `gubbins_skipped_few_genomes`. Validators accept a skipped summary without other files.
- **Outgroup:** names must be cohort samples (run start error otherwise). If the ordering
  pass predicts Gubbins will drop an outgroup taxon, Gubbins runs without `--outgroup` and
  warning `outgroup_unusable`.
- **summary.json** (schema_version 2): tool version, full command, input path + sha256,
  params (seed, iterations, converge method, filter percentage, tree builder, outgroup,
  threads, invariant site correction), `algorithm_version`, `ran`, `skip_reason`,
  `converged`, `iterations_run`, `max_iterations`, `convergence_reason`, diversity
  (mean, max, variable sites, thresholds, warned), `longest_branch_share` (node, share),
  `samples_left_out`, block count, SNPs inside/outside recombination, r/m (SNP ratio),
  per-tip recombination fraction (`Bases in Recombinations Excluding Gaps` / `Genome
  Length`) with the top 5 outliers, `warnings`, outputs.
- **Block to gene table:** `gubbins/recombination_blocks_genes.tsv` (block start/end, taxa
  count, genes overlapped, first/last gene ordering-genome coordinates, any gene
  `in_panaroo_filtered`) from the recombination GFF and `gene_map.tsv`.

### 3.3 Core tree with IQ-TREE (new step `core_tree`)

New module `src/pipeline/CoreTreeWorkflow.py`; `_run_raxml_tree` stays only for Snippy;
`_resolve_core_tree_alignment` and the line 1699 log message are deleted (P09).

- **Input:** `gubbins/core_gene.filtered_polymorphic_sites.fasta`, used as is.
- **Invariant-site correction** `phylo.iqtree_invariant_correction`:
  - `fconst` (default, Gubbins' own IQ-TREE route): `mask_gubbins_aln.py --aln <ordered
    aln> --gff gubbins/core_gene.recombination_predictions.gff --out <tmp>` (gubbins env),
    count A/C/G/T of the columns monomorphic among called bases, delete the masked file,
    pass `-fconst a,c,g,t` with `-m <phylo.iqtree_model>`. Branch lengths are
    substitutions per core site (comparable with Snippy). Check logged: polymorphic +
    constant + uncalled columns = alignment length.
  - `asc`: columns with fewer than two states among called bases are dropped and counted,
    written to `trees/core_recomb_iqtree/input.variable_sites.fasta`, `-m
    <iqtree_model>+ASC`. Branch lengths are still per site: the +ASC model implies the
    invariant sites (checked on a 16-genome subset; "per variable site" was wrong).
- **Command:** `iqtree2 -s <input> [-fconst ...] -m MFP -B 1000 --bnni --seed <phylo.seed>
  -T <threads> --prefix trees/core_recomb_iqtree/core_recomb`. Threads from the thread
  plan (new key `core_tree`: min(total, 16)), never `-T AUTO`.
- **Fast mode** above `phylo.iqtree_fast_above_genomes` (default 1000): `--fast --alrt
  <iqtree_bootstraps>`, no `-B`/`--bnni` (IQ-TREE refuses UFBoot with `--fast`); support
  type `sh_alrt`.
- **Few genomes:** 3 genomes: Gubbins was skipped (n < 4), IQ-TREE on the ordered
  alignment with `-m MFP`, no `-B`, id `core`, `recombination_filtered: false`, warning
  `core_tree_no_support_few_genomes`. 1 or 2 genomes: trivial tree
  (`_write_trivial_tree`), same id and warning. 4 or more with Gubbins skipped cannot
  happen with defaults; if configured, same `core` tree with bootstrap.
- **No polymorphic sites:** star tree, note `no_variable_sites`, no IQ-TREE call.
- **Failure:** the step fails and the run fails.
- **Kept outputs** in `trees/core_recomb_iqtree/`: `.iqtree`, `.log`, `.treefile`
  (unrooted, support labels), `.contree`, `.model.gz`, `run.json` (stamp, 3.7);
  checkpoint and `.splits.nex` removed.
- **Parsed** from `.iqtree`/`.log`: best model ("Best-fit model: ... chosen according to
  BIC"), log-likelihood, UFBoot convergence (warning `core_tree_ufboot_not_converged` on
  "bootstrap analysis did not converge"), sequences failing the composition chi2 test
  (count and names), identical sequences collapsed (count; support among them is
  meaningless, noted). RF and normalized RF to the Gubbins final tree are put in the
  entry's notes (recombination was reconstructed on the Gubbins tree, not this one).
- `iqtree2` and `mask_gubbins_aln.py` added to `_tool_env_key` (`gubbins` env) and
  `iqtree2` to `_collect_versions`.

### 3.4 fastGEAR input

- Reads `core_alignment/core_gene_alignment.ordered.aln` directly (the run folder is
  bind-mounted at its own path); the `fastgear/input/` copy is deleted.
- Implementation checks: how fastGEAR's MATLAB reader treats `N` versus `-` (if `N` is
  not read as missing, stage a copy with `N` written as `-` and keep it in
  `fastgear/input/`); that fastGEAR writes nothing beside its input.
- `fastgear/input.json`: input path + sha256, fastGEAR settings, image digest, specs
  file sha256, `algorithm_version`. Reuse needs all to match plus valid outputs; runs
  without the stamp rerun.
- Recombination coordinates refer to the ordered alignment; `map_alignment_interval`
  maps them to genes. Expect lineages and clusters nearly unchanged (BAPS on SNPs is
  order-independent) and recombination segments to change.

### 3.5 Rooting (cheap post step, from raw unrooted trees)

- Rule: `phylo.outgroup` set: root on the edge above the outgroup's MRCA (single sample or
  comma list); if not monophyletic, or an outgroup taxon is not in the tree, midpoint
  with warning `outgroup_not_monophyletic` / `outgroup_unusable`. Otherwise midpoint
  (longest tip to tip path).
- Raw trees are kept unrooted: IQ-TREE `trees/core_recomb_iqtree/core_recomb.treefile`
  (or `core.treefile`), Gubbins `gubbins/core_gene.final_tree.tre`, Snippy
  `trees/raw/snippy_core.raxml.nwk` (RAxML-NG best tree; an old run's
  `trees/snippy_raxml.tree` is copied there and left in place for the scripts). `write_trees_index(run_dir, cfg)` roots each
  raw tree by the rule, renames the Snippy `Reference` tip to the reference stem
  (`DataInput` forces it to be a cohort genome), writes `trees/<id>.nwk` and rebuilds
  `trees/trees.json` from the stamps (3.6). It runs at the end of `core_tree` and of
  `phylogeny` (and on every resume). Changing `phylo.outgroup` reruns Gubbins (its
  `--outgroup` param) and then the core tree; the rooting itself is this cheap step.
- Code: Qt-free `src/pipeline/Newick.py` takes `TreeNode`, `NewickParser`, `prune_to` and
  `midpoint_root` from `phylo.py` (moved; `phylo.py` imports them back), adds a writer,
  `root_on_outgroup`, and support handling: numeric internal labels become
  `TreeNode.support`; rerooting moves each value with its edge. Check on real data: the
  (split, support) set is identical before and after rooting.
- **Names:** `DataInput` rejects at run start sample stems containing whitespace or any of
  `( ) [ ] , : ; ' "` (IQ-TREE and Newick rewrite them). Before a tree is listed, its tip
  set must equal the expected samples minus `samples_missing`; otherwise it is left out
  with warning `tree_tip_mismatch`.

### 3.6 Tree files, stamps and `trees/trees.json`

| id | file | built from | default |
|---|---|---|---|
| `core_recomb` | `trees/core_recomb.nwk` | IQ-TREE on Gubbins filtered polymorphic sites | yes |
| `core` | `trees/core.nwk` | IQ-TREE on the ordered alignment, only when Gubbins was skipped | yes, when present |
| `gubbins_final` | `trees/gubbins_final.nwk` | Gubbins final tree (RAxML inside Gubbins, no support), re-rooted by the rule | no |
| `snippy_core` | `trees/snippy_core.nwk` | Snippy core.full.aln, RAxML-NG GTR+G, `Reference` tip renamed | no |

- `trees/core_gene_tree.nwk` is no longer produced (moved to the backup on the reference
  run, section 5).
- **Stamps:** each tree has `trees/<id>.json` written by its producer: `id`, `input_path`,
  `input_sha256`, `params`, `tool_version`, `threads`, `seed`, `algorithm_version`,
  `raw_file`, and the producer-known fields of its `trees.json` entry. `trees.json` is
  built only from stamps whose `input_sha256` matches the current upstream file
  (`core_recomb`/`core`: the IQ-TREE input; `gubbins_final`: the ordered alignment
  sha in `gubbins/summary.json`; `snippy_core`: `snippy/core/core.full.aln`); a
  non-matching stamp is left out with warning `tree_stale`. Old runs reused without a
  Snippy stamp get one written from the existing files (no Snippy call).
- **Invalidation at rerun start:** a `core_alignment` rerun deletes `trees/core_recomb.*`,
  `trees/core.*`, `trees/gubbins_final.*`; a Gubbins rerun deletes `trees/core_recomb.*`,
  `trees/gubbins_final.*`. So a downstream failure never leaves an old tree beside new
  upstream results.

### 3.7 Steps, freshness and run_state

Step order: panaroo, **core_alignment**, gubbins, **core_tree**, phylogeny (Snippy only),
fastgear_clades. `core_alignment` and `core_tree` are in `RESULT_OPTIONAL_STEPS`.

Every stamp holds an `algorithm_version` constant per module, bumped on any logic change.
Inputs are fingerprinted by sha256.

| Step | Reused when (all must match, plus outputs validate) | Otherwise |
|---|---|---|
| core_alignment | sha256 of source alignment, EMBL header (and filtered header), `gene_presence_absence.csv`, ordering genome GFF3, external genome; `phylo.core_alignment`, `phylo.order_reference`, `phylo.ordering_min_placed`; resolved ordering genome; `algorithm_version`; output sha256 equals the file | rebuild (about 1 min) |
| gubbins | input sha256 == `ordering.json` output sha256; params incl. threads; Gubbins version; `algorithm_version` (old summaries without sha: rerun) | rerun |
| core_tree | `trees/core_recomb_iqtree/run.json`: input sha256 (Gubbins filtered sites, or ordered alignment for `core`), recombination GFF sha256 (fconst), seed, model, bootstraps, bnni, fast mode, invariant correction, threads, IQ-TREE version, `algorithm_version` | rerun |
| phylogeny (Snippy) | Snippy tree validates (raw new or legacy name); no Gubbins gate | rerun Snippy |
| fastgear_clades | `fastgear/input.json` (3.4) | rerun |

- `gubbins_regenerated` gating and the line 1699 message are deleted; each step logs the
  real reason ("input alignment changed", "parameters changed", "no input stamp").
- `run_state.json` / `run_report.json`: new top-level `warnings` list of
  `{step, code, message}`. Codes: `ordering_many_unplaced`, `gubbins_high_diversity`,
  `gubbins_lineage_structure`, `gubbins_not_converged`, `gubbins_samples_left_out`,
  `gubbins_skipped_few_genomes`, `outgroup_unusable`, `outgroup_not_monophyletic`,
  `core_tree_no_support_few_genomes`, `core_tree_ufboot_not_converged`,
  `tree_tip_mismatch`, `tree_stale`.
- Step info: core_alignment (ordering_method, ordering_genome, genes_total, genes_placed,
  genes_appended, genes_reversed, multi_copy, alignment, sha256, source_alignment,
  diversity), gubbins (converged, iterations_run, max_iterations, seed, threads,
  converge_method, outgroup, diversity, longest_branch_share, skipped, samples_left_out,
  input_sha256, warnings), core_tree (tree_id, file, model, support_type, replicates,
  fast_mode, invariant_correction, constant_counts, seed, threads, sites, ufboot_converged,
  composition_failed, identical_collapsed, rf_to_gubbins), rooting (method, outgroup).
- `run_config.json` `inputs.phylo` snapshots the effective settings; `thread_plan.core_tree`
  added. Validators: `validate_phylogeny_outputs` checks the Snippy tree only (new or
  legacy name); new `validate_core_alignment_outputs`, `validate_core_tree_outputs`;
  `validate_gubbins_outputs` accepts a skipped summary. `expected_outputs.json`,
  `run_expectations.py` and `classify_expected_path` follow the new names.

### 3.8 Config keys (`bactolomew_config.json` section `phylo`, `ToolEnv.configured_phylo_settings`)

| Key | Default | Meaning |
|---|---|---|
| `core_alignment` | `"full"` | Panaroo alignment used by every core consumer (`full` or `filtered`) |
| `order_reference` | `""` | empty = auto; a cohort sample name; or a genome FASTA path |
| `ordering_min_placed` | `0.95` | auto ordering genome must place this share of core genes |
| `outgroup` | `""` | sample name or comma list; empty = midpoint |
| `seed` | `1` | Gubbins and IQ-TREE seed |
| `gubbins_iterations` | `20` | Gubbins maximum iterations |
| `gubbins_converge_method` | `"weighted_robinson_foulds"` | Gubbins `--converge-method` |
| `gubbins_filter_percentage` | `25` | Gubbins `--filter-percentage` |
| `gubbins_tree_builder` | `"raxml"` | Gubbins `--tree-builder` |
| `gubbins_min_genomes` | `4` | skip Gubbins below this |
| `gubbins_diversity_warn` | `0.01` | warn above this mean pairwise core p-distance |
| `gubbins_max_diversity_warn` | `0.02` | warn above this max pairwise core p-distance |
| `gubbins_longest_branch_warn` | `0.5` | warn above this share of clonal-frame SNPs on one branch |
| `iqtree_model` | `"MFP"` | ModelFinder (`+ASC` appended only in `asc` mode) |
| `iqtree_invariant_correction` | `"fconst"` | `fconst` or `asc` |
| `iqtree_bootstraps` | `1000` | UFBoot replicates (SH-aLRT replicates in fast mode) |
| `iqtree_bnni` | `true` | `--bnni` with UFBoot |
| `iqtree_fast_above_genomes` | `1000` | `--fast` with SH-aLRT above this |

### 3.9 Edge cases

- 1 or 2 genomes: no Gubbins, trivial `core` tree, notes. 3 genomes: no Gubbins, IQ-TREE
  on the ordered alignment without bootstrap (`core`).
- High diversity or one dominant branch: warnings only, never a refusal.
- Ordering genome missing core genes: placeability rule, appended and counted, warning.
- Gubbins drops gappy taxa: predicted in the ordering pass, `samples_left_out` +
  warning, `trees.json` `samples_missing`; the GUI header says "N genomes not in this
  tree" and greys or omits them as now.
- Huge cohorts: streaming ordering, numpy diversity on variable columns, IQ-TREE fast mode
  above 1000; Gubbins runtime logged with the genome count (settings not changed).
- Old runs (no `core_alignment/`, no `trees.json`): validate unchanged; a resume reruns
  core_alignment, Gubbins, core_tree and fastGEAR, never Snippy.

### 3.10 Tests (very few, critical concepts only)

1. Ordering: synthetic alignment, EMBL header, presence table (with a `;` cell) and GFF3
   with a minus-strand gene; ordered columns follow genome order, forward slices equal the
   Panaroo block, the reversed slice equals its reverse complement,
   `map_alignment_interval` returns the right gene and coordinates, unplaced genes are
   appended and counted.
2. Tree choice: with `trees.json` the GUI default is the recombination-filtered tree and
   a linked Gene presence/absence follows the Phylogeny source; without it (old run) the
   Gubbins final tree wins over Snippy and `core_gene_tree.nwk`.
3. Non-convergence: Gubbins stdout with "Maximum number of iterations (20) reached." gives
   `converged: false` and a `gubbins_not_converged` warning.
Existing tests touching old names (`tests/_pan_fixture.py`, `test_pan_panels.py` 146,
`test_row_link.py` 105/197/287, `test_tree_view.py` 836-877) are adjusted, not extended.
Everything else is checked on the reference run (section 6).

## 4. Contract between the pipeline lane and the GUI lane

Paths (relative to the run folder):

| Path | Producer | GUI reads |
|---|---|---|
| `trees/trees.json` | pipeline (`write_trees_index`) | yes: the only source of tree choice when present |
| `trees/<id>.nwk` | pipeline | yes, via `trees.json` `file` (relative to `trees/`) |
| `trees/<id>.json` (stamps), `trees/raw/`, `trees/core_recomb_iqtree/` | pipeline | no |
| `core_alignment/ordering.json`, `gene_map.tsv` | pipeline | no (not in this rework) |
| `gubbins/summary.json` | pipeline | no (its fields are copied into `trees.json`) |
| `run_state.json` / `run_report.json` `warnings` | pipeline | optional |

`trees.json`, `schema_version` 1. Every key below is always present; unknown values are
`null`, empty lists are `[]`. The GUI ignores unknown keys.

```json
{
  "schema_version": 1,
  "generated_at": "2026-10-10T18:00:00+03:00",
  "default_tree": "core_recomb",
  "rooting": {"method": "midpoint", "outgroup": null, "note": null},
  "ordering_genome": "GCA_003931755.1",
  "gubbins": {
    "ran": true,
    "converged": true,
    "iterations_run": 7,
    "max_iterations": 20,
    "converge_method": "weighted_robinson_foulds",
    "samples_left_out": []
  },
  "warnings": [
    {"step": "gubbins", "code": "gubbins_lineage_structure", "message": "..."}
  ],
  "trees": [
    {
      "id": "core_recomb",
      "label": "Core genome, recombination-filtered",
      "file": "core_recomb.nwk",
      "method": "IQ-TREE 2.4.0, ModelFinder (BIC), ultrafast bootstrap",
      "tool_version": "2.4.0",
      "alignment": "gubbins/core_gene.filtered_polymorphic_sites.fasta",
      "alignment_kind": "core_variable_sites",
      "sites": 28672,
      "invariant_correction": "fconst",
      "model": "GTR+F+R4",
      "log_likelihood": -123456.7,
      "support": {"type": "ufboot", "replicates": 1000, "threshold": 95},
      "rooting": {"method": "midpoint", "outgroup": null, "note": null},
      "recombination_filtered": true,
      "ordering_genome": "GCA_003931755.1",
      "branch_length_unit": "substitutions/site",
      "samples": 153,
      "samples_missing": [],
      "reference_sample": null,
      "seed": 1,
      "threads": 16,
      "fast_mode": false,
      "notes": []
    }
  ]
}
```

Field values
- `id`: `core_recomb` | `core` | `gubbins_final` | `snippy_core`. At most one entry per id;
  `trees` is in display order (default first).
- `label` (UI text, US spelling; the GUI shows it verbatim and never invents labels for
  `trees.json` runs):
  - `core_recomb`: "Core genome, recombination-filtered"
  - `core`: "Core genome (Gubbins skipped)"
  - `gubbins_final`: "Gubbins final tree, recombination-filtered (no support)"
  - `snippy_core`: "Whole-genome SNPs (Snippy), not recombination-corrected"
- `alignment_kind`: `core_variable_sites` (core_recomb) | `core_alignment` (core) |
  `core_snps` (gubbins_final) | `snippy_core_full` (snippy_core).
- `invariant_correction`: `fconst` | `asc` | `null` (trees not built by IQ-TREE on
  variable sites).
- `support.type`: `ufboot` (threshold 95) | `sh_alrt` (threshold 80) | `none` (threshold
  `null`, replicates 0). Support values are numeric internal node labels in the `.nwk`
  (percent, one number per node); leaves are sample names; no `Reference` tip.
- `rooting.method`: `midpoint` | `outgroup`; `rooting.outgroup`: the configured outgroup as a
  list of sample names, or `null` (also kept when it fell back to midpoint); the GUI joins it
  with ", "; `note` says why an outgroup fell back.
- `branch_length_unit`: `"substitutions/site"` (core_recomb with fconst, core,
  snippy_core) | `"substitutions/site (+ASC model)"` (core_recomb with asc) | `"SNPs
  (reconstructed point mutations)"` (gubbins_final). The GUI scale bar shows it verbatim.
- `samples`: tip count; `samples_missing`: cohort samples absent from the tree.
- `reference_sample`: the Snippy reference stem for `snippy_core`, else `null`.
- top-level `gubbins`: `null` when there is no Gubbins summary; `ran: false` when skipped.
  `converged: false` or `null` with `ran: true` means not converged.
- top-level `warnings`: the run warnings whose step is core_alignment, gubbins, core_tree
  or phylogeny.

GUI rules
- Default tree: `default_tree` if it names a listed entry whose file exists; else the
  first entry with `recombination_filtered: true`; else the first entry.
- Every listed file is rooted as its `rooting` says; for `trees.json` trees the GUI root
  mode defaults to "As in file". A small Root dropdown also offers "Midpoint" and
  "fastGEAR clades" (`ROOT_GROUPS`, so `group_root` keeps a caller).
- Header: the source label, "Root: midpoint" / "Root: outgroup X" / "Root: fastGEAR
  clades", "N genomes not in this tree" when `samples_missing` is non-empty, and a "Gubbins
  did not converge" warning chip (tooltip with iterations) when `gubbins.ran` and
  `gubbins.converged` is not `true` and the shown tree is `recombination_filtered`.
- Support: only in the hover text of an internal node or branch ("UFBoot 97, well
  supported at >= 95"), only when `support.type` is not `none`; never drawn by default;
  numeric internal labels are never drawn as names.
- Source selection: the Phylogeny panel has a source dropdown with the `trees.json`
  labels; the choice is saved in the panel state (`save_state` / `restore_state`), not in
  layouts; a restored id that is not listed falls back to the default. A Gene
  presence/absence panel linked through `TreeRowLink` uses its leading tree's source; an
  unlinked one uses the default. `_reference_leaf_order` memo keys include the source id.
  Every other tree consumer (sample tree, exports) uses the default.
- Old runs without `trees.json` (GUI fallback, in this order): `gubbins/core_gene.final_tree.tre`
  ("Gubbins final tree, recombination-filtered", unit "SNPs (reconstructed point
  mutations)"), `trees/snippy_raxml.tree` ("Snippy SNPs, not recombination-corrected",
  `Reference` leaf mapped as today), `trees/core_gene_tree.nwk` ("Core genome, not
  recombination-corrected"), Panaroo's tree; unit "substitutions/site" for the RAxML
  files. Legacy files are unrooted: the GUI roots them at midpoint and says "Root:
  midpoint".
- GUI lane files: `model/pan_data.py` (one `tree_sources()` from `trees.json` or the
  fallback, replacing `resolve_raxml_tree_path` / `raxml_tree_source` /
  `tree_source_candidates`); `panaroo_core_accessory.py` (dropdown, header, chip, unit,
  memo key); `panels_controller.py` (default root "As in file", Root dropdown);
  `panels/row_link.py` (linked heatmap follows the tree's source); `tree_view.py`
  (support hover, no label drawing). Avoid the names retired in `test_retired_shims.py`
  (`set_tree_source`, `tree_source_key`, `tree_source_label`, `tree_source_path`).
  Update the CLAUDE.md Phase 9 sentence about `ROOT_GROUPS`.

## 5. Rerun plan for the reference run

Run only after the code, the review fixes and the targeted tests are done.

1. **Guards:** no other pipeline driver running on this run; the leftover cgMLST schema
   build (pid 3042661, `DB/cgmlst`) is untouched and chewBBACA must be "reuse" in the plan.
2. **Backup** to `/home/bals/Documents/dev/BACTolomew_backups/2026-10-10_phylogeny/`
   (about 1 GB): `gubbins/` (547 MB, including `input/` and the stale `tmp1tj7krj1/`),
   `fastgear/` (417 MB), `trees/`, `snippy/reference.txt`, `run_config.json`,
   `run_state.json`, `run_report.json`, `expected_outputs.json`, and
   `cache/run_78eaa9126850/{tables_meta.json,ui_cache.json}`. Plus a manifest of every
   file in the run folder (path, size, mtime_ns) for the after check.
3. **Prepare:** move `trees/core_gene_tree.nwk` into the backup (no longer produced).
   `trees/snippy_raxml.tree` is copied to `trees/raw/` by the code (3.5) and stays in place. The missing input
   stamps make Gubbins and fastGEAR rerun.
4. **Driver** `<scratch>/phylo/rerun/driver.py plan|run`, reusing `build_cfg()` from the
   P01 to P06 driver (`<scratch>/rerun/rerun_driver.py`: GUI Resume `DataInput`, COG
   dropped, 96 cores). `plan` prints each step from the Runner's own reuse predicates and
   refuses unless exactly: core_alignment RERUN, gubbins RERUN, core_tree RERUN,
   fastgear_clades RERUN; assembly_qc, bakta, mlst, amrfinder_mutations, mobsuite,
   panaroo, phylogeny (Snippy: stamp and rooting only), chewbbaca, diamond_markers,
   diamond_reference_loci reuse; collect_annotations re-merges in memory without a
   DIAMOND search.
5. **Run** in the background with a log. Expected wall time: core_alignment about 1 to
   2 min; Gubbins 25 to 85 min; mask + IQ-TREE 15 to 60 min (28,672 sites + fconst, 153
   taxa, 16 threads, UFBoot with bnni); fastGEAR about 2 h. Total about 3 to 4.5 h.
6. **After:**
   - Manifest diff: only `core_alignment/`, `gubbins/`, `trees/`, `fastgear/`,
     `run_config.json`, `run_state.json`, `run_report.json`, `expected_outputs.json` may
     change; `snippy/`, `panaroo/`, `bakta/`, `chewbbaca/`, `annotations/`, `dmnd/` must
     be byte-identical (size and mtime_ns).
   - `validate_completed_run` 0 errors; `run_state.warnings` reviewed (expected:
     `gubbins_lineage_structure`; `gubbins_not_converged` only if it did not converge).
   - GUI cache rebuild (`<scratch>/rerun/rebuild_cache.py`): `run_config.json` is part
     of the cache fingerprint.
   - App check (`.claude/skills/bactolomew-app-check`): Phylogeny and Gene
     presence/absence screenshots with the new default tree, header label, root label,
     source dropdown, support hover, a switch to the Snippy tree with the linked heatmap
     following.

## 6. Verification and before/after metrics

Baseline values come from the backup; script in `<scratch>/phylo/metrics/`.

- **Ordering:** ordering genome and candidate table, genes placed / appended / reversed /
  multi-copy / fragmented (expected GCA_003931755.1, about 2610 placed of 2676); a spot
  check of 5 genes (forward and reversed): ordered slice equals the Panaroo block (or its
  reverse complement), ordering-genome coordinates match the GFF3.
- **Gubbins:** converged and iterations; recombination blocks (4129 before), SNPs
  inside/outside and r/m (2.70 before), per-branch r/m and rho/theta for the main clades;
  longest-branch share of clonal-frame SNPs (about 0.81 before) and the filtered sites
  variable within each side of that branch; per-tip recombination fraction outliers;
  blocks spanning more than one gene and blocks spanning genes not adjacent in the
  ordering genome (the P39 artefact: should fall sharply); hotspot genes, with and without
  `in_panaroo_filtered` genes (r/m both ways); diversity (0.0030 / 0.0153 expected).
- **Trees:** RF and normalized RF between the new `core_recomb` and: old
  `core_gene_tree.nwk`, old and new Gubbins final trees, Snippy tree; splits with UFBoot
  >= 95 that conflict with the old trees; UFBoot distribution (median, share >= 95);
  UFBoot convergence, composition test failures, identical sequences collapsed; model
  and log-likelihood; root position versus the old clade root.
- **Rooting/support integrity:** the (split, support) set of the `.treefile` equals that
  of `trees/core_recomb.nwk`.
- **fastGEAR:** lineage and cluster membership before/after (contingency tables, adjusted
  Rand index for lineages and clusters separately; before: 2 lineages of 6 and 147, 4
  clusters), then counts of recent and ancestral recombination events.
- **GUI clade grouping:** whether each fastGEAR lineage and cluster is monophyletic in the
  new tree under midpoint rooting; split clade bands; Gene presence/absence row order
  change (rank correlation); screenshots before (backup state) and after.
- **Reproducibility:** seeds and thread counts recorded. Multithreaded IQ-TREE with a fixed
  seed is not bit-reproducible (scratch check: same topology, log-likelihood differing in
  the third decimal, UFBoot median 96 vs 97), so compare topologies (RF), not files.

## 7. Manuscript scripts in `src/scripts` that read changed outputs (list only)

- `plot_st492_raxml_tree_metadata.R`: `gubbins/core_gene.final_tree.tre`,
  `fastgear/sample_clades.tsv`, `fastgear/output/recombinations_recent.txt`; roots on
  fastGEAR lineage 1.
- `plot_st492_tree_pangenome_phandango.R`: Gubbins final tree, `fastgear/sample_clades.tsv`.
- `plot_st492_tree_resistance_phandango.R`: Gubbins final tree, `fastgear/sample_clades.tsv`.
- `plot_st492_tree_virulence_phandango.R`: Gubbins final tree, `fastgear/sample_clades.tsv`.
- `build_supplementary_raw_outputs.py`: Gubbins files, `trees/core_gene_tree.nwk`,
  `trees/snippy_raxml.tree`, fastGEAR outputs (old names disappear).
- `build_st492_cgmlst_mst.py`: `fastgear/sample_clades.tsv`.
- `tune_cgmlst_force.py`: fastGEAR clades through the GUI loaders.
- Not affected (hard-coded paths outside the run): `correct_nwk.py`, `Pangenome.R`.

## Rejected critique points

- Deepest-split distance before Gubbins: replaced by the post-Gubbins longest-branch share
  (3.2), which measures the same structure without building a tree first; per-lineage
  site counts are a section 6 metric, not a pipeline output.
- `--alrt 1000` beside UFBoot: it would make internal labels two numbers ("80/95") and
  complicate the contract; UFBoot with `--bnni` is enough. SH-aLRT is used in fast mode.
- Config key for Gubbins `--invariant-site-correction`: kept at Gubbins' default and
  recorded; the `gubbins_final` tree is labeled with its SNP unit instead.
- Rotating the largest contig to start at dnaA: optional, affects one boundary per
  genome; not now.

## Implementation notes

### Pipeline lane fixes after review (2026-10-10)

- `rooting.outgroup` in `trees.json` is a list of names or null (contract updated); the GUI
  joins it with ", " (`model/pan_data._names_text`). `phylo.outgroup` given as a JSON list
  is accepted and joined with commas (`ToolEnv.configured_phylo_settings`).
- Sample names may only hold letters, digits and `_ . - |` (`DataInput.NEWICK_UNSAFE_NAME`):
  IQ-TREE rewrites every other character, which used to drop the default tree silently.
- Freshness: `ordering.json` records run-relative paths and the reuse check compares sha256
  only, so a moved run folder is reused. Trees are deleted after a `core_alignment` rerun
  only when the ordered alignment's sha256 changed (checked on a scratch run: identical
  output kept the trees), and after Gubbins only when Gubbins fails; otherwise the core
  tree and `trees.json` follow the Gubbins output's sha256. Thread counts are recorded
  but no longer part of the Gubbins and core tree reuse keys.
- Multithreaded IQ-TREE with a fixed seed is not bit-reproducible (same topology, last
  digits of the likelihood and support differ); documented in `CoreTreeWorkflow` and 6.
- Gubbins (`ALGORITHM_VERSION` 2): the lineage-structure share merges the root's two
  branches (from `node_labelled.final_tree.tre`) into one edge; predicted dropouts are
  recomputed from `ordering.json` `missing_fraction` (every sample) with the current
  `gubbins_filter_percentage`; blocks touching genes appended at the end are flagged
  (`touches_appended_gene`) and left out of `blocks_spanning_several_genes`, counted in
  `blocks_touching_appended_genes`; `genome_span` shows how far apart a block's first and
  last genes sit in the ordering genome. Placing the appended genes with the next
  candidate genome is not done.
- A star tree (no variable sites) is stamped `support.type: none`, method "star tree".
- Snippy with fewer than 4 sequences (a 3-genome cohort) writes the only possible tree
  (3 taxa: branch lengths from pairwise p-distances) instead of calling RAxML-NG, with
  warning `snippy_tree_few_genomes`; checked with the real Runner on the 3-genome scratch run.
- `trees.json` warnings are replaced on each write, so a tree stale before the Snippy step
  and rebuilt by it leaves no warning behind.
- The legacy `trees/snippy_raxml.tree` is copied to `trees/raw/`, not moved, so the
  manuscript scripts keep their input. The `TreesIndex` docstring and 3.5 now say that
  changing `phylo.outgroup` reruns Gubbins and the core tree.
- Dead helpers removed: `column_base_counts`, `gubbins_ran`, `read_trees_index`,
  `is_numeric_label`.

### GUI lane fixes after review (2026-10-10)
- Outgroup display: `rooting.outgroup` is a list in `trees.json`; `pan_data._names_text`
  joins it with ", " ("Root: outgroup A, B"). Already in place when reviewed.
- Layouts: the Phylogeny panel's `source` and `root` params pick run data, so
  `LayoutStore` leaves them out of the default layout, user layouts and the last window
  (`layout_store.RUN_DATA_PARAMS`, `_arrangement_panels`). The run's window slots (the
  session, duplicated windows) keep them, so a reopened or duplicated window draws the same
  tree. Applying a layout never changes the drawn tree.
- One render per source switch: `RaxmlTreeView.load_tree_file(path, root_mode=,
  branch_length_unit=, support_kind=, support_threshold=)` sets everything before the one
  load and render (`set_branch_length_unit` / `set_support_kind` removed).
  `set_tree_display_options` records a restored root before switching the source, so a
  restore renders once too (checked: 1 render each on the synthetic run).
- The source dropdown elides in the middle (`_MiddleElideComboBox`), so a narrow Classic
  panel keeps the end that tells trees apart ("Whole-geno...-corrected").
- Export: the default file name is `tree_<source id>.pdf`, and the figure has a caption
  above the tree with the source label and the root ("Core genome, recombination-filtered.
  Root: midpoint").
- Heatmap axis title shortened to "Samples (N not in tree)"; the tooltip keeps the detail.
- Hover-budget test: `HoverBudgetTests` (bus only) failed after `test_pan_panels` because
  that module's windows, left to `deleteLater`, were deleted inside the first timed event
  loop. Its `setUp` now flushes deferred deletes first. The real-panel hover test
  (`HeavyPanelHoverBudgetTests`) passes alone every time; it once read 31 ms in a long
  combined run (borderline machine load), and passed in the final combined run.
- Not changed: combined "SH-aLRT/UFBoot" labels. The contract gives one number per node,
  and `CoreTreeWorkflow` passes either `-B` or `--alrt`, never both, so guessing which half
  to show could mislabel support. The reference run's current (non-converged, legacy)
  Gubbins tree shows no warning chip until the rerun writes `trees.json`, as the contract
  says (the old `gubbins/summary.json` has no convergence field).

## Reference run rerun (2026-10-10)

Run `runs/sample_78eaa9126850/run_78eaa9126850`, through the pipeline's own Runner
(`<scratch>/phylo/rerun/driver.py plan|run`, `build_cfg()` of the P01 to P06 driver: GUI
Resume `DataInput`, COG dropped, 96 cores). Default `phylo` settings (no config
overrides): full Panaroo alignment, auto ordering genome, seed 1, 20 Gubbins iterations,
`weighted_robinson_foulds`, IQ-TREE `MFP` + `fconst` + UFBoot 1000 with `--bnni`, midpoint
root.

**Backup** (copies, verified by file count, byte size and `diff -rq`):
`/home/bals/Documents/dev/BACTolomew_backups/2026-10-10_phylogeny/`: `run/{gubbins,trees,fastgear}`
(15, 2 and 290 files), `run/{run_config,run_state,run_report,expected_outputs}.json`,
`run/snippy/reference.txt`, the whole GUI cache `cache/run_78eaa9126850` (5 files,
469 MB), and `manifest_before.tsv` / `manifest_after.tsv` (path, size, mtime of every file
in the run folder). `trees/core_gene_tree.nwk` was removed from the run after the copy
was verified (no longer produced).

**Plan check** (Runner's own reuse predicates): assembly_qc, bakta, mlst,
amrfinder_mutations, mobsuite, panaroo, phylogeny (Snippy), chewbbaca, diamond_markers,
diamond_reference_loci and the marker hit search reused; core_alignment ("no ordering
record"), gubbins ("no input stamp"), core_tree and fastgear_clades rerun.

**Bug found and fixed during the rerun:** the first attempt stopped in Gubbins' first
iteration with `Bus error (core dumped)` in pyjar's joint ancestral reconstruction.
Docker gives a container a private 64 MB `/dev/shm`, and Gubbins' multithreaded pyjar
puts the alignment in shared memory. `ToolEnv` now runs every docker tool command with
`--ipc host` (`DOCKER_IPC_ARGS`), so tools see the host's `/dev/shm` as they do natively and
under Singularity. The old Gubbins outputs were untouched by the failure (temp folder);
the second attempt reused the ordered alignment written by the first (same sha256).

**What ran, wall times:**

| Step | Result | Wall time |
|---|---|---|
| core_alignment | ran (attempt 1), reused in attempt 2 | 20 s |
| gubbins | ran, 64 threads | 42.4 min |
| core_tree (mask + IQ-TREE 2.4.0, 16 threads) | ran | 18.6 min (IQ-TREE 18.2 min) |
| phylogeny (Snippy) | reused; stamp written, legacy tree copied to `trees/raw/` | 2 s |
| fastgear_clades | ran | 40.7 min |
| collect_annotations | in-memory re-merge, marker hits reused | 2.1 min |
| everything else | reused | |
| run_workflow (attempt 2) | success | 104 min |
| GUI cache rebuild (`rerun/rebuild_cache.py`) | needed (`run_config.json` is in the fingerprint) | 3.9 min |

Manifest diff: only `core_alignment/`, `gubbins/`, `trees/`, `fastgear/`, the four run
JSON files and the container home `.bactolomew_home/` (Gubbins' numba cache) changed;
`snippy/`, `panaroo/`, `bakta/`, `chewbbaca/`, `annotations/`, `dmnd/` and the rest are
identical in size and mtime. `validate_completed_run`: 0 errors; run_state and run_report
status success, outputs validated. Rebuilt GUI tables are equal to the backup's
(`polars.DataFrame.equals` on all three parquet files).

**Ordering:** `auto_most_contiguous`, ordering genome GCA_003931755.1 (1 contig, N50
3,943,275; the only candidate evaluated, as it passed `ordering_min_placed`). 2610 of 2676
core genes placed, 66 appended (65 absent from the genome, 1 refound only), 1305
reverse-complemented, 0 multi-copy, 6 fragmented. Diversity: mean pairwise p-distance
0.0032, max 0.0154 (GCA_023499085.1 vs GCA_958320995.1), 67,806 variable columns; no
diversity warning; no predicted Gubbins dropouts (max missing fraction 0.12).

**Gubbins:** did NOT converge: stopped at the maximum of 20 iterations
(`gubbins_not_converged` in run_state, run_report and trees.json, so the GUI chip shows).
1246 recombination blocks (4129 before), SNPs inside 115,524 / outside 31,533, r/m (SNP
ratio) 3.66 (2.70 before). Longest-branch share of clonal-frame SNPs 0.92 on the root
edge (Node_5 + Node_151), so `gubbins_lineage_structure` warns. No samples left out.
Top per-tip recombination fractions: GCA_002017515.1 0.146, GCA_034330085.1 0.123,
GCA_014267385.1 0.118. Blocks over more than one gene: 406 (plus 6 touching appended
genes, flagged); blocks whose genes are not adjacent in the ordering genome (the P39
artefact): 680 of 708 multi-gene blocks before, 0 of 412 after.

**Core tree (`core_recomb`):** IQ-TREE 2.4.0, best-fit model K3Pu+F+I (BIC), log-likelihood
-3,764,523.29; 26,861 filtered polymorphic sites plus `-fconst` 773093,519148,524816,771182
(polymorphic + constant + uncalled = 2,627,569 columns, checked). UFBoot converged; 150
internal splits, median support 86, 61 (41%) >= 95, 51 < 70. Six sequences fail the
composition chi2 test: the six genomes of fastGEAR lineage 1 (GCA_0234990xx), expected for
the distant lineage. Two pairs of identical sequences. (split, support) set of the
`.treefile` equals that of the rooted `trees/core_recomb.nwk`.

| core_recomb vs | RF | normalized | its UFBoot >= 95 splits absent there |
|---|---|---|---|
| new Gubbins final tree | 40 | 0.133 | 0 |
| old Gubbins final tree | 128 | 0.427 | 10 |
| old `core_gene_tree.nwk` (RAxML-NG, unfiltered) | 158 | 0.527 | 16 |
| Snippy tree | 168 | 0.560 | 15 |

Old vs new Gubbins final tree: RF 116 (0.387).

**fastGEAR** (on the ordered alignment, `fastgear/input.json` stamp written): lineages and
clusters unchanged (2 lineages of 6 and 147, clusters of 4, 6, 60 and 83; adjusted Rand
index 1.000 for both). Recent recombination rows 5670 (6336 before), ancestral 7 (3
before): the segments follow the new column order, as expected. Under the midpoint root
of the new default tree both lineages and clusters 1, 2 and 4 are monophyletic; cluster 3
is not (it is not in the old trees either).

**trees/trees.json:** default `core_recomb`, rooting midpoint, ordering genome
GCA_003931755.1, Gubbins `converged: false` (20 of 20), warnings `gubbins_not_converged`
and `gubbins_lineage_structure`. Trees: `core_recomb` ("Core genome, recombination-filtered",
UFBoot 1000, threshold 95), `gubbins_final` ("Gubbins final tree, recombination-filtered
(no support)", SNP branch lengths), `snippy_core` ("Whole-genome SNPs (Snippy), not
recombination-corrected", `Reference` tip renamed to GCA_028622775.1). All 153 samples in
each tree. The GUI's `tree_sources()` reads the index with `core_recomb` as default.

**Open points:** Gubbins non-convergence at 20 iterations on this mixed cohort (the root
edge holds 92% of the clonal-frame SNPs); a per-lineage Gubbins run (lineage 2, 147
genomes) or more iterations are the options, not done. `trees/core_recomb_iqtree/` was
created with mode 0700 (temporary folder renamed in place); harmless for the owner.

## Before and after (rechecked)

Adversarial recheck of the impact report (2026-10-10). Every number was recomputed from the
run files and the pre-rerun backup (`BACTolomew_backups/2026-10-10_phylogeny/run`) with
independent scripts (Biopython trees, own GFF and per-branch parsing; scratch
`phylo/impact_check`). Numbers that held are listed briefly; corrections are marked
**Corrected**.

**Ordering genome and gene placement.** Confirmed: 2,676 core genes, 2,610 placed (97.5%),
1,305 reversed, 66 appended (65 absent, 1 refind only), 6 fragmented, length 2,627,569 bp,
67,806 variable columns before and after. Every gene block of the ordered alignment equals
its Panaroo block (reverse complemented where flagged; 2,676 of 2,676 checked), placed genes
follow the chromosome with no position inversion, appended genes are all at the end.
**Corrected:** GCA_003931755.1 is not the only single-contig assembly: three cohort genomes
are complete one-contig chromosomes (GCA_003931755.1 3,943,275 bp, GCA_003931775.1
3,943,105 bp, GCA_014267385.1 3,909,403 bp); GCA_003931755.1 won the N50 tie-break.
The 66 appended genes are 2.5% of genes but 78,469 bp, 3.0% of the alignment; 71 of the
1,246 Gubbins blocks touch that appended tail, where neighbors are not physical.

**Gubbins.** Confirmed: before 4,129 blocks (median 616 bp, max 15,863), union 1,382,872 bp
(52.6%), SNPs in/out 88,918 / 32,894, r/m 2.70, 208 of 305 branches with blocks, per-branch
r/m median 5.0 (IQR 2.0 to 13.4), rho/theta median 0.50, Node_151 holds 80.6% of clonal
SNPs, 28,672 filtered sites, random seed (`-p 1459`), stopped at 5 iterations. After 1,246
blocks (median 229 bp, max 67,061), union 1,197,444 bp (45.6%), 115,524 / 31,533, r/m 3.66,
198 branches with blocks, r/m median 5.7 (IQR 1.8 to 30.9), rho/theta median 0.37, root edge
28,935 clonal SNPs (91.8%), 26,861 filtered sites, seed 1, stopped at 20 iterations
(not converged). Multi-gene blocks 708 to 412, of which joining non-adjacent placed genes
680 to 0.
**Corrected:**
- "Fewer, longer blocks" holds for the mean (1,187 to 5,363 bp) and the maximum, but the
  median block got shorter (616 to 229 bp). The new block set is a few very long tracts plus
  many short ones.
- The 0 non-adjacent multi-gene blocks is guaranteed by construction (neighbors in the ordered
  alignment are neighbors on the one chromosome); it shows the artefact cannot occur, not
  that the calls are biologically better. `gubbins/summary.json` says 406 multi-gene blocks
  while `recombination_blocks_genes.tsv` has 412 (different counting rule; minor).
- Recombination is not "mostly on the edge between the two lineages": root-edge blocks cover
  only 61,536 bp (2.3% of the alignment; 215,408 bp before). Almost all block coverage is
  within lineage 2 (tip blocks alone cover 824,296 bp).
- The cohort-wide r/m of 3.7 mixes two regimes. The root edge carries 28,935 clonal and
  8,512 recombinant SNPs (r/m 0.29). Without it, r/m is 41 (107,012 / 2,598) after versus
  12.6 (80,199 / 6,365) before: within lineage 2 the clonal frame lost 59% of its SNPs.
- Per-tip recombination: the median cumulative fraction is 0.21 after versus 0.22 before
  (root to tip, own genome length). The top tips are not lineage 1 genomes. By own terminal
  branch the top is GCA_002017515.1 (0.146 of the alignment, 0.15 of its genome), then
  GCA_034330085.1, GCA_014267385.1, GCA_030550955.1, GCA_026558295.1 (before:
  GCA_034330085.1 at 0.090).

**Core tree (IQ-TREE 2.4.0).** Confirmed: K3Pu+F+I by BIC (w-BIC 0.999), A-G = C-T 5.61,
A-T = C-G 0.79, p-inv 0.880, logL -3,764,523.29, 1,000 UFBoot converged, seed 1, `--bnni`,
fconst 773,093 / 519,148 / 524,816 / 771,182 (2,615,100 sites with the 26,861 variable
ones). Six lineage 1 genomes fail the composition test (p about 1%); 2 identical pairs.
**Corrected or qualified:** by AIC the best model is TVM+F+I (w-AIC 0.68) and GTR+F+I+R5 gets
0.19, so "no rate heterogeneity beyond +I" is a BIC result. The fconst counts include the
roughly 41,000 variable columns Gubbins removed as recombinant (counted as constant). The
composition failure is expected when the test runs on SNP columns where these six carry the
minor allele at most sites; it says nothing about recombination history. 38 of 151 internal
branches are shorter than 0.5 expected substitutions over the whole alignment
(effectively polytomies; 4 of them still have UFBoot 95 or more).

**Tree comparison.** Confirmed exactly: new core vs new Gubbins nRF 0.133 (0 strong
conflicts), vs old Gubbins 0.427 (10), old core 0.527 (16), Snippy 0.560 (15); ST492
0.122 / 0.344 / 0.478 / 0.478; ST2 0.140 / 0.561 / 0.596 / 0.702; old core vs Snippy 0.460;
old vs new Gubbins 0.387. Strong conflicts within ST492: 0 / 4 / 8 / 5; within ST2:
0 / 7 / 8 / 11.
**Qualified:** the old trees already disagreed with each other by nRF 0.45 to 0.47 (old core vs
old Gubbins 0.453, old Gubbins vs Snippy 0.473). The new tree is about as far from the old
trees as the old trees were from each other; the fine topology was never stable. What
matters are the 10 to 16 strongly supported (UFBoot 95 or more) new splits the old trees
contradict.

**Patristic correlation.** Confirmed: new core vs new Gubbins 1.000 / 0.999, old Gubbins
1.000 / 0.846, old core 0.965 / 0.818, Snippy 0.909 / 0.792; old trees among themselves
Spearman 0.958 to 0.973; within ST492 0.82 / 0.74 / 0.72; within ST2 0.63 / 0.55 / 0.51.
**Corrected:** "the relative depth of lineage 1 vs lineage 2 is unchanged" is false. The mean
lineage 1 to lineage 2 distance divided by the mean within lineage 2 distance is 370 in the
new core tree, 283 in the old Gubbins tree, but 12.9 in the old core tree and 6.8 in the
Snippy tree: removing recombination collapses lineage 2 relative to the lineage split.
Within lineage 2 alone, Pearson against the old Gubbins, old core and Snippy trees is only
0.75, 0.75 and 0.70.

**Topology and rooting.** Confirmed: midpoint root on the lineage 1 edge in all five trees
(6 versus 147); cluster 4 (83) a clade with UFBoot 100; enclosing clades 84 (+GCA_008987025.1,
UFBoot 98), 90 (+3 more ST2 and the 3 cluster 3 ST492, 70), 94 (8 ST2 and those 3, 97),
140 (55); cluster 1 branches first in lineage 2 (cluster 1 UFBoot 100, the remaining 143
UFBoot 91); cluster 3 is not monophyletic in any tree.
**Corrected:** GCA_008987025.1 was the sister of cluster 4 in the old Gubbins and Snippy trees,
but not in the old unfiltered core tree, where cluster 4's sister was a group of 26 ST2.
The move of GCA_026558295.1, GCA_026563615.1 and GCA_026566215.1 next to cluster 4 is new
(before they joined cluster 4 only at the 134 to 142 genome level). It rests on a UFBoot 70
node inside a UFBoot 97 node, and GCA_026558295.1 is among the five tips with the most own
branch recombination, so the placement depends on the non-converged recombination calls.
"Intermediate forms" is an overstatement: they are ST492 genomes outside cluster 4, in the
same well-supported 94 genome clade.

**Bootstrap support.** Confirmed: 150 splits, median 86, 61 (41%) at 95 or more, 99 (66%) at 70
or more, 51 below 70; clades of only ST492 n=88 (median 81.5, 38% / 61%), only ST2 n=50
(89.5, 46% / 68%); clusters 1, 2, 4 and the lineage split at 100; the rooted file keeps the
same splits and supports as IQ-TREE's unrooted treefile.

**fastGEAR.** Confirmed: lineages and clusters identical (`sample_clades.tsv` byte-identical);
recent rows 6,335 to 5,669, median 638 to 789 bp, mean 824 to 1,832, max 10,993 to 31,277,
total 5.22 to 10.38 Mb, median per strain 6,795 to 13,744 bp (147 lineage 2 strains; lineage 1
receives none).
**Corrected:** ancestral events are 6 after and 2 before (not 7 and 3; the earlier count, also
in the rerun section above, included a header line). New logBF 4.2 to 272.0; old 2.7 and
89.6. The recent rows are per-strain records: they collapse to 1,690 distinct segments after
and 2,119 before (fastGEAR's own header says 339 and 812). The doubled total is not explained
by joining tracts that Panaroo's order split (joining keeps the total); fastGEAR calls more
sequence recombinant when the signal is contiguous.

**GUI.** Default source `core_recomb` and the "Gubbins did not converge" chip confirmed in the
code and in `phylo/final/02_phylogeny_panel_default_tree.png`. **Corrected:** the root shows
as the dropdown value "As in file (midpoint)", and at default panel width the source label is
elided ("Core genome...on-filtered"). Screenshots exist in `phylo/final/`. The old default
(Snippy tree first) is confirmed from the 2026-10-09 code backup. Row order figures (Spearman
0.926, cluster 3 runs 52 + 8 to 60) were not independently recomputed.

**Reproducibility.** The seed is fixed (Gubbins `-p 1`, IQ-TREE `--seed 1`), but no byte-identical
rerun of the 153-genome Gubbins step was done.

## Rework 2 (2026-10-12): Gubbins per lineage

User decisions of 2026-10-12, planned here; nothing below is implemented yet. Gubbins runs
once per fastGEAR lineage, each lineage's recombination calls mask only that lineage's
genomes in the full ordered core alignment, and one IQ-TREE tree is built for the whole
cohort on that masked alignment. The same automatic rule applies to every dataset.

### R2.1 Step order and dependencies (Runner)

New order: panaroo -> core_alignment -> **fastgear_clades** -> gubbins -> core_tree ->
phylogeny (Snippy only) -> chewbbaca -> ... (unchanged after that).

| Step | Reads | Writes |
|---|---|---|
| core_alignment | Panaroo alignment, header, presence table, GFF3 | `core_alignment/core_gene_alignment.ordered.aln`, `gene_map.tsv`, `ordering.json` (unchanged) |
| fastgear_clades | the ordered alignment | `fastgear/` (unchanged, now before Gubbins) |
| gubbins | ordered alignment, `fastgear/sample_clades.tsv` (`FastGearLineage`), `ordering.json` | `gubbins/` (R2.6), including the masked full alignment |
| core_tree | `gubbins/core_gene_alignment.masked.aln` | `trees/core_recomb_iqtree/`, `trees/core_recomb.*`, `trees.json` |
| phylogeny (Snippy) | Snippy inputs | unchanged; still independent of Gubbins |

The fastGEAR block moves up in `Runner._run_workflow_impl`, right after core_alignment.
fastGEAR, Gubbins and IQ-TREE ran serially before too, so the order change adds no wall
time by itself.

Inside the gubbins step:
1. Resolve the mode and the groups (R2.2). A group is a lineage, or `all` in single mode.
2. Decide per group: run, reuse (R2.3) or skip (too small, too few variable sites).
3. Stage each group's input `gubbins/lineage_<id>/input.aln`: the group's rows of the
   ordered alignment, all columns, in ordered-alignment row order. Columns are not removed,
   so Gubbins coordinates equal the ordered alignment's coordinates for every group.
4. Run the groups concurrently (`ThreadPoolExecutor`, one worker per group) within the
   `gubbins` thread budget (`thread_plan["gubbins"]`, 64 now). Split: threads_i =
   max(1, floor(T * v_i / sum v)), where v_i is the group's variable column count; the
   largest group gets the remainder. Reference run: lineage 1 gets 1 thread, lineage 2
   gets 63. Every call writes into a temp folder inside `gubbins/` that replaces
   `lineage_<id>/` only on success.
5. Parse each group (convergence, blocks, r/m, per-branch statistics, left-out taxa), as
   `GubbinsWorkflow` does now, but per group.
6. Mask each group with Gubbins' own `mask_gubbins_aln.py`, then merge (R2.5).
7. Write the combined GFF and gene table, `gubbins/summary.json` (schema 3) and the
   run warnings.
Any group failure (non-zero exit, missing outputs) fails the step and the run, as now:
there is no silent fallback to unmasked genomes for a lineage that should have run.

### R2.2 Mode and group rule (dataset agnostic)

- `phylo.gubbins_mode`: `auto` (default) | `single` | `per_lineage`.
- `auto`: per lineage when `fastgear/sample_clades.tsv` names 2 or more lineages, else one
  run on all genomes (yesterday's behavior). `per_lineage` with only one lineage behaves
  like `single` (recorded as requested `per_lineage`, resolved `single`, note).
- `single`: one group `all` with every genome, outputs at the top of `gubbins/` exactly as
  now (`core_gene.*`), so old runs, validators and the manuscript scripts that read
  `gubbins/core_gene.final_tree.tre` keep working. `phylo.gubbins_min_genomes` (4) still
  skips Gubbins entirely.
- A lineage is not run when it has fewer than `phylo.gubbins_min_lineage_size` genomes
  (default 5), or fewer than 3 variable columns among its called bases (Gubbins'
  `--min-snps` default is 3: no block can be called). Its genomes stay unmasked; the group
  entry records `ran: false` and `skip_reason` (`below_min_lineage_size` |
  `too_few_variable_sites`); warning `gubbins_lineage_skipped`.
- Why 5: Gubbins 3.4.3 exits below 3 sequences (`common.py` line 193, "Three or more
  sequences are required for a meaningful phylogenetic analysis."). With 3 or 4 taxa the
  unrooted tree has 0 or 1 internal branches, so the reconstruction carries almost no
  information. Configured values below 3 are refused at run start.
- Genomes without a lineage in `sample_clades.tsv` (should not happen; the fastGEAR
  validator checks the sample set): unmasked, reason `no_lineage`, same warning.
- fastGEAR unavailable or failed: in `auto` the step falls back to `single` with warning
  `gubbins_lineages_unavailable` (message includes the fastGEAR error). fastGEAR stays a
  required step: the Runner catches its exception, runs gubbins and core_tree in single
  mode, runs Snippy, then ends the run as failed with the fastGEAR error (before chewBBACA,
  so nothing else runs on a failed run). A resume reruns fastGEAR; if it then reports 2 or
  more lineages, Gubbins and the core tree rerun per lineage (R2.3). In `per_lineage`
  (explicit) a missing fastGEAR result fails the gubbins step.
- Per-lineage runs use no `--outgroup`: a cohort outgroup usually sits in another lineage,
  and Gubbins' root only affects its own reconstruction (midpoint by default). The cohort
  tree is still rooted by `phylo.outgroup` in the rooting step (3.5). Single mode passes
  `--outgroup` as now. Consequence: in per-lineage mode a change of `phylo.outgroup`
  reruns only the rooting step.

### R2.3 Reuse and freshness

Per-group stamp `gubbins/lineage_<id>/run.json` (single mode: inside `summary.json`):
sorted member names, `input_sha256` (sha256 of the staged input, computed while it is
streamed; deterministic because the row order and line format are fixed), params (seed,
iterations, converge method, filter percentage, tree builder, outgroup in single mode
only), Gubbins version, `algorithm_version` (GubbinsWorkflow `ALGORITHM_VERSION` 3),
threads (recorded, not part of the key, as decided on 2026-10-10).

| Change | What reruns |
|---|---|
| ordered alignment sha256 changes | every group (input sha differs), mask, core tree |
| fastGEAR reruns, same lineage memberships | nothing in Gubbins: groups are matched by member set, not by lineage id; a renumbered lineage's folder is renamed to the new id |
| one lineage gains or loses a genome | that lineage only; the others are reused; then mask and core tree |
| a lineage crosses `gubbins_min_lineage_size` | that lineage runs or is dropped; mask and core tree |
| resolved mode changes (single <-> per lineage) | every group; the old layout is removed (top-level `core_gene.*` or the `lineage_*/` folders) |
| Gubbins params, version or `algorithm_version` | every group |
| `phylo.outgroup` (per lineage) | rooting only (3.5) |

- The mask and merge (R2.5) is redone whenever any group ran or the group set changed
  (about 2 to 4 min on the reference run). The masked alignment is reused when
  `summary.json` `masked_alignment` records the same ordered alignment sha, the same per-
  group GFF sha256 list, the same unmasked list and an output sha256 equal to the file.
- core_tree: its stamp's `input_sha256` is the masked alignment's sha256; the recombination
  GFF sha and the fconst counts leave the stamp. `CoreTreeWorkflow.ALGORITHM_VERSION` 2.
- `trees.json`: `core_recomb` is fresh when its stamp matches the masked alignment sha in
  `gubbins/summary.json`.
- The 2026-10-10 Gubbins output (schema 2, no `mode`) is treated as stale: Gubbins and the
  core tree rerun on resume; fastGEAR and Snippy are reused.

### R2.4 Gubbins 3.4.3 options (checked in `bactolomew-tools:0.1.0`, env `gubbins`)

- `run_gubbins.py --version`: 3.4.3.
- `--converge-method {weighted_robinson_foulds,robinson_foulds,recombination}`, default
  `weighted_robinson_foulds`. Code (`gubbins/common.py` 574-590, 1373-1412):
  `weighted_robinson_foulds` needs a weighted RF of exactly 0 to any earlier tree, so the
  branch lengths must repeat too (very strict; the reason the cohort run never converged);
  `robinson_foulds` needs an identical topology; `recombination` stops when the current
  iteration's per-taxon block coordinates (`.tab` file: start..end per taxon) equal those
  of any earlier iteration, checked from iteration 2. Stdout lines, already parsed:
  "Convergence after N iterations: Recombinations observed before." /
  "...: Tree observed before." / "Maximum number of iterations (N) reached.".
  New default `phylo.gubbins_converge_method` = `recombination` (already an allowed value
  in `_PHYLO_CHOICES`). Convergence is recorded per group.
- Masking: `mask_gubbins_aln.py --aln <fasta> --gff <recombination_predictions.gff> --out
  <fasta> [--out-fmt fasta] [--missing-char C]` (default `-`). It reads the whole alignment
  with Biopython, and for each GFF line masks the listed taxa over start..end (1-based,
  inclusive) only when every listed taxon is in the alignment (`taxon_set.issubset`).
  Already mapped to the `gubbins` env in `_tool_env_key`. We pass `--missing-char N`, so
  masked bases stay distinguishable from Panaroo's gaps (`-`, absent gene or indel);
  IQ-TREE treats both as missing.
- Minimum taxa: 3 after Gubbins' own `--filter-percentage` filtering (exit otherwise).
- Outgroup: one name or a comma list forming a clade; Gubbins exits if the outgroup was
  filtered out. Not used per lineage (R2.2).
- Seeds: `--seed` as now (same seed for every group).

### R2.5 Masking merge algorithm

1. For each group that ran: `mask_gubbins_aln.py --aln lineage_<id>/input.aln --gff
   lineage_<id>/core_gene.recombination_predictions.gff --out <tmp>/lineage_<id>.masked.aln
   --missing-char N` (gubbins env, run concurrently, one per group). Every block's taxa
   are members of that group by construction, so no block is skipped by the subset check;
   any skipped GFF line (taxa not in the input) is counted and fails the step.
2. Stream the ordered alignment once, in its row order. For each genome: if its group ran
   and Gubbins kept it, write its row from that group's masked file; otherwise write the
   original row (unmasked) and record why (`lineage_skipped`, `dropped_by_gubbins`,
   `no_lineage`).
3. Checks during the merge (cheap, numpy per row): the row length equals the alignment
   length; a masked row differs from the original only where it now holds `N`; the number
   of newly masked positions equals the union of that genome's own blocks in its group's
   GFF minus positions already `N`. A mismatch fails the step. This is the property "each
   genome masked only by its own lineage's calls".
4. Write `gubbins/core_gene_alignment.masked.aln` through a temp file, record its sha256,
   per-genome masked bases and fraction (`gubbins/masked_bases.tsv`: sample, group,
   masked_bases (A/C/G/T), masked_positions, masked_fraction, unmasked_reason), delete the staged inputs and masked
   per-group files (sha256 of the staged inputs is in the stamps; they are rebuilt from the
   ordered alignment when needed).
- Memory: `mask_gubbins_aln.py` holds one group in memory (lineage 2: about 390 MB); the
   merge streams one row at a time.
- Location: under `gubbins/`, not `core_alignment/`: the gubbins step writes it, and its
  freshness follows the Gubbins calls; `core_alignment/` stays owned by the core_alignment
  step and its reuse check.
- Single mode uses the same code (one group `all`), so core_tree always reads the masked
  alignment. The fconst and asc routes of 3.3 are deleted with
  `phylo.iqtree_invariant_correction` (an old config value is ignored, as any unknown key).
- **Limitation (recorded in `summary.json` notes, `trees.json` notes of `core_recomb` and
  here):** recombination between lineages is not masked. Imports from one lineage into the
  other that happened on the long inter-lineage branch, or before the lineages split, are
  invisible to per-lineage runs; they inflate the inter-lineage branch, not the
  within-lineage branching. Recombination within a skipped lineage is not masked either.

### R2.6 Output layout

Per-lineage mode:
- `gubbins/lineage_<id>/`: the standard Gubbins outputs with prefix `core_gene`
  (`.final_tree.tre`, `.node_labelled.final_tree.tre`, `.recombination_predictions.gff`
  and `.embl`, `.per_branch_statistics.csv`, `.filtered_polymorphic_sites.fasta` and
  `.phylip`, `.branch_base_reconstruction.embl`, `.summary_of_snp_distribution.vcf`,
  `.log`), `command.stdout.log`, `command.stderr.log`, `run.json` (stamp and parsed
  values). Skipped groups get only `run.json` with `ran: false`.
- `gubbins/recombination_predictions.gff`: all groups' blocks, column 9 extended with
  `lineage="<id>"`.
- `gubbins/recombination_blocks_genes.tsv`: the current table plus a `lineage` column;
  genes from `gene_map.tsv` (cluster name, ordering-genome locus tag and coordinates).
- `gubbins/core_gene_alignment.masked.aln`, `gubbins/masked_bases.tsv`.
- `gubbins/summary.json`, schema_version 3: `mode` (requested), `resolved_mode`,
  `mode_reason` (`two_or_more_lineages` | `one_lineage` | `fastgear_unavailable` |
  `configured`), tool version, params, `min_lineage_size`, `algorithm_version`, `groups`:
  one entry per group with `id`, `members` (names), `genomes`, `ran`, `skip_reason`,
  `variable_sites`, `diversity` (mean, max), `threads`, `command`, `input_sha256`,
  `converged`, `iterations_run`, `max_iterations`, `convergence_reason`,
  `recombination_blocks`, `snps_inside_recombination`, `snps_outside_recombination`,
  `r_m_snp_ratio`, `longest_branch_share`, `tip_recombination_fraction_top5`,
  `samples_left_out`, `wall_time_s`, `warnings`; plus `masked_alignment` (path, sha256,
  inputs, `samples_unmasked` with reasons), aggregated `converged` (true only when every
  group that ran converged), `not_converged` (group ids), `warnings`, `notes`.
- Single mode: top-level `core_gene.*` as now plus the masked alignment and
  `masked_bases.tsv`; `summary.json` schema 3 with one group `all`.
- Diversity and lineage-structure warnings are computed per group (the per-group variable
  columns come from the staging pass). `gubbins_lineage_structure` now flags a dominant
  branch inside a group.
- Validators: `validate_gubbins_outputs` accepts schema 2 (old runs, single layout) and
  schema 3 (either layout: the files of every group with `ran: true`, the masked
  alignment, the combined GFF); a skipped run keeps working as now.
  `expected_outputs.json` and `run_expectations.py` list the per-lineage files from the
  resolved groups (only known after fastGEAR: the expected list is rewritten when the
  gubbins step starts).

### R2.7 Core tree and trees.json

- IQ-TREE on `gubbins/core_gene_alignment.masked.aln`: `iqtree2 -s <masked> -m MFP -B 1000
  --bnni --seed <seed> -T <threads> --prefix trees/core_recomb_iqtree/core_recomb`, no
  `-fconst`, no `+ASC`; fast mode and few-genome rules of 3.3 unchanged. Branch lengths in
  substitutions per site of the full core alignment.
- `core_recomb` entry: `alignment` `gubbins/core_gene_alignment.masked.aln`,
  `alignment_kind` `core_alignment_masked`, `sites` = alignment length, `invariant_correction`
  `null`, `recombination_filtered: true` when at least one group ran, label unchanged
  ("Core genome, recombination-filtered"), new key `samples_unmasked` (list of names, `[]`
  when none), notes: the inter-lineage limitation (per-lineage mode) and the unmasked
  genomes with reasons.
- When no group ran (all lineages skipped): the tree is `core` ("Core genome (Gubbins
  skipped)") on the unmasked ordered alignment, `recombination_filtered: false`.
- `gubbins_final`: listed only in single mode. Per-lineage Gubbins trees stay in
  `gubbins/lineage_<id>/` and are never listed as cohort trees.
- `trees.json` `schema_version` 2. Top-level `gubbins` keeps the schema 1 keys with
  aggregated values (`ran`: any group ran; `converged`: all ran groups converged;
  `iterations_run`: the largest; `max_iterations`; `converge_method`; `samples_left_out`:
  now `[]` in per-lineage mode, since dropped genomes stay in the tree unmasked) and adds
  `mode`, `resolved_mode` and `lineages`: per group `id`, `genomes`, `ran`, `skip_reason`,
  `converged`, `iterations_run`, `max_iterations`, `convergence_reason`,
  `recombination_blocks`, `r_m_snp_ratio`. Every key present, `null` when unknown.
- GUI (`model/pan_data.py`, `panaroo_core_accessory.py`): reads schema 1 and 2. The chip
  shows when the shown tree is `recombination_filtered` and a lineage that ran did not
  converge: "Gubbins did not converge (lineage 2)" (several: "lineages 1, 2"); tooltip
  with each lineage's iterations ("Lineage 2 stopped after 20 of 20 iterations."). Schema
  1 files and single mode: the current chip text and rule. `samples_unmasked` non-empty:
  "N genomes not recombination-masked" in the header notes, names in the tooltip. Old
  runs without `trees.json` keep the 2026-10-10 fallback.

### R2.8 Docker shared memory

- Replace `DOCKER_IPC_ARGS = ("--ipc", "host")` with `--shm-size <size>` in both docker
  command builders of `ToolEnv` (`_container_prefix` and the second builder near line
  380). Config `container.shm_size` (string, docker syntax), default `"8g"`.
- Sizing: pyjar puts the SNP reconstruction array (polymorphic sites x internal nodes,
  int8) and the base patterns in `SharedMemoryManager` blocks; on the cohort that is about
  10 to 30 MB per block, yet the default 64 MB failed with a bus error, so the real peak
  (several blocks, the manager's copies) is measured, not derived. A tmpfs `--shm-size` is
  a ceiling, not an allocation: memory is used only when written, so a generous default
  costs nothing. When `resources.max_ram_gb` is set, shm pages count against the
  container's `--memory` cgroup, so `shm_size` is capped at half of it.
- Check during the rerun: sample the container's `/dev/shm` use every 10 s
  (`docker exec <id> df -B1 /dev/shm`) for the lineage 2 run; record the peak in the
  progress notes; the default must be at least 4 times that peak. Success criterion: no
  "Bus error" in either lineage's logs. Singularity is unchanged (host `/dev/shm`).

### R2.9 GUI: long-branch break on the Distance scale

- Rule (dataset agnostic, on the rendered rooted tree, Distance layout only): sort edges
  by length L1 >= L2 >= ...; for m = 1..min(4, edge count) let B_m be the m longest edges
  and H_m the 95th percentile root-to-tip depth with the edges of B_m counted as 0. Take
  the smallest m with L_m > k * H_m and L_(m+1) <= k * H_m; break the edges of B_m. No such
  m: no break. k = 10 (constant `LONG_BRANCH_BREAK_FACTOR` in `tree_view.py`). A midpoint
  root on a long edge splits it into two long root edges, which this rule finds as m = 2.
  Reference tree: lineage 1 to lineage 2 mean distance about 370 times the within-lineage
  2 distance, so m = 2 breaks the two root edges and nothing else.
- Drawing: a broken edge is drawn with length 0.5 * H_m, with a break mark (two short
  parallel slanted strokes across the edge at its middle, `COLORS` text color) and a small
  label above it with the true length and unit ("0.0123 subst/site"). Hover on the edge
  says "Branch shortened for display; true length 0.0123 substitutions/site."
- The scale bar is drawn from the unbroken scale (edges outside B_m keep their scale), so
  it stays correct for the unbroken part; its caption gains "broken branches not to
  scale" while a break is drawn.
- Toggle: a "True branch lengths" checkbox beside the Distance/Cladogram choice, default
  off, disabled with a tooltip when the rule breaks nothing; saved in the panel state
  (`save_state` / `restore_state`) like the root choice.
- Leaf order: only x positions change; y (leaf order) and therefore the linked Gene
  presence/absence row order (`TreeRowLink`) and `_reference_leaf_order` are unaffected.
- Exports: `export_figure` renders the scene, so exports show exactly what is drawn
  (breaks and labels included, or true lengths when toggled); the caption adds "Long
  branches shortened (marked)." when breaks are drawn.
- Fit to window (commits c4b33c4..d11b88d) uses the shortened width.

### R2.10 Config keys (added or changed)

| Key | Default | Meaning |
|---|---|---|
| `phylo.gubbins_mode` | `"auto"` | `auto` (per lineage when fastGEAR has 2 or more lineages) \| `single` \| `per_lineage` |
| `phylo.gubbins_min_lineage_size` | `5` | lineages smaller than this are not run (values below 3 refused) |
| `phylo.gubbins_converge_method` | `"recombination"` (was `weighted_robinson_foulds`) | Gubbins `--converge-method` |
| `phylo.iqtree_invariant_correction` | removed | the core tree reads the masked full alignment |
| `container.shm_size` | `"8g"` | docker `--shm-size` for every tool container |

`run_config.json` `inputs.phylo` snapshots the new keys; step info of gubbins gains
`mode`, `resolved_mode`, `groups` (id, genomes, ran, converged, iterations_run).

### R2.11 Edge cases

- One lineage (or fastGEAR finds one): `auto` resolves to `single`, yesterday's behavior
  apart from the converge method and the masked-alignment core tree.
- A lineage below the minimum size, or with fewer than 3 variable columns: not run,
  genomes unmasked, `gubbins_lineage_skipped`, listed in `samples_unmasked`.
- All lineages skipped: no Gubbins call, `core` tree on the unmasked ordered alignment.
- fastGEAR failed or unavailable: `single` with `gubbins_lineages_unavailable` (R2.2).
- A genome Gubbins drops (above `--filter-percentage` within its lineage): kept in the
  cohort tree unmasked, reason `dropped_by_gubbins`, warning `gubbins_samples_unmasked`;
  predicted from `ordering.json` `missing_fraction` as now. It is no longer missing from
  the default tree (previously `samples_missing`).
- A lineage whose staged alignment is identical for all members (no variable columns):
  skipped as `too_few_variable_sites`.
- Reference run's lineage 1: 6 genomes with only 13 variable columns (checked on the
  ordered alignment); it runs (6 >= 5, 13 >= 3) and is expected to call few or no blocks.
  Lineage 2: 147 genomes, 41,227 variable columns (67,806 cohort).

### R2.12 Tests (very few)

1. Masking merge: synthetic ordered alignment (two lineages and one skipped singleton),
   two per-lineage GFFs including a block at the same columns in the other lineage; each
   genome is masked only by its own lineage's blocks, skipped and dropped genomes keep
   their original row, the masked-count check fails on a corrupted row.
2. Mode resolution: lineage tables with 1 lineage (single), 2 lineages (per lineage), a
   lineage below `gubbins_min_lineage_size` (skipped, reason), missing fastGEAR (single
   plus `gubbins_lineages_unavailable`).
3. Long-branch break: a tree with one edge 300 times the rest: the leaf order equals the
   unbroken layout's, the edge is drawn shorter with its true-length label, the toggle
   restores true lengths; a balanced tree breaks nothing.
Existing tests that hard-code the fconst route or `gubbins_final` in `trees.json` are
adjusted, not extended. Everything else is checked on the reference run.

### R2.13 Reference run rerun plan

1. Guards: no other driver on the run; the cgMLST schema build stays untouched; chewBBACA
   must be "reuse" in the plan.
2. Backup to `/home/bals/Documents/dev/BACTolomew_backups/2026-10-12_gubbins_per_lineage/`
   (verified by file count, size and `diff -rq`): `run/gubbins/` (yesterday's cohort run),
   `run/trees/`, `run/{run_config,run_state,run_report,expected_outputs}.json`, the GUI
   cache `cache/run_78eaa9126850/`, and `manifest_before.tsv` / `manifest_after.tsv`. The
   2026-10-10 backup keeps the original state.
3. Driver `<scratch>/phylo2/rerun/driver.py plan|run` reusing yesterday's `build_cfg()`.
   `plan` refuses unless exactly: core_alignment reuse (ordering unchanged, same sha256),
   fastgear_clades reuse (`input.json` stamp matches), gubbins RERUN (schema 2 summary,
   algorithm version), core_tree RERUN, phylogeny (Snippy) reuse, everything else reuse,
   collect_annotations in-memory re-merge.
4. Expected: resolved mode per lineage, lineage 1 (6) and lineage 2 (147) both run,
   concurrently (1 and 63 threads). Wall time: lineage 1 seconds; lineage 2 about 15 to
   45 min (41,227 variable columns versus 67,806; no root edge; convergence on the calls
   may stop well before 20 iterations); masking and merge about 3 to 5 min; IQ-TREE about
   15 to 30 min: the masked full alignment compresses to about 11,500 distinct column
   patterns (measured with yesterday's cohort calls as a proxy; 25,218 unmasked), against
   yesterday's 26,861 sites plus fconst; GUI cache rebuild about 4 min. Total about 40 to
   85 min, fastGEAR reused. A fresh run adds fastGEAR (about 41 min, as yesterday), so a
   full run is about as long as yesterday's 104 min, or shorter if lineage 2 converges
   early. Disk: the masked alignment adds about 400 MB.
5. Docker shm check (R2.8) during lineage 2.
6. After: manifest diff (only `gubbins/`, `trees/`, the four run JSONs and the container
   home may change; `core_alignment/`, `fastgear/`, `snippy/`, `panaroo/`, `bakta/`,
   `chewbbaca/`, `annotations/`, `dmnd/` identical); `validate_completed_run` 0 errors;
   warnings reviewed; GUI cache rebuild; app check: Phylogeny with the long-branch break,
   the toggle, the chip (only if a lineage did not converge), an export, the linked Gene
   presence/absence row order.

### R2.14 Metrics against yesterday

Baselines: the current outputs (copied to the 2026-10-12 backup: cohort Gubbins with
`weighted_robinson_foulds`, 20 iterations, not converged; `core_recomb` with fconst) and
the original state (`BACTolomew_backups/2026-10-10_phylogeny`). Scripts in
`<scratch>/phylo2/metrics/`.

- Gubbins per lineage: converged, iterations, convergence reason, wall time, threads.
  Blocks, SNPs inside/outside, r/m. Lineage 2 against yesterday's cohort run restricted to
  lineage 2 (r/m 41: 107,012 / 2,598 clonal SNPs off the root edge; 1,246 cohort blocks,
  root edge 61,536 bp). Lineage 1 blocks (expected about 0).
- Masking: per-genome masked bases and fraction against yesterday's per-genome masked
  union (cohort GFF); Jaccard of masked positions per genome; genomes whose masking
  changed most; columns masked in all genomes.
- Core tree: model and log-likelihood, distinct patterns, UFBoot convergence and
  distribution (median, share >= 95, < 70), composition-test failures (yesterday: the 6
  lineage 1 genomes), identical sequences.
- Topology: RF and normalized RF of the new `core_recomb` against yesterday's
  `core_recomb` (current), yesterday's Gubbins final tree, the old Gubbins, old core and
  Snippy trees; whole tree, within ST492, within ST2; strongly supported (UFBoot >= 95)
  conflicts each way.
- Distances: between-lineage over within-lineage-2 mean patristic ratio (370 yesterday);
  Pearson and Spearman of within-lineage-2 patristic distances against yesterday's
  `core_recomb`; pairwise SNP distances within ST492 from the masked alignment against the
  unmasked ordered alignment (share of pairs crossing the usual cluster thresholds).
- Clades: lineages and clusters 1, 2, 4 monophyletic under midpoint rooting with support;
  cluster 4 inside ST2; placement of GCA_008987025.1 and of GCA_026558295.1,
  GCA_026563615.1, GCA_026566215.1 (the UFBoot 70 node of yesterday).
- GUI: Gene presence/absence row order change (Spearman against yesterday), long-branch
  break screenshots, shm peak.

### R2 GUI lane notes (2026-10-12)

- Long-branch rule: `phylo.long_branch_breaks(root, factor)` (Qt free, `BranchBreaks`), called by `RaxmlTreeView` on the drawn rooted tree with `LONG_BRANCH_BREAK_FACTOR = 10` (`tree_view.py`). It replaces the old neoTree rule in `phylo.build_layout` (`LONG_BRANCH_SHARE`, `broken_node`, `shorten_long_branch`), which no app code used. A broken edge is drawn at 0.5 * H with a screen-stable break mark (`_BranchBreakItem`: a gap and two slanted strokes) and its true length ("0.00513 subst/site") above it; for the lower of two broken sibling edges (a midpoint root on the long branch) the label goes below, so the two labels never overlap.
- Toggle: the checkbox reads "True lengths" (not "True branch lengths": the Classic panel's toolbar is too narrow); its tooltip says "True branch lengths ...". Saved as the panel param `true_lengths: "on"` (absent when off, so existing saved states are unchanged); it is a display option, so layouts keep it.
- Reference run, current `core_recomb` (midpoint): m = 2, H = 3.81e-5 subst/site, the root edges to lineage 1 (6 genomes, 0.00513) and lineage 2 (147, 0.00508) are broken. `gubbins_final` also breaks its two root edges; `snippy_core` breaks nothing (its longest edge is about 1 times H).
- Header (trees.json schema 2): the chip shows when `gubbins.resolved_mode == "per_lineage"` and a lineage with `ran: true` has `converged` not true ("Gubbins did not converge (lineage 2)"; tooltip "Lineage 2 stopped after 20 of 20 iterations."); otherwise the schema 1 rule and text. `samples_unmasked` gives "N genomes not recombination-masked" with the names in the tooltip.

### R2 Pipeline lane notes (2026-10-12)

Implemented as planned in R2.1 to R2.11, except where noted. Nothing was rerun on the reference run.

- **Code:** `GubbinsWorkflow.py` (rewritten: `resolve_gubbins_groups`, `stage_group_alignments`, `split_threads`, parallel runs with `_Processes` (a failing group terminates the others), `merge_masked_alignment`, `write_combined_gff`, schema 3 summary, `gubbins_reuse_issue`, `expected_gubbins_files`), `CoreTreeWorkflow.py` (`ALGORITHM_VERSION` 2: IQ-TREE on `gubbins/core_gene_alignment.masked.aln`; the fconst and +ASC code is deleted), `TreesIndex.py` (schema 2, `gubbins_index_block`, `samples_unmasked`, `gubbins_final` only in single mode), `RunValidation.py` (schema 2 and 3 layouts; the unused `gubbins_samples_left_out` is deleted), `Runner.py` (fastGEAR before Gubbins; a fastGEAR failure lets gubbins, core_tree and Snippy finish and then fails the run; `expected_outputs.json` rewritten with the resolved groups when the gubbins step starts), `ToolEnv.py` (`--shm-size`, the new `phylo` keys), `CoreAlignment.fasta_record`.
- **Reuse:** a group's stamp is keyed by its member set plus the ordered alignment's sha256 (the staged input is a deterministic function of both, so the plan check needs no staging). The staged input's sha256 is still recorded and compared when the step runs. A step that reruns reuses every unchanged group; masking and the merge are always redone (13 s on the cohort).
- **Merge check:** stricter than R2.5 step 3: the changed columns of each masked row must equal its own blocks (minus columns already N) column by column, not only in count.
- **Gubbins' missing-data filter:** Gubbins 3.4.3 logs "Excluded sequence X ... missing data" but X stays in the final tree and the GFF (checked on a scratch subset with `gubbins_filter_percentage` 10). The pipeline follows the outputs: a genome is `dropped_by_gubbins` only when it is absent from its group's final tree, so this path is covered by the unit test only.
- **Single mode:** Gubbins reads the ordered alignment itself (no staged copy); `gubbins/core_gene.*`, `gubbins_final` and `samples_left_out` stay as before.
- **Scratch checks** (`<scratch>/phylo2/pipeline/`, subsets of the reference run's ordered alignment and fastGEAR table, the pipeline's own functions; `verify_mask.py` re-derives each genome's masked columns from the GFFs independently):
  - 20 genomes (6 lineage 1, 7 ST492 and 7 ST2 of lineage 2): per lineage, 1 and 15 threads, both converged (4 and 6 iterations, `recombination`), 2 and 264 blocks; masked rows match their own blocks exactly, none masked by the other lineage's 1,926 or 783,191 columns. IQ-TREE 4.8 min. A rerun reuses everything. Swapping the lineage ids reuses both runs (folders renamed), gives a byte-identical masked alignment and reuses the core tree. Without fastGEAR, auto falls back to single with `gubbins_lineages_unavailable` and removes the lineage folders; explicit `per_lineage` without fastGEAR is refused.
  - 12 lineage 2 genomes: auto resolves single (`one_lineage`), top-level layout, `gubbins_final` listed.
  - 4 lineage 1 plus 12 lineage 2 genomes: lineage 1 skipped (`below_min_lineage_size`, `gubbins_lineage_skipped`), its 4 genomes unmasked and listed in `samples_unmasked` and the tree notes.
  - Full cohort (a byte-identical scratch copy of the ordered alignment, 64 Gubbins threads, 16 IQ-TREE threads): lineage 1 (1 thread) 5 s, converged, 2 blocks (13 variable columns). Lineage 2 (63 threads) 22.2 min, **did not converge** at 20 iterations even with `recombination`; 946 blocks, SNPs inside/outside 109,417 / 2,619, r/m 41.8. No "Bus error"; the peak `/dev/shm` use sampled every 5 s was 149 MB (default 8g). Masking and merge 13 s; mean masked fraction 0.20. IQ-TREE 28.1 min: K3Pu+F+I+G4, logL -3,832,198.17, 39,717 variable columns, UFBoot converged, median 88, 40% at 95 or more. Between-lineage over within-lineage-2 mean patristic distance 487 (370 yesterday); within lineage 2 Pearson and Spearman 0.995 against yesterday's `core_recomb`; RF 74 (normalized 0.247) to it. IQ-TREE's composition test now fails for 140 of 153 sequences; the stamp keeps the full list, the tree notes show 20 names.
- **Docs:** CLAUDE.md phylogeny line updated.

### R2 Pipeline lane fixes after review (2026-10-12)

- **expected_outputs.json after the gubbins step:** the Runner rewrites it from the step's summary (run and reuse branches): `expected_gubbins_files(summary)` and the core tree id from the groups that ran, so a lineage skipped later for too few variable columns (or every group skipped, tree `core`) no longer leaves expected files that never appear (a finished run read as partial).
- **Thread budget:** `split_threads` gives 1 thread each and splits the rest of the budget by variable columns, so the sum is the budget (the old floor of 1 could exceed it). `_run_parallel` starts the jobs largest first, each once its threads fit in the budget, so many small lineages queue instead of starting together; `mask_gubbins_aln.py` jobs (each loads its group in memory) run at most 4 at once. The per-container `--memory` and `--shm-size` are ceilings, not usage; lineage inputs are disjoint row subsets of the cohort alignment, so running lineages together needs about the memory of one cohort run.
- **Stopping the other lineages:** a docker command gets `--cidfile`; on a failure the others' containers are killed with `docker kill` (terminating the `docker run` client alone left the container running). Checked: a 120 s container stopped 8 s in when its sibling failed, no container left.
- **Stamps:** single mode also writes its run's stamp to `gubbins/run.json` as soon as Gubbins finishes, and per lineage each `lineage_<id>/run.json` is now written when that lineage finishes (not after all of them), so a step that fails later (mask, merge, tables) or a sibling failure keeps the finished runs for the resume. Checked: single mode with the merge forced to fail after Gubbins (86 s), the resume reused the run (8 s). `planned_gubbins_files` lists `run.json` in both modes.
- **Lineage renumbering:** reused entries get their warnings rebuilt for the current id and settings (`group_warnings`). Rejected part: a pure relabel still reruns the step (staging, masking, merge; Gubbins itself is reused), because the combined GFF, the block table, `masked_bases.tsv` and the summary carry the lineage ids. Checked: ids swapped, both runs reused, masked alignment byte-identical, warnings name the new ids.
- **Core tree notes on reuse:** the core tree stamp keeps its IQ-TREE notes apart (`iqtree_notes`); when the tree is reused, `refresh_gubbins_notes` rewrites the RF notes, masking notes and `samples_unmasked` from the current Gubbins summary (no IQ-TREE rerun), so a relabel or a changed unmasked set with a byte-identical masked alignment no longer leaves stale ids or genomes in trees.json. Checked on the swap above (RF notes now name the swapped lineages).
- **masked_bases.tsv:** `masked_bases` counts masked A/C/G/T only, the new `masked_positions` also counts gaps turned into N, `masked_fraction` is over the genome's called bases. On the 20-genome check gaps were up to 80% of a lineage 1 genome's masked positions. A block on an internal Gubbins branch masks every genome below it, so a high fraction is not "recombined more".
- **Removed:** the dead `invariant_site_correction` Gubbins parameter (stamps from before this fix rerun once), the top-level `longest_branch_share` in the gubbins step info (now per group).
- **IQ-TREE reproducibility (not a bug):** same seed and input at 16 threads can change poorly supported splits (seen inside lineage 1, UFBoot 2 to 44) and so the leaf order; supported clades are stable. Documented in the CoreTreeWorkflow docstring; the core tree is reused whenever its input is unchanged, so this shows only after a real rerun.
- **Checks:** `tests.test_pipeline_phylogeny`, `test_pipeline_genome_set`, `test_tree_sources` pass; 20-genome scratch run (`scratchpad/phylo2/fix_pipeline/runs/A`): per lineage rerun gives a masked alignment byte-identical to the verified one, IQ-TREE rebuilt (median UFBoot 96), validators clean, switch to single and back works; the reference run still validates with 0 errors (read only).

### R2 GUI lane fixes after review (2026-10-12)

- **Toolbar overlap (fixed):** the tree toolbar (`_TreeToolbar` in `panaroo_core_accessory.py`) hides the captions that do not fit, the "Tree" title before "Branch scale" (the scale dropdown keeps its accessible name and tooltip). Classic's 440 px Phylogeny panel: no widget overlaps, "True lengths" and "Export..." whole (before: checkbox over the dropdown arrow); Pangenome (517 px) keeps "Tree". The toggle stays beside the Distance choice as R2.9 says.
- **Notes row clipped (fixed):** the notes under the tree source (`_NotesRow`) stand one under the other when they do not fit side by side, so "Gubbins did not converge (lineage 2)" and "1 genome not recombination-masked" are shown whole in Classic.
- **Break labels on one path and merged marks (fixed):** broken edges are taken top to bottom, then left to right; their labels go above and below in turn and their marks sit at 35% and 65% of the edge in turn (one break: middle, label above). Root "fastGEAR clades" on the reference tree (two broken edges on one row) and the two root edges in Classic's short panel now show two separate marks and labels. The C2 clade chip near the root is the existing chip placement, not changed here.
- **SNP labels (fixed):** a true length of 100 or more is shown to one decimal with thousands separators ("36,991.5 SNPs"), below 100 at 3 significant digits as before (substitutions/site).
- **Unmasked tooltip (fixed):** lists the core tree's notes that give the reason per genome group ("6 genome(s) not recombination-masked (lineage not run by Gubbins): ..."); without them it names the genomes without a reason.
- **Convergence not reported (fixed):** a lineage that ran with `converged` null still counts (the pipeline counts it as not converged too), but the tooltip says "Lineage 1: convergence not reported (4 of 20 iterations)." instead of "stopped after".
- **Dead layout code (fixed):** `phylo.build_layout`, `TreeLayout`, `_depths` and `group_runs` had no app caller (the Phylogeny panel lays out in `tree_view`); deleted with `test_mrca_and_group_runs`. `test_layout_matches_the_tree_view_order` now compares `prepare_tree`'s leaves with the view. Module docstring corrected.
- **Rejected: per-lineage fallback when `lineages` is missing.** `TreesIndex.gubbins_index_block` always writes `lineages` as a list, and in per-lineage mode the ids are fastGEAR lineage ids (never null). An empty list falls back to the top-level `converged`, which the pipeline aggregates over the lineages that ran, so the chip is still right.
- **Checks:** targeted tests (`test_tree_view`, `test_tree_sources`, `test_pan_panels`, `test_layout_migration`, `test_retired_shims`) 126 OK. App check on the reference run (Classic with Root midpoint, fastGEAR clades, as in file, `gubbins_final`; Pangenome): 0 uncaught exceptions, no toolbar overlaps. Synthetic per-lineage trees.json (lineage 2 not converged; lineage 1 skipped with a reason note; converged null): chip, tooltips and note layout as above. Scripts and screenshots: `scratchpad/phylo2/fix_gui/`.

## Rework 2 reference run rerun (2026-10-12)

Run `runs/sample_78eaa9126850/run_78eaa9126850`, through the pipeline's own Runner
(`<scratch>/phylo2/rerun/driver.py plan|run`, yesterday's `build_cfg()`: GUI Resume
`DataInput`, COG dropped, 96 cores). Default `phylo` settings, no config overrides:
`gubbins_mode` auto, `gubbins_min_lineage_size` 5, `gubbins_converge_method`
`recombination`, 20 iterations, seed 1; IQ-TREE `MFP`, UFBoot 1000 with `--bnni`; docker
`--shm-size 8g` (the `container.shm_size` default; `resources.max_ram_gb` unset).

**Backup** (copies, verified by file count, byte total of the files and `diff -rq`):
`/home/bals/Documents/dev/BACTolomew_backups/2026-10-12_gubbins_per_lineage/`:
`run/{gubbins,trees,core_alignment}` (15, 15 and 3 files; 159,991,756, 530,969 and
409,025,762 bytes), `run/{run_config,run_state,run_report,expected_outputs}.json` (`cmp`
equal), the GUI cache `cache/run_78eaa9126850` (5 files, 469,414,689 bytes) and
`manifest_before.tsv` / `manifest_after.tsv` (path, size, mtime of all 12,479 files of the
run folder before). The 2026-10-10 backup keeps the original state.

**Plan check** (Runner's own reuse predicates, before the run): genome set unchanged;
assembly_qc, bakta, mlst, amrfinder_mutations, mobsuite, panaroo, core_alignment,
fastgear_clades, phylogeny (Snippy), chewbbaca, diamond_markers, diamond_reference_loci and
the marker hit search reused; gubbins ("the Gubbins output predates per-lineage runs and
the masked core alignment") and core_tree ("the core tree step changed") rerun. Resolved
plan: per lineage (`two_or_more_lineages`), lineage 1 (6 genomes) and lineage 2 (147), both
run. `run_state.json` afterwards agrees: every other step `reused_existing_output`,
collect_annotations an in-memory re-merge with the marker hit TSVs reused.

**What ran, wall times:**

| Step | Result | Wall time |
|---|---|---|
| steps before gubbins | reused | 11 s |
| gubbins: lineage 1 (1 thread) | ran, in parallel with lineage 2 | 6.3 s |
| gubbins: lineage 2 (63 threads) | ran | 1,325.6 s (22.1 min) |
| gubbins: `mask_gubbins_aln.py` per lineage + merge | ran | 14 s |
| gubbins step (staging, runs, masking, tables) | success | 1,352 s (22.5 min) |
| core_tree (IQ-TREE 2.4.0, 16 threads) | ran | 1,634 s (27.2 min); IQ-TREE itself 1,611 s |
| phylogeny (Snippy), chewbbaca, DIAMOND | reused | |
| collect_annotations | in-memory re-merge, marker hits reused | 126 s |
| run_workflow | success | 3,143 s (52.4 min) |
| GUI cache rebuild (`rebuild_cache.py`) | needed (`run_config.json` is in the fingerprint) | 241 s |

Against the plan's estimate (R2.13: 40 to 85 min) the run took 52 min; yesterday's rerun
took 104 min, of which 41 min was fastGEAR (reused today). A fresh run would add fastGEAR,
so about 93 min in total, a little less than yesterday's 104 min.

**Docker shared memory:** no "Bus error" in either lineage's logs
(`gubbins/lineage_*/command.*.log`). The peak `/dev/shm` use in the lineage 2 container,
sampled every 10 s (131 samples), was 149 MB (148,725,760 bytes) against the 8 GB ceiling
(about 58 times the peak; the R2.8 criterion asked for at least 4 times). IQ-TREE used none.

**Gubbins per lineage** (`gubbins/summary.json` schema 3, `resolved_mode` per_lineage,
Gubbins 3.4.3, converge method `recombination`):

| | Lineage 1 | Lineage 2 |
|---|---|---|
| Genomes | 6 | 147 |
| Variable columns in its rows | 13 | 41,227 |
| Iterations | 4 of 20 | 20 of 20 |
| Converged | yes ("Recombinations observed before") | **no** (maximum number of iterations reached) |
| Recombination blocks | 2 | 946 |
| SNPs inside / outside recombination | 10 / 3 | 109,417 / 2,619 |
| r/m (SNP ratio) | 3.33 | 41.78 |
| Longest-branch share of clonal SNPs | 0.33 (1 SNP) | 0.044 (GCA_002017515.1) |
| Masked A/C/G/T bases (all genomes) | 1,845 | 76,755,677 |
| Masked positions incl. gaps set to N | 8,737 | 78,575,812 |
| Masked fraction per genome (of its called bases) | mean 0.00012, max 0.00019 | mean 0.204, median 0.241, max 0.249 |
| Top tip recombination fractions | | GCA_002017515.1 0.147, GCA_034330085.1 0.123, GCA_014267385.1 0.118 |

Lineage 2 again did not converge at 20 iterations, as in the full-cohort scratch run, with
the same numbers (946 blocks, r/m 41.8): Gubbins is deterministic for the same input and seed.
Run warnings: only `gubbins_not_converged` (lineage 2). `gubbins_lineage_structure` is gone,
because no branch within a lineage holds most of the clonal SNPs. No genome was left out or
left unmasked (`samples_unmasked` empty).

**Masking checked independently** (`<scratch>/phylo2/pipeline/verify_mask.py` on the run):
153 genomes, 2,627,569 columns, 78,584,549 positions newly set to N; every genome's newly
masked columns equal the union of its own lineage's GFF blocks exactly (0 mismatches). The
1,166,115 columns blocked only in lineage 2 are masked in no lineage 1 genome, and the 197
columns blocked only in lineage 1 in no lineage 2 genome.

**Core tree (`core_recomb`):** IQ-TREE 2.4.0 on `gubbins/core_gene_alignment.masked.aln`
(2,627,569 sites, of them about 2.588 million constant (98.49%), 38,320 parsimony-informative, 11,427
distinct patterns; no `-fconst`, branch lengths in substitutions/site). Best-fit model
K3Pu+F+I+G4 (BIC), log-likelihood -3,832,198.17 (identical to the full-cohort scratch run).
UFBoot converged; 150 internal splits, median support 88, 59 (39%) at 95 or more, 48 below
70. 140 of 153 sequences fail IQ-TREE's composition chi2 test, compared with yesterday's 6
lineage 1 genomes: the masked rows have about 20% N, so the test is not informative here,
and the notes keep the list. One pair of identical sequences (GCA_031156395.1 and
GCA_029679835.1). RF to the lineage Gubbins trees: lineage 1 4 (0.667, 6 genomes),
lineage 2 60 (0.208).

| new `core_recomb` vs | RF (normalized) | within ST492 | within ST2 | new UFBoot >= 95 splits absent there | its UFBoot >= 95 splits absent in new |
|---|---|---|---|---|---|
| yesterday's `core_recomb` (cohort Gubbins, fconst) | 74 (0.247) | 56 (0.311) | 18 (0.158) | 2 | 1 |
| yesterday's cohort Gubbins final tree | 82 (0.273) | 54 (0.300) | 26 (0.228) | 1 | 0 |
| old Gubbins final tree (2026-10-10 backup) | 138 (0.460) | 68 (0.378) | 68 (0.596) | 10 | 0 |
| old `core_gene_tree.nwk` (RAxML-NG, unfiltered) | 162 (0.540) | 88 (0.489) | 70 (0.614) | 15 | 0 |
| Snippy tree | 172 (0.573) | 90 (0.500) | 80 (0.702) | 17 | 0 |

**Clades** (midpoint root, as in `trees/core_recomb.nwk`): lineages 1 and 2 and fastGEAR
clusters 1, 2 and 4 are monophyletic with UFBoot 100; cluster 3 is not (as before).
Cluster 4 (83 ST492 genomes) sits inside lineage 2 in a 94-genome clade (UFBoot 94) with 8
ST2 and 3 ST492 genomes of cluster 3; a 46-genome ST2 clade (UFBoot 100) is its sister, so
"cluster 4 inside ST2" holds. GCA_008987025.1 is again sister to cluster 4 (UFBoot 97;
98 yesterday). GCA_026558295.1, GCA_026563615.1 and GCA_026566215.1 are in the same 4-genome
clade as yesterday (UFBoot 96).

**Distances:** the ratio of the mean between-lineage patristic distance to the mean within
lineage 2 is 487 (370 yesterday); within lineage 2 the mean is 3.31e-5 subst/site
(2.8e-5 yesterday). Pearson and Spearman of the within-lineage-2 patristic distances
against yesterday's `core_recomb`: 0.995 and 0.995. Pairwise SNP distances
(`<scratch>/phylo2/metrics/snp_dist.py`: columns where both genomes have A/C/G/T and
differ), masked against the unmasked ordered alignment:

| Set | Pairs | Unmasked median (p5 to p95) | Masked median (p5 to p95) | Pairs <= 10 / 20 / 50 SNPs, unmasked | masked |
|---|---|---|---|---|---|
| ST492 in lineage 2 (87) | 3,741 | 2,704 (60 to 8,932) | 30 (10 to 99) | 13 / 32 / 127 | 195 / 828 / 3,309 |
| ST2 (60, all lineage 2) | 1,770 | 2,834 (362 to 8,593) | 69 (30 to 144) | 2 / 4 / 16 | 8 / 29 / 339 |
| ST492 all (93, with the 6 of lineage 1) | 4,278 | 2,758 | 32 (max 34,718) | 27 / 47 / 142 | 210 / 843 / 3,324 |

The maximum of the third row comes from pairs across lineages: recombination between
lineages is not masked (trees.json note). Any SNP threshold on this cohort must use the
masked distances, and only within a lineage.

**trees/trees.json** (schema 2): default `core_recomb`, rooting midpoint, ordering genome
GCA_003931755.1; `gubbins` block with `mode` auto, `resolved_mode` per_lineage,
`converged` false, `lineages` [1: ran, converged, 4 of 20, 2 blocks, r/m 3.33; 2: ran, not
converged, 20 of 20, 946 blocks, r/m 41.78]; warning `gubbins_not_converged`. Trees:
`core_recomb` (`alignment` `gubbins/core_gene_alignment.masked.aln`, `alignment_kind`
`core_alignment_masked`, `sites` 2,627,569, `invariant_correction` null, `samples_unmasked`
empty, notes: composition test, identical sequences, RF to each lineage's Gubbins tree, the
inter-lineage limitation) and `snippy_core`. `gubbins_final` is no longer listed, and its
`trees/gubbins_final.{nwk,json}` were removed; the per-lineage Gubbins trees stay in
`gubbins/lineage_<id>/`.

**Validation:** `validate_completed_run` 0 errors (15 steps); run_state and run_report
status success, outputs validated, warning `gubbins_not_converged` only;
`expected_outputs.json` names the per-lineage files and the masked alignment, and every
`gubbins/` and `trees/` path in it exists. Manifest diff: only `gubbins/` (29 added, 3
changed, 12 removed: yesterday's top-level `core_gene.*` and logs), `trees/` (10 changed,
2 removed) and the four run JSON files changed; `core_alignment/`, `fastgear/`,
`snippy/`, `panaroo/`, `bakta/`, `chewbbaca/`, `annotations/`, `dmnd/` and the rest are
identical in size and mtime. The rebuilt GUI tables equal the backup's
(`polars.DataFrame.equals` on all three parquet files).

**App check** (`<scratch>/phylo2/rerun/app_check.py`, offscreen, harness QSettings): 0
uncaught exceptions. Classic Phylogeny on `core_recomb`, root "As in file (midpoint)": chip
"Gubbins did not converge (lineage 2)" (tooltip "Lineage 2 stopped after 20 of 20
iterations."), no unmasked note; both root edges broken (lineage 1, 6 genomes, "0.00809
subst/site"; lineage 2, 147 genomes, "0.00804 subst/site"); "True lengths" draws them at
full length (the tree is then a thin column at the right) with the leaf order unchanged; a
PNG export keeps the break marks and labels, and its scale bar reads "(broken branches not
to scale)". The Gene presence/absence rows are linked and equal the tree's 153 leaves. Row
order against yesterday's tree (same `prepare_tree` and groups): Spearman 0.992, 84 of 153
rows in the same place.

**Open points:** lineage 2 still does not converge at 20 iterations with the
`recombination` method; more iterations (`phylo.gubbins_iterations`) is the remaining
option, not tried. Scripts and logs: `<scratch>/phylo2/rerun/` (driver, run.log,
shm_samples.txt, cache.log, app check screenshots) and `<scratch>/phylo2/metrics/`.

## Rework 2 before and after (rechecked)

Adversarial recheck (2026-10-12) of the per-lineage impact report. Every number was
recomputed from the current run, the 2026-10-12 backup (cohort Gubbins, `core_recomb` with
fconst) and the 2026-10-10 backup (original state) with independent scripts (own FASTA and
GFF parsing, numpy pair distances, Biopython trees; the repo's `phylo.prepare_tree` and
`long_branch_breaks` only for the GUI row order, cross-checked by an own implementation of
the break rule). Scratch: `<scratch>/phylo2/impact_check/` (`gub.py`, `snp.py`, `trees.py`,
`conf2.py`, `gui.py`, `c4.py`, `depth.py` and their `.txt` outputs). Numbers that held are
listed briefly; corrections are marked **Corrected**.

**Gubbins.** Confirmed from `per_branch_statistics.csv` and the GFFs: original 4,129 blocks,
88,918 / 32,894 SNPs in / out, r/m 2.70, 5 iterations (random seed); cohort 1,246 blocks,
115,524 / 31,533, r/m 3.66, 20 of 20 iterations, not converged; root edge 8,512 / 28,935
(r/m 0.29), rest 107,012 / 2,598 (r/m 41.2). Per lineage: lineage 1 converged after 4
iterations ("Recombinations observed before"), 2 blocks, 10 / 3, r/m 3.33; lineage 2 20 of
20, not converged, 946 blocks, 109,417 / 2,619, r/m 41.78; total 948 blocks (the combined
GFF has the same 948 blocks; no block names a genome outside its own lineage).
**Qualified:** lineage 2 per lineage and the cohort run's within-lineage-2 part agree in
counts (946 against 950 blocks off the cohort root edge, which alone held 296 blocks; SNPs
in / out 109,417 / 2,619 against 107,012 / 2,598), but not in where the blocks lie: the
per-genome masked share has Spearman 0.917 only, and lineage 2 genomes are more masked per
lineage (mean 0.204 against 0.176) although the cohort run also masked its 296 root-edge
blocks in every lineage 2 genome.

**Masking.** Confirmed: rebuilding the per-lineage mask from the two lineage GFFs reproduces
`gubbins/core_gene_alignment.masked.aln` exactly (0 differing cells; 78,584,549 cells set to
N, 76,757,522 of them A/C/G/T). Masked share of called bases: original mean 0.190, median
0.218, max 0.258; cohort 0.169 / 0.202 / 0.216; per lineage 0.196 / 0.241 / 0.249; ST492 in
lineage 2 median 0.242 / 0.212 / 0.243; ST2 median 0.141 / 0.120 / 0.156; per-genome
Spearman per lineage against cohort 0.917, against original 0.818. Called bases per genome
are identical in Panaroo's and the ordered alignment.
**Qualified:** the 1,845 masked lineage 1 bases are the total over the six genomes (340 to
485 in five, 0 in GCA_023499005.1), mean share 0.00012. "Most of it is old imports shared by
cluster 4" holds: in each cluster 4 genome 97.6% to 100% (median 99.1%) of its masked bases
come from blocks whose taxa include all 83 cluster 4 genomes (98 such blocks, 17 of them on
the cluster 4 stem itself).

**Pairwise SNP distances** (both genomes A/C/G/T after masking, and different). Confirmed
exactly for ST492 in lineage 2 (87, 3,741 pairs), ST2 (60, 1,770) and cluster 4 (83, 3,403):
medians, p5 / p95, maxima, pair counts at 10 / 25 / 50 SNPs, single-linkage groups,
crossing pairs and Spearman, as in the impact report (for example ST492 in lineage 2 per
lineage median 30, pairs 195 / 1,383 / 3,309, groups at 10 SNPs 10 / 63 / 19; crossing
against cohort 44 / 174 / 47, Spearman 0.988, median absolute difference 1).
**Corrected:**
- "Uncorrected distances inflate by about 70 times": in ST492 of lineage 2 the per-pair
  median ratio is 74 and the ratio of medians 90 (cluster 4: 69 and 95); in ST2 it is 44
  and 41. Say "about 70 to 90 times in ST492, about 40 times in ST2".
- "The original masking roughly doubled the distances": 2.2 times in ST492 of lineage 2
  (median 67 against 30), 1.6 times in ST2 (111 against 69).
- Comparable-site share, all pairs instead of 300 random ones: cluster 4 median 75.6% (74.6%
  to 76.2%), ST2 83.8% (76.5% to 99.4%).

**Middle East trio.** Confirmed: masked, GCA_026558295.1 is 29 SNPs from GCA_016536005.2
(ST2) and 39 (median 50) from cluster 4; the Saudi pair is 3 SNPs apart and 86 / 89 from
the nearest cluster 4 genome. Tree placement (4-genome clade UFBoot 96, Saudi pair 100,
joining cluster 4 at the 90-genome node, UFBoot 54 against 70 in the cohort tree) confirmed.
**Corrected:**
- "7,648 to 10,677 SNPs from cluster 4" mixed a minimum and a median: unmasked, the nearest
  cluster 4 genome is 7,648 to 9,880 SNPs away (medians 8,747 to 10,677).
- Masked, the Saudi pair's nearest genome outside the pair is also ST2 (GCA_016536005.2, 81
  and 84 SNPs), closer than any cluster 4 genome (86 and 89). By median, though, the trio is
  closer to cluster 4 (50 to 103) than to ST2 (70 to 124), because ST2 is diverse.
- "Different capsule loci, KL22 and KL234, from cluster 4's KL2" is wrong: cluster 4 is
  KL104 (45 of 83), KL2 (18), KL4 (8), KL120 (8) and KL9 (4). The trio's loci (KL22, KL234)
  occur in no cluster 4 genome, so the capsule argument holds, but it must not name KL2 as
  cluster 4's locus.
- ST492 is not monophyletic, and besides lineage 1 and the trio a fourth lineage 2 ST492
  genome sits apart: SAMN04040115 (fastGEAR cluster 1, USA, KL22). The ST492 clades
  83 / 6 / 2 / 1 / 1 are cluster 4, lineage 1, the Saudi pair, GCA_026558295.1 and
  SAMN04040115.

**Trees.** Confirmed exactly: normalized RF of every pair of the six trees, whole tree and
within ST492 / ST2 (per lineage against cohort 0.247 / 0.311 / 0.158, against cohort
Gubbins final 0.273, original Gubbins 0.460, original core 0.540, Snippy 0.573); strong
conflict totals (per lineage against cohort 2, cohort against per lineage 1, against
original Gubbins 10, original core 15, Snippy 17; cohort against them 10, 16, 15); support
counts (150 splits, 59 at 95 or more, 43 at 70 to 94, 48 below 70; cohort 61 / 38 / 51,
median 86); clades, monophyly of lineages and clusters 1, 2, 4, enclosing clades around
cluster 4 (84 at 97, 90 at 54, 94 at 94), ST2 in 10 clades (largest 46), 16 / 19 / 25 in the
original trees.
**Corrected:**
- The per-lineage median UFBoot is 87.5, not 88.
- The within-ST492 split of the strong conflicts depends on which support a pruned split
  keeps when several original splits collapse into it (minimum or maximum): per lineage
  against original Gubbins 2 to 4 (ST2 7), original core 5 to 7 (ST2 8), Snippy 5 to 7 (ST2
  11); cohort against the same trees 2 to 4, 6 to 8, 3 to 5 (ST2 7, 8, 11). The impact
  report's 3 / 6 / 6 and 3 / 7 / 4 fall inside these ranges, and the earlier recheck's 4 / 8
  / 5 is the maximum rule. Report the totals; the within-ST492 parts only as ranges.

**Lineage depth.** Confirmed: mean between-lineage over mean within-lineage-2 patristic
distance 487 (cohort 370, original Gubbins 283, original core 12.9, Snippy 6.8);
within-lineage-2 Spearman against per lineage: cohort 0.995, original Gubbins 0.80,
original core 0.77, Snippy 0.73.
**Corrected:** "lineage 2 shrinks further" is wrong in absolute terms. The mean within
lineage 2 distance grew from 2.77e-5 to 3.31e-5 subst/site (+20%); the ratio rose because the
between-lineage distance grew more, from 1.02e-2 to 1.62e-2 (+58%; root edges 0.0051 to
0.0081). This fits two changes of this round: the cohort run's 296 root-edge blocks are no
longer masked (inter-lineage recombination is not masked by design), and IQ-TREE now sees
the masked columns as missing data instead of the fconst counts that treated them as
constant. The two effects were not separated (that needs a rerun).

**GUI.** Confirmed with the repo's `prepare_tree` (lineage groups) under the file root,
midpoint and "fastGEAR clades" roots: lineages and clusters 1 to 4 are one run each in both
corrected `core_recomb` trees, cluster 3 is split into 2 runs in the original Gubbins, core
and Snippy trees and in the cohort Gubbins tree; row order against per lineage: cohort
`core_recomb` Spearman 0.992 (84 of 153 rows in place; 82 under "fastGEAR clades"),
original trees 0.931 / 0.931 / 0.932. Long-branch rule (factor 10, own implementation and
`long_branch_breaks` agree): the 2 root edges, 0.00809 and 0.00804 subst/site, reference
depth 4.33e-5; the deepest tip is 187 times that depth. Cluster 4 sits next to cluster 3 in
the row order (83 rows, then 60 of cluster 3), which is why cluster 3 stays one band.
`trees.json` lists only `core_recomb` and `snippy_core` and records lineage 1 converged (4)
and lineage 2 not (20).

**Manuscript scripts.** Confirmed: the four tree R scripts read the missing
`gubbins/core_gene.final_tree.tre`. **Corrected:** `build_supplementary_raw_outputs.py`
does not just skip the missing files: `copy` returns None and `note(None, ...)` raises
AttributeError, so the script stops at the first missing file in section 07
(`gubbins/core_gene.final_tree.tre`). `gubbins/summary.json` exists and would be copied, but
its description ("per-iteration summary") no longer fits schema 3.

**Not rechecked:** the GUI chip text and screenshots (taken from the rerun's app check, not
rerun here). Nothing in the run was rerun or edited.

## Branch lengths including recombination (2026-10-12)

A second set of branch lengths for the corrected tree, same topology: not a new tree search
and not the uncorrected tree.

**Pipeline (`CoreTreeWorkflow.run_recombination_refit`, end of the `core_tree` step).**
IQ-TREE `-te trees/core_recomb_iqtree/core_recomb.treefile` (the raw unrooted tree) on the
unmasked `core_alignment/core_gene_alignment.ordered.aln`, `-m <the core tree's model>`
(no ModelFinder, no bootstrap), `-keep-ident`, `--seed <phylo seed>`, `-T <core_tree
threads>`. Outputs in `trees/core_recomb_refit/` (`core_recomb_with_recombination.treefile`,
`.iqtree`, `.log`, `run.json`). The result must have the core tree's leaves and splits, or it
is refused. `run.json` is its own reuse stamp (`key`: sha256 of the input tree and of the
unmasked alignment, model, seed, IQ-TREE version, algorithm version; an unreadable IQ-TREE
version matches the rest of the key) with `status` built, skipped (`skip_reason`: `disabled`,
`gubbins_skipped`, `no_iqtree_tree`, `nothing_masked`) or failed; the core tree's own reuse
check ignores it. Any failure is the warning `core_tree_refit_failed` (the core_tree step
still succeeds) and keeps a refit built before; a later run retries. Config:
`phylo.recombination_branch_lengths` (default true). Step info:
`recombination_branch_lengths`.

**trees.json (schema_version stays 2).** When the refit is current (built from the current
raw tree and unmasked alignment, `CoreTreeWorkflow.current_refit_stamp`), `write_trees_index`
writes `trees/core_recomb_with_recombination.nwk`: rooted on the same bipartition as
`core_recomb.nwk` at the same share of the root edge (`Newick.root_on_split`), its supports
copied by split, its children ordered like `core_recomb.nwk`. The `core_recomb` entry then
gains:

```
"branch_length_sets": [
  {"id": "mutations", "label": "Mutations only", "file": "core_recomb.nwk", "default": true,
   "tree_length": ...},
  {"id": "with_recombination", "label": "Including recombination",
   "file": "core_recomb_with_recombination.nwk", "default": false,
   "alignment": "core_alignment/core_gene_alignment.ordered.aln", "model": "...",
   "branch_length_unit": "substitutions/site", "method": "...", "log_likelihood": ...,
   "tree_length": ...}
]
```

Absent key or a missing file: no choice (old runs unchanged). Without a current refit the
key and the file are removed. The validator checks a current refit's files and listing.

**GUI (Phylogeny toolbar).** "Branch lengths: Mutations only | Including recombination"
(`panaroo_core_accessory._LengthsSwitch`), shown when the drawn tree lists two or more sets
whose files exist, disabled on the Topology scale. Only the lengths change
(`RaxmlTreeView.set_branch_lengths`): the root (midpoint too), the Order and the clusters are
computed on the mutations only tree, so rows never move; each drawn edge takes the other
fit's length of the same split, the root edge split in the default tree's proportion. The
long-branch break and the scale bar follow the shown lengths; the scale bar label (kept at
its screen size) names them. A file that is not the same tree is refused: the mutations only
lengths are drawn and named, and the switch is disabled with the reason. The choice is saved
in the panel state (`branch_lengths`); exports name the lengths in the caption ("Branch
lengths: including recombination; root and leaf order from the mutations only lengths") and
the default file name.

### Branch length toggle: reference run rerun (2026-10-12)

Backup first: `BACTolomew_backups/2026-10-12_bl_toggle/` (`run/trees`, the four run JSONs,
the GUI cache, `manifest_before.tsv` / `manifest_after.tsv`; `diff -rq` and `cmp` equal).
Plan check with the Runner's own reuse predicates (`<scratch>/bltoggle/rerun/driver.py`):
every step reused (core_alignment, fastgear_clades, gubbins and core_tree included), only
the refit sub-step to build. Run through `Runner.run_workflow`: 549 s in all, of which the
refit (IQ-TREE 2.4.0, K3Pu+F+I+G4, `-te`, 16 threads) 382 s and collect_annotations (in-memory
re-merge) 128 s; GUI cache rebuild 238 s (`run_config.json` is in its fingerprint).
`validate_completed_run`: 0 errors; warnings: only `gubbins_not_converged` (lineage 2, as
before). Manifest: new files only under `trees/core_recomb_refit/` and
`trees/core_recomb_with_recombination.nwk`; `trees.json` gained the `branch_length_sets`;
`core_recomb.nwk`, `snippy_core.nwk` and the core tree's IQ-TREE files are byte identical.
A second plan afterwards reuses the refit.

Same leaves (153), splits, supports, leaf order in the file and root bipartition (lineage 1
| lineage 2, root share 0.5016 in both). Branch lengths (substitutions/site, mutations only ->
including recombination):

| Edges | Mutations only | Including recombination | Ratio |
|---|---|---|---|
| Whole tree (303 edges) | 0.017119 | 0.240703 | 14.1 |
| Lineage split (root edge, lineage 1 \| lineage 2) | 0.016124 | 0.143695 | 8.9 |
| Within lineage 2 (147 genomes, 292 edges) | 0.000994 | 0.096998 | 97.6 |
| Within cluster 4 (83 genomes, 164 edges) | 0.000291 | 0.017199 | 59.2 |
| Branch leading to cluster 4 | 0.000004 | 0.004046 | about 1,090 |
| Within lineage 1 (6 genomes, 10 edges) | 0.000001 | 0.000010 | 6.8 |

Median per-edge ratio 35. The masked recombination sits mostly inside lineage 2 (Gubbins
r/m 41.8 there), so the within-lineage branches grow far more than the lineage split; with
recombination included the root edge no longer needs the long-branch break.

GUI (app check, Pangenome preset, core_recomb, midpoint root, By clade size, Distance):
the switch shows under the toolbar; the leaf order is identical in both states (153 rows);
the scale bar reads "0.00001 substitutions/site, mutations only (broken branches not to
scale)" and "0.01 substitutions/site, including recombination"; the panel state saves
`branch_lengths`; on Topology the switch is disabled with "Branch lengths are drawn on the
Distance scale only.".
