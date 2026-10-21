# Pipeline fixes 2026-10-09 (P01 to P06)

> A dated log, kept as history. Current state (2026-10-21): BACTolomew no longer ships the
> CARD and AMRF marker databases that P02 and P04 filtered. Gene categories come from the
> `bakta_evidence` step (Bakta, AMRFinderPlus, MOB-suite), and marker databases are only the
> ones users add (Settings > Databases, `bactolomew db`). `pipeline/MarkerDbBuild.py` and the
> runtime guard still exist for databases built that way. The cgMLST (P01, P03) and genome-set
> (P05, P06) fixes still describe the current code.

## P02: CARD variant/knockout/overexpression models and AMRFinderPlus POINT references screened as acquired AMR genes

### What changed (code)
- New one-time DB build step `src/pipeline/MarkerDbBuild.py` (`PYTHONPATH=src .venv/bin/python -m pipeline.MarkerDbBuild [--db-root DB] [--dry-run]`). There was no script that built `DB/seq/markers` before; the FASTAs had been dropped in by hand. The step:
  - CARD_AMR: keeps protein homolog models (`HM|`), drops variant (`VM|`), knockout (`KM|`) and overexpression (`OM|`) models (as RGI does for presence calls). Records without a prefix (raw CARD broadstreet FASTA) are checked by ARO against `DB/metadata/CARD/card.json` model types.
  - AMRF_*: drops references whose `ReferenceGeneCatalog.txt` rows are all `POINT`, `POINT_DISRUPT` or `AMR-SUSCEPTIBLE`. Accessions that also have an acquired-gene row (AMR, BIOCIDE, METAL, ...) are kept; accessions missing from the catalog are kept and listed.
  - Archives the unfiltered FASTA once as `DB/seq/marker_sources/<stem>.unfiltered.fasta.gz` (deterministic gzip) and always filters from that archive, so a rebuild reproduces the same file byte for byte (verified: second run writes nothing, mtimes unchanged). Writes `DB/seq/marker_sources/reports/<stem>.removed.tsv` (header, accession, reason) and `DB/seq/marker_sources/marker_build_report.json` (rules, catalog version, per DB counts, sha256 of source and output). The folder is outside `DB/seq/markers` because `collect_fastas` globs that folder recursively.
- `src/pipeline/ResistanceMetadata.py`: new guard `drop_resistance_target_hits()` (CARD hits whose header prefix is VM/KM/OM or whose ARO only has non-homolog models in card.json; AMRF hits whose accession only has POINT/POINT_DISRUPT/AMR-SUSCEPTIBLE rows). `enrich_marker_hits_with_resistance()` calls it first, so such a hit can never be reported, even with an old DB or old DIAMOND TSVs; it logs one warning with counts per DB.
- `src/pipeline/Runner.py` (`collect_annotations`) and `src/annotation_explorer_polars/data.py` (GUI cache build) call the guard before `collapse_diamond_annotations`, so a valid lower-ranked homolog hit on the same CDS is not lost behind a dropped target hit.
- Stale `.dmnd` / hit TSVs (needed so the rebuilt FASTAs reach a run): `src/pipeline/DiamondWorkflow.py` gets `dmnd_is_current()` (sidecar `<stem>.dmnd.source.json` with FASTA size, mtime and sha256, written by `build_diamond_db_for_files`; legacy .dmnd without the sidecar is current only if not older than its FASTA), `hit_tsv_is_current()` (TSV not older than its .dmnd), and `build_diamond_db_for_files` / `run_diamond_annotations` rebuild instead of skipping stale files. `Runner.py`: the `diamond_markers` step reuses the .dmnd files only when `_marker_dmnds_current(CFG)`, and `collect_annotations` reuses the hit TSVs only when `_marker_hits_current(CFG)`. Before this, both steps reused anything that existed, so a rebuilt marker FASTA would never have reached an existing run.
- `src/annotation_explorer_polars/cache.py`: `TABLE_CACHE_SCHEMA_VERSION` 13 -> 14, so existing GUI caches are rebuilt with the guard and the new indexes on next open.

### Commands and timings
- `PYTHONPATH=src .venv/bin/python -m pipeline.MarkerDbBuild --dry-run` (1.4 s), then the real build (1.34 s wall, 242 MB RSS), then a second build to confirm idempotence (all `written=False`).
- Archived sources were byte-compared (`zcat | cmp`) with the backup `BACTolomew_backups/2026-10-09_p01-p06/DB/seq/markers/*.fasta`: all 6 identical.
- sha256 CARD_AMR.fasta 61bf7d53... -> ad0a336c...; AMRF_core-AMR.fasta 0dc58348... -> 38c5a754...; other marker FASTAs unchanged.

### Per database: entries before -> after, dropped categories
| DB | before | after | dropped |
|---|---|---|---|
| CARD_AMR | 4984 | 4773 | 211: VM 178, KM 18, OM 13, HM records whose ARO is overexpression-only in card.json 2 (baeS, baeR; added in Review fixes) |
| AMRF_core-AMR | 8263 | 7995 | 268 sequences: POINT only 220, POINT+POINT_DISRUPT 22, AMR-SUSCEPTIBLE 26 |
| AMRF_plus-AMR | 95 | 95 | 0 (4 accessions have AMR+POINT rows: kept) |
| AMRF_core-BMRG | 39 | 39 | 0 |
| AMRF_plus-BMRG | 218 | 218 | 0 |
| AMRF_plus-VIRULENCE | 968 | 968 | 0 |

The review's 1394 POINT / 21 POINT_DISRUPT / 26 AMR-SUSCEPTIBLE are catalog rows; one reference sequence carries many POINT rows, so 268 sequences go. Kept but not in the 2026-03-24.1 catalog (no metadata, listed in the report): 10 in AMRF_core-AMR (e.g. Bla1, DfrA51, Sat1, ACC-8, ACT-9, IMI-28, LEN-171/172, MprF), 1 in AMRF_plus-AMR (CrpP), 3 in AMRF_plus-VIRULENCE. A. baumannii models removed from CARD: Abau_gyrA_FLO, Abau_parC_FLO, LpxA, LpxC (VM), Abau_OprD_IPM, carO, Abau_AbuO (KM). The CARD_AMR FASTA is an older CARD release (4775 HM) than the 4.0.1 metadata (6052 HM); the build only filters, it does not change the release.

### .dmnd files
- Per-run, under `<run>/dmnd/queries/`; the reference run's CARD_AMR.dmnd and AMRF_core-AMR.dmnd (built 2026-05-29) are now older than their FASTA, so the next pipeline run on it rebuilds those two .dmnd files. Hit TSVs: a run that has a matching `annotations/bakta_vs_queries/query_fasta.json` re-searches only those two DBs; the reference run predates that stamp, so its next `collect_annotations` re-searches all 13 marker DBs once (corrected in Review fixes below; this line used to say only two TSVs). Not rerun here, as instructed.
- Built in scratch for the check (docker `bactolomew-tools:0.1.0`, diamond 2.2.3): `diamond dbinfo` CARD_AMR 4984 -> 4775 sequences.

### Direct DIAMOND check (scratch `p01p06/amr/dmnd_check`)
24 Bakta CDS from the reference run (3 genomes x gyrA, parC, oprD, pmrB, lpxC, blaOXA-23, blaOXA-66, adeS), blastx with the pipeline's options against the old and the new CARD_AMR and AMRF_core-AMR, best hit per DB and enrichment with the old and the new code:
| gene | old | new |
|---|---|---|
| gyrA | AMRF gyrA_G79C QUINOLONE; CARD VM Abau_gyrA_FLO fluoroquinolone | no AMR hit |
| parC | AMRF parC_S84F QUINOLONE; CARD VM Abau_parC_FLO | no AMR hit |
| oprD | CARD KM Abau_OprD_IPM carbapenem;... (82.2% id) | no AMR hit |
| lpxC | AMRF lpxC_P30L COLISTIN; CARD VM LpxC | no AMR hit |
| pmrB | AMRF pmrB_S17R COLISTIN | no AMR hit |
| adeS | AMRF adeS_F12S TETRACYCLINE; CARD HM adeS | CARD HM adeS only (glycylcycline;tetracycline) |
| blaOXA-23 | AMRF blaOXA-23 / blaOXA / BETA-LACTAM; CARD HM OXA-23 carbapenem;penicillin beta-lactam | identical |
| blaOXA-66 | AMRF blaOXA-66 BETA-LACTAM; CARD HM OXA-66 | identical |

### Data before -> after on the reference run (existing DIAMOND TSVs, old vs new code, no DIAMOND rerun; `p01p06/amr/before_after.py`, 30 s)
- Raw hits: CARD_AMR 135770 -> 134404, AMRF_core-AMR 208214 -> 205589 (3991 dropped by the guard).
- Best hit per CDS and DB (what the GUI's marker_hits holds): CARD_AMR 6766 -> 5586 (-1180), AMRF_core-AMR 5130 -> 2505 (-2625); AMRF_plus-AMR, core-BMRG, plus-BMRG unchanged.
- CDS with any AMR-family marker (AMR, core-AMR, plus-AMR; the "A" tag and the sample card's "Resistance genes" count): 8551 -> 5929. Per genome median 54 -> 36, mean 55.9 -> 38.8. Of the 2622 CDS that lose AMR, 1993 have no marker left, 306 become VIRULENCE, 170 MOBILE, 153 BMRG.
- Genes that lose the AMR call (genomes): gyrA, parC, lpxA, lpxC, rpsL, rpsJ, rpoB, ftsI, baeS, tolC 153; oprD 152; baeR 147; tetR(B) 133 (CARD OM tetR); opgE 126; pmrA, pmrB 98; tuf 76; ompR 8. Removed hit rows per model are in `p01p06/amr/removed_hits_reference_run.tsv`.
- CDS with resistance metadata: 8738 -> 6116.

### How the data is seen now
- Files: `DB/seq/markers/CARD_AMR.fasta` and `AMRF_core-AMR.fasta` no longer hold target models; `DB/seq/marker_sources/` holds the unfiltered sources and the removed-entry reports. A run's `annotations/bakta_vs_queries/*-vs-CARD_AMR.tsv` and `*-vs-AMRF_core-AMR.tsv` will lose those hits once the DIAMOND step reruns; until then the guard removes them when annotations are collected.
- GUI: on next open the cache (schema 14) is rebuilt; gyrA, parC, oprD, lpxA/C, pmrA/B, rpoB, ... are no longer tagged A, are absent from the Predictors CARD_AMR/AMRF_core-AMR filter, the sample card's resistance gene list and drug class cards, and the gene card's resistance fields. Acquired genes (blaOXA-23/66/72, blaADC, blaNDM-1, blaTEM-1, aac/aph, qacEdelta1, sul, tet(B) ...) are unchanged.

## P04: AMRF and CARD metadata indexes kept the first catalog row per accession

### What changed (code)
- `src/pipeline/ResistanceMetadata.py`: `_load_amrf_catalog()` (replaces `_load_amrf_accession_index`) indexes only reportable rows (every subtype except POINT, POINT_DISRUPT, AMR-SUSCEPTIBLE). Accessions with only target rows go into `target_accessions` (used by the guard). When an accession has several distinct reportable rows, the fields all rows agree on are kept (for example the gene family) and the allele and every disputed field are left empty; such accessions are listed in `ambiguous_accessions`. No first-row-wins anywhere.
- `_load_card_index()` returns `CardIndex`: rows whose ARO only has non-homolog models in card.json are excluded (390 AROs), `by_aro` and `by_protein` use the same no-guessing merge; lookup is by the hit's own ARO (from the header) first, the protein accession only as fallback.

### Data before -> after
- AMRF index keys 36949 -> 36085; 864 accessions are target-only (not indexed); 417 ambiguous keys, all nucleotide accessions (contigs carrying several genes), 0 protein accessions.
- CARD index by_aro 6440 -> 6088, by_protein 6324 -> 6086; 3 ambiguous protein keys (e.g. NP_214519.2 gyrB, allele left empty).
- Examples: WP_031379793.1 (adeS) old adeS_F12S / TETRACYCLINE, new: not indexed (target). WP_000116450.1 (gyrA) old gyrA_G79C, new: not indexed. WP_001046004.1 old and new blaOXA-23 / blaOXA / BETA-LACTAM / CARBAPENEM. WP_015585966.1 (fexA, AMR + 2 POINT rows): new payload from the AMR row (family fexA, no allele). CARD AAC75291.1 (E. coli gyrA, two variant AROs): old "E. coli gyrA conferring resistance to fluoroquinolones", new not indexed.
- Reference run, best hits (old vs new code on the existing TSVs): hits carrying a point-mutation allele 2625 -> 0. Top AMRF alleles before: pmrC_R125P 218, adeL_G27R 182, blaOXA-66 181, then 153 each for ftsI_A515V, baeS_D89V, baeR_S104N, parC_S84F, rpsJ_V57L, rpoB_S583L, lpxC_P30L, pmrA_E8D, pmrB_S17R, lpxA_G68D, gyrA_G79C, rpsL_K43R. After: blaOXA-66 181, blaADC-30 149, qacEdelta1 124, blaOXA-72 108, blaOXA-23 70, blaTEM-1 46, blaADC-73 26, blaNDM-1 24, aacA16 12. Of the hits present before and after, 0 changed allele, gene family, drug class or antibiotic (the bad payloads were all on target-only accessions, which P02 removes).
- Note: "pmrA_M12I" in `src/scripts/build_st492_res_panel.py` is a Bakta gene name (53 genomes carry it in the Bakta `Gene` column), not a value from these indexes; it is not changed here.

## P02 + P04: tests and downstream consumers
- New `tests/test_pipeline_amr_markers.py` (Qt-free, tiny fixtures; 20 tests at first, trimmed to 4 concept tests in Review fixes): AMRF index (target-only, acquired allele, mixed AMR+POINT accession, ambiguous accession), CARD index (variant ARO excluded, ambiguous protein, ARO first), guard (drops VM/KM/OM and AMRF target accessions, works without card.json, no point allele ever reported, guard before collapse keeps a lower-ranked HM hit), MarkerDbBuild (counts, archive, reports, idempotent rebuild with unchanged mtime, dry run, unprefixed CARD records by ARO), .dmnd freshness (legacy mtime rule, stamp content, stale TSV, rebuild with stamp, Runner reuse checks).
- `src/scripts` that read outputs that change once the cache is rebuilt and/or the DIAMOND step reruns (not modified): build_st492_res_panel.py (marker_hits; its Colistin row lists pmrA, pmrB, lpxA, lpxC), build_supplementary_annotation_workbook.py, build_supplementary_gene_calls.py, build_supplementary_raw_outputs.py (bakta_vs_queries TSVs), build_st492_synteny_contexts.py, build_st492_synteny_contexts_allctx.py, verify_st492_synteny_allctx.py (ARG coloring from AMR DB hits; their FASTA lookups `|qacEdelta1 ` in CARD_AMR and WP_004201167.1 in AMRF_core-AMR still resolve), list_blaoxa72.py (ConsensusCategory). build_st492_vf_panel.py (VFDB only) and build_st492_phage_tables.py (external XLSX) are not affected.

### Review fixes
- Freshness, one stale DB forced a re-search of every DB (confirmed): `collect_annotations` passed `force=not (reuse_marker_hits and query_matches)`, and `reuse_marker_hits` is False as soon as one TSV is older than its .dmnd, so a rebuilt CARD_AMR also re-searched Phaster, VFDB and the rest. Now `force=not query_matches` (`src/pipeline/Runner.py`): only a changed annotation query invalidates every TSV; otherwise `run_diamond_annotations` re-searches each DB whose TSV is older than its .dmnd and keeps the rest.
- A failed DIAMOND search left the old TSV, which the new `query_fasta.json` could then certify (confirmed): `run_diamond_annotations` (`src/pipeline/DiamondWorkflow.py`) now deletes the stale TSV when its search fails and logs it. The TSV is then missing, so the stamp is not written (`_all_files_exist` is false), `validate_marker_hit_outputs` reports the step's missing file, and the next run searches that DB again. The em-dash in that log message was replaced.
- baeS/baeR kept as HM by the DB build but always dropped by the guard (confirmed): `MarkerDbBuild.card_classifier` now checks every record's ARO against card.json, also HM-prefixed ones. Rebuilt `DB/seq/markers/CARD_AMR.fasta` 4775 -> 4773 (dropped `HM|gb|BAA15934.1|ARO:3000829|baeS` and `HM|gb|BAA15935.1|ARO:3000828|baeR`, E. coli K-12, overexpression-only in CARD 4.0.1; reason "CARD non-homolog model in card.json (by ARO)" in `DB/seq/marker_sources/reports/CARD_AMR.removed.tsv`); AMRF files unchanged (`written=False`). Previous file kept in scratch `amr_fix/CARD_AMR.before_baeSR.fasta`. The reference run has 0 hits to these two entries, so its data does not change; the build and the guard now agree, so the "rebuild the marker databases" warning can be cleared by a rebuild. The reference run's CARD_AMR.dmnd was already stale, so no extra rerun is needed.
- CARD FASTA vs metadata release (note, no code change): 12 HM records in CARD_AMR.fasta carry AROs absent from card.json / aro_index 4.0.1 (for example ACT-9 ARO:3001826, Lhon_ampC_BLA ARO:3004292, RSA-2 ARO:3004445, HERA-1 to HERA-8); hits to them get empty Resistance* fields. The reference run has 0 such hits. On the AMRF side LEN-171/LEN-172 (WP_430689812.1/WP_430689813.1) are kept but missing from the 2026-03-24.1 catalog and do appear among the reference run's AMRF_core-AMR subjects without metadata (listed under `kept_not_in_catalog`).
- Not changed: AMRF_plus-VIRULENCE has no `resistance_metadata` entry in run_config, so the runtime guard does not cover it; its DB filter removes 0 records, so nothing reaches it. Verified OK (guard and DB filter drop exactly the same references on the reference run's TSVs; P04 ambiguity only on nucleotide accessions): no change.
- Tests trimmed to the user's few-tests policy: `tests/test_pipeline_amr_markers.py` 20 -> 4 tests, one per concept: the DB build drops target models (including baeS-like HM records), the guard drops target hits from an old DB without losing a lower-ranked homolog hit, the metadata index never guesses (P04), and a rebuilt DB re-searches only its own TSV and a failed search leaves no stale TSV. The Runner `force` expression itself is not unit-tested (it needs a full `collect_annotations`); it is a one-line change. All `tests/test_pipeline_*.py` pass (14 tests).

## P05: genome sample names must be unique

What changed (code):
- `src/pipeline/DataInput.py`: new module functions `genome_sample_key` (file name minus `.gz` and one FASTA suffix: `.fasta`/`.fa`/`.fna`, case-insensitive), `find_genome_name_collisions` (groups of two or more different files sharing either `Path.stem`, the sample name the pipeline uses, or `genome_sample_key`; the same file given twice is not a collision) and `format_genome_name_collisions`. `DataInput.__post_init__` raises `ValueError` listing every colliding group with full paths, right after the genome list is normalized and before the output folder is created.
- `src/annotation_explorer_polars/startup_dialog.py` (`NewRunDialog`): "Add Genomes" now goes through `_add_genome_paths`, which adds non-colliding files and refuses the rest with one warning ("Duplicate genome names") naming each refused file and the already listed file it collides with; `_validate_and_accept` blocks OK if a collision is present anyway.
- `BACTOLOMEW_PIPELINE_GRAPH.txt`: the genome node now describes the real check.
- CLI: unchanged (`collect_fastas` still uses `rglob`); collisions found across subfolders now stop the run in DataInput.

Reference run check (run_78eaa9126850): `inputs.genomes` holds 153 files, all `/Data/ST492/final_fasta/*.fna`; 153 distinct stems and 153 distinct normalized names, 0 collisions. `collect_fastas(/Data/ST492/final_fasta)` (read only) also finds 153 files, 0 collisions. The reference run's data is unchanged by P05.

Demo (scratch `p01p06/inputs/p05/genomes`: `batchA/S1.fna`, `batchB/S1.fna`, `S2.fasta`, `S2.fna`, `S3.fna`):
- Before (code from the backup tarball): DataInput accepted 5 files as 3 sample names (S1 x2, S2 x2, S3), so two isolates would have been merged into others.
- After: `python src/main.py --genomes .../p05/genomes ...` exits 1 in 0.7 s with `ValueError: genomes: sample names (...) must be unique; ... S1: .../batchA/S1.fna, .../batchB/S1.fna; S2: .../S2.fasta, .../S2.fna`. No output folder is created and no tool starts.

How the data is seen now: no change for existing runs (outputs, cache and GUI identical). A new run with colliding names cannot start (CLI error, GUI picker refuses the file).

Tests: `tests/test_pipeline_genome_set.py` `GenomeNamesTest` (duplicate names refused by DataInput), `tests/test_startup_genome_names.py` (one Qt test: the picker refuses a colliding file; isolated QSettings, offscreen).

## P06: changed genome set in an existing output folder

What changed (code):
- `src/pipeline/RunValidation.py`: `recorded_genome_samples(output_dir)` reads `inputs.genomes` from `run_config.json`, then `run_state.json`, then `run_report.json`, and falls back to the Panaroo `gene_presence_absence.csv` sample columns for a folder with none of them; `genome_set_change_error(output_dir, genomes)` compares sample names (stems; moving the genome files to another folder is fine) and returns a message with the recorded and given counts and the added and removed samples (first 10 each).
- `src/pipeline/Runner.py`: new `GenomeSetChangedError`; the first thing `_run_workflow_impl` does is call `genome_set_change_error` and raise, before `run_config.json`, `expected_outputs.json`, `run.lock.json` or `run_state.json` are written and before any tool runs. Policy: refuse (no flag); use a new output folder for a new genome set (the GUI already does this for every new run; Resume passes the recorded genomes so it is unaffected).
- Defense in depth, data level: `validate_panaroo_outputs`, `validate_gubbins_outputs`, `validate_phylogeny_outputs` and `validate_chewbbaca_outputs` take an optional `genomes` list and then check that the Panaroo sample columns, the Gubbins final tree leaves, both tree files' leaves (the Snippy tree's `Reference` leaf is mapped to the stem in `snippy/reference.txt`) and the chewBBACA `results_alleles.tsv` FILE column equal the run's genome stems. The Runner passes `CFG.genomes`, so a stale cohort output is not reused (Panaroo reuse now also requires the sample set to match; otherwise Panaroo reruns, and fastGEAR follows), and `validate_step_outputs` / `validate_completed_run` (used by the GUI loader in deep mode) pass `inputs.genomes`. Overhead: under 3 ms per validator on the reference run; `validate_completed_run` on the reference run 10.84 s before, 10.98 s after (noise), 0 errors before and after.

Demo (scratch `p01p06/inputs/p06`): a finished v1 run was simulated by writing `out/run_config.json` with the Runner's own `_build_run_config` for `genomes_v1` (S1, S2; sample_set_id `sample_bd070cbe18f9`) and a Panaroo header with S1, S2.
- Before (backup code, stopped at the first tool step by a patched path helper): the run went ahead with v2 (S1, S2, S3), kept `sample_set_id` `sample_bd070cbe18f9` and the old `run_id` (a fresh id for v2 would be `sample_45b7b3dc0135`), rewrote `run_config.json` with 3 genomes and would have reused the 2-sample Panaroo table.
- After: `python src/main.py --genomes .../genomes_v2 --output .../p06/out ...` exits 1 in 0.5 s with `GenomeSetChangedError: ... (2 samples recorded in run_config.json, 3 given). Added: S3. ...`; `run_config.json` sha1 unchanged (04c4730b...), no other file written. With `genomes_v3` (S1 only): `... 1 given). Removed: S2.` With `run_config.json` deleted, the Panaroo header fallback still refuses.

Reference run: the guard passes for its own 153 genomes; adding one genome is refused. The four cohort validators report 0 errors with the genome list (Panaroo 153 columns, Gubbins tree 153 leaves, core tree 153 leaves, Snippy tree 152 + `Reference` = GCA_028622775.1, chewBBACA 153 rows). No data changed.

How the data is seen now: unchanged for existing runs (no output, cache or GUI change). A CLI run into a folder built from another genome set stops at once with the message above instead of producing trees, cgMLST and marker hits that leave out or misattribute isolates.

Tests: `tests/test_pipeline_genome_set.py` `GenomeSetTest` (the Runner refuses a changed genome set and leaves the folder untouched; a stale cohort output is caught while a genome the tool left out is accepted; marker hits are tied to the exact query FASTA).

Manuscript scripts in `src/scripts`: none consume outputs changed by P05 or P06 (no output changed).

Test runs (P05/P06): `tests/test_pipeline_genome_set.py` 19 tests OK (0.04 s); `tests/test_startup_genome_names.py` 3 tests OK. Full suite (`unittest discover`, offscreen, 1229 s): 707 tests, 706 OK; the one error was `test_runner_refuses_before_writing_anything` failing to import `pipeline.Runner` (`dmnd_is_current` missing from DiamondWorkflow) while another item's edit to DiamondWorkflow/Runner was in flight; rerun right after: OK, and all `test_pipeline_*` modules plus the startup test (60 tests) pass together.
- P02/P04 test run: full suite (`unittest discover -s tests`, 727 tests, 1335 s, machine shared with other agents) passed except the timing test `test_pan_panels.HeavyPanelHoverBudgetTests.test_a_hover_reaches_four_real_heavy_panels_within_30_ms`, which passed on rerun with `BACTOLOMEW_TEST_TIMING_FACTOR=2` as CLAUDE.md prescribes for a loaded machine (2 tests OK).

### Review fixes

Review of P05/P06 (six findings), checked against the code:
- Medium, Gubbins/chewBBACA leaving a genome out made a run fail: real; fixed. The exact sample-set check only applies when nothing shows the tool was given this genome set. `Runner` writes `gubbins/genome_set.json` and `chewbbaca/genome_set.json` after the tool returns, logs and records `samples_left_out` in the step summary, and the validators then accept missing samples (a stamp for another set, or any sample outside the set, is still stale). New today: folders written before the stamps existed (every finished run so far) fall back to the tool's own input record, `gubbins/input/core_gene_alignment.aln` headers equal to the genome set, or chewBBACA `out/logging_info.txt` "Number of inputs" equal to the genome count (`RunValidation._stamp_errors(tool_input_matches=...)`). Without this an older run where Gubbins filtered a gappy genome could not be opened in the GUI any more: the cache bump to schema 14 forces a rebuild, which runs the deep validation, and any error refuses the run. Checked on real data (scratch `p06_review/legacy_check.py`, copies of the reference run's outputs, 1.2 s): reference run 0 errors for both validators (its staged Gubbins input holds the 153 samples, chewBBACA logged 153 inputs); a copy of the Gubbins tree without GCA_023499075.1 and the 153-sample input: 0 errors; the same tree with a 152-sample input (another set): "missing 1"; chewBBACA without one row and 153 logged inputs: 0 errors, with 152 logged inputs: "missing 1". The phylogeny trees keep the exact check: the core tree is built from the full Panaroo alignment and a failed Snippy sample stops the step.
- Low, marker-hit TSVs reused in folders already rewritten by the old code: real; already fixed during P02/P06 by `annotations/bakta_vs_queries/query_fasta.json` (sha256 of the annotation FASTA the hits were searched from; a folder without it counts as stale, so `collect_annotations` re-searches). Consequence for the reference run: all 13 marker TSVs are recomputed once, not 2.
- Low, external Snippy reference leaf reported as unexpected: real; already fixed (`tree_sample_names` drops `Reference` when the reference is not one of the run's genomes).
- Low, `tests/test_pipeline_genome_set.py` failing to import alone: real; already fixed (the module adds `src` to `sys.path`; checked with `cd tests && python -m unittest test_pipeline_genome_set`).
- Low, GUI picker checking unresolved paths, O(n^2): real; already fixed (`DataInput.genome_name_keys` names a file by its resolved target in both places; the picker builds its name map once per add).
- Low, the guard compares sample names only: not fixed, out of scope. A genome file replaced by another assembly under the same name is not detected, and when Panaroo reruns for the same set the downstream steps follow their own validators. Same behavior as before P06; use a new output folder after replacing assemblies.

Tests trimmed to the concepts (user request): `tests/test_pipeline_genome_set.py` 4 tests (was 31), `tests/test_startup_genome_names.py` 1 test (was 4). All `tests/test_pipeline_*.py` plus the startup test: 31 tests OK (3.7 s).

## P01 and P03: cgMLST schemas adapted with chewBBACA PrepExternalSchema; self-score guard

### What changed (code)
- `src/pipeline/ChewbbacaSchemaBuild.py`: the hand-written schema writer (which only renamed headers, took the first allele of each locus as its only representative and wrote `.schema_config` / `.genes_list` itself, so the old "must be adapted with PrepExternalSchema" check always passed) is replaced. Loci are staged as one FASTA per locus (headers `<locus>_<allele>`) and adapted with `chewBBACA.py PrepExternalSchema` (chewBBACA 3.5.4 in `bactolomew-tools:0.1.0`, `--bsr 0.6`, `--t 11`, `--t 4` for Mycoplasma). PrepExternalSchema writes every allele in coding orientation, drops alleles that are not a complete CDS and picks the representatives in `short/`. `adapted_schema_issues()` then requires: `.genes_list` and `.schema_config` present, one `short/` file per locus, `.genes_list` = locus FASTAs, no representative or allele that is a CDS only on the reverse strand, no invalid representative. A schema is installed (the old folder is renamed aside, the new one renamed in, and the old one restored if that second rename fails; two renames, not one atomic step) only after that; a failed build leaves the existing schema untouched. A prebuilt chewBBACA folder is copied only when its representatives are forward CDS (runtime caches `pre_computed`, `loci_modes`, `short/self_scores` are not copied), otherwise its loci are re-adapted.
- Reports next to each schema: `DB/cgmlst/<schema>/schema_build_report/` holds `build_info.json` (builder, chewBBACA version, BSR, genetic code, source, expected/source/kept loci and alleles, representatives, rejected counts), `invalid_alleles.txt` and `invalid_loci.txt` (PrepExternalSchema's per-allele and per-locus reasons), `summary_stats.tsv` and `prep_external_schema.log`. chewBBACA ignores the folder (it reads only `*.fasta` at the schema root; the smoke AlleleCalls ran with it present).
- `src/pipeline/ChewBBACAWorkflow.py`:
  - `validate_chewbbaca_settings` refuses a schema whose representatives are CDS only on the reverse strand (`schema_orientation_issue`), so an unadapted schema can no longer produce a run with half the loci LNF.
  - Before AlleleCall, `check_self_scores_before_call`: when `short/self_scores` exists but some representative has no entry, the cache is deleted (warning logged) so AlleleCall recomputes it for every representative.
  - After AlleleCall, `check_self_scores_after_call` fails the step (RuntimeError, cache deleted) when `short/self_scores` is missing, when any representative has no self-score, or when stdout has any "Could not get the self-score" line (counted in the message).
  - Rule used: every representative must have a self-score (missing ones are what makes chewBBACA drop BLAST matches). Strict count equality, as the review suggested, is enforced only for a freshly built schema (build script). At run time it would fail valid schemas: chewBBACA 3.5.4 `allele_call.py` keeps self-scores of representative candidates that were not added ("This includes self-score for candidates that are not added", and entries of candidates converted to NIPH are not deleted), so after an inferred-mode AlleleCall the cache can legitimately hold more entries than representatives.
- `src/annotation_explorer_polars/cgmlst_graph_data.py`: `INF-<id>` (a novel allele chewBBACA inferred) is now read as a called allele; before it was parsed as missing. LNF, NIPH, PLOT3/5, ASM, ALM stay missing. The old run had 0 INF calls (P03), so this had no visible effect until now.
- `src/scripts/build_chewbbaca_schemes_from_ncs.py`: now a thin driver over `build_schema` for the cgMLST.org NCS downloads in `/Data/cgmlst/download_tmp/<slug>/{targets.csv,alleles.zip}`: build into `--staging-dir`, validate, smoke AlleleCall with `--no-inferred` (the schema gains no alleles; this computes `short/self_scores`) on `--smoke-genome` files, require 0 "Could not get the self-score" lines and self-scores = representatives, then with `--swap` move the old schema to `--old-dir` and the new one into place. Added in this continuation:
  - Idempotent: a schema whose `build_info.json` says PrepExternalSchema built it from the same `alleles.zip` (not modified since `built_at`) and that still passes `adapted_schema_issues` is skipped (`--force` rebuilds). Checked read only on the installed A. baumannii schema: skipped for its own zip, not for another zip (4.7 min, almost all of it the orientation scan of 1.58 M alleles, with the leftover build using the CPUs).
  - A backup already in `--old-dir` (the original schema) is never overwritten; a later one gets a timestamp suffix.
  - At the end the driver checks `DB/cgmlst/scheme_species_catalog.json` against the schema folders. Now: 34 entries, 34 folders, 0 differences. The rebuild keeps every folder name, so no catalog entry changes.

### Commands and timings
- A. baumannii (done by the previous session, report `p01p06/cgmlst/report_abaumannii.json`, log `ab_build.log`):
  `PYTHONPATH=src .venv/bin/python src/scripts/build_chewbbaca_schemes_from_ncs.py --output-dir DB/cgmlst --staging-dir DB/cgmlst/.p01_staging --old-dir <scratch>/p01p06/cgmlst/old --smoke-genome <scratch>/p01p06/smoke/SAMN04040115.fna --smoke-genome <scratch>/p01p06/smoke/20038ENE2.fna --smoke-dir <scratch>/p01p06/cgmlst/smoke_runs --only-existing --match abaumannii --cpu 64 --swap --report .../report_abaumannii.json`: 631 s wall (build 336 s, smoke AlleleCall 96 s), 61 MB RSS for the driver.
- Other 33 schemas: the same command with `--skip acinetobacter_baumannii --cpu 16 --jobs 6`, started 2026-10-09 10:58 (pid 3042661, log `p01p06/cgmlst/rest_build.log`, report `report_rest.json` written at the end). Still running when this was written (11:32): 28 of 33 swapped in (all `[ok]`, 64 s to 1007 s each), 5 still building (mycobacterium_tuberculosis_bovis_africanum_canettii, pseudomonas_aeruginosa, serratia_marcescens, staphylococcus_capitis, streptococcus_pyogenes). providencia_stuartii and staphylococcus_argenteus finished after the table below was taken. Their counts must be read from `report_rest.json` / each `schema_build_report/build_info.json` when it finishes.
- Original schemas: `<scratch>/p01p06/cgmlst/old/<schema>` (moved there by `--swap`); before stats of all 34: `p01p06/cgmlst/before_stats.json` (`measure.py`).

### Per schema: before -> after (schemas swapped in so far)
Before: every schema had 1 representative per locus and about half of them (column "reverse reps") were CDS only on the reverse strand; only A. baumannii had a `self_scores` file, with 1 entry. After: 0 reverse-strand alleles or representatives, self-scores = representatives (computed by the smoke AlleleCall), 0 "Could not get the self-score" lines.

| schema | loci before | alleles before | reverse reps before | loci after | alleles after | reps after | self-scores | rejected alleles | rejected loci | code |
|---|---|---|---|---|---|---|---|---|---|---|
| acinetobacter_baumannii | 2390 | 1611136 | 1167 | 2390 | 1580762 | 2585 | 2585 | 30374 | 0 | 11 |
| bacillus_anthracis | 3803 | 285856 | 2006 | 3803 | 285274 | 4113 | 4113 | 582 | 0 | 11 |
| bordetella_pertussis | 2983 | 42649 | 1467 | 2983 | 42649 | 3097 | 3097 | 0 | 0 | 11 |
| brucella_melitensis | 2704 | 119768 | 1321 | 2704 | 113826 | 2829 | 2829 | 5942 | 0 | 11 |
| brucella_spp | 1764 | 71630 | 851 | 1764 | 71606 | 1830 | 1830 | 24 | 0 | 11 |
| burkholderia_mallei_rki | 3328 | 242749 | 1834 | 3328 | 242747 | 3478 | 3478 | 2 | 0 | 11 |
| burkholderia_pseudomallei | 4221 | 730422 | 2207 | 4221 | 729220 | 4451 | 4451 | 1202 | 0 | 11 |
| campylobacter_jejuni_coli | 637 | 830906 | 319 | 637 | 808548 | 750 | 750 | 22358 | 0 | 11 |
| citrobacter_freundii | 3250 | 628242 | 1602 | 3250 | 628242 | 3324 | 3324 | 0 | 0 | 11 |
| citrobacter_freundii_portucalensis_braakii_europaeus | 2307 | 339754 | 1128 | 2307 | 339754 | 2359 | 2359 | 0 | 0 | 11 |
| clostridioides_difficile | 2147 | 472389 | 1077 | 2147 | 472246 | 2229 | 2229 | 143 | 0 | 11 |
| clostridium_perfringens | 1431 | 198784 | 926 | 1431 | 198784 | 1449 | 1449 | 0 | 0 | 11 |
| corynebacterium_diphtheriae | 1553 | 139376 | 764 | 1553 | 139376 | 1564 | 1564 | 0 | 0 | 11 |
| corynebacterium_pseudotuberculosis | 1577 | 14612 | 852 | 1577 | 14612 | 1579 | 1579 | 0 | 0 | 11 |
| cronobacter_sakazakii_malonaticus | 2831 | 344630 | 1359 | 2831 | 344160 | 2861 | 2861 | 470 | 0 | 11 |
| enterobacter_hormaechei | 2178 | 421393 | 1164 | 2178 | 421393 | 2230 | 2230 | 0 | 0 | 11 |
| enterococcus_faecalis | 1972 | 530236 | 1013 | 1972 | 521028 | 2007 | 2007 | 9208 | 0 | 11 |
| enterococcus_faecium | 1423 | 666591 | 871 | 1423 | 563071 | 1625 | 1625 | 103520 | 0 | 11 |
| francisella_tularensis | 1147 | 20120 | 614 | 1147 | 19291 | 1153 | 1153 | 829 | 0 | 11 |
| klebsiella_oxytoca_grimontii_michiganensis_pasteurii | 2536 | 519896 | 1320 | 2536 | 519896 | 2606 | 2606 | 0 | 0 | 11 |
| legionella_pneumophila | 1521 | 143380 | 754 | 1521 | 139152 | 1532 | 1532 | 4228 | 0 | 11 |
| listeria_monocytogenes | 1701 | 1121939 | 908 | 1701 | 1102346 | 1791 | 1791 | 19593 | 0 | 11 |
| morganella_morganii | 2462 | 92729 | 1190 | 2462 | 92729 | 2467 | 2467 | 0 | 0 | 11 |
| mycobacteroides_abscessus | 2904 | 274983 | 1514 | 2904 | 272506 | 2936 | 2936 | 2477 | 0 | 11 |
| mycoplasma_gallisepticum | 425 | 10223 | 53 | 425 | 9859 | 425 | 425 | 364 | 0 | 4 |
| paenibacillus_larvae | 2419 | 26255 | 1164 | 2418 | 26252 | 2420 | 2420 | 3 | 1 (ERICI_RS18885) | 11 |
| proteus_mirabilis | 2574 | 128477 | 1294 | 2574 | 128477 | 2582 | 2582 | 0 | 0 | 11 |

Rejected alleles are Ridom alleles with an internal stop codon or no start codon in either orientation (reasons per allele in `invalid_alleles.txt`); chewBBACA cannot call them. The largest share is E. faecium (103520 of 666591, 15.5%); it keeps all 1423 loci. Mycoplasma gallisepticum was built with genetic code 4 (TGA = Trp).

### A. baumannii smoke AlleleCall (2 ST492 genomes, `--no-inferred`, old schema vs new)
- 20038ENE2 (in the reference run): old run EXC 1155, INF 0, LNF 1235; new EXC 2218, INF 104, NIPH 21, PLOT 3, ALM 3, ASM 4, LNF 37. Called loci as the GUI's cgMLST loader counts them (`_load_chewbbaca_profiles_cached`): 1155 -> 2322 of 2390.
- SAMN04040115 (also in the reference run): old run EXC 1166, INF 0, LNF 1224; new EXC 2302, INF 47, LNF 13 (2349 of 2390 called).
- Per locus, old reference-run calls vs new smoke calls: where both schemas called a locus the allele ID is the same (1141 of 1141 for 20038ENE2, 1154 of 1154 for SAMN04040115). 14 (20038ENE2: NIPH 14) and 12 (SAMN04040115: NIPHEM 8, NIPH 4) loci called with the old schema are now NIPH/NIPHEM: real protein-level paralogs (e.g. ACICU_RS07090/RS07095 have two extra copies on 20038ENE2_contig00016) that the old schema missed because half its representatives were reverse-only. Correct behavior, not a regression.
- Reference run, old schema (unchanged, not rerun): EXC 174318, INF 0, LNF 187555 (1235 to 1236 LNF per genome in the rows checked), NIPHEM 3797, 32461 "Could not get the self-score" lines.

### Tests
- `tests/test_pipeline_cgmlst_schema.py` trimmed to the concepts (17 -> 6 tests): a Ridom zip is adapted with reverse-strand alleles oriented, non-CDS alleles and empty loci reported in `schema_build_report/`; an unoriented prebuilt schema is refused by the orientation check and re-adapted; `run_chewbbaca` drops a degenerate cache before AlleleCall and accepts fresh scores; the step fails when a representative has no self-score and when chewBBACA logs "Could not get the self-score"; the real PrepExternalSchema (tools image, skipped without docker) orients reverse alleles. Run: 6 tests OK in 6.2 s.

### How the data is seen now
- Files: every rebuilt `DB/cgmlst/<schema>` has oriented alleles, several representatives per locus where alleles diverge, `short/self_scores` with one entry per representative, and `schema_build_report/`. The A. baumannii schema is final. A run's `chewbbaca/out/*` changes only when chewBBACA reruns.
- Reference run: NOT rerun here. Its `chewbbaca/out/` still holds the old calls (about 1235 LNF per genome), and the cgMLST step now reruns by itself: its calls have no `chewbbaca/schema_stamp.json` (see Review fixes below), so the pipeline does not reuse them. In inferred mode (the pipeline default) AlleleCall adds the novel alleles and their self-scores to `DB/cgmlst/acinetobacter_baumannii`, so its allele and representative counts will grow past the numbers above after that run.
- GUI (after the rerun and the cache rebuild): the cgMLST force graph uses about 2300 called loci per genome instead of about 1155, novel alleles (INF) count as called, distances and the called-loci % change, and duplicate (identical) profiles become rarer.
- GUI, per locus: a few loci per genome (14 and 12 in the two smoke genomes) move from called to NIPH/NIPHEM and count as missing in the force graph; every locus called by both schemas keeps its allele ID.
- Allele numbers: PrepExternalSchema keeps the Ridom allele IDs, and chewBBACA numbers a novel allele as the highest kept ID plus 1. In 2 A. baumannii loci the highest Ridom IDs were rejected (18 IDs in all; e.g. ACICU_RS01015 keeps 8 of 86 alleles, highest kept ID 69, and the smoke run calls INF-70 and INF-71, which are Ridom numbers of other sequences). No effect inside BACTolomew (distances are computed within one schema, IDs in it stay unique); only matters when comparing allele numbers with external cgmlst.org profiles.

### Review fixes
- Medium, chewBBACA reuse ignored the schema (old calls from the unoriented schema passed as reusable): fixed. `run_chewbbaca` now writes `chewbbaca/schema_stamp.json` (schema folder plus the sha256 of its `schema_build_report/build_info.json`; only the build record is hashed because inferred-mode AlleleCalls grow the schema). `Runner._validate_chewbbaca_reuse_outputs` reruns chewBBACA when the stamp is missing or the schema's build record changed (`ChewBBACAWorkflow.chewbbaca_schema_stamp_issue`). Checked read only on the reference run: "chewBBACA calls have no schema stamp", so its chewBBACA step reruns without moving `chewbbaca/` aside. The GUI's deep validation (`validate_chewbbaca_outputs`) does not check the stamp on purpose: existing runs keep opening, with their old calls, until they are rerun.
- Low, schema install not atomic and no rollback: fixed. `ChewbbacaSchemaBuild._replace_output_dir` (also the GUI "Add schema" path) no longer deletes the old schema before the new one is in place: it renames the old one aside, renames the new one in, restores the old one if that fails, and deletes the old one only after success (checked in scratch with a simulated failed rename: the old schema stays in place). The driver's `--swap` restores the backup into `DB/cgmlst/<schema>` when `staging.replace(target)` fails. The "atomic rename" wording above is corrected.
- Low, SAMN04040115 "not run before": corrected above (old EXC 1166, LNF 1224 from the reference run's `results_statistics.tsv`).
- Low, loci now NIPH: confirmed (14 and 12, 0 allele ID differences on shared calls) and noted above; no code change.
- Low, Ridom allele numbers reused by novel alleles in 2 loci: noted above; no code change.
- Leftover schema build (pid 3042661): finished, `done: 33 ok, 0 failed` in `rest_build.log`. The table above still lists only the schemas swapped in when it was written; the rest are in `report_rest.json`.
- Test: the existing run_chewbbaca concept test also checks that the calls are stamped and that a rebuilt schema makes them stale (no new test). All `tests/test_pipeline_*.py`: 30 tests OK.

## Reference run rerun (run_78eaa9126850, 2026-10-09 12:00 to 12:06)

Only the steps the fixes affect were rerun, through the pipeline's own `Runner.run_workflow`, from a scratch driver (`<scratch>/rerun/rerun_driver.py`) that builds `DataInput` the way the GUI's Resume does (`startup_dialog._data_input_from_run_config`: same 153 genomes, Snippy reference, run_id, Bakta DB and tool settings), with 96 cores. The driver checks the step plan with the Runner's own reuse predicates first and refuses to run when any unrelated step would rerun.

### Preparation
- COG (pre-existing blocker, not caused by P01 to P06): `run_config.json` still listed `DB/seq/markers/COGorg24_COG.fasta`, which was moved to `DB/_removed` on 2026-08-11 and is not in the marker catalog, so neither the pipeline (`DataInput`) nor the GUI cache build could use the run. The driver drops COG from the marker list; the run's COG hit TSV (86 MB) and COG .dmnd (3.2 GB) were moved (not deleted) to `BACTolomew_backups/2026-10-09_rerun_moved_out/`. The run now has 13 marker DBs, and COG is no longer a Predictor in the GUI for this run (cache marker_hits 535240 -> 53758 rows, 477677 of the old rows were COG).
- Query stamp: the run predates `annotations/bakta_vs_queries/query_fasta.json`, so `collect_annotations` would have re-searched all 13 marker DBs. The driver regenerated the annotation query FASTA exactly as `collect_annotations` does (574416 CDS, sha256 6c243444...), checked every existing hit TSV against it (every qseqid present and qlen equal to the sequence length, 0 mismatches in 607200 rows), and only then wrote the stamp. The Runner then re-searched only the TSVs older than their rebuilt .dmnd.
- Plan check (all from the Runner's predicates): genome set unchanged; assembly_qc, Bakta, MLST, MOB-suite, Panaroo, Gubbins, phylogeny, fastGEAR and reference loci reused; chewBBACA reruns (no `schema_stamp.json`); stale .dmnd exactly AMRF_core-AMR and CARD_AMR. No validator or the P06 guard blocked the unchanged genome set, so no code fix was needed.
- A. baumannii schema before the run: 2390 loci, 1580762 alleles, 2585 representatives, 2585 self-scores.

### What ran (wall times, run_workflow 321.8 s in all)
| step | result | time |
|---|---|---|
| assembly_qc, bakta, mlst, mobsuite, panaroo, gubbins, phylogeny, fastgear_clades, diamond_reference_loci | reused_existing_output | reuse checks 12 s |
| chewbbaca | AlleleCall, 153 genomes, 32 threads, inferred mode, rebuilt A. baumannii schema | 116 s |
| diamond_markers | rebuilt AMRF_core-AMR.dmnd and CARD_AMR.dmnd only (with `.source.json` stamps); the other 11 kept (mtimes 2026-05-29) | 4 s |
| collect_annotations | re-searched AMRF_core-AMR (11 s) and CARD_AMR (11 s); the other 11 hit TSVs skipped (mtimes still 2026-06-03, byte-identical to the backup) | 184 s |

`run_state.json` and `run_report.json`: status success, outputs_validated true, every step success, chewBBACA samples_left_out []. `validate_completed_run`: 0 errors. New files: `chewbbaca/schema_stamp.json`, `chewbbaca/genome_set.json`, `dmnd/queries/{AMRF_core-AMR,CARD_AMR}.dmnd.source.json`, `annotations/bakta_vs_queries/query_fasta.json`. `run_config.json` and `expected_outputs.json` were rewritten by the Runner (13 markers, 96 cores).

### P01 / P03 on the run (old calls from the backup vs new)
| | before | after |
|---|---|---|
| EXC (sum over 153 genomes) | 174318 | 348648 |
| INF | 0 | 858 |
| LNF | 187555 (median 1220 per genome) | 3481 (median 13) |
| NIPH / NIPHEM | 0 / 3797 | 337 / 10544 |
| PLOT3 / PLOT5 / LOTSC / ALM / ASM | 0 | 129 / 231 / 1 / 439 / 1002 |
| called loci per genome (GUI rule: number or INF) | median 1155, min 871 | median 2322, min 1916 |
| loci called in at least 95% of genomes | 925 | 1896 |
| distinct allele profiles | 143 | 153 |
| "Could not get the self-score" lines | 32461 | 0 |

The schema grew as expected in inferred mode: 1581620 alleles (+858 novel), 2589 representatives, 2589 self-scores, 0 reverse-only; the P03 guard passed.

### P02 / P04 on the run
- CARD_AMR hits: 135770 rows (HM 134404, VM 891, KM 305, OM 170) -> 134404 rows, all HM; CDS with a hit 6766 -> 5586. AMRF_core-AMR: 208214 -> 205589 rows, CDS 5130 -> 2505. No new hit appeared; every hit kept has the same bitscore (e-values are slightly smaller with the smaller DBs, none crossed the cutoff). The runtime guard dropped nothing (the DBs are clean now).
- GUI cache `cache/run_78eaa9126850` rebuilt like the app loader does (`validate_output_dir(deep=True)` 0 errors in 11 s, `load_output_data` 240 s, 267 s wall, 6.6 GB RSS): schema 14, annotations 574416 rows (unchanged), CDS with the A tag 8551 -> 5929, CDS with resistance drug classes 8738 -> 6116, A-tagged gyrA and parC 153 -> 0, blaOXA-23 66 and blaOXA-66 180 unchanged. These match the numbers measured earlier from the old TSVs with the guard.
- Backups of everything replaced: `BACTolomew_backups/2026-10-09_p01-p06/{run,cache_run_78eaa9126850}`.

## Data impact (final, rechecked)

Every number below was recomputed from the files (backup `BACTolomew_backups/2026-10-09_p01-p06/{run,cache_run_78eaa9126850,DB}` vs the live run, cache and `DB/`), with independent scripts in `<scratch>/impact_check/` (`cg.py`, `mst.py`, `schemas.py`, `amr.py`, `idx.py`, `atb.py`). Corrections to the first draft of this report are marked "Corrected".

### P01: cgMLST schemas rebuilt with PrepExternalSchema
- Code: `ChewbbacaSchemaBuild.py` builds each schema with chewBBACA PrepExternalSchema (forward-oriented CDS alleles, several representatives per locus, invalid alleles listed), installs by rename-aside with rollback and writes `schema_build_report/build_info.json` (Corrected: the file is inside `schema_build_report/`, not at the schema root). `ChewBBACAWorkflow` stamps the calls (`chewbbaca/schema_stamp.json`) so calls from another schema build are not reused. `src/scripts/build_chewbbaca_schemes_from_ncs.py` drives the one-time rebuild.
- DB/cgmlst (checked by re-reading all 34 schemas, each with its own translation table): 0 forward-invalid and 0 reverse-only alleles everywhere; `.genes_list` equals the locus files. Loci unchanged except paenibacillus_larvae 2419 -> 2418 (ERICI_RS18885). A. baumannii before: 1611136 alleles, 818561 not a forward CDS (788187 of them reverse-only, 30374 invalid either way). Largest drops of invalid alleles: E. faecium 103520, M. tuberculosis 40875, A. baumannii 30374. mycoplasma_gallisepticum uses table 4: 364 alleles dropped (the 8417 "invalid" in the before column of `cgmlst_schemas_after.tsv` were judged with table 11).
- A. baumannii schema: 2390 loci; alleles 1611136 -> 1580762 (build) -> 1581620 (after the run, +858 INF); representatives 2390 -> 2585 -> 2589. P. aeruginosa representatives 3867 -> 7155.
- Reference run `chewbbaca/out/`: Corrected: `results_alleles.tsv`, `results_statistics.tsv`, `results_contigsInfo.tsv`, `loci_summary_stats.tsv` and `logging_info.txt` were rewritten, not new (they existed before); `cds_coordinates.tsv`, `invalid_cds.txt` and `paralogous_*.tsv` are byte-identical. New files: `chewbbaca/schema_stamp.json`, `chewbbaca/genome_set.json`. Same 153 FILE rows, same order.

| reference run, 153 genomes | before | after |
|---|---|---|
| LNF per genome min/median/max (total) | 1207/1220/1519 (187555) | 2/13/394 (3481) |
| called loci per genome (number or INF) min/median/max | 871/1155/1183 | 1916/2322/2359 |
| loci called in >=95% / >=99% / 100% of genomes | 925 / 680 / 439 | 1896 / 1502 / 1052 |
| loci never called | 1199 | 5 |
| EXC / INF / NIPH / NIPHEM | 174318 / 0 / 0 / 3797 | 348648 / 858 / 337 / 10544 |
| ASM / ALM / PLOT3 / PLOT5 / LOTSC / PAMA | all 0 | 1002 / 439 / 129 / 231 / 1 / 0 |
| comparable loci per pair median (min) | 1135 (698) | 2283 (1566) |
| GUI distance, all 11628 pairs: median, IQR, max | 149, 59-196, 898 | 294, 118-380, 1801 |
| GUI distance median within ST492 (93) / ST2 (60) / ST492 vs ST2 (min, median) | 58 / 82 / 103, 181 | 112 / 175 / 208, 339.5 |
| strict core distance (439 -> 1052 loci) median all / ST492 / ST2; ST492 vs ST2 min | 73 / 11 / 27.5; 43 | 159 / 37 / 79; 104 |
| zero-distance pairs (GUI / strict core) | 83 / 118 | 2 / 15 |
| distinct profiles (missing kept as missing / strict core) | 143 / 109 | 152 / 145 |
| GUI MST: total, zero-length edges, median edge, max edge | 3187, 31, 5, 649 | 8179, 2, 17.5, 1454 |
| mean pairwise distance in the force graph info line | 182.67 | 365.81 |
| single linkage clusters at 5/10/20 alleles, GUI distance (largest) | 74/49/34 (37/45/45) | 122/101/69 (10/10/37) |
| single linkage clusters at 5/10/20, strict core | 33/26/18 | 88/59/36 |

- Distance per comparable locus is about 13% in both (149/1135, 294/2283); absolute distances roughly doubled. Spearman of all pair distances before vs after: 0.985 (GUI), 0.984 (strict core); against Gubbins recombination-filtered SNP distances: 0.940 -> 0.954.
- Selected ST492 pairs (GUI distance): GCA_029679835.1/GCA_031156395.1 0 -> 0; GCA_958318645.1/GCA_958321385.1 0 -> 0; GCA_016474395.1/GCA_022017345.1 0 -> 1; GCA_958321045.1/SAMN09202764 16 -> 32; GCA_031875465.1/GCA_048215925.1 62 -> 114; GCA_031156315.1/GCA_031156365.1 6 -> 92; GCA_023499075.1/GCA_964331735.1 893 -> 1801.
- Corrected: the distinct-profile count after is 152 (missing kept as missing), not 153 as in the rerun table above.
- Corrected (MST topology): the global ranking of distances is preserved (Spearman 0.985), but the MST itself changed a lot: only 70 of its 152 edges are shared before and after. What stays: fastGEAR lineage 1 (6 genomes) is still one subtree attached to the rest by a single edge (649 before, 1454 after, to GCA_023499055.1 in both). Corrected: 2 zero-length MST edges remain (4 genomes, the two 0 -> 0 pairs above), not "2 nodes".
- GUI: the cgMLST force graph reads `results_alleles.tsv` directly (no cache), so it shows the new calls now. The Sample tree, the Annotations table and the cache are not affected by P01. Any fixed allele threshold used before (for example <=10 alleles) must be re-derived before interpreting clusters.
- Affected outside the app: `src/scripts/build_st492_cgmlst_mst.py` and `build_st_cgmlst_matrix.py` (read this run's `results_alleles.tsv`) and their R plotters `plot_st492_cgmlst_mst.R`, `plot_st_cgmlst_matrix.R`; `build_supplementary_raw_outputs.py`; `tune_cgmlst_force.py` (dev). Corrected: `build_supplementary_annotation_workbook.py` is not affected by P01 (it only carries the metadata column `cgmlst_complex_type`).

### P03: self-score guard
- Before: only A. baumannii had `short/self_scores`, with 1 entry that matched none of its 2390 representatives; the other 33 schemas had none. After: all 34 schemas have exactly one self-score per representative (0 missing, 0 extra); A. baumannii 2589 = 2589 after the run.
- "Could not get the self-score" lines in the run's chewBBACA stdout: 32461 -> 0.
- No separate GUI change: the missing self-scores were a main cause of the false LNF calls counted under P01. A run now fails at the chewBBACA step instead of producing such calls.

### P02: target models removed from the AMR marker DBs
- Marker FASTAs: CARD_AMR 4984 -> 4773 entries (VM 178, KM 18, OM 13 and 2 HM-prefixed baeS/baeR records removed); AMRF_core-AMR 8263 -> 7995. Unfiltered sources and removal reports in `DB/seq/marker_sources/`.
- Run: only `CARD_AMR.dmnd`, `AMRF_core-AMR.dmnd` and their two hit TSVs were rebuilt; the other 11 hit TSVs are byte-identical. CARD_AMR TSV rows 135770 (HM 134404, VM 891, KM 305, OM 170) -> 134404 (all HM), CDS 6766 -> 5586; AMRF_core-AMR rows 208214 -> 205589, CDS 5130 -> 2505; AMRF_plus-AMR 598 rows / 598 CDS unchanged. 0 new query/subject pairs, 0 bitscore changes, 0 best-hit switches (raw TSVs and cache).
- Cache (schema 13 -> 14, 574416 annotations, same AnnotationIDs and 153 samples): A-tagged CDS 8551 -> 5929 (per genome 34/54/80 -> 20/36/64); distinct A-tagged genes per genome 33/51/58 -> 20/35/41; Sample tree ATB 28/35/39 -> 21/30/32. 0 CDS gained A.
- Tag transitions: A -> none 1841; A -> V 306 (lpxA, lpxC); AI -> MI 170 (tetR(B) 169 + tetR 1); A -> B 153 (tolC); AI -> I 152 (oprD).
- CDS losing A (CDS in genomes): opgE 218 in 126, baeS 207 in 153, tetR(B) 169 in 133, gyrA, parC, ftsI, lpxA, lpxC, tolC, rpsL, rpsJ, rpoB 153 each, oprD 152, baeR 147, pmrA 98, pmrB 98, tuf 93 in 76, pmrA_M12I 53, ompR 8, tetR 1 and 1 CDS without a gene name. Corrected: 20 gene names plus one unnamed CDS, and not "all chromosomal": tetR(B) is the Tn10 repressor on a mobile element (CARD OM model).
- Models gone: CARD VM Abau_gyrA_FLO, Abau_parC_FLO, LpxA, LpxC (153 each), Ecol_EFTu_KIR 93; KM Abau_AbuO 153, Abau_OprD_IPM 152; OM tetR 170 (134 genomes). AMRF: pmrC_R125P 218, adeL_G27R 182, adeR_D20N 149, adeS_F12S 147, tufA_L121Q 93, and 153 each for gyrA_G79C, parC_S84F, lpxA_G68D, lpxC_P30L, pmrA_E8D, pmrB_S17R, ftsI_A515V, rpoB_S583L, rpsL_K43R, rpsJ_V57L, baeS_D89V, baeR_S104N. Corrected: pmrB_S17R sat on 98 pmrB + 54 baeS + 1 unnamed CDS (not 55 baeS), which is why the GUI search "pmrB" returned 153 rows before and 98 after.
- Per gene: gyrA and parC A with mutation labels (153) -> no tag, predictor or allele; oprD AI -> I (152); adeR/adeS/adeL keep A from CARD HM with alleles adeR (149), adeS (147), adeL (181). CDS with a point-mutation allele label: 2625 -> 0.
- Acquired genes unchanged (CDS / genomes, tags identical): blaOXA-23 66/43 (64 AI, 2 A), blaOXA-72 108/86, blaOXA-24 4/2, blaOXA-235 8/2, blaNDM-1 24/19, blaOXA-66 180/138, blaOXA 10/10; no OXA-58-like hit before or after. The GUI search "blaOXA-23" still returns 78 rows in 49 samples (a substring search: 70 rows with allele blaOXA-23 plus 8 blaOXA-235).
- Sample 20038ENE2: card Resistance genes 53 -> 35, ATB classes 32 -> 26 (Colistin, Quinolone, Elfamycin, Monobactam, Nybomycin-like, Polyketide gone; Fluoroquinolone stays through efflux genes); Sample tree AMRF_core-AMR 32 -> 14, CARD_AMR 41 -> 33, ATB 38 -> 31.
- Affected outside the app: `build_st492_res_panel.py` (its Colistin row pmrA, pmrB, pmrA_M12I, lpxA, lpxC, lpsB) and `plot_st492_tree_resistance_phandango.R` through its CSV, `build_supplementary_annotation_workbook.py`, `build_supplementary_gene_calls.py`, `build_supplementary_raw_outputs.py`, the synteny scripts. Not affected: `build_st492_vf_panel.py`, `build_st492_phage_tables.py`.

### P04: no first-row-wins in the metadata indexes
- Indexes (old vs new code on the same catalogs): AMRF keys 36949 -> 36085, 864 target-only, 417 ambiguous (0 protein-style accessions); CARD by_aro 6440 -> 6088, by_protein 6324 -> 6086, 390 target AROs, 3 ambiguous protein keys (AEZ36150.1, NP_214519.2, NP_214520.1).
- Run: 8689 AMR hit rows kept with the same best subject; 0 changes in any of the six Resistance* fields. No visible change on this run beyond P02.

### P05: unique genome sample names
- Run: 0 collisions among 153 genomes; adding `/other/20038ENE2.fasta` is reported as a collision with `/Data/ST492/final_fasta/20038ENE2.fna`. Cache samples and chewBBACA FILE rows: same 153 before and after. No data or GUI change for this run.

### P06: genome set guard
- Run: the guard returns no error for the run's 153 genomes; +NEW1 is refused (153 recorded, 154 given) and -SAMN04040115 is refused (152 given). `validate_completed_run` after the rerun: 0 errors (the before state was not re-validated here). `gubbins/genome_set.json` does not exist (Gubbins was reused). No data or GUI change for this run.

### Caveats
- The COG change (marker_hits 535240 -> 53758 rows, 477677 of them COG; sample card Marker hits 3386 -> 307 and Marker DBs 11 -> 10; db_total_hits; COG Predictor) is not P02. It happened in this rerun: COG was dropped from `run_config.json` (its DB was removed on 2026-08-11) and its TSV and .dmnd were moved to `BACTolomew_backups/2026-10-09_rerun_moved_out/`. `build_supplementary_raw_outputs.py` will no longer find the COG TSV.
- The "before" GUI screenshots used the current code on the backup data, not the old code.
- Schema checks for the 33 non-A. baumannii schemas prove orientation, validity and self-score counts; their smoke run used A. baumannii genomes and only proves that they run.
- Point mutations are no longer reported at all; the old labels were reference names of best BLASTp hits, not called mutations.
- Manuscript figures built from `results_alleles.tsv` or from AMR hits are out of date until regenerated; no script was run.

## P20: AMR point mutations (AMRFinderPlus --organism)

### Problem
- Since P02 the DIAMOND screen no longer reports point-mutation references (AMRF POINT, CARD VM) as acquired genes, so known resistance mutations (gyrA S81L, parC S84L in ST2/ST492) were not reported at all. Best BLASTp hits cannot call a mutation; only a residue check against the reference can.

### Fix
- New pipeline step `amrfinder_mutations` (`src/pipeline/AmrFinderWorkflow.py`), run by `Runner` right after `mlst` (it needs the Bakta outputs and the MLST species). Per genome, in the container's `bakta` env (`ToolEnv`: tool `amrfinder` -> env `bakta`): `amrfinder --protein <bakta>.faa --nucleotide <bakta>.fna --gff <bakta>.gff3 --annotation_format bakta --organism <org> --plus --database <bakta_db>/amrfinderplus-db/latest --threads 4 --name <sample> --output amrfinder/<sample>.tsv --quiet`; genomes run in parallel (thread plan `amrfinder`, min(cores, 32), 4 threads each).
- Organism (`resolve_amrfinder_organism`): the config override `amrfinder.organism` if set, else the majority MLST scheme mapped to an AMRFinderPlus taxgroup (`MLST_SCHEME_ORGANISMS`, checked against the DB's `taxgroup.tsv`; a tie or an ambiguous scheme such as `neisseria` gives none), else the cgMLST schema's species. Without an organism the step succeeds with no calls and says why in `amrfinder/summary.json` (`skipped`).
- `amrfinder/summary.json`: organism and how it was chosen, AMRFinderPlus version, DB dir and version, samples, total rows and POINT rows. The step's record in `run_report.json` carries organism, DB version and POINT rows; `amrfinder` is in `tools.versions`.
- Reuse: `validate_amrfinder_outputs` (summary, sample set, one report per genome with the expected header) plus `amrfinder_reuse_issue` (organism or DB version changed, or a report older than its Bakta proteins). The step is in `RunValidation.RESULT_OPTIONAL_STEPS`: validated when a run recorded it, not required, so runs made before P20 still open (they simply have no mutation calls).
- Data: only POINT and POINT_DISRUPT rows are read; acquired-gene rows stay in the reports but are not used (gene calls stay with DIAMOND). `attach_resistance_mutations` sets `ResistanceMutations` (';'-joined element symbols) and `ResistanceMutationClasses` (';'-joined lower-case subclass and class) on the CDS, matched by Bakta locus tag (AMRFinderPlus `Protein id`), else by contig and coordinates. Used by both `Runner.collect_annotations` and the GUI build (`data.build_annotations_from_output`).
- GUI: tag `X` (resistance mutation) for every CDS with `ResistanceMutations`; cache schema 14 -> 15 and the cache fingerprint now includes `amrfinder/*.tsv` and `summary.json`. X is drawn in the AMR red (`annotation_table._tag_color`, gene card `_igv_like_color`, `plasmid_model._gene_color`); the Annotations table draws X as a red letter followed by the shared mutation badge (`assets/mutation_badge.svg`) with the tooltip from `resistance_mutation.mutation_note` ("Usually not a resistance gene; carries gyrA_S81L (quinolone), a known resistance mutation (AMRFinderPlus)"). Tags filter: "X: Resistance mutation". Gene card: a "Resistance mutation" chip and a "Resistance mutations" section (badge, note, mutations, drug classes). Sample card: a separate "Resistance mutations" row (count and symbols); the mutations' drug classes join the ATB classes with the source "AMRFinderPlus mutation". The neighbourhood and genome track badges are in the views lane.

### Verification (not on the reference run)
- 3 genomes copied to scratch (20038ENE2 and 21001CNE5, ST492; GCA_000804885.1, ST2), run through `run_amrfinder_mutations` with the pipeline's container wrapper: organism Acinetobacter_baumannii from MLST scheme `abaumannii_2` (3/3), AMRFinderPlus 4.0.23, DB 2025-03-25.1, about 19 s per genome. 60 report rows, 6 POINT rows: gyrA_S81L and parC_S84L (QUINOLONE) in each genome. Each landed on the Bakta gyrA / parC CDS by locus tag (e.g. 20038ENE2 PIDAKJ_00437 gyrA, PIDAKJ_03473 parC); those CDS get tag X and no A (P02 removed their A).
- Reuse check: unchanged inputs reuse; a newer Bakta `.faa` forces a rerun; an extra genome fails validation.
- Expected on the reference run (153 genomes, not run here): gyrA and parC tagged X in most or all genomes; runtime about 6 to 7 minutes at 8 genomes in parallel. The GUI cache rebuilds once (schema 15).
- Tests: `tests/test_amr_point_mutations.py` (POINT rows tag the matching CDS by locus tag and by coordinates, acquired rows and STRESS point rows are ignored, the badge note; organism from MLST and skip without one).

### Review follow-up
- Non-fatal: a genome AMRFinderPlus fails on no longer stops the run. It is listed under `failed` in `summary.json` (no report, so no calls for it), the validator accepts its missing report, and `amrfinder_reuse_issue` reruns the step next time. Any other error in the step records it as `skipped` (with the error as the reason) and the run goes on; `validate_completed_run` ignores a skipped optional step.
- Only `Type` AMR point rows count (`_amr_point_rows`), so STRESS or VIRULENCE point mutations reported with `--plus` never get X. A POINT_DISRUPT call keeps its symbol with the suffix " (gene disrupted)" so the note and cards say the gene is disrupted.
- Calls that match no Bakta CDS (e.g. promoter mutations from the AMR_DNA tables) are logged as a warning by `attach_resistance_mutations` (none in the 153 genome scratch run).
- When AMRFinderPlus assessed the CDS as another gene than Bakta named it (pmrB_I164F on Bakta `baeS`), the note adds "which reads this CDS as pmrB".
- A reused step records the same fields as a fresh one (organism source and POINT rows read back from `summary.json`).

### Reference run (run_78eaa9126850, 153 genomes, 2026-10-09)
- Driver: the pipeline's own `Runner.run_workflow` with the GUI Resume `DataInput` (96 cores; thread plan `amrfinder` 32, so 8 genomes at a time with 4 threads each). Plan check first with the Runner's reuse checks: assembly_qc, Bakta, MLST, MOB-suite, Panaroo, Gubbins, phylogeny, fastGEAR, chewBBACA, diamond_markers and diamond_reference_loci all reused, the marker hit TSVs current (collect_annotations re-merged in memory only, no DIAMOND search); only `amrfinder_mutations` ran.
- File check (mtime and size of every file in the run folder, before and after): the only additions are `amrfinder/` (153 reports plus `summary.json`); the only other changes are `run_config.json` (`amrfinder_threads` and thread plan entry, `created_at`), `expected_outputs.json` (`amrfinder/summary.json`), `run_state.json` and `run_report.json`. Nothing else was rewritten.
- Wall times: `amrfinder_mutations` 678 s (11.3 min, so the estimate of 6 to 7 minutes above was low), collect_annotations 127 s, whole `run_workflow` 831 s. GUI cache rebuild the way the app does it (`validate_output_dir(deep=True)` then `load_output_data` then `write_cached_tables`): 238 s (deep validation 10.5 s with 0 errors, load 217 s), schema 15, 574416 annotation rows.
- Run state: status success, every step `outputs_validated`, `validate_completed_run` 0 errors; the step record has organism Acinetobacter_baumannii (MLST scheme `abaumannii_2`, 153/153 genomes), DB 2025-03-25.1, 321 POINT rows, no failed genomes; `tools.versions.amrfinder` 4.0.23.
- Reports: 3861 rows in all (AMR/AMR 3160 acquired genes, kept but not used; AMR/POINT 321; STRESS/BIOCIDE 138 and STRESS/METAL 242, ignored). No POINT_DISRUPT rows. All 321 POINT rows landed on a Bakta CDS by locus tag (0 unmatched, 0 coordinate fallbacks needed).
- MLST: 93 genomes ST492, 60 ST2.

| Mutation | Genomes (= CDS) | ST492 | ST2 | Drug class (subclass) | Bakta gene of the CDS |
|---|---|---|---|---|---|
| gyrA_S81L | 147 | 87 | 60 | quinolone | gyrA |
| parC_S84L | 141 | 87 | 54 | quinolone | parC |
| ftsI_A515V | 19 | 0 | 19 | beta-lactam (carbapenem, sulbactam-durlobactam) | ftsI |
| adeS_G336S | 6 | 6 | 0 | efflux | adeS |
| parC_E88K | 4 | 0 | 4 | quinolone | parC |
| pmrB_I164F | 2 | 2 | 0 | colistin | baeS (AMRFinderPlus reads it as pmrB) |
| pmrA_M12I | 1 | 0 | 1 | colistin | ompR (AMRFinderPlus reads it as pmrA) |
| pmrB_P233S | 1 | 0 | 1 | colistin | no gene name (product "ATP-binding protein", GCA_029713925.1 LOCDHM_01311) |

- Per gene (CDS with X): gyrA 147, parC 145 (141 S84L plus 4 E88K, never both), ftsI 19, adeS 6, baeS 2, ompR 1, unnamed 1; 321 CDS with X in total, one mutation per CDS. Every genome carries at least one X.
- Per genome: 120 gyrA_S81L + parC_S84L (85 ST492, 35 ST2); 17 ftsI_A515V + gyrA_S81L + parC_S84L (ST2); 6 adeS_G336S only (ST492, the GCA_0234990xx series, wild type gyrA and parC); 4 gyrA_S81L + parC_E88K (ST2); 2 gyrA_S81L only (ST2); 2 gyrA_S81L + parC_S84L + pmrB_I164F (ST492); 1 ftsI_A515V + gyrA_S81L + parC_S84L + pmrA_M12I and 1 ftsI_A515V + gyrA_S81L + parC_S84L + pmrB_P233S (both ST2).
- Drug classes on the X CDS: quinolone 292, carbapenem;sulbactam-durlobactam;beta-lactam 19, efflux 6, colistin 4.
- CDS with both A and X: 6, all the adeS_G336S calls (adeS keeps its DIAMOND A tag as a regulator of the AdeABC pump), so they draw in the A color with the badge. No other X CDS has A.
- Spot checks against the Bakta proteins in the cache: 20038ENE2 (ST492, gyrA PIDAKJ_00437 / parC PIDAKJ_03473), GCA_000804885.1 (ST2, GOFMDK_02370 / GOFMDK_00973) and GCA_000810385.3 (ST2, OOILJH_01810 / OOILJH_00715) all read HGDL(81)AVY in gyrA and HGDL(84)ACY in parC, i.e. L at 81 and 84 as called. Over all 153 genomes the residue agrees with the call every time: gyrA 81 is L in the 147 called and S in the 6 others; parC 84 is L in the 141 S84L calls and S in the other 12; parC 88 is K in the 4 E88K calls and E or G elsewhere.
