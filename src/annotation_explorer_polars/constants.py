from __future__ import annotations

from typing import List, Tuple

# Expected cache filenames
ANNOTATIONS_FILENAME = "annotations.parquet"
REFERENCE_LOCUS_FILENAME = "reference_locus_hits.parquet"
MARKER_HITS_FILENAME = "marker_hits.parquet"
RUN_CONFIG_FILENAME = "run_config.json"
RUN_REPORT_FILENAME = "run_report.json"
RUN_STATE_FILENAME = "run_state.json"
EXPECTED_OUTPUTS_FILENAME = "expected_outputs.json"
RUN_LOCK_FILENAME = "run.lock.json"
RUN_STOP_FILENAME = "run.stop"
MARKER_SEARCH_TEXT_COLUMN = "MarkerSearchText"

# Expected columns
EXPECTED_ANNOTATION_COLUMNS: List[str] = [
    "Sample",
    "Contig",
    "Start",
    "Stop",
    "Strand",
    "Gene",
    "Product",
    "Tags",
    "Identity",
    "Evalue",
    "Bitscore",
    "Method",
    "Predictors",
    "MarkerCategory",
    MARKER_SEARCH_TEXT_COLUMN,
    "NCBI_NRP_ID",
    "Refseq_ID",
    "Id",
    "Locus",
    "StartType",
    "RBS",
    "Protein_MW",
    "pI",
    "Hypothetical",
    "Nucleotide_Sequence",
    "Amino_Acid_Sequence",
    "AnnotationID",
    "ConsensusCategory",
    "ConsensusSource",
    "ConsensusDonorProduct",
    "ResistanceMetadataSource",
    "ResistanceAlleles",
    "ResistanceGeneFamilies",
    "ResistanceDrugClasses",
    "ResistanceAntibiotics",
    "ResistanceMechanisms",
    "ResistanceMutations",
    "ResistanceMutationClasses",
    "MobSuitePrediction",
    "MobSuiteIsPlasmid",
    "MobSuiteSourceContig",
    "MobSuitePlasmidID",
    "MobSuiteMoleculeType",
    "MobSuitePrimaryClusterID",
    "MobSuiteSecondaryClusterID",
    "MobSuiteLength",
    "MobSuiteRepliconTypes",
    "MobSuiteRelaxaseTypes",
    "MobSuiteMPFType",
    "MobSuiteOriTTypes",
    "MobSuitePredictedMobility",
    "MobSuiteMashNearestNeighbor",
    "MobSuiteMashNeighborDistance",
    "GeneID",
]

EXPECTED_REFERENCE_LOCUS_COLUMNS: List[str] = [
    "Sample",
    "Contig",
    "Start",
    "Stop",
    "Strand",
    "Identity",
    "Evalue",
    "Bitscore",
    "Method",
    "MarkerCategory",
    "Product",
]

# Metadata filtering sentinel for missing values
MISSING_FILTER_TOKEN = "__MISSING__"

# Normalized truthy values used for MOB-suite plasmid calls.
PLASMID_TRUE_VALUES = frozenset({"true", "1", "yes", "y", "p", "plasmid"})

# Annotation columns shown in the table (excluding Identity/Evalue/Bitscore/Refseq)
ANNOTATION_COLUMNS_TO_SHOW: List[Tuple[str, str]] = [
    ("Start", "Start"),
    ("Stop", "Stop"),
    ("Strand", "Strand"),
    ("Tags", "Tags"),
    ("Predictors", "Predictors"),
    ("Gene", "Gene"),
    ("Product", "Product"),
    ("MarkerDisplay", "Marker"),
    ("Details", ""),
]

# Country aliases used by Geo map normalization.
# Keys should be lowercase human-entered variants; values are canonical labels.
COUNTRY_ALIASES = {
    # United States
    "usa": "United States",
    "u s a": "United States",
    "u.s.a": "United States",
    "u.s.a.": "United States",
    "us": "United States",
    "u.s.": "United States",
    "united states": "United States",
    "united states of america": "United States",
    "america": "United States",
    # United Kingdom
    "uk": "United Kingdom",
    "u.k.": "United Kingdom",
    "gb": "United Kingdom",
    "g.b.": "United Kingdom",
    "great britain": "United Kingdom",
    "britain": "United Kingdom",
    "england": "United Kingdom",
    "scotland": "United Kingdom",
    "wales": "United Kingdom",
    "northern ireland": "United Kingdom",
    # Korea
    "korea, south": "South Korea",
    "south korea": "South Korea",
    "korea south": "South Korea",
    "republic of korea": "South Korea",
    "rok": "South Korea",
    "korea, north": "North Korea",
    "north korea": "North Korea",
    "korea north": "North Korea",
    "democratic people's republic of korea": "North Korea",
    "dprk": "North Korea",
    # China / Taiwan
    "pr china": "China",
    "p.r. china": "China",
    "people's republic of china": "China",
    "peoples republic of china": "China",
    "mainland china": "China",
    "republic of china": "Taiwan",
    "taiwan, roc": "Taiwan",
    "taiwan roc": "Taiwan",
    # Russia / Czechia / common updates
    "russia": "Russian Federation",
    "czech republic": "Czechia",
    "slovak republic": "Slovakia",
    # Cote d'Ivoire
    "cote d ivoire": "Cote d'Ivoire",
    "cote d'ivoire": "Cote d'Ivoire",
    "ivory coast": "Cote d'Ivoire",
    # Lao / Vietnam / Iran / Syria / Tanzania
    "laos": "Lao People's Democratic Republic",
    "lao pdr": "Lao People's Democratic Republic",
    "viet nam": "Vietnam",
    "iran": "Iran (Islamic Republic of)",
    "syria": "Syrian Arab Republic",
    "tanzania": "United Republic of Tanzania",
    # Long official names often shortened in metadata
    "venezuela": "Venezuela (Bolivarian Republic of)",
    "bolivia": "Bolivia (Plurinational State of)",
    "moldova": "Republic of Moldova",
    "brunei": "Brunei Darussalam",
    "micronesia": "Micronesia (Federated States of)",
    "cape verde": "Cabo Verde",
    "swaziland": "Eswatini",
    "macedonia": "North Macedonia",
    "burma": "Myanmar",
    "east timor": "Timor-Leste",
    "vatican": "Holy See",
    "vatican city": "Holy See",
    # Congo disambiguation
    "drc": "Democratic Republic of the Congo",
    "dr congo": "Democratic Republic of the Congo",
    "drc congo": "Democratic Republic of the Congo",
    "congo kinshasa": "Democratic Republic of the Congo",
    "democratic republic of congo": "Democratic Republic of the Congo",
    "democratic republic of the congo": "Democratic Republic of the Congo",
    "congo brazzaville": "Republic of the Congo",
    "republic of congo": "Republic of the Congo",
    "republic of the congo": "Republic of the Congo",
    # Gulf / Middle East variants
    "uae": "United Arab Emirates",
    "u.a.e.": "United Arab Emirates",
    "palestine": "State of Palestine",
    "palestinian territories": "State of Palestine",
    "palestinian territory": "State of Palestine",
    "occupied palestinian territory": "State of Palestine",
    "west bank and gaza": "State of Palestine",
    "west bank": "State of Palestine",
    "gaza strip": "State of Palestine",
    # Caribbean and small-state abbreviations
    "bahamas, the": "Bahamas",
    "the bahamas": "Bahamas",
    "gambia, the": "Gambia",
    "the gambia": "Gambia",
    "antigua & barbuda": "Antigua and Barbuda",
    "trinidad & tobago": "Trinidad and Tobago",
    "st kitts and nevis": "Saint Kitts and Nevis",
    "st. kitts and nevis": "Saint Kitts and Nevis",
    "st lucia": "Saint Lucia",
    "st. lucia": "Saint Lucia",
    "st vincent and the grenadines": "Saint Vincent and the Grenadines",
    "st. vincent and the grenadines": "Saint Vincent and the Grenadines",
    "sao tome and principe": "Sao Tome and Principe",
    "sao tome & principe": "Sao Tome and Principe",
    # Europe and nearby
    "bosnia": "Bosnia and Herzegovina",
    "bosnia herzegovina": "Bosnia and Herzegovina",
    "republic of ireland": "Ireland",
    "the netherlands": "Netherlands",
    "holland": "Netherlands",
    "turkiye": "Turkey",
    # Alternate punctuation/spelling
    "cabo verde": "Cabo Verde",
    "eswatini": "Eswatini",
    "myanmar (burma)": "Myanmar",
    "timor leste": "Timor-Leste",
    "syria arab republic": "Syrian Arab Republic",
    "iran islamic republic of": "Iran (Islamic Republic of)",
    "venezuela bolivarian republic of": "Venezuela (Bolivarian Republic of)",
    "bolivia plurinational state of": "Bolivia (Plurinational State of)",
    "micronesia federated states of": "Micronesia (Federated States of)",
    "lao peoples democratic republic": "Lao People's Democratic Republic",
}
