# UCT Applied Proteomics

This repository contains a Python-based workflow for reference-based taxonomic assignment of archaeological bone samples using LC-MS/MS peptide data. The project is built around a set of reusable analysis scripts for loading peptide data, cleaning it, assigning biological groups, creating marker libraries, and comparing unknown samples against reference groups.

## Default workflow

The workflow implemented by the scripts follows this pattern:

1. Load peptide/protein data from Excel files
2. Preprocess and clean the peptide tables
3. Assign sample groups and taxonomy labels
4. Build a reference marker library
5. Compare unknown samples against reference markers
6. Export summary reports and Excel outputs

## Typical analysis pipeline

```python
from Scripts.loaders import load_peptides
from Scripts.preprocessing import adjust_peptides
from Scripts.identifiers import get_group_names
from Scripts.reference_analysis import reference_markers_robust
from Scripts.unknown_analysis import analyse_unknown, analyse_taxonomy, analyse_species
from Scripts.exporters import export_results, export_taxonomy

# 1. Load raw peptide data from Excel sheets in a folder
peptides = load_peptides("data/raw")

# 2. Remove false annotations, duplicates, and unnecessary columns
peptides = adjust_peptides(peptides, duplicates=False, cut_int=None)

# 3. Assign group names and optional taxonomy labels
# Example taxonomy map: {"cow": "Bovidae", "pig": "Suidae"}
peptides = get_group_names(
    peptides,
    taxonomy_map={"cow": "Bovidae", "pig": "Suidae"},
    merge_map=None
)

# 4. Build reference marker libraries
reference_groups = sorted(peptides["Group"].unique())
taxonomy_groups = sorted(peptides["Taxonomy"].unique())

species_markers = reference_markers_robust(
    peptides,
    reference_groups,
    rep_in=1,
    rep_out=1
)

taxonomy_markers = reference_markers_robust(
    peptides,
    reference_groups,
    taxonomy_groups=taxonomy_groups,
    rep_in=1,
    rep_out=1
)

# 5. Analyse unknown samples
unknown_groups = ["Unknown_1", "Unknown_2"]
unknown_results = analyse_unknown(peptides, species_markers, unknown_groups, replicates=1, keys=True)

taxonomy_results, best_clusters = analyse_taxonomy(peptides, taxonomy_markers, unknown_groups, replicates=1)

species_results = analyse_species(
    best_clusters,
    species_markers,
    {
        "Bovidae": ["cow", "sheep"],
        "Suidae": ["pig"]
    },
    peptides,
    unknown_groups,
    replicates=1,
)

# 6. Export Excel reports
export_results(
    unknown_results,
    species_markers,
    peptides,
    output_folder="results",
    total_rep=10,
    ref_replicates=2,
    un_replicates=1,
    keys=True,
)

export_taxonomy(
    taxonomy_results,
    species_results,
    taxonomy_markers,
    species_markers,
    peptides,
    output_folder="results",
    total_rep=10,
    ref_replicates=2,
    un_replicates=1,
)
```

## Script overview

### 1. loaders.py
Responsible for reading Excel files and loading peptide/protein tables into pandas DataFrames.

- `load_peptides_proteins(folder, ignored=None)`
- `load_peptides(folder, ignored=None)`
- `load_marker_library(folder, file, key=True)`

### 2. preprocessing.py
Cleans peptide and protein data before further analysis.

- `adjust_peptides(peptides, duplicates=False, cut_int=None)`
- `remove_protein_duplicates(proteins, cut_int=None)`

This step removes untrusted peptide annotations, drops duplicates, removes unnecessary columns, and can optionally filter low-intensity values.

### 3. identifiers.py
Assigns sample groups from sample IDs, and optionally applies taxonomy grouping.

- `get_group_names(peptides, taxonomy_map=None, merge_map=None)`

The function derives a biological group label from the SampleID and can map groups into higher taxonomic categories.

### 4. reference_analysis.py
Builds a marker library of peptide markers that are characteristic for each reference group.

- `reference_markers(peptides, reference_groups, min_replicates=1)`
- `reference_markers_robust(peptides, reference_groups, taxonomy_groups=None, rep_in=1, rep_out=1)`

These functions find peptides that are robustly present in a target group and absent or rare in comparison groups.

### 5. unknown_analysis.py
Compares unknown samples against the reference marker library to estimate the best matches.

- `analyse_unknown(...)`
- `analyse_taxonomy(...)`
- `analyse_species(...)`

This produces match scores and confidence values for each unknown sample against reference species or taxonomy groups.

### 6. exporters.py
Exports analysis outputs to CSV and Excel files.

- `export_marker_library(...)`
- `export_results(...)`
- `export_taxonomy(...)`

Outputs typically include summary tables, detailed marker information, and reference libraries for downstream interpretation.

## Recommended project layout

```text
.
├── data/
│   ├── raw/
│   └── results/
├── Scripts/
│   ├── loaders.py
│   ├── preprocessing.py
│   ├── identifiers.py
│   ├── reference_analysis.py
│   ├── unknown_analysis.py
│   └── exporters.py
├── README.md
└── requirements.txt
```

## Notes

- Input data are assumed to be Excel files with one or more sheets containing peptide data.
- The workflow is designed around a reference-based classification approach where unknown samples are scored against known reference groups.
- Marker robustness is controlled by replicate thresholds (`rep_in`, `rep_out`, `replicates`) used during library construction and sample matching.

## Typical output

The workflow generates files such as:

- `Analysis_report[...].xlsx`
- `Analysis_report[...].xlsx_taxonomy`
- CSV marker libraries
- Final summary tables for unknown sample assignment

## License

This project does not currently specify a license in the repository files reviewed here. If needed, add one explicitly before public distribution.
