A Python workflow for reference-based taxonomic assignment of archaeological bone samples using LC-MS/MS peptide data. This project provides a modular analysis pipeline for comparing unknown proteomic samples against reference marker libraries to identify taxonomic origin.

## Workflow Overview

1. **Full Pipeline** – Load data in .xlsx format → preprocess → build reference markers → analyze unknown samples
2. **Data Loading & Processing** – Load and clean LC-MS/MS data
3. **Library Creation** – Build reference marker libraries from processed data (single-level or hierarchical)
4. **Unknown Sample Analysis** – Compare unknown samples against existing reference libraries

## Running the Analysis

Execute the main script:

```bash
python main.py
```

### Configuration

Edit `main.py` to customize:
- **Reference groups:** Animals/species with known origin
- **Unknown samples:** Samples/species to be identified
- **Taxonomy mapping:** Grouping of species into higher taxonomic categories
- **Merger mapping:** Grouping of species together into a larger cluster (for example, chicken and hen)
- **Processed data:** CSV file of loaded and preprocessed data

Example configuration:
```python
REFERENCE_GROUPS = ["cat", "cattle", "deer", "dog", "duck", ...]
UNKNOWN_SAMPLES = ["1", "2", "3", ...]
taxonomy_groups = {
    "carnivores": ["dog", "cat"],
    "birds": ["duck", "hen", "quail", "turkey"],
    ...
}
```

You will be prompted to select an analysis type:
- `full` – Complete pipeline from raw data to results
- `load` – Data loading and preprocessing only
- `lib` – Create reference marker libraries
- `unk` – Analyze unknown samples using existing libraries

Once data have been loaded and preprocessed, they are saved into Dataframe folder, from which it may easily be loaded for further exploration

## Project Structure

```
.
├── Data/                   # Raw LC-MS/MS Excel files
├── Dataframe/              # Processed peptide DataFrames (CSV)
├── Libraries/              # Reference marker libraries (CSV)
├── Results/                # Analysis output Excel files
├── Logs/                   # Analysis logs with timestamps
├── Scripts/
│   ├── loaders.py          # Load Excel data into DataFrames
│   ├── preprocessing.py    # Clean and filter peptide data
│   ├── identifiers.py      # Assign sample groups and taxonomy
│   ├── reference_analysis.py  # Build marker libraries
│   ├── unknown_analysis.py    # Compare unknown samples
│   └── exporters.py        # Export results to Excel/CSV
└── main.py                 # Main entry point
```

## Key Functions

### Data Loading (`loaders.py`)
- `load_peptides(folder, ignored=None)` – Load peptide data from Excel files
- `load_marker_library(folder, file, key=True)` – Load reference marker library from CSV

### Preprocessing (`preprocessing.py`)
- `adjust_peptides(peptides, duplicates=False, cut_int=None)` – Remove false annotations, duplicates, and low-intensity peptides

### Sample Grouping (`identifiers.py`)
- `get_group_names(peptides, taxonomy_map=None, merge_map=None)` – Extract group names from sample IDs and apply optional taxonomy mapping

### Marker Library Creation (`reference_analysis.py`)
- `reference_markers_robust(peptides, reference_groups, taxonomy_groups=None, rep_in=1, rep_out=1)` – Build robust marker libraries with configurable replicate thresholds

### Unknown Sample Analysis (`unknown_analysis.py`)
- `analyse_unknown(peptides, markers, unknown_groups, replicates=1, keys=True)` – Single-level species matching
- `analyse_taxonomy(peptides, markers, unknown_groups, replicates=1)` – Higher-level taxonomy assignment
- `analyse_species(cluster, markers, taxonomy, peptides, unknown_groups, replicates=1)` – Species-level matching within assigned taxonomy

### Results Export (`exporters.py`)
- `export_marker_library(output_folder, markers, peptides, ...)` – Save marker library to CSV
- `export_results(analysed_samples, markers, peptides, ...)` – Export single-level analysis to Excel
- `export_taxonomy(analysed_taxa, analysed_species, ...)` – Export hierarchical analysis to Excel

## Analysis Modes

### Species-Level Analysis
Compares unknown samples directly against all reference species in a single matching step.

### Hierarchical Taxonomy Analysis
Two-step approach: unknown samples are first matched to higher taxonomy groups (e.g., carnivores, birds), then matched against species within the best-matching group.

## Output Files

Analysis generates:
- **Analysis_report[*].xlsx** – Species-level results with summary, detailed markers, and reference library
- **Analysis_report[*]_taxonomy.xlsx** – Hierarchical results with taxonomy assignments and species-level matches
- **Marker_library[*].csv** – Reference peptide marker library
- **Analysis_log_[timestamp].txt** – Complete analysis log with all output and user inputs

## Interpretation
Resulting Excel file contains the following columns:
- **Markers Found** - Number of reference markers found for selected reference species within the unknown sample
- **Confidence** - Percentage of reference markers for selected reference species out of ALL reference markers found
- **Reference Markers** - Number of reference markers available in the reference library for selected species
- **Match Score** - Equal to marker recovery: percentage of how many reference markers from the reference library were found within the unknown sample
- **Replicates** - Number of both technical and biological replicates present for the unknown sample  
   
Reference species with the highest confidence for selected unknown sample is highlighted by green colour and bold text  
Confidence should always be interpreted in the context of the reference species library, as it is computed only from all present species

## Citation
To be added
