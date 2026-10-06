"""
This main script calls upon all other modules and utilises various functions to perform full analysis
of LC-MS/MS proteomic data
"""
import re
import datetime
import builtins

from Scripts.loaders import *
from Scripts.unknown_analysis import *
from Scripts.exporters import *
from Scripts.preprocessing import *
from Scripts.identifiers import *
from Scripts.reference_analysis import *

class DualLogger:
    """
    Intercepts a console stream (stdout or stderr) and writes it to both
    the terminal and a shared log file.
    """
    def __init__(self, terminal_stream, log_file):
        self.terminal = terminal_stream
        self.log = log_file

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

if __name__ == "__main__":
    FOLDER_PATH = Path(__file__).resolve().parent
    DATA_FOLDER = FOLDER_PATH / "Data"
    DATAFRAME_FOLDER = FOLDER_PATH / "Dataframe"
    LIBRARY_FOLDER = FOLDER_PATH / "Libraries"
    RESULTS_FOLDER = FOLDER_PATH / "Results"
    LOGS_FOLDER = FOLDER_PATH / "Logs"

    if not LOGS_FOLDER.exists():
        LOGS_FOLDER.mkdir(parents=True)

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file_path = LOGS_FOLDER / f"Analysis_log_{timestamp}.txt"

    # 1. Open the file once in append mode
    log_file = open(log_file_path, "a", encoding="utf-8")

    # 2. Redirect both stdout (normal text) and stderr (errors/tracebacks)
    sys.stdout = DualLogger(sys.stdout, log_file)
    sys.stderr = DualLogger(sys.stderr, log_file)

    original_input = builtins.input
    def logged_input(prompt=""):
        # original_input handles printing the prompt (which gets caught by DualLogger)
        user_response = original_input(prompt)

        # Write exactly what the user typed into the log file, plus a newline
        log_file.write(f"{user_response}\n")
        log_file.flush()

        return user_response

    # Replace the standard input function globally
    builtins.input = logged_input

    REFERENCE_GROUPS = ["cat", "cattle", "deer", "dog", "duck", "goat", "hen", "horse", "pig", "quail", "rabbit", "sheep", "turkey"]
    UNKNOWN_SAMPLES = ["1", "2", "3", "4", "5", "n", "ko"]

    LIBRARY_NAME_TAX = "Marker_library_taxonomy[15-15].csv"
    LIBRARY_NAME_SPEC = "Marker_library_species[15-15].csv"

    PROCESSED_DATA = ""

    IGNORED = None

    MERGER = None


    taxonomy_groups = {
        "carnivores": ["dog", "cat"],
        "non-carnivores": ["cattle", "deer", "goat", "horse", "pig", "rabbit", "sheep"],
        "birds": ["duck", "hen", "quail", "turkey"]
    }


    TAXONOMY = {
        animal : group
        for group, animals in taxonomy_groups.items()
        for animal in animals
    }


    process = input("-"*10 + "LC-MS/MS module" + "-"*10 + "\n"
                    "What do you want to do?\n"
                    "Perform the whole pipeline - from raw data to results? - Type 'full'\n"
                    "Load and process LC-MS/MS data? - Type 'load'\n"
                    "Load processed data and find reference markers? - Type 'lib'\n"
                    "Load processed data and analyse unknown samples? - Type 'unk'\n")

    if process == "full":
        print("You chose to run your data through the whole pipeline - loading, preprocessing and finding"
              "reference markers\n"
              "-----Check list:-----\n"
              f"I have all data files in {DATA_FOLDER}\n"
              "I have defined both reference sample and unknown sample names\n")
        input("Press Enter to continue...")
        # create a folder if it doesn't already exist
        dataframe_path = Path(DATAFRAME_FOLDER)
        if not dataframe_path.exists():
            dataframe_path.mkdir(parents=True)
            print(f"\nCreated a folder for processed datasets: {dataframe_path}")
        pep_df = load_peptides(DATA_FOLDER, ignored=None)
        pep_adj = adjust_peptides(pep_df, duplicates=False)
        atype = input(
            "Do you wish to perform species-level analysis or use higher taxonomy? - type 'taxonomy' or 'species'\n")
        if atype == "taxonomy":
            print("You chose higher taxonomy analysis")
            pep_ide = get_group_names(pep_adj, taxonomy_map=TAXONOMY, merge_map=MERGER)
            pep_ide.to_csv(DATAFRAME_FOLDER / PROCESSED_DATA)
            unk_df = pep_ide[pep_ide["Group"].isin(UNKNOWN_SAMPLES)]
            ref_df = pep_ide[pep_ide["Group"].isin(REFERENCE_GROUPS)]
            TOTAL_REPLICATES = int(ref_df.groupby("Group")["SampleID"].nunique().min())
            REP_IN = int(input(f"Please specify the number of target group replicates (range 1-{TOTAL_REPLICATES}): "))
            print("First step - taxonomy level markers:")
            library_tax = reference_markers_robust(ref_df, REFERENCE_GROUPS, list(taxonomy_groups.keys()), REP_IN)
            LIBRARY_TAX = input("\nSet file name of the taxonomy library:\n"
                                 f"Leave empty to use the example name: 'Marker_library_taxonomy_[{TOTAL_REPLICATES}-{REP_IN}]'\n")
            if LIBRARY_TAX == "":
                LIBRARY_TAX = f'Marker_library_taxonomy'
            export_marker_library(LIBRARY_FOLDER, library_tax, ref_df, REP_IN, TOTAL_REPLICATES, filename=LIBRARY_TAX)

            print("Second step - species-specific markers:")
            library_spec = reference_markers_robust(ref_df, REFERENCE_GROUPS,rep_in=REP_IN)
            LIBRARY_SPECIES = input("\nSet file name of the species library:\n"
                                 f"Leave empty to use the example name: 'Marker_library_species'\n")
            if LIBRARY_SPECIES == "":
                LIBRARY_SPECIES = f"Marker_library_species"
            export_marker_library(LIBRARY_FOLDER, library_spec, ref_df, REP_IN, TOTAL_REPLICATES, filename=LIBRARY_SPECIES)
            print(f"Reference marker libraries are now in {LIBRARY_FOLDER}.")
            results, taxonomy = analyse_taxonomy(unk_df, library_tax, UNKNOWN_SAMPLES)
            final_results = analyse_species(taxonomy, library_spec, taxonomy_groups, unk_df, UNKNOWN_SAMPLES)

            export_taxonomy(results, final_results, library_tax, library_spec, pep_ide, RESULTS_FOLDER,
                            TOTAL_REPLICATES, ref_replicates=REP_IN)
            print(f"Hierarchy analysis results are now in {RESULTS_FOLDER}.")

        elif atype == "species":
            print("You chose species-level analysis")
            pep_ide = get_group_names(pep_adj, merge_map=MERGER)
            pep_ide.to_csv(DATAFRAME_FOLDER / PROCESSED_DATA)
            unk_df = pep_ide[pep_ide["Group"].isin(UNKNOWN_SAMPLES)]
            ref_df = pep_ide[pep_ide["Group"].isin(REFERENCE_GROUPS)]
            TOTAL_REPLICATES = int(ref_df.groupby("Group")["SampleID"].nunique().min())
            REP_IN = int(input(f"Please specify the number of target group replicates (range 1-{TOTAL_REPLICATES}): "))
            library_spec = reference_markers_robust(ref_df, REFERENCE_GROUPS, rep_in=REP_IN)
            LIBRARY_SPECIES = input("\nSet file name of the library:\n"
                                 f"Leave empty to use the example name: 'Marker_library_species_[{TOTAL_REPLICATES}-{REP_IN}]'\n")
            if LIBRARY_SPECIES == "":
                LIBRARY_SPECIES = f"Marker_library_species[{TOTAL_REPLICATES}-{REP_IN}]"
            export_marker_library(LIBRARY_FOLDER, library_spec, ref_df, REP_IN, TOTAL_REPLICATES, filename=LIBRARY_SPECIES)
            print(f"Reference marker library is now in {LIBRARY_FOLDER}.")
            results = analyse_unknown(unk_df, library_spec, UNKNOWN_SAMPLES)
            export_results(results, library_spec, pep_ide, RESULTS_FOLDER, TOTAL_REPLICATES, REP_IN)
            print(f"Species-level analysis results are now in {RESULTS_FOLDER}.")
        else:
            print("Invalid input!")
            sys.exit(-1)

    elif process == "load":
        print("You chose to analyse LC-MS data - All available data will be concatenated and preprocessed\n"
              "-----Check list:-----\n"
              f"I have all data files in {DATA_FOLDER}\n"
              "I have set up taxonomy map and/or merger map\n"
              f"I have defined reference sample names\n")
        input("Press Enter to continue...")
        pep_df = load_peptides(DATA_FOLDER, ignored=None)
        pep_adj = adjust_peptides(pep_df, duplicates=False)
        pep_ide = get_group_names(pep_adj, taxonomy_map=TAXONOMY)

        # create a folder if it doesn't already exist
        dataframe_path = Path(DATAFRAME_FOLDER)
        if not dataframe_path.exists():
            dataframe_path.mkdir(parents=True)
            print(f"\nCreated a folder for processed datasets: {dataframe_path}")
        pep_ide.to_csv(dataframe_path / PROCESSED_DATA)
        print(f"\nSaved processed data to {dataframe_path / PROCESSED_DATA}.")
    elif process == "lib":
        print("You chose to create a library of reference markers - Preprocessed data will be loaded and filtered "
              "accordingly to obtain group-specific markers\n")
        analysis = input("What kind of analysis are you looking for?\n"
                        "Two-step hierarchical analysis? - Type 'hierarchy'\n"
                        "Single step species-level analysis? - Type 'default'\n")

        # create a folder if it doesn't already exist
        library_path = Path(LIBRARY_FOLDER)
        if not library_path.exists():
            library_path.mkdir(parents=True)
            print(f"\nCreated a folder for processed datasets: {library_path}")

        if analysis == "default":
            print("You chose the default species-level analysis\n")
            print("-----Check list:-----\n"
                  f"I have processed dataframe in {DATAFRAME_FOLDER}\n"
                  f"I have defined reference sample names\n")
            input("Press Enter to continue...")

            print(f"\nLoading data from {PROCESSED_DATA}...")
            pep_df = pd.read_csv(DATAFRAME_FOLDER / PROCESSED_DATA)
            ref_df = pep_df[pep_df["Group"].isin(REFERENCE_GROUPS)]
            print("\nData loaded successfully!")
            TOTAL_REPLICATES = int(ref_df.groupby("Group")["SampleID"].nunique().min())
            REP_IN = int(input(f"Please specify the number of target group replicates (range 1-{TOTAL_REPLICATES}): "))

            library = reference_markers_robust(ref_df, REFERENCE_GROUPS, rep_in=REP_IN)
            LIBRARY_NAME = input("\nSet file name of the library - for inspiration see below:\n"
                                 f"Marker_library_species\n"
                                 "Leave empty to use the example name:\n")
            if LIBRARY_NAME == "":
                LIBRARY_NAME = f"Marker_library_species"
            library_path = export_marker_library(LIBRARY_FOLDER, library, pep_df, REP_IN,
                                                 TOTAL_REPLICATES, filename=LIBRARY_NAME)
            print(f"Reference marker library is now in {library_path}.")

        elif analysis == "hierarchy":
            print("You chose to create a hierarchy of reference markers - first, markers of taxonomy groups will be created"
                  " and then, species-level analysis will be performed\n")
            print("-----Check list:-----\n"
                  f"I have processed dataframe in {DATAFRAME_FOLDER}\n"
                  f"I have defined reference sample names\n"
                  f"I have provided hierarchical grouping of reference samples\n")
            input("Press Enter to continue...")

            print(f"\nLoading data from {DATAFRAME_FOLDER}...")
            pep_df = pd.read_csv(DATAFRAME_FOLDER / PROCESSED_DATA)
            pep_ref = pep_df[pep_df["Group"].isin(REFERENCE_GROUPS)]
            print("\nData loaded successfully!")
            TOTAL_REPLICATES = int(pep_ref.groupby("Group")["SampleID"].nunique().min())
            REP_IN = int(input(f"Please specify the number of target group replicates (range 1-{TOTAL_REPLICATES}): "))

            print("First step - taxonomy level markers:")
            library = reference_markers_robust(pep_ref, REFERENCE_GROUPS, list(taxonomy_groups.keys()), REP_IN)
            LIBRARY_NAME = input("\nSet file name of the library:\n"
                                 f"Leave empty to use the example name: 'Marker_library_taxonomy'\n")
            if LIBRARY_NAME == "":
                LIBRARY_NAME = f'Marker_library_taxonomy'
            export_marker_library(LIBRARY_FOLDER, library, pep_ref, REP_IN, TOTAL_REPLICATES, filename=LIBRARY_NAME)

            print("Second step - species-specific markers:")
            library = reference_markers_robust(pep_ref, REFERENCE_GROUPS, rep_in=REP_IN)
            LIBRARY_NAME = input("\nSet file name of the library:\n"
                                 f"Leave empty to use the example name: 'Marker_library_species'\n")
            if LIBRARY_NAME == "":
                LIBRARY_NAME = f"Marker_library_species"
            export_marker_library(LIBRARY_FOLDER, library, pep_ref, REP_IN, TOTAL_REPLICATES, filename=LIBRARY_NAME)

            print(f"Reference marker libraries are now in {LIBRARY_FOLDER}.")
        else:
            print("Invalid input!")
            sys.exit(-1)
    elif process == "unk":
        print("You chose to analyse unknown samples by comparing them to previously found reference markers\n"
            "-----Check list:-----\n"
            f"I have all processed unknown sample data {DATAFRAME_FOLDER}\n"
            "I have defined unknown sample names\n"
            "I have created reference markers in library\n")
        input("Press Enter to continue...")

        pep_df = pd.read_csv(DATAFRAME_FOLDER / PROCESSED_DATA)
        unk_df = pep_df[pep_df["Group"].isin(UNKNOWN_SAMPLES)]
        ref_df = pep_df[pep_df["Group"].isin(REFERENCE_GROUPS)]
        TOTAL_REPLICATES = int(ref_df.groupby("Group")["SampleID"].nunique().min())
        analysis = input("What kind of analysis are you looking for?\n"
                        "Two-step hierarchical analysis? - Type 'hierarchy'\n"
                        "Single step species-level analysis? - Type 'default'\n")
        # create a folder if it doesn't already exist
        results_path = Path(RESULTS_FOLDER)
        if not results_path.exists():
            results_path.mkdir(parents=True)
            print(f"\nCreated a folder for processed datasets: {results_path}")

        if analysis == "default":
            print("You chose default species-level analysis\n"
                  "-----Check list:-----\n"
                  "I have defined unknown sample names\n"
                  "I have loaded unknown sample data\n"
                  f"I have available species_level marker library in {LIBRARY_FOLDER}\n")
            input("Press Enter to continue...")

            match = re.search(r"\[\d+-(\d+)\]", LIBRARY_NAME_SPEC)
            used_reps = int(match.group(1))
            library = load_marker_library(LIBRARY_FOLDER, LIBRARY_NAME_SPEC)
            results = analyse_unknown(unk_df, library, UNKNOWN_SAMPLES)
            export_results(results, library, pep_df, RESULTS_FOLDER, TOTAL_REPLICATES, used_reps)

        elif analysis == "hierarchy":
            print("You chose two-step taxonomy based analysis - unknown samples are first analysed against a higher"
                  " taxonomy cluster and based on the result analysed against individual species.\n"
                  "-----Check list:-----\n"
                  "I have defined unknown sample names\n"
                  f"I have available both taxon-level and species-level marker libraries in {LIBRARY_FOLDER}\n"
                  f"I have defined taxonomy of species:\n{taxonomy_groups}\n")
            input("Press Enter to continue...")

            match = re.search(r"\[\d+-(\d+)\]", LIBRARY_NAME_TAX)
            used_reps = int(match.group(1))

            library_tax = load_marker_library(LIBRARY_FOLDER, LIBRARY_NAME_TAX)
            results, taxonomy = analyse_taxonomy(unk_df, library_tax, UNKNOWN_SAMPLES)
            library_spec = load_marker_library(LIBRARY_FOLDER, LIBRARY_NAME_SPEC)
            final_results = analyse_species(taxonomy, library_spec, taxonomy_groups, unk_df, UNKNOWN_SAMPLES)
            export_taxonomy(results, final_results, library_tax, library_spec, pep_df, RESULTS_FOLDER,
                            TOTAL_REPLICATES,
                            ref_replicates=used_reps)
            print(f"Hierarchy analysis results are now in {RESULTS_FOLDER}.")
        else:
            print("Invalid input!")
            sys.exit(-1)
    else:
        print("Invalid input!")
        sys.exit(-1)
