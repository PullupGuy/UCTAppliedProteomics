"""
This module provides means to load all kinds of proteomic data into dataframe structure

load_peptides_proteins extracts Excel file into separate protein and peptide dataframes, providing
the option to ignore certain sampleIDs

load_peptides extracts Excel file into a single peptides dataframe, providing
the option to ignore certain sampleIDs

load_marker_library_csv reads CSV file and creates a dictionary of all PeptideKeys/Sequence with
keys as Species
"""
import sys

import pandas as pd
from pathlib import Path

def load_peptides_proteins(folder, ignored=None):
    """
    Loads all .xlsx files in folder and returns every sheet as a pandas dataframe
    separates every file into Protein and Peptide dataframes based on _pep suffix
    Automatically parses sheet name to get sample id information
    :param folder: directory of proteomic data in .xlsx format
    :param ignored: specific sampleID to ignore when loading the data
    :return: all recorded proteins and peptides as a separate dataframes
    """
    # set lists for proteins and peptides
    protein_dfs = []
    peptide_dfs = []

    # find all Excel files in specified folder
    files = list(Path(folder).glob("*.xlsx"))

    if len(files) == 0:
        print("!Error: no .xlsx files found, unable to continue!")
        return pd.DataFrame()

    print(f"Found {len(files)} .xlsx files. Starting to parse...")

    for file_path in files:
        # extract file name from files Path objects
        filename = file_path.name
        print(f"Parsing {filename}...")

        try:
            # read every list in Excel file (arg sheet_name=None) with openpyxl
            # result is a dictionary of dataframes
            xls_dict = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")

            for sheet_name, df in xls_dict.items():
                # if dataframe is empty skip to the next sheet_name
                if df.empty:
                    continue

                # set is_peptide to True or False
                is_peptide = "pep" in sheet_name

                # extract sample ID from sheet_name - remove _pep and any left whitespaces around sample ID
                sample_id = sheet_name.replace("_pep", "").strip()

                # if ignored group is in sample_id continue to the next sheet
                if ignored and str(ignored) in sample_id:
                    continue

                # create new dataframe columns
                df["SampleID"] = sample_id

                # based on is_peptide, decide whether data are peptides or proteins
                if is_peptide:
                    peptide_dfs.append(df)
                else:
                    protein_dfs.append(df)

            print(f"Parsed {len(xls_dict)} sheets from {filename}.")

        except Exception as e:
            print(f"!Error parsing {filename}: {e}!")

    print("Joining data...")

    # concatenate all proteins/peptides into a single dataframe while resetting indexes
    # if proteins_dfs is empty, concatenate empty dataframes -> prevents errors
    all_proteins = pd.concat(protein_dfs, ignore_index=True) if protein_dfs else pd.DataFrame()
    all_peptides = pd.concat(peptide_dfs, ignore_index=True) if peptide_dfs else pd.DataFrame()

    print(f"DONE. Total proteins recorded: {len(all_proteins)}")
    print(f"DONE. Total peptides recorded: {len(all_peptides)}\n")

    return all_proteins, all_peptides

def load_peptides(folder, ignored=None):
    """
    Loads all .xlsx files in folder and returns every sheet as a pandas dataframe.
    Assumes ALL data are peptides.
    Automatically parses sheet name to get sample id information
    :param folder: directory of proteomic data in .xlsx format
    :param ignored: specific sampleID to ignore when loading the data
    :return: single dataframe containing all recorded peptides
    """
    peptide_dfs = []

    # create a list of all .xlsx files in the folder
    files = list(Path(folder).glob("*.xlsx"))

    if len(files) == 0:
        print("\n!Error: no .xlsx files found, unable to continue with the analysis!")
        sys.exit(-1)

    print(f"\nFound {len(files)} .xlsx files. Starting to parse...")

    for file_path in files:
        # extract filename from Path object
        filename = file_path.name
        print(f"Parsing {filename}...")

        try:
            # read EVERY Excel sheet in file_path with openpyxl engine
            xls_dict = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")

            for sheet_name, df in xls_dict.items():
                # if there are no data skip to the next sheet
                if df.empty:
                    continue

                # extract sample_id by stripping all white spaces from sheet_name
                sample_id = sheet_name.strip()

                # if ignored group is in sample_id continue to the next sheet
                if ignored and str(ignored) in sample_id:
                    continue

                # add all required columns so that no information is lost
                df["SampleID"] = sample_id

                peptide_dfs.append(df)

            print(f"Parsed {len(xls_dict)} sheets from {filename}.")

        except Exception as e:
            print(f"!Error parsing {filename}: {e}")

    print("Joining data...")

    # concatenate all data into a single dataframe and reset indexes
    if peptide_dfs:
        all_peptides = pd.concat(peptide_dfs, ignore_index=True)

    else:
        all_peptides = pd.DataFrame()

    print(f"DONE! Total peptides recorded: {len(all_peptides)}")

    return all_peptides

def load_marker_library(folder, file, key=True):
    """
    Scans folder and reads CSV file
    creates a dictionary of all found peptide markers where key is species
    ignores Protein column as it is unnecessary
    :param folder: folder where desired CSV file is located
    :param file: desired CSV file
    :param key: whether to use PeptideKey or Sequence as a unique identifier
    :return: dictionary of found peptide markers which species as keys
    """
    # create file_path object
    file_path = Path(folder) / Path(file)

    # read cvs file and save into dataframe
    df = pd.read_csv(file_path)
    if key:
        # set which columns are to be loaded
        required_columns = ['Species', 'Peptide Marker Key']

        # check whether all required columns are present
        if not all(col in df.columns for col in required_columns):
            print(f"!CSV does not contain all required columns: {required_columns}!")
            return {}   # empty dictionary

        # create values (set of Peptide Marker Keys) and add Species as keys
        marker_sets = df.groupby("Species")["Peptide Marker Key"].agg(set).to_dict()
        print(f"\nLoaded marker library from {file}")
    else:
        required_columns = ['Species', 'Peptide Sequence']

        # check whether all required columns are present
        if not all(col in df.columns for col in required_columns):
            print(f"!CSV does not contain all required columns: {required_columns}!")
            return {}  # empty dictionary

        # create values (set of Peptide Marker Keys) and add Species as keys
        marker_sets = df.groupby("Species")["Peptide Sequence"].agg(set).to_dict()
    return marker_sets