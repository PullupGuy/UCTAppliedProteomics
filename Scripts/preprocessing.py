"""
This module provides several functions for data preprocessing

adjust_peptides creates PeptideKey column so that it takes modifications into account and removes duplicate peptides
 based on SampleID and unique identifier (Sequence or PeptideKey) and also removes redundant information columns;
 moreover, low-intensity peptides might be removed as well as all falsely annotated peptides

remove_protein_duplicates removes duplicate proteins based on SampleID and protein identifier, removes redundant
 information columns and optionally cuts off low-intensity proteins
"""
import sys

def adjust_peptides(peptides, duplicates=False, cut_int=None):
    """
    analyses peptide dataframe and adjust it to become more opt for further analysis
    PeptideKey column is added so that modifications are also taken into account
    all duplicate peptides are removed and columns with useless information are dropped
    all falsely annotated peptides are also removed
    optionally low-intensity peptides might be removed
    :param peptides: dataframe of all parsed peptides
    :param duplicates: if True, only remove duplicate peptides based on peptide Sequence - ignores modifications
    :param cut_int: intensity quantile based on which low-intensity peptides are removed
    :return: dataframe of adjusted peptides - without duplicates etc.
    """
    print("\nPreprocessing peptides...")
    # check if SampleID is present
    if "SampleID" not in peptides.columns:
        print("!Error: SampleID column not found, unable to process!")
        sys.exit(-1)

    # make a copy of peptides
    adj_pep = peptides.copy()
    initial_count = len(adj_pep)

    adj_pep = adj_pep[adj_pep["OK"] == True]
    final_count = len(adj_pep)
    print(f"Removed {initial_count-final_count} falsely annotated peptides.")

    initial_count = len(adj_pep)
    if not duplicates:
        # create a new column PeptideKey as a combination of Modifications and Sequence
        adj_pep["Modifications"] = adj_pep["Modifications"].fillna("Unmodified")
        adj_pep["PeptideKey"] = adj_pep["Sequence"] + " | " + adj_pep["Modifications"]

        # remove all duplicate values based on PeptideKey and SampleID column
        adj_pep = adj_pep.drop_duplicates(subset=["PeptideKey", "SampleID"])
    else:
        adj_pep = adj_pep.drop_duplicates(subset=["Sequence", "SampleID"])

    cols_to_remove = ["User Flag", "Ratio [%]", "Peptide", "Cmpd.", "Row", "#Cmpds."]
    adj_pep = adj_pep.drop(columns=cols_to_remove, errors="ignore")
    print("Dropped columns with useless information.")

    final_count = len(adj_pep)
    print(f"Removed {initial_count-final_count} duplicate peptides.")

    if cut_int:
        limit = adj_pep["Int."].quantile(cut_int)
        adj_pep = adj_pep[adj_pep["Int."] >= limit]
        print(f"Removed low-intensity peptides (based on {cut_int} quantile).")
        final_count = len(adj_pep)

    print(f"Remaining unique peptide entries: {final_count}")
    print("Preprocessing completed!")

    return adj_pep

def remove_protein_duplicates(proteins, cut_int=None):
    """
    analyses peptide dataframe and adjust it to become more opt for further analysis
    all duplicate proteins are removed, columns with redundant information are dropped
    optionally low-intensity proteins might be removed
    :param proteins: dataframe of all parsed proteins
    :param cut_int: intensity quantile based on which low-intensity proteins are removed
    :return: dataframe of adjusted peptides - without duplicates etc.
    """
    # check if SampleID is present
    if "SampleID" not in proteins.columns:
        print("!Error: SampleID column not found, unable to process!")
        return proteins

    # make a copy of protein dataframe
    adj_prot = proteins.copy()

    # remove all duplicate values based on Protein and SampleID column
    initial_count = len(adj_prot)
    adj_prot = adj_prot.drop_duplicates(subset=["Protein", "SampleID"])

    adj_prot = adj_prot.drop(columns=["User Flag", "Ratio [%]", "Peptide", "Cmpd.", "Row"], ignore_errors = True)
    print("Dropped columns with useless information.")
    final_count = len(adj_prot)

    print(f"Removed {initial_count-final_count} duplicate proteins.")

    if cut_int:
        limit = adj_prot["Int."].quantile(cut_int)
        adj_pep = adj_prot[adj_prot["Int."] >= limit]
        print(f"Removed low-intensity peptides (based on {cut_int} quantile).")
        final_count = len(adj_pep)

    print(f"Remaining unique protein entries: {final_count}\n")

    return adj_prot