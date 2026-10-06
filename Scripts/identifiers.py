"""
This module provides means to extract Group (Species) name from SampleID:

get_group_names - extracts group from SampleID based on these conditions
if SampleID contains "_", split it and take the first part of the split
if SampleID contains "_" 2 or more times, split it and take the second part of the split
if SampleID does not contain "_", ignore the last letter
allows user to optionally merge several groups into a single new Group and/or assign
higher taxonomy to samples
"""

import sys
import numpy as np

def get_group_names(peptides, taxonomy_map=None, merge_map=None):
    """
    extracts Group name from SampleID:
    if SampleID contains "_", split it and take the first part
    if SampleID contains "_" 2 or more times, split it and take the second part
    if SampleID does not contain "_", ignore the last two letters
    option to add pre-defined higher taxonomy grouping and/or to merge several groups into one
    :param peptides: dataframe of all available peptides
    :param taxonomy_map: dict, optional map for joining Groups together according to taxonomy
    :param merge_map: dict, optional map for joining Groups together as desired
    :return: dataframe with Group and/or Taxonomy column added
    """
    print("Assigning biological groups to analysed data...")
    if "SampleID" not in peptides:
        print("Error: SampleID not found in peptides dataframe!")
        sys.exit(-1)

    peptide_df = peptides.copy()

    # convert SampleID into string
    sample_ids = peptide_df["SampleID"].astype(str)

    # count how made underscores does sample id contain
    underscore_count = sample_ids.str.count("_")

    # conditions
    one_underscore = underscore_count == 1
    more_underscore = underscore_count >= 2

    # choices for how to extract Group
    none_choice = sample_ids.str[:-1]   # ignore the last letter
    one_choice = sample_ids.str.split("_").str[0]   # split at "_" and take the first part of the split
    more_choice = sample_ids.str.split("_").str[1]  # split at "_" and take the second part of the split

    # create Group column based on conditions and choices, while default is to ignore the last letter
    peptide_df["Group"] = np.select([one_underscore, more_underscore], [one_choice, more_choice],
                                    default=none_choice)

    peptide_df["Group"] = peptide_df["Group"].replace("cow", "cattle")

    # create a Taxonomy column from Groups based on taxonomy_map
    if taxonomy_map is not None:
        # replace() finds dictionary keys in Group column and replaces them with dictionary values
        peptide_df["Taxonomy"] = peptide_df["Group"].replace(taxonomy_map)

    # merge Groups according to merge_map
    if merge_map is not None:
        # replace() finds dictionary keys in Group column and replaces them with dictionary values
        peptide_df["Group"] = peptide_df["Group"].replace(merge_map)
    print("DONE! Biological groups and taxonomy successfully added to your data")
    return peptide_df