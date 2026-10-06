"""
This module provides the means to analyse data and created a reference marker
library of reference samples, with option to set a minimum number of replicates required
for marker to be counted (robust)

reference_markers filters peptides based on Sequence; allows minimum replicates required
for marker to be counted together for target and background sets; this represents historic simple evaluation

reference_markers_robust allows the user to specify minimum number of replicates separately for
target group (which markers are being searched for) and background groups (against which we are filtering)
this allows the user to create significantly more robust markers that are present only in a single group
additionally, you may find markers on species-level or higher taxonomy-level
"""

def reference_markers(peptides, reference_groups, min_replicates = 1):
    """
    Analyses all found peptides and filters them according to Sequence and Group
    all peptides for each sample are found and then compared to all other samples
    what is left are characteristic peptides for the specific sample
    :param peptides: dataframe of all discovered peptides
    :param reference_groups: list of reference samples (Group)
    :param min_replicates: minimum number of replicates in which peptide must be present to provide robust markers
    :return: a dictionary of unique markers for all reference samples
    """
    # check is Group is present and if not exit the function
    if "Group" not in peptides.columns:
        print("No Group column in peptide dataframe, unable to perform analysis.")
        return {}   # empty dictionary

    # dictionary where all peptide sets are stored
    peptide_sets = {}

    # find all peptides associated with reference sample and add them to peptide_sets
    for sample in reference_groups:
        sub_ref = peptides[peptides["Group"] == sample]
        sample_ids = sub_ref["SampleID"].unique()
        num_replicates = len(sample_ids)

        # number of technical replicates cannot be lower than minimum required replicates
        if num_replicates < min_replicates:
            print(
                f"Error for reference group {sample}: fewer replicates ({num_replicates}) than"
                f" required ({min_replicates}), unable to perform analysis.")
            peptide_sets[sample] = set()
            continue

        # for every Sequence count in how many unique SampleIDs is present and name the column num_samples_with_marker
        presence_counts = (sub_ref.groupby('Sequence')['SampleID'].nunique()
                           .reset_index(name='num_samples_with_marker'))

        # filter peptides that appear in at least min_replicates samples
        robust_peptides = presence_counts[presence_counts['num_samples_with_marker'] >= min_replicates]['Sequence']

        # append found robust peptides into peptide dictionary
        current_peptides = set(robust_peptides)
        peptide_sets[sample] = current_peptides

    # dictionary where all unique markers are stored
    marker_sets = {}

    # filter unique markers for every reference sample
    for target_sample in reference_groups:
        target_set = peptide_sets[target_sample]

        # create a set to which all other species will be appended
        background_set = set()

        # target_set is compared to all other reference species which are aggregated in background_set
        # .update() is ideal, as already present peptides are not added again
        for other_sample in reference_groups:
            if other_sample != target_sample:
                background_set.update(peptide_sets[other_sample])

        # markers are the difference between target_set and background_set
        # peptides that are present ONLY in target_set are saved
        unique_markers = target_set - background_set

        # unique_markers for target_sample are added to marker library
        marker_sets[target_sample] = unique_markers

        # if there are some unique markers, create a list a show how many markers there are
        if len(unique_markers) > 0:
            marker_list = list(unique_markers)
            print(f"Found {len(marker_list)} characteristic peptides for sample {target_sample}.")
        else:
            print(f"Found no characteristic peptides for sample {target_sample}.")
    print("Finished constructing peptide marker library\n")
    return marker_sets

def reference_markers_robust(peptides, reference_groups, taxonomy_groups=None, rep_in = 1, rep_out = 1):
    """
    filters peptides (PeptideKey) among groups to acquire group-specific markers
    allows the user to specify minimum number of replicates separately for
    target group (which markers are being searched for) and background groups (against which we are filtering)
    this allows the user to create significantly more robust markers that are present only in a single group
    additionally, you may find markers on species-level or higher taxonomy-level
    :param peptides: dataframe of all discovered peptides
    :param reference_groups: list of reference samples (Groups)
    :param taxonomy_groups: list of pre-defined higher taxonomy groups
    :param rep_in: minimum number of replicates of target group to provide robust markers
    :param rep_out: minimum number of replicates of background group to provide robust markers
    :return: a dictionary of unique markers for all reference samples/taxonomy groups
    """
    print("\nStarting to construct marker library...")

    if taxonomy_groups:
        # check is Group is present and if not exit the function
        if "Group" not in peptides.columns:
            print("No Group column in peptide dataframe, unable to perform analysis.")
            return {}   # empty dictionary

        # check is Group is present and if not exit the function
        if "Taxonomy" not in peptides.columns:
            print("No Taxonomy column in peptide dataframe, unable to perform analysis.")
            return {}   # empty dictionary

        cols_to_str = ["Taxonomy", "SampleID", "PeptideKey", "Group"]
        for col in cols_to_str:
            if col in peptides.columns:
                peptides[col] = peptides[col].astype(str)

        # dictionary where all peptide sets are stored
        peptide_sets_in = {}
        peptide_sets_out = {}

        # find all peptides associated with reference sample and add them to peptide_sets
        for sample in list(taxonomy_groups):
            sub_ref = peptides[peptides["Taxonomy"] == sample]

            # for every PeptideKey count in how many unique SampleIDs is present and name the column num_samples_with_marker
            presence_counts = (sub_ref.groupby('PeptideKey')['SampleID'].nunique()
                                .reset_index(name='num_samples_with_marker'))

            # filter peptides that appear in at least rep_in samples
            robust_peptides_in = presence_counts[presence_counts['num_samples_with_marker'] >= rep_in]['PeptideKey']

            # append robust peptides to the peptide_in dictionary
            peptide_sets_in[sample] = set(robust_peptides_in)

            # filter peptides that appear in at least rep_out samples
            robust_peptides_out = presence_counts[presence_counts['num_samples_with_marker'] >= rep_out]['PeptideKey']

            # append robust peptides to the peptide_out dictionary
            peptide_sets_out[sample] = set(robust_peptides_out)

        if peptide_sets_in and peptide_sets_out:
            print("Successfully assessed robustness of potential markers...")
        else:
            print("Error: Unable to assess robustness of potential markers! Quitting...")
            return {}

        # dictionary where all unique markers are stored
        marker_sets = {}

        # filter unique markers for every reference sample
        for target_sample in taxonomy_groups:
            target_set = peptide_sets_in[target_sample]

            # create a set into which all other species will be appended
            background_set = set()

            # target_set is compared to all other reference species aggregated in background_set
            # .update() is ideal, as already present peptides are not added again
            for other_sample in taxonomy_groups:
                if other_sample != target_sample:
                    background_set.update(peptide_sets_out[other_sample])

            # markers are the difference between target_set and background_set
            # peptides that are present ONLY in target_set are saved
            unique_markers = target_set - background_set

            # unique_markers for target_sample are added to marker library dictionary
            marker_sets[target_sample] = unique_markers

            # if there are some unique markers, create a list a show how many markers there are
            if len(unique_markers) > 0:
                marker_list = list(unique_markers)
                print(f"Found {len(marker_list)} characteristic peptides for taxonomic cluster {target_sample}.")
            else:
                print(f"Found no characteristic peptides for taxonomic cluster {target_sample}.")

        print("Successfully finished constructing peptide marker library based on taxonomic similarity.")
        return marker_sets

    else:
        cols_to_str = ["Group", "SampleID", "PeptideKey"]
        for col in cols_to_str:
            if col in peptides.columns:
                peptides[col] = peptides[col].astype(str)

        # dictionary where all peptide sets are stored
        peptide_sets_in = {}
        peptide_sets_out = {}

        # find all peptides associated with reference sample and add them to peptide_sets
        for sample in list(reference_groups):
            sub_ref = peptides[peptides["Group"] == sample]

            # for every PeptideKey count in how many unique SampleIDs is present and name the column num_samples_with_marker
            presence_counts = (sub_ref.groupby('PeptideKey')['SampleID'].nunique()
                                .reset_index(name='num_samples_with_marker'))

            # filter peptides that appear in at least rep_in samples
            robust_peptides_in = presence_counts[presence_counts['num_samples_with_marker'] >= rep_in]['PeptideKey']

            # append robust peptides to the peptide_in dictionary
            peptide_sets_in[sample] = set(robust_peptides_in)

            # filter peptides that appear in at least rep_out samples
            robust_peptides_out = presence_counts[presence_counts['num_samples_with_marker'] >= rep_out]['PeptideKey']

            # append robust peptides to the peptide_out dictionary
            peptide_sets_out[sample] = set(robust_peptides_out)

        if peptide_sets_in and peptide_sets_out:
            print("Successfully assessed robustness of potential markers...")
        else:
            print("Error: Unable to assess robustness of potential markers! Quitting...")

        # dictionary where all unique markers are stored
        marker_sets = {}

        # filter unique markers for every reference sample
        for target_sample in reference_groups:
            target_set = peptide_sets_in[target_sample]

            # create a set into which all other species will be appended
            background_set = set()

            # target_set is compared to all other reference species aggregated in background_set
            # .update() is ideal, as already present peptides are not added again
            for other_sample in reference_groups:
                if other_sample != target_sample:
                    background_set.update(peptide_sets_out[other_sample])

            # markers are the difference between target_set and background_set
            # peptides that are present ONLY in target_set are saved
            unique_markers = target_set - background_set

            # unique_markers for target_sample are added to marker library dictionary
            marker_sets[target_sample] = unique_markers

            # if there are some unique markers, create a list a show how many markers there are
            if len(unique_markers) > 0:
                marker_list = list(unique_markers)
                print(f"Found {len(marker_list)} characteristic peptides for sample {target_sample}.")
            else:
                print(f"Found no characteristic peptides for sample {target_sample}.")

        print("Successfully finished constructing peptide marker library.")
        return marker_sets