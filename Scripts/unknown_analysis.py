"""
This module provides the means for user to analyse unknown samples and try matching
them to previously created reference marker library
allows to specify how made replicates are required for unknown sample identifier to be counted as robust

analyse_unknown analyses all unknown samples by comparing their unique peptide identifiers (PeptideKey/Sequence) to
a predefined peptide marker library. For every reference species, a match_score and confidence are computed based
on the similarity between the unknown sample and the reference species

analyse_taxonomy analyses all unknown samples by comparing their unique peptide identifiers (PeptideKey) to
a predefined peptide marker library. For every reference taxonomy group, a match_score and confidence are computed based
on the similarity between the unknown sample and the taxonomy group. Part of the analysis is identification of
taxonomy group with the highest confidence.

analyse_species analyses all unknown samples by comparing their unique peptide identifiers (PeptideKey) to
a predefined peptide marker library. Every sample is compared only to reference species within previously identified
taxonomy group with the highest confidence, a match_score and confidence are computed based on the similarity
between the unknown sample and the select reference species.
"""

import pandas as pd

def analyse_unknown(peptides, markers, unknown_groups, replicates=1, keys=True):
    """
    analyses all found peptides and filters them according to PeptideKey/Sequence and Group
    the only considered Groups are those found in unknown_groups
    peptides are associated with unknown samples and then compared to marker library
    match_score is calculated based on the number of proteins found both in the library and in the unknown sample
    what is left are unique peptides for the unknown sample
    presence of peptides in different technical replicates is also taken into account with default setting being 1
    :param peptides: dataframe of all discovered and adjusted peptides
    :param markers: library of all identified marker peptides
    :param unknown_groups: list of unknown groups to be analysed (IDs)
    :param replicates: minimum number of replicates in which peptide must be present to provide robust identifiers
    :param keys: whether to use PeptideKey (True) or Sequence (False) as an identifier
    :return: dataframe containing important information about unknown samples
    """
    print("\nBeginning analysis of unknown samples...")
    # if there are no unknown samples
    if len(unknown_groups) == 0:
        print("Error: no unknown samples found for analysis!")
        return pd.DataFrame()  # empty dataframe

    # take only peptides in unknown_groups
    unknown_df = peptides[peptides["Group"].isin(unknown_groups)]

    # empty list of results where data will be appended
    results = []

    # loop over every unknown sample and extract the number of replicates based on unique SampleIDs
    for un_sample in unknown_groups:
        print(f"Analyzing {un_sample}...")
        sub_unknown = unknown_df[unknown_df["Group"] == un_sample]
        sample_ids = sub_unknown["SampleID"].unique()
        num_replicates = len(sample_ids)

        # number of technical replicates cannot be lower than minimum required replicates
        if num_replicates < replicates:
            print(f"Error for {un_sample} group: fewer replicates than required {replicates},"
                  f" analysis cannot be performed.")

            # for every reference species (extracted from markers keys)
            for ref in markers.keys():
                # calculate the number of all reference markers
                total_markers = len(markers.get(ref, set()))

                if keys:
                    # append 0% match scores to unknown sample for all references
                    results.append({
                        'Unknown Sample': un_sample,
                        "Reference Group": ref,
                        'Markers Found': 0,
                        "Confidence [%]": 0.0,
                        "Reference Markers": total_markers,
                        'Match Score [%]': 0.0,
                        'Found Marker Keys': [],
                        'Replicates': num_replicates
                    })
                else:
                    # append 0% match scores to unknown sample for all references
                    results.append({
                        'Unknown Sample': un_sample,
                        "Reference Group": ref,
                        'Markers Found': 0,
                        "Confidence [%]": 0.0,
                        "Reference Markers": total_markers,
                        'Match Score [%]': 0.0,
                        'Found Marker Sequence': [],
                        'Replicates': num_replicates
                    })
            continue

        if keys:
            # for every PeptideKey count in how many unique SampleIDs is present and name the column num_sample_with_marker
            presence_counts = (sub_unknown.groupby('PeptideKey')['SampleID'].nunique()
                               .reset_index(name='num_samples_with_marker'))

            # set index so that column num_samples_with_marker is accessible with PeptideKey
            presence_counts = presence_counts.set_index('PeptideKey')
        else:
            # for every PeptideKey count in how many unique SampleIDs is present and name the column num_sample_with_marker
            presence_counts = (sub_unknown.groupby('Sequence')['SampleID'].nunique()
                               .reset_index(name='num_samples_with_marker'))

            # set index so that column num_samples_with_marker is accessible with PeptideKey
            presence_counts = presence_counts.set_index('Sequence')

        total_marker_count = 0
        temp_results = []

        # loop thorough every set of markers in markers
        for ref, marker_set in markers.items():

            # convert a set to a list that is iterable
            marker_list = list(marker_set)

            # calculate the number of all markers within reference set
            total_markers = len(marker_list)

            # if there are no total markers, append a special case of nulls
            if total_markers == 0:
                if keys:
                    temp_results.append({
                        "Unknown Sample": un_sample,
                        "Reference Group": ref,
                        "Markers Found": 0,
                        "Reference Markers": total_markers,
                        "Match Score [%]": 0.0,
                        "Found Marker Keys": [],
                        "Replicates": num_replicates})
                else:
                    temp_results.append({
                        "Unknown Sample": un_sample,
                        "Reference Group": ref,
                        "Markers Found": 0,
                        "Reference Markers": total_markers,
                        "Match Score [%]": 0.0,
                        "Found Marker Sequence": [],
                        "Replicates": num_replicates})
                continue

            # create dataframe with rows corresponding to marker_list
            # for markers are not present in presence_counts, assign NaN values
            marker_presence = presence_counts.reindex(marker_list)

            # boolean mask that checks if markers are defined AND are longer than MIN_REPLICATES
            # result is markers that are frequent enough
            robust_mask = (marker_presence["num_samples_with_marker"].notna() &
                           (marker_presence["num_samples_with_marker"] >= replicates))

            # count the number of True values in robust_mask -> number of present markers
            count_markers = robust_mask.sum()

            # create a list of markers which fulfill robust_mask and are in marker_list
            found_markers = [m for m in marker_list if robust_mask.get(m, False)]

            # calculate the match score (as percentage)
            match_score = (count_markers / total_markers) * 100
            total_marker_count += count_markers

            if keys:
                # append computed information into temp_results list
                temp_results.append({
                    "Unknown Sample": un_sample,
                    "Reference Group": ref,
                    "Markers Found": count_markers,
                    "Reference Markers": total_markers,
                    "Match Score [%]": match_score,
                    "Found Marker Keys": found_markers,
                    "Replicates": num_replicates})
            else:
                # append computed information into temp_results list
                temp_results.append({
                    "Unknown Sample": un_sample,
                    "Reference Group": ref,
                    "Markers Found": count_markers,
                    "Reference Markers": total_markers,
                    "Match Score [%]": match_score,
                    "Found Marker Sequence": found_markers,
                    "Replicates": num_replicates})

        # Calculate confidence and append to actual results
        for res in temp_results:
            if total_marker_count > 0:
                confidence = (res["Markers Found"] / total_marker_count) * 100
            else:
                confidence = 0.0

            res["Confidence [%]"] = confidence
            results.append(res)

    print("Completed unknown sample analysis!")

    # force dataframe type on a list of results
    return pd.DataFrame(results)

def analyse_taxonomy(peptides, markers, unknown_groups, replicates=1):
    """
    analyses all found peptides and filters them according to PeptideKey and Group
    the only considered Groups are those found in unknown_groups
    peptides are associated with unknown samples and then compared to marker library
    match_score is calculated based on the number of proteins found both in the library and in the unknown sample
    what is left are unique peptides for the unknown sample
    presence of peptides in different technical replicates is also taken into account with default setting being 1
    :param peptides: dataframe of all discovered and adjusted peptides
    :param markers: library of all identified marker peptides
    :param unknown_groups: list of unknown groups to be analysed (ID)
    :param replicates: minimum number of replicates in which peptide must be present to provide robust markers
    :return: dataframe containing important information about unknown samples and
             a dictionary containing taxonomy cluster with the best match
    """
    print("\nBeginning analysis of unknown samples against higher taxonomy groups...")
    # if there are no unknown samples
    if len(unknown_groups) == 0:
        print("Error: no unknown samples found for analysis!")
        return pd.DataFrame()  # empty dataframe

    # empty list of results where data will be appended
    results = []
    best_clusters = {}

    # loop over every unknown sample and extract the number of replicates based on unique SampleIDs
    for un_sample in unknown_groups:
        print(f"Analyzing {un_sample}...")
        sub_unknown = peptides[peptides["Group"] == un_sample]
        sample_ids = sub_unknown["SampleID"].unique()
        num_replicates = len(sample_ids)

        # number of technical replicates cannot be lower than minimum required replicates
        if num_replicates < replicates:
            print(f"Error for {un_sample} group: fewer replicates than required {replicates},"
                    f" analysis cannot be performed.")

            # for every reference species (extracted from markers keys)
            for ref in markers.keys():
                # calculate the number of all reference markers
                total_markers = len(markers.get(ref, set()))

                # append 0% match scores to unknown sample for all references
                results.append({
                    'Unknown Sample': un_sample,
                    "Taxonomy Group": ref,
                    'Markers Found': 0,
                    "Confidence [%]": 0.0,
                    "Reference Markers": total_markers,
                    'Match Score [%]': 0.0,
                    'Found Marker Keys': [],
                    'Replicates': num_replicates
                })
            continue

        # for every PeptideKey count in how many unique SampleIDs is present and name the column num_sample_with_marker
        presence_counts = (sub_unknown.groupby('PeptideKey')['SampleID'].nunique()
                           .reset_index(name='num_samples_with_marker'))

        # set index so that column num_samples_with_marker is accessible with PeptideKey
        presence_counts = presence_counts.set_index('PeptideKey')

        best_match = 0
        total_marker_count = 0
        temp_results = []
        # loop thorough every set of markers in markers
        for ref, marker_set in markers.items():

            # convert a set to a list that is iterable
            marker_list = list(marker_set)

            # calculate the number of all markers within reference set
            total_markers = len(marker_list)

            # if there are no total markers, append a special case of nulls
            if total_markers == 0:
                temp_results.append({
                    "Unknown Sample": un_sample,
                    "Taxonomy Group": ref,
                    "Markers Found": 0,
                    "Reference Markers": total_markers,
                    "Match Score [%]": 0.0,
                    "Found Marker Keys": [],
                    "Replicates": num_replicates,
                    "Best Cluster": None})
                continue

            # create dataframe with rows corresponding to marker_list
            # for markers are not present in presence_counts, assign NaN values
            marker_presence = presence_counts.reindex(marker_list)

            # boolean mask that checks if markers are defined AND are longer than MIN_REPLICATES
            # result is markers that are frequent enough
            robust_mask = (marker_presence["num_samples_with_marker"].notna() &
                           (marker_presence["num_samples_with_marker"] >= replicates))

            # count the number of True values in robust_mask -> number of present markers
            count_markers = robust_mask.sum()

            # create a list of markers which fulfill robust_mask and are in marker_list
            found_markers = [m for m in marker_list if robust_mask.get(m, False)]

            # calculate the match score (as percentage)
            match_score = (count_markers / total_markers) * 100

            if count_markers > best_match:
                best_match = count_markers
                best_cluster = ref

            total_marker_count += count_markers

            # append computed information into results list
            temp_results.append({
                "Unknown Sample": un_sample,
                "Taxonomy Group": ref,
                "Markers Found": count_markers,
                "Reference Markers": total_markers,
                "Match Score [%]": match_score,
                "Found Marker Keys": found_markers,
                "Replicates": num_replicates})

        best_clusters.update({un_sample: best_cluster})

        for res in temp_results:
            if total_marker_count > 0:
                confidence = (res["Markers Found"] / total_marker_count) * 100
            else:
                confidence = 0.0

            res["Confidence [%]"] = confidence
            results.append(res)

    print("Completed taxon_level unknown sample analysis!")

    # force dataframe type on a list of results
    return pd.DataFrame(results), best_clusters

def analyse_species(cluster, markers, taxonomy, peptides, unknown_groups, replicates=1):
    """
    Follow-up unknown sample analysis which take as an input higher taxonomy information and limits the search
    to species within that taxonomy cluster and perform similar comparison to reference marker library
    :param cluster: dictionary of sample:best cluster information that limits the search
    :param markers: library of species-level markers
    :param taxonomy: dictionary of assigned taxonomy order (e.g. {"Birds": ["Hen", "Duck"]})
    :param peptides: dataframe of all unknown sample preprocessed peptides
    :param unknown_groups: list of unknown sample groups
    :param replicates: minimum number of replicates required for analysis
    :return: dataframe of analysis species-level results
    """
    print("\nBeginning analysis of unknown samples based on acquired taxonomy information...")

    if len(unknown_groups) == 0:
        print("Error: no unknown samples found for analysis!")
        return pd.DataFrame()

    if "PeptideKey" not in peptides.columns:
        print("Error: Peptide key column missing!")
        return pd.DataFrame()

    results = []

    # loop through all unknown samples and associated taxa
    for un_sample, taxon in cluster.items():
        print(f"Analysing {un_sample} for species within taxon {taxon}...")
        # get every candidate species from previous analysis for a concrete unknown sample
        candidate_species = taxonomy.get(taxon, [])

        if not candidate_species:
            print(f"Warning: No species found for taxon '{taxon}'. Skipping.")
            continue

        # filter peptides only for our chosen sample
        sub_unknown = peptides[peptides["Group"] == un_sample]
        sample_ids = sub_unknown["SampleID"].unique()
        num_replicates = len(sample_ids)

        # replicate control
        if num_replicates < replicates:
            print(f"{un_sample} group: fewer replicates than required ({replicates}), assigning null values.")

            # assign zero values if there are not enough replicates
            for species in candidate_species:
                total_markers = len(markers.get(species, set()))
                results.append({
                    'Unknown Sample': un_sample,
                    'Assigned Taxon': taxon,  # Save assigned taxon for future reference
                    "Reference Group": species,
                    'Markers Found': 0,
                    "Confidence [%]": 0.0,
                    "Reference Markers": total_markers,
                    'Match Score [%]': 0.0,
                    'Found Marker Keys': [],
                    'Replicates': num_replicates
                })
            continue

        # asses in how many samples are unique PeptideKeys present
        presence_counts = (sub_unknown.groupby('PeptideKey')['SampleID'].nunique()
                           .reset_index(name='num_samples_with_marker'))
        presence_counts = presence_counts.set_index('PeptideKey')

        total_marker_count = 0
        temp_results = []

        # the analysis itself - only through candidate species
        for species in candidate_species:

            # get all markers for select species
            marker_set = markers.get(species, set())
            marker_list = list(marker_set)
            total_markers = len(marker_list)

            if total_markers == 0:
                print(f"Found zero markers for {species}, assigning null values!")
                temp_results.append({
                    "Unknown Sample": un_sample,
                    'Assigned Taxon': taxon,
                    "Reference Group": species,
                    "Markers Found": 0,
                    "Reference Markers": total_markers,
                    "Match Score [%]": 0.0,
                    "Found Marker Keys": [],
                    "Replicates": num_replicates
                })
                continue

            marker_presence = presence_counts.reindex(marker_list)
            robust_mask = (marker_presence["num_samples_with_marker"].notna() &
                           (marker_presence["num_samples_with_marker"] >= replicates))

            count_markers = robust_mask.sum()
            found_markers = [m for m in marker_list if robust_mask.get(m, False)]
            match_score = (count_markers / total_markers) * 100

            total_marker_count += count_markers

            temp_results.append({
                "Unknown Sample": un_sample,
                'Assigned Taxon': taxon,
                "Reference Group": species,
                "Markers Found": count_markers,
                "Reference Markers": total_markers,
                "Match Score [%]": match_score,
                "Found Marker Keys": found_markers,
                "Replicates": num_replicates
            })

        for res in temp_results:
            if total_marker_count > 0:
                confidence = (res["Markers Found"] / total_marker_count) * 100
            else:
                confidence = 0.0

            res["Confidence [%]"] = confidence
            results.append(res)

    results_df = pd.DataFrame(results)
    print("Completed species-level unknown sample analysis!")

    return results_df
