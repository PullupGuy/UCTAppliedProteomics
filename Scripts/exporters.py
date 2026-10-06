"""
This module provides several exportation options for LC-MS/MS analysis

export_marker_library enables the user to conveniently export created marker library
to CSV file, from which the data can be easily loaded

export_results creates an Excel file with analysis results in user-friendly format containing several sheets:
        first sheet - list of assessment of each sample compared to every reference group
        second sheet - list of all peptide markers associated with each unknown sample
        third sheet - library of all identified reference peptide markers

export_taxonomy creates an Excel file with analysis results in user-friendly format containing several sheets:
    first sheet - list of percentual match of each sample to each taxonomy group
    second sheet - list of percentual match of each sample to species within the select taxonomy group
    third sheet - detailed information on all found taxonomy group markers for every unknown sample
    fourth sheet - detailed information on all found species-level markers for every unknown sample
    fifth sheet - library of all found taxonomy group markers
    sixth sheet - library of all found species-level group markers
"""

from pathlib import Path
import pandas as pd

def export_marker_library(output_folder, markers, peptides, ref_replicates, total_rep, keys=True, filename=None):
    """
    Export a dictionary of peptide marker identifiers (PeptideKey or Sequence) into CSV file
    :param output_folder: defines into which folder should the library be exported
    :param markers: dictionary of found reference markers
    :param ref_replicates: number of replicates required
    :param peptides: dataframe of all found peptides
    :param total_rep: total number of replicates available
    :param keys: whether to use PeptideKey (True) or Sequence (False) as an identifier
    :param filename: name of the file to be used
    """
    # create a folder if it doesn't already exist
    folder_path = Path(output_folder)
    if not folder_path.exists():
        folder_path.mkdir(parents=True)
        print(f"Created a folder for reference libraries: {folder_path}")

    if not filename:
        # set a filename
        filename = f"Marker_library[{total_rep}-{ref_replicates}]"

    else:
        filename = filename + f"[{total_rep}-{ref_replicates}]"

    # create a full path by joining folder and filename paths for CSV
    path_csv = folder_path / Path(f"{filename}.csv")
    if not len(markers) == 0:
        if keys:
            # create a 1:1 map of PeptideKey-Protein
            prot_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["Protein"]

            # using a list comprehension on all species and associated markers
            lib_data = [
                {'Species': sp, 'Peptide Marker Key': m}
                for sp, m_list in markers.items()
                for m in m_list]

            # force a dataframe data type on the list
            lib_df = pd.DataFrame(lib_data)

            # associate every Peptide Key Marker with Protein from prot_map
            lib_df["Protein"] = lib_df["Peptide Marker Key"].map(prot_map).fillna("Unknown Protein")

            # set the order of columns in the dataframe
            lib_df = lib_df[["Species", "Peptide Marker Key", "Protein"]]
        else:
            # create a 1:1 map of PeptideKey-Protein
            prot_map = peptides.drop_duplicates("Sequence").set_index("Sequence")["Protein"]

            # using a list comprehension on all species and associated markers
            lib_data = [
                {'Species': sp, 'Peptide Sequence': m}
                for sp, m_list in markers.items()
                for m in m_list]

            # force a dataframe data type on the list
            lib_df = pd.DataFrame(lib_data)

            # associate every Peptide Key Marker with Protein from prot_map
            lib_df["Protein"] = lib_df["Peptide Sequence"].map(prot_map).fillna("Unknown Protein")

            # set the order of columns in the dataframe
            lib_df = lib_df[["Species", "Peptide Sequence", "Protein"]]

        # export to csv
        lib_df.to_csv(path_csv, index=False)

        print(f"Saved library {filename} to CSV.\n")
    else:
        print("Unable to create reference marker library as no markers have been found.")

def export_results(analysed_samples, markers, peptides, output_folder, total_rep, ref_replicates,
                             un_replicates=1, keys=True):
    """
        Creates an Excel file with analysis results in user-friendly format containing several sheets:
        first sheet - list of assessment of each sample compared to every reference group
        second sheet - list of all peptide markers associated with each unknown sample
        third sheet - library of all identified reference peptide markers
        :param analysed_samples: dataframe of all analysed unknown samples
        :param markers: dictionary of found reference markers
        :param peptides: dataframe of all discovered and adjusted peptides
        :param output_folder: name of folder to be saved into
        :param ref_replicates: number of replicates used for peptide markers
        :param total_rep: total number of reference group replicates available
        :param un_replicates: number of unknown sample replicates used
        :param keys: whether to use PeptideKey (True) or Sequence (False) as an identifier
        :return: excel workbook with 3 sheets containing different data
        """
    print("\nInitiating data exportation to Excel ...")
    # create a folder if it doesn't already exist
    folder_path = Path(output_folder)
    if not folder_path.exists():
        folder_path.mkdir(parents=True)
        print(f"Created a folder for saved results: {folder_path}")

    # set filename
    filename = f"Analysis_report[{total_rep}-{ref_replicates}-{un_replicates}]"

    # create a full path by joining folder and filename paths
    full_path = folder_path / Path(f"{filename}.xlsx")

    # helper function for nice formatting in Excel
    def highlight_max_confidence(df):
        """
        Applies CSS styling to highlight the best-matching reference group for each unknown sample.

        This function is designed to be used with the Pandas Styler object via the .apply() method
        with axis=None. It creates a parallel DataFrame of the same shape as the input,
        filled with CSS strings, which Pandas then maps 1:1 to the exported Excel cells.

        :param df: dataframe containing the analysed unknown samples and their match scores
        :return: A DataFrame of the exact same shape as 'df', containing CSS style strings
        """

        # Create an empty DataFrame with the exact same dimensions (index and columns) as the input 'df'.
        # This acts as a blank canvas. Pandas will overlay this on the actual data during Excel export.
        styles = pd.DataFrame('', index=df.index, columns=df.columns)

        # Iterate through every unique unknown sample present in the dataset.
        # Since one unknown sample is compared against multiple reference groups,
        # it spans multiple rows. We need to find the best row for each sample independently.
        for sample in df['Unknown Sample'].unique():

            # Create a boolean mask (True/False series) to isolate only the rows
            # that belong to the current unknown sample in the loop.
            mask = df['Unknown Sample'] == sample

            # Safety check: Ensure that the 'Confidence [%]' column for this specific sample
            # actually contains numerical data. If all values are NaN (Not a Number),
            # the idxmax() function below would crash.
            if not df.loc[mask, 'Confidence [%]'].isna().all():
                # Find the row index where the 'Confidence [%]' is at its absolute maximum
                # for the current sample. (This identifies the most probable reference match).
                max_idx = df.loc[mask, 'Confidence [%]'].idxmax()

                # Apply the CSS styling string to the entire row (all columns, denoted by ':')
                # at the identified index in our parallel 'styles' DataFrame.
                # #c6efce is a standard light green hex color.
                styles.loc[max_idx, :] = 'background-color: #c6efce; font-weight: bold;'

        # Return the fully populated styling DataFrame back to the Pandas Styler engine
        return styles

    # ExcelWriter is used to create multiple sheets
    with pd.ExcelWriter(full_path, engine="openpyxl") as writer:

        analysed_samples = analysed_samples.rename(columns={
            "Reference Group": "Reference Species",
            "Found Marker Keys": "Marker Key",
            "Found Marker Sequence": "Marker Sequence"})

        # the first list - table of percentual match of unknown samples to reference library
        print("Creating the first Excel sheet of unknown sample assessment information...")
        if not (analysed_samples.empty and not markers):
            # remove IDs to make the table more readable
            if keys:
                summary_df_final = analysed_samples.drop(columns=["Marker Key"], errors="ignore")
            else:
                summary_df_final = analysed_samples.drop(columns=["Marker Sequence"], errors="ignore")

            cols = ["Unknown Sample", "Reference Species", "Markers Found", "Confidence [%]",
                    "Reference Markers", "Match Score [%]", "Replicates"]
            summary_df_final = summary_df_final[cols]

            # round the percentual match to first decimal digit
            summary_df_final["Match Score [%]"] = summary_df_final["Match Score [%]"].round(1)
            summary_df_final["Confidence [%]"] = summary_df_final["Confidence [%]"].round(1)

            # export to Excel sheet without any indexes
            styled_final = summary_df_final.style.apply(highlight_max_confidence, axis=None)
            styled_final.to_excel(writer, sheet_name="Unknown_sample_assessment", index=False)
        else:
            print("Unable to create Unknown sample assignment sheet because there are either no analysed samples "
                  "or no markers")

        # second sheet - table of concrete peptides associated with samples
        print("Creating second Excel sheet of detailed marker information...")
        if not analysed_samples.empty and ("Marker Key" or "Marker Sequence") in analysed_samples.columns:

            if keys:
                # write every element in the dataframe column to different rows in Excel
                detail_df = analysed_samples.explode("Marker Key")

                # create a Pandas Series of PeptideKey: Protein pairs based on peptides dataframe
                # drop_duplicates is necessary as two identical peptides can be found in 2 samples and are not unique
                prot_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["Protein"]
                mz_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["m/z meas."]
                charge_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["z"]

                # for every Marker Key find a corresponding protein and add it to Protein column
                detail_df["Protein"] = detail_df["Marker Key"].map(prot_map).fillna("Unknown Protein")
                detail_df["Charge"] = detail_df["Marker Key"].map(charge_map).fillna("NA")
                detail_df["m/z"] = detail_df["Marker Key"].map(mz_map).fillna("NA")

                detail_df = detail_df.dropna(subset=["Marker Key"])

                # define which columns are wanted in the Excel sheet
                cols_to_keep = ["Unknown Sample", "Reference Species", "Marker Key", "m/z", "Charge", "Protein"]
            else:
                # write every element in the dataframe column to different rows in Excel
                detail_df = analysed_samples.explode("Marker Sequence")

                # create a Pandas Series of PeptideKey: Protein pairs based on peptides dataframe
                # drop_duplicates is necessary as two identical peptides can be found in 2 samples and are not unique
                prot_map = peptides.drop_duplicates("Sequence").set_index("Sequence")["Protein"]
                mz_map = peptides.drop_duplicates("Sequence").set_index("Sequence")["m/z meas."]
                charge_map = peptides.drop_duplicates("Sequence").set_index("Sequence")["z"]

                # for every Marker Key find a corresponding protein and add it to Protein column
                detail_df["Protein"] = detail_df["Marker Sequence"].map(prot_map).fillna("Unknown Protein")
                detail_df["Charge"] = detail_df["Marker Sequence"].map(charge_map).fillna("NA")
                detail_df["m/z"] = detail_df["Marker Sequence"].map(mz_map).fillna("NA")

                detail_df = detail_df.dropna(subset=["Marker Sequence"])

                # define which columns are wanted in the Excel sheet
                cols_to_keep = ["Unknown Sample", "Reference Species", "Marker Sequence", "m/z", "Charge", "Protein"]

            # write desired columns and all associated data to Excel sheet
            detail_df[cols_to_keep].to_excel(writer, sheet_name='Markers_Detail', index=False)
        else:
            print("Unable to create second Excel sheet (Detailed Markers).")

        # third sheet - creation of reference library
        print("Creating third Excel sheet listing the whole reference marker library...")

        if not len(markers) == 0:
            if keys:
                # using a list comprehension on all species and associated markers
                lib_data = [
                    {'Species': sp, 'Peptide Marker Key': m}
                    for sp, m_list in markers.items()
                    for m in m_list]

                # force a dataframe data type on the list
                lib_df = pd.DataFrame(lib_data)

                # create Pandas Series maps based on the original peptides dataframe
                prot_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["Protein"]
                mz_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["m/z meas."]
                charge_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["z"]

                # map the values to the Reference Library dataframe
                lib_df["m/z"] = lib_df["Peptide Marker Key"].map(mz_map).fillna("NA")
                lib_df["Charge"] = lib_df["Peptide Marker Key"].map(charge_map).fillna("NA")
                lib_df["Protein"] = lib_df["Peptide Marker Key"].map(prot_map).fillna("Unknown Protein")

                # optionally reorder columns for a cleaner look in Excel
                lib_df = lib_df[["Species", "Peptide Marker Key", "m/z", "Charge", "Protein"]]
            else:
                # using a list comprehension on all species and associated markers
                lib_data = [
                    {'Species': sp, 'Peptide Sequence': m}
                    for sp, m_list in markers.items()
                    for m in m_list]

                # force a dataframe data type on the list
                lib_df = pd.DataFrame(lib_data)

                # create Pandas Series maps based on the original peptides dataframe
                prot_map = peptides.drop_duplicates("Sequence").set_index("Sequence")["Protein"]
                mz_map = peptides.drop_duplicates("Sequence").set_index("Sequence")["m/z meas."]
                charge_map = peptides.drop_duplicates("Sequence").set_index("Sequence")["z"]

                # map the values to the Reference Library dataframe
                lib_df["m/z"] = lib_df["Peptide Sequence"].map(mz_map).fillna("NA")
                lib_df["Charge"] = lib_df["Peptide Sequence"].map(charge_map).fillna("NA")
                lib_df["Protein"] = lib_df["Peptide Sequence"].map(prot_map).fillna("Unknown Protein")

                # optionally reorder columns for a cleaner look in Excel
                lib_df = lib_df[["Species", "Peptide Sequence", "m/z", "Charge", "Protein"]]

            # write to Excel
            lib_df.to_excel(writer, sheet_name='Reference_Library', index=False)
        else:
            print("Unable to create third Excel sheet (reference marker library) as no markers have been found.")

        # loop through all created Excel sheets
        for sheet_name in writer.sheets:
            worksheet = writer.sheets[sheet_name]

            # loop through all columns
            for column in worksheet.columns:
                max_length = 0

                # get column letter
                column_letter = column[0].column_letter

                # loop for every cell and find the longest text
                for cell in column:
                    max_length = max(max_length, len(str(cell.value)))

                # set column width plus a little reserve
                adjusted_width = max_length + 1
                worksheet.column_dimensions[column_letter].width = adjusted_width

    print(f"Saved results {filename} to Excel.")

def export_taxonomy(analysed_taxa, analysed_species, taxonomy_markers, species_markers, peptides,
                             output_folder, total_rep, ref_replicates, un_replicates=1):
    """
    Create an Excel file with analysis results in user-friendly format containing several sheets:
    first sheet - list of percentual match of each sample to each taxonomy group
    second sheet - list of percentual match of each sample to species within the select taxonomy group
    third sheet - detailed information on all found taxonomy group markers for every unknown sample
    fourth sheet - detailed information on all found species-level markers for every unknown sample
    fifth sheet - library of all found taxonomy group markers
    sixth sheet - library of all found species-level group markers
    :param analysed_taxa: dataframe of all unknown samples analysed against higher taxonomy
    :param analysed_species: dataframe of all unknown samples analysed on species-level
    :param output_folder: name of folder to be saved into
    :param peptides: dataframe of all discovered and adjusted peptides
    :param taxonomy_markers: library of found peptide markers for higher taxonomy groups
    :param species_markers: library of found species-specific markers
    :param ref_replicates: number of reference replicates used for peptide markers
    :param total_rep: total number of replicates available
    :param un_replicates: number of unknown replicates used
    :return: Excel workbook with 6 sheets containing various information from the analysis
    """
    print("\nInitiating data exportation to Excel ...")
    # create a folder if it doesn't already exist
    folder_path = Path(output_folder)
    if not folder_path.exists():
        folder_path.mkdir(parents=True)
        print(f"Created a folder for saved results: {folder_path}")

    # set filename
    filename = f"Analysis_report[{total_rep}-{ref_replicates}-{un_replicates}]_taxonomy"

    # create a full path by joining folder and filename paths
    full_path = folder_path / Path(f"{filename}.xlsx")

    # helper function for nice formatting in Excel
    def highlight_max_confidence(df):
        """
        For every unknown sample find the row with the highest confidence
        and add formatting style to that whole row
        :param df: dataframe of analysed unknown samples
        """
        # empty dataframe for CSS styles
        styles = pd.DataFrame('', index=df.index, columns=df.columns)

        # loop over every unknown sample
        for sample in df['Unknown Sample'].unique():
            # mask to separate current sample
            mask = df['Unknown Sample'] == sample

            # safety check - do we have data?
            if not df.loc[mask, 'Confidence [%]'].isna().all():
                # find the number of row at which there is maximum confidence
                max_idx = df.loc[mask, 'Confidence [%]'].idxmax()

                # add CSS style to the whole row with maximum confidence: here its light green and bold text
                styles.loc[max_idx, :] = 'background-color: #c6efce; font-weight: bold;'

        return styles

    # ExcelWriter is used to create multiple sheets
    with pd.ExcelWriter(full_path, engine="openpyxl") as writer:

        analysed_species = analysed_species.rename(columns={
            "Reference Group": "Reference Species",
            "Found Marker Keys": "Marker Key"})

        analysed_taxa = analysed_taxa.rename(columns={
            "Found Marker Keys": "Marker Key"})

        # the first list - table of percentual match of higher taxonomy
        print("Creating the first Excel sheet of higher taxonomy information...")

        if not (analysed_taxa.empty and not taxonomy_markers):
            # remove IDs to make the table more readable
            summary_df = analysed_taxa.drop(columns=["Found Marker Keys"], errors="ignore")

            cols = ["Unknown Sample", "Taxonomy Group", "Markers Found", "Confidence [%]",
                    "Reference Markers", "Match Score [%]", "Replicates"]
            summary_df = summary_df[cols]
            # round the percentual match to first decimal digit
            summary_df["Match Score [%]"] = summary_df["Match Score [%]"].round(1)
            summary_df["Confidence [%]"] = summary_df["Confidence [%]"].round(1)

            # export to Excel sheet without any indexes
            styled_summary = summary_df.style.apply(highlight_max_confidence, axis=None)
            styled_summary.to_excel(writer, sheet_name="Higher_Taxonomy_Assessment", index=False)
        else:
            print("Unable to create the first sheet of higher taxonomy assessment.")

        # the second list - table of percentual match on species-level
        print("Creating the second Excel sheet of species-level information...")
        if not (analysed_species.empty and not species_markers):

            cols = ["Unknown Sample", "Assigned Taxon", "Reference Species", "Markers Found", "Confidence [%]", "Reference Markers", "Match Score [%]", "Replicates"]
            summary_df_final = analysed_species[cols]

            # round the percentual match to first decimal digit
            summary_df_final["Match Score [%]"] = summary_df_final["Match Score [%]"].round(1)
            summary_df_final["Confidence [%]"] = summary_df_final["Confidence [%]"].round(1)

            # export to Excel sheet without any indexes
            styled_final = summary_df_final.style.apply(highlight_max_confidence, axis=None)
            styled_final.to_excel(writer, sheet_name="Species-level_Assessment", index=False)
        else:
            print("Unable to create the second sheet of species-level assessment.")

        # third sheet - table of concrete higher taxonomy marker peptides associated with unknown samples
        print("Creating the third Excel sheet of higher taxonomy detailed marker information...")

        if not analysed_taxa.empty and "Marker Key" in analysed_taxa.columns:

            # write every element in the dataframe column to different rows in Excel
            detail_df = analysed_taxa.explode("Marker Key")

            # create a Pandas Series of PeptideKey: Protein pairs based on peptides dataframe
            # drop_duplicates is necessary as two identical peptides can be found in 2 samples and are not unique
            prot_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["Protein"]
            mz_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["m/z meas."]
            charge_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["z"]

            # for every Marker Key find a corresponding protein and add it to Protein column
            detail_df["Protein"] = detail_df["Marker Key"].map(prot_map).fillna("Unknown Protein")
            detail_df["Charge"] = detail_df["Marker Key"].map(charge_map).fillna("NA")
            detail_df["m/z"] = detail_df["Marker Key"].map(mz_map).fillna("NA")

            detail_df = detail_df.dropna(subset=["Marker Key"])

            # define which columns are wanted in the Excel sheet
            cols_to_keep = ["Unknown Sample", "Taxonomy Group", "Marker Key", "m/z", "Charge", "Protein"]

            # write desired columns and all associated data to Excel sheet
            detail_df[cols_to_keep].to_excel(writer, sheet_name='Detailed_taxonomy_markers', index=False)
        else:
            print("Unable to create the third Excel sheet of detailed higher taxonomy markers.")

        # fourth sheet - table of concrete species-level marker peptides associated with unknown samples
        print("Creating the fourth Excel sheet of species-level detailed marker information...")
        if not analysed_species.empty and "Marker Key" in analysed_species.columns:

            # write every element in the dataframe column to different rows in Excel
            detail_df = analysed_species.explode("Marker Key")

            # create a Pandas Series of PeptideKey: Protein pairs based on peptides dataframe
            # drop_duplicates is necessary as two identical peptides can be found in 2 samples and are not unique
            prot_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["Protein"]
            mz_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["m/z meas."]
            charge_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["z"]

            # for every Marker Key find a corresponding protein and add it to Protein column
            detail_df["Protein"] = detail_df["Marker Key"].map(prot_map).fillna("Unknown Protein")
            detail_df["Charge"] = detail_df["Marker Key"].map(charge_map).fillna("NA")
            detail_df["m/z"] = detail_df["Marker Key"].map(mz_map).fillna("NA")

            detail_df = detail_df.dropna(subset=["Marker Key"])

            # define which columns are wanted in the Excel sheet
            cols_to_keep = ["Unknown Sample", "Reference Species", "Marker Key", "m/z", "Charge", "Protein"]

            # write desired columns and all associated data to Excel sheet
            detail_df[cols_to_keep].to_excel(writer, sheet_name='Detailed_species_markers', index=False)
        else:
            print("Unable to create the fourth Excel sheet of detailed species markers.")

        # fifth sheet - creation of reference library for higher taxonomy
        print("Creating the fifth Excel sheet listing the whole higher taxonomy reference marker library...")

        if not len(taxonomy_markers) == 0:
            # using a list comprehension on all species and associated markers
            lib_data = [
                {'Taxon': sp, 'Peptide Marker Key': m}
                for sp, m_list in taxonomy_markers.items()
                for m in m_list]

            # force a dataframe data type on the list
            lib_df = pd.DataFrame(lib_data)

            # create Pandas Series maps based on the original peptides dataframe
            prot_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["Protein"]
            mz_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["m/z meas."]
            charge_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["z"]

            # map the values to the Reference Library dataframe
            lib_df["m/z"] = lib_df["Peptide Marker Key"].map(mz_map).fillna("NA")
            lib_df["Charge"] = lib_df["Peptide Marker Key"].map(charge_map).fillna("NA")
            lib_df["Protein"] = lib_df["Peptide Marker Key"].map(prot_map).fillna("Unknown Protein")

            # optionally reorder columns for a cleaner look in Excel
            lib_df = lib_df[["Taxon", "Peptide Marker Key", "m/z", "Charge", "Protein"]]

            # write to Excel
            lib_df.to_excel(writer, sheet_name='Reference_Taxonomy_Library', index=False)

        else:
            print("Unable to create the fifth Excel sheet of higher taxonomy reference marker library")

        # sixth sheet - creation of reference library on species-level
        print("Creating the sixth Excel sheet listing the whole species-level reference marker library...")

        if not len(species_markers) == 0:
            # using a list comprehension on all species and associated markers
            lib_data = [
                {'Species': sp, 'Peptide Marker Key': m}
                for sp, m_list in species_markers.items()
                for m in m_list]

            # force a dataframe data type on the list
            lib_df = pd.DataFrame(lib_data)

            # create Pandas Series maps based on the original peptides dataframe
            prot_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["Protein"]
            mz_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["m/z meas."]
            charge_map = peptides.drop_duplicates("PeptideKey").set_index("PeptideKey")["z"]

            # map the values to the Reference Library dataframe
            lib_df["m/z"] = lib_df["Peptide Marker Key"].map(mz_map).fillna("NA")
            lib_df["Charge"] = lib_df["Peptide Marker Key"].map(charge_map).fillna("NA")
            lib_df["Protein"] = lib_df["Peptide Marker Key"].map(prot_map).fillna("Unknown Protein")

            # optionally reorder columns for a cleaner look in Excel
            lib_df = lib_df[["Species", "Peptide Marker Key", "m/z", "Charge", "Protein"]]

            # write to Excel
            lib_df.to_excel(writer, sheet_name='Reference_Species_Library', index=False)

        else:
            print("Unable to create the sixth Excel sheet of species-level reference marker library.")

            # loop through all created Excel sheets
        for sheet_name in writer.sheets:
            worksheet = writer.sheets[sheet_name]

            # loop through all columns
            for column in worksheet.columns:
                max_length = 0

                # get column letter
                column_letter = column[0].column_letter

                # loop for every cell and find the longest text
                for cell in column:
                    max_length = max(max_length, len(str(cell.value)))

                # set column width plus a little reserve
                adjusted_width = max_length + 1
                worksheet.column_dimensions[column_letter].width = adjusted_width

    print(f"Saved results {filename} to Excel.")