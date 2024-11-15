import bisect
import os
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pyteomics import mass
from brainpy import isotopic_variants
import h5py
from tqdm import tqdm
from fastDatabaseSearchFunctions import FastDbSearchFunctions


# to change; change to index search with nesteted B tree

class DbSearchFunctions():
    def write_read_proteoform_file(
        self, h5_file_path, database_file,
        min_mw, max_mw, ion_type, disulfide
    ):
        # if h5 proteoform file already exists, use already created file
        if os.path.exists(h5_file_path):
            print(f"Using {h5_file_path} as proteoform file.")
        else:
            print("Generating proteoform file. Please be patient.")

            # read in all key proteoforms from xml file
            (
                proteoform_ids,
                proteoform_sequences,
                proteoform_n_terminal_modifications,
                proteoform_disulfide_positions
            ) = self.load_database_xml(database_file, disulfide)


            #apply mass filter based on user input
            (
                filtered_ids,
                filtered_sequences,
                filtered_n_term_mods,
                filtered_mws,
                filtered_disulfide_positions
            ) = self.filter_proteoforms_by_mass(
                proteoform_ids,
                proteoform_sequences,
                proteoform_n_terminal_modifications,
                min_mw,
                max_mw,
                proteoform_disulfide_positions
            )

            # generate theoretical ions for each proteoform and store in a list of lists
            theo_b_ions, theo_y_ions = self.theoretical_frag_generator_multiple_search(
                filtered_sequences,
                filtered_n_term_mods,
                filtered_disulfide_positions,
                ion_type
            )

            # zero pad each array to allow h5 file to be saved
            b_max_length = max(len(arr) for arr in theo_b_ions)
            y_max_length = max(len(arr) for arr in theo_y_ions)
            padded_theo_b_ions = [
                np.pad(arr, (b_max_length - len(arr), 0), "constant")
                for arr in theo_b_ions
            ]
            padded_theo_y_ions = [
                np.pad(arr, (y_max_length - len(arr), 0), "constant")
                for arr in theo_y_ions
            ]
            theo_b_ions_array = np.array(padded_theo_b_ions)
            theo_y_ions_array = np.array(padded_theo_y_ions)

            max_disulfide_length = max(len(arr) for arr in filtered_disulfide_positions)
            padded_disulfide_positions = [
                np.pad(
                    arr,
                    (max_disulfide_length - len(arr), 0),
                    "constant",
                    constant_values = "999.9999"
                )
                for arr in filtered_disulfide_positions
            ]
            disulfide_position_array = np.array(padded_disulfide_positions)

            # write h5 file in the same path as the directory
            with h5py.File(h5_file_path, "w") as file:
                dt = h5py.special_dtype(vlen=str)
                file.create_dataset("id", data=filtered_ids, dtype=dt)
                file.create_dataset("sequence", data=filtered_sequences, dtype=dt)
                file.create_dataset("n_term_mod", data=filtered_n_term_mods)
                file.create_dataset("mw", data=filtered_mws)
                file.create_dataset("b_ions", data=theo_b_ions_array)
                file.create_dataset("y_ions", data=theo_y_ions_array)
                file.create_dataset("disulfide_positions", data = disulfide_position_array)

        # read h5 file and save lists
        with h5py.File(h5_file_path, "r") as file:
            proteoform_ids = list(file["id"])
            proteoform_ids = [proteoform_id.decode("utf-8") for proteoform_id in proteoform_ids]
            proteoform_sequences = list(file["sequence"])
            proteoform_sequences = [sequence.decode("utf-8") for sequence in proteoform_sequences]
            proteoform_n_term_mods = list(file["n_term_mod"])
            proteoform_mws = list(file["mw"])
            theo_b_ions = list(file["b_ions"])
            theo_y_ions = list(file["y_ions"])
            proteoform_disulfide_positions = []
            # remove padding values
            temp_disulfide_positions = list(file["disulfide_positions"])
            for pos_list in temp_disulfide_positions:
                new_pos_list = [x for x in pos_list if x != 999.9999]
                proteoform_disulfide_positions.append(new_pos_list)

        return (
            proteoform_ids,
            proteoform_sequences,
            proteoform_n_term_mods,
            proteoform_mws,
            proteoform_disulfide_positions,
            theo_b_ions,
            theo_y_ions
        )


    def load_database_xml(self, input_file, disulfide):
        n_terminal_modifications = {
            "none": 0.000000,
            "acetylation": 42.010565,
            #"formylation": 27.994915,
        }

        proteoform_ids = []
        proteoform_sequences = []
        proteoform_n_terminal_modifications = []
        proteoform_disulfide_positions = []

        with open(input_file, "r") as file:
            for record in SeqIO.parse(file, "uniprot-xml"):
                protein_name = str(record.name)

                protein_full_name = str(record.description)
                if len(protein_full_name) >= 25:
                    protein_full_name = protein_full_name[:24] + "..."

                protein_sequence = str(record.seq)
                # some protein sequences contain X?
                protein_sequence = protein_sequence.replace("X", "G")

                # generate different truncated products using the different forms
                truncation_ids = []
                truncation_sequences = []
                truncation_disulfide_positions = []

                # determine index of disulfide bonded cysteines, store in disulfide_res
                disulfide_res = []
                if disulfide == "ox":
                    for feature in record.features:
                        if feature.type == "disulfide bond":
                            if isinstance(feature.location.start, int):
                                if isinstance(feature.location.end, int):
                                    disulfide_res.append(feature.location.start)
                                    disulfide_res.append(feature.location.end)


                for feature in record.features:
                    if feature.type in [
                    "chain",
                    "peptide"
                    ]:
                        cut_site_start = feature.location.start
                        cut_site_end = feature.location.end

                        if isinstance(cut_site_start, int) and isinstance(cut_site_end, int):
                            truncation_ids.append(
                                f"{protein_name} {protein_full_name} "
                                f"({cut_site_start + 1}-{cut_site_end})"
                            )
                            truncation_sequences.append(
                                protein_sequence[cut_site_start:cut_site_end]
                            )

                            # make a temp list for disulfides within this chain
                            # then write corrected positions to list
                            disulfide_temp = np.zeros(cut_site_end - cut_site_start)
                            for res in disulfide_res:
                                if cut_site_start <= res <= cut_site_end:
                                    if res - cut_site_start <= len(disulfide_temp)-1:
                                        # won't work if final residue (to do)
                                        disulfide_temp[res - cut_site_start] = -1.007825
                            truncation_disulfide_positions.append(disulfide_temp)

                            if cut_site_start == 0: # remove n terminal residue
                                truncation_ids.append(
                                    f"{protein_name} {protein_full_name} "
                                    f"({cut_site_start + 2}-{cut_site_end})"
                                )
                                truncation_sequences.append(
                                    protein_sequence[cut_site_start+1:cut_site_end]
                                )

                                disulfide_temp = np.zeros(cut_site_end - cut_site_start - 1)
                                for res in disulfide_res:
                                    if cut_site_start <= res <= cut_site_end:
                                        if res - cut_site_start - 1 <= len(disulfide_temp)-1:
                                            disulfide_temp[res - cut_site_start - 1] = -1.007825
                                truncation_disulfide_positions.append(disulfide_temp)

                #chain_count = 0
                #for feature in record.features:
                #    if feature.type in ["chain", "peptide"]:
                #        chain_count += 1

                #if chain_count == 0: # when using processed xml files
                cut_site_start = 0
                cut_site_end = len(protein_sequence)
                if isinstance(cut_site_start, int) and isinstance(cut_site_end, int):
                    truncation_ids.append(
                        f"{protein_name} {protein_full_name} "
                        f"({cut_site_start + 1}-{cut_site_end})"
                    )
                    truncation_sequences.append(
                        protein_sequence[cut_site_start:cut_site_end]
                    )

                    # make a temp list for disulfides within this chain
                    # then write corrected positions to list
                    disulfide_temp = np.zeros(cut_site_end - cut_site_start)
                    for res in disulfide_res:
                        if cut_site_start <= res <= cut_site_end:
                            if res - cut_site_start <= len(disulfide_temp)-1:
                                disulfide_temp[res - cut_site_start] = -1.007825
                    truncation_disulfide_positions.append(disulfide_temp)
                    if cut_site_start == 0: # remove n terminal residue
                        truncation_ids.append(
                            f"{protein_name} {protein_full_name} "
                            f"({cut_site_start + 2}-{cut_site_end})"
                        )
                        truncation_sequences.append(
                            protein_sequence[cut_site_start+1:cut_site_end]
                        )
                        disulfide_temp = np.zeros(cut_site_end - cut_site_start - 1)
                        for res in disulfide_res:
                            if cut_site_start <= res <= cut_site_end:
                                if res - cut_site_start - 1 <= len(disulfide_temp)-1:
                                    disulfide_temp[res - cut_site_start - 1] = -1.007825
                        truncation_disulfide_positions.append(disulfide_temp)

                # generate all proteoforms (3 for each truncation)
                for i, truncation in enumerate(truncation_sequences):
                    for n_mod in n_terminal_modifications:
                        truncation_id = truncation_ids[i]
                        proteoform_id = f"{truncation_id} N-term: {n_mod}"
                        proteoform_ids.append(proteoform_id)
                        proteoform_sequences.append(truncation)
                        proteoform_n_terminal_modifications.append(
                            n_terminal_modifications[n_mod]
                        )
                        proteoform_disulfide_positions.append(
                            truncation_disulfide_positions[i]
                        )

        return (
            proteoform_ids, proteoform_sequences,
            proteoform_n_terminal_modifications, proteoform_disulfide_positions
        )


    def filter_proteoforms_by_mass(
        self, proteoform_ids, proteoform_sequences, proteoform_n_terminal_modifications,
        lower_limit, upper_limit, proteoform_disulfide_positions
    ):
        # lists after filtereing
        filtered_proteoform_ids = []
        filtered_proteoform_sequences = []
        filtered_proteoform_n_terminal_modifications = []
        filtered_proteoform_mws = []
        filtered_disulfide_positions = []

        # calculate molecular weight and apply filter
        for i, sequence in enumerate(proteoform_sequences):
            try:
                mw = (
                    ProteinAnalysis(sequence).molecular_weight() +
                    proteoform_n_terminal_modifications[i]
                )

                if lower_limit <= mw <= upper_limit:
                    filtered_proteoform_ids.append(proteoform_ids[i])
                    filtered_proteoform_sequences.append(sequence)
                    filtered_proteoform_n_terminal_modifications.append(
                        proteoform_n_terminal_modifications[i]
                    )
                    filtered_proteoform_mws.append(mw)
                    filtered_disulfide_positions.append(proteoform_disulfide_positions[i])
            except:
                pass

        return (
            filtered_proteoform_ids,
            filtered_proteoform_sequences,
            filtered_proteoform_n_terminal_modifications,
            filtered_proteoform_mws,
            filtered_disulfide_positions
        )


    def count_matches_for_proteins(
        self, theoretical_ions_for_each_proteoform,
        observed_ions, ppm_threshold
    ):
        # takes input of list of lists (each list is a list of theoretical ions for each proteoform)
        # binary search w/ ppm tolerance helper function
        def binary_search(arr, target, ppm_threshold):
            low = bisect.bisect_left(arr, target * (1 - ppm_threshold / 1e6))
            high = bisect.bisect_right(arr, target * (1 + ppm_threshold / 1e6))
            return high - low

        # sort for binary search
        observed_ions.sort()

        # list to store the number of matches for each protein
        matches_per_proteoform = []

        # iterate through each proteoform and count matches
        for theoretical_ions in tqdm(
            theoretical_ions_for_each_proteoform,
            total=len(theoretical_ions_for_each_proteoform)
        ):
            theoretical_ions.sort() # sort for binary search
            matches = 0 # initialize count
            for ion in theoretical_ions:
                match_test = binary_search(observed_ions, ion, ppm_threshold)
                if match_test >= 1:
                    matches += 1

            matches_per_proteoform.append(matches)

        return matches_per_proteoform


    def report_matches_for_protein(
        theoretical_ions, observed_ions,
        ppm_threshold, all_charge=False
    ):
        # binary search w/ ppm tolerance helper function
        def binary_search(arr, target, ppm_threshold):
            low = bisect.bisect_left(arr, target * (1 - ppm_threshold / 1e6))
            high = bisect.bisect_right(arr, target * (1 + ppm_threshold / 1e6))
            return high - low, low, range(low, high)

        matched_theoretical_ion_indices = []
        matched_observed_ion_indices = []

        # sort both lists for binary search; nearly always sorted but this is still there
        theoretical_ions.sort()
        observed_ions.sort()

        for i, theo_ion in enumerate(theoretical_ions):
            num_matches, match, all_matches = binary_search(observed_ions, theo_ion, ppm_threshold)
            if num_matches >= 1:
                if all_charge:
                    for matched_index in all_matches:
                        matched_theoretical_ion_indices.append(i)
                        matched_observed_ion_indices.append(matched_index)
                else:
                    matched_theoretical_ion_indices.append(i)
                    matched_observed_ion_indices.append(match)

        return matched_theoretical_ion_indices, matched_observed_ion_indices


    def theoretical_frag_generator_single_search(
        sequence,
        n_terminal_modification_mass,
        disulfide_positions,
        ion_type = "b/y"
    ):
        if ion_type == "c/z•":
            cz_ions = True
        elif ion_type == "b/y":
            cz_ions = False

        # lists of b and y ions
        theo_b_ions_list = []
        theo_y_ions_list = []

        sequence_mass = mass.fast_mass(sequence, ion_type="M") # intact mass
        aa_masses = np.array(
        [mass.fast_mass(aa, ion_type="b") for aa in sequence] # individual residue masses
    )

        cumulative_disulfides = np.cumsum(disulfide_positions)

        if cz_ions:
            b_ion_masses = (
                np.cumsum(aa_masses) +
                np.cumsum(disulfide_positions) +
                n_terminal_modification_mass +
                17.026549 # + nh3
            )
            y_ion_masses = (
                sequence_mass +
                np.sum(disulfide_positions) -
                np.cumsum(aa_masses) -
                #corrects positioning
                np.insert(cumulative_disulfides, len(cumulative_disulfides), 0)[1:] -
                17.026549 + 1.00783 # - nh3 + extra proton to balance electron charge + electron
            )
        else:
            b_ion_masses = (
                np.cumsum(aa_masses) +
                np.cumsum(disulfide_positions) +
                n_terminal_modification_mass
            )
            y_ion_masses = (
                sequence_mass +
                np.sum(disulfide_positions) -
                np.cumsum(aa_masses) -
                np.insert(cumulative_disulfides, len(cumulative_disulfides), 0)[1:]
            )

        theo_b_ions_list.append(b_ion_masses[:-1]) # exclude the dehydrated full protein
        theo_y_ions_list.append(y_ion_masses[:-1])

        # sort for binary search
        theo_y_ions_list.sort()

        return theo_b_ions_list, theo_y_ions_list


    def theoretical_frag_generator_multiple_search(
        self,
        sequence_list,
        n_terminal_modification_mass_list,
        disulfide_positions,
        ion_type = "b/y"
    ):
        if ion_type == "c/z•":
            cz_ions = True
        elif ion_type == "b/y":
            cz_ions = False

        # lists of b and y ions
        theo_b_ions_list = []
        theo_y_ions_list = []

        for i, (sequence, n_terminal_modification_mass) in tqdm(
            enumerate(zip(sequence_list, n_terminal_modification_mass_list)),
            desc='Processing sequences',
            total=len(sequence_list)
        ):
            sequence_mass = mass.fast_mass(sequence, ion_type="M")
            aa_masses = np.array([mass.fast_mass(aa, ion_type="b") for aa in sequence])

            cumulative_disulfides = np.cumsum(disulfide_positions[i])

            if cz_ions:
                b_ion_masses = (
                    np.cumsum(aa_masses) +
                    np.cumsum(disulfide_positions[i]) +
                    n_terminal_modification_mass +
                    17.026549 # + nh3
                )
                y_ion_masses = (
                    sequence_mass +
                    np.sum(disulfide_positions[i]) -
                    np.cumsum(aa_masses) -
                    np.insert(cumulative_disulfides, len(cumulative_disulfides), 0)[1:] -
                    17.026549 + 1.00783
                )
            else:
                b_ion_masses = (
                    np.cumsum(aa_masses) +
                    np.cumsum(disulfide_positions[i]) +
                    n_terminal_modification_mass
                )
                y_ion_masses = (
                    sequence_mass +
                    np.sum(disulfide_positions[i]) -
                    np.cumsum(aa_masses) -
                    np.insert(cumulative_disulfides, len(cumulative_disulfides), 0)[1:]
                )

            theo_b_ions_list.append(b_ion_masses[:-1])
            theo_y_ions_list.append(y_ion_masses[:-1])

        return theo_b_ions_list, theo_y_ions_list


    def averagine(mono_mass, charge, abundance):
        averagine_units = mono_mass / 111.0543
        num_c = int(np.round(4.9384 * averagine_units))
        num_h = int(np.round(7.7583 * averagine_units))
        num_n = int(np.round(1.3577 * averagine_units))
        num_o = int(np.round(1.4773 * averagine_units))
        num_s = int(np.round(0.0417 * averagine_units))

        # molecule with closest integer number of averagines
        averagine_molecule_mono_mass = 12 * num_c + 1.007825035 * num_h + 15.99491463 * num_o\
                                       + 14.003074 * num_n + 31.9720707 * num_s

        # make up the difference with protons and recalculate mass
        diff_h = np.round((mono_mass - averagine_molecule_mono_mass) / 1.007825035)
        num_h = int(num_h + diff_h)
        averagine_molecule_mono_mass = 12 * num_c + 1.007825035 * num_h + 15.99491463 * num_o\
                                       + 14.003074 * num_n + 31.9720707 * num_s
        averagine_molecule = {
            "H": num_h,
            "C": num_c,
            "O": num_o,
            "N": num_n,
            "S": num_s
        }

        # calculate isotopic distribution using brainpy
        theoretical_isotopic_cluster = isotopic_variants(
            averagine_molecule,
            npeaks=20,
            charge=charge
        )

        # store envelope in x and y lists
        x=[]
        y=[]
        for peak in theoretical_isotopic_cluster:
            x.append(peak.mz)
            y.append(peak.intensity)

        # shift by small remaining distance
        mono_theo_mass = x[0]
        diff = mono_theo_mass - (mono_mass + charge * 1.007825035) / charge
        x2 = []
        for peak in x:
            x2.append(peak - diff)

        y = [inten / max(y) * abundance for inten in y]
        return(x2,y)


    def database_search(
        self, observed_ions, theo_b_ions, theo_y_ions, accuracy,
        proteoform_ids, proteoform_sequences, proteoform_n_term_mods, proteoform_mws,
        proteoform_disulfide_positions
    ):
        # count ion matches for each proteoform
        print('Searching b or c-type ions...')
        fast_search = FastDbSearchFunctions(accuracy)
        b_match_list = fast_search.count_matches_for_proteins(theo_b_ions, observed_ions)
        #b_match_list = self.count_matches_for_proteins(theo_b_ions, observed_ions, accuracy)
        print('Searching y or z-type ions...')
        fast_search = FastDbSearchFunctions(accuracy)
        y_match_list = fast_search.count_matches_for_proteins(theo_y_ions, observed_ions)
        #y_match_list = self.count_matches_for_proteins(theo_y_ions, observed_ions, accuracy)
        total_match_list = [b + y for b, y in zip(b_match_list, y_match_list)]

        # sort and find the top 100 matches
        num_matches = 100
        indexed_match_list = list(enumerate(total_match_list))
        sorted_match_list = sorted(indexed_match_list, key=lambda x: x[1], reverse=True)
        num_matches = min(len(proteoform_sequences), num_matches)
        top_matches = [index for index, _ in sorted_match_list[:num_matches]]

        matched_ids = [proteoform_ids[i] for i in top_matches]
        matched_sequences = [proteoform_sequences[i] for i in  top_matches]
        matched_n_term_mods = [proteoform_n_term_mods[i] for i in  top_matches]
        matched_mws = [round(proteoform_mws[i]) for i in  top_matches]
        matched_disulfide_positions = [proteoform_disulfide_positions[i] for i in top_matches]
        num_matched_b_ions = [b_match_list[i] for i in top_matches]
        num_matched_y_ions = [y_match_list[i] for i in top_matches]

        return (
            matched_ids,
            matched_sequences,
            matched_n_term_mods,
            matched_mws,
            matched_disulfide_positions,
            num_matched_b_ions,
            num_matched_y_ions
        )
