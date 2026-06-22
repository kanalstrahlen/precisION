import bisect
import os
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pyteomics import mass
from brainpy import isotopic_variants
import h5py
from tqdm import tqdm

from logger import info, warn, error, success


class DbSearchFunctions():
    def write_read_proteoform_file(
        self, h5_file_path, database_file,
        min_mw, max_mw, ion_type, disulfide
    ):
        info("[b]Generating/loading proteoform file...")
        # if h5 proteoform file already exists, use already created file
        if os.path.exists(h5_file_path):
            info(f"Using {h5_file_path} as proteoform file.")
        else:
            info("Generating proteoform file. Please be patient.")

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

            # build fragment-mass index for searching
            b_index_masses, b_index_owners = self.build_fragment_index(theo_b_ions)
            y_index_masses, y_index_owners = self.build_fragment_index(theo_y_ions)

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
                file.create_dataset("b_index_masses", data=b_index_masses)
                file.create_dataset("b_index_owners", data=b_index_owners)
                file.create_dataset("y_index_masses", data=y_index_masses)
                file.create_dataset("y_index_owners", data=y_index_owners)
                file.create_dataset("disulfide_positions", data = disulfide_position_array)

        # read h5 file and save lists
        with h5py.File(h5_file_path, "r") as file:
            # asstr decodes all strings in one go
            proteoform_ids = list(file["id"].asstr()[:])
            proteoform_sequences = list(file["sequence"].asstr()[:])
            proteoform_n_term_mods = list(file["n_term_mod"])
            proteoform_mws = list(file["mw"])
            # load the fragment index; rebuild from the padded ion arrays for older files
            if "b_index_masses" in file:
                b_index_masses = np.array(file["b_index_masses"])
                b_index_owners = np.array(file["b_index_owners"])
                y_index_masses = np.array(file["y_index_masses"])
                y_index_owners = np.array(file["y_index_owners"])
            else:
                b_index_masses, b_index_owners = self.build_fragment_index(np.array(file["b_ions"]))
                y_index_masses, y_index_owners = self.build_fragment_index(np.array(file["y_ions"]))
            # remove padding values
            temp_disulfide_positions = np.array(file["disulfide_positions"])
            proteoform_disulfide_positions = [
                row[row != 999.9999].tolist() for row in temp_disulfide_positions
            ]

        return (
            proteoform_ids,
            proteoform_sequences,
            proteoform_n_term_mods,
            proteoform_mws,
            proteoform_disulfide_positions,
            b_index_masses,
            b_index_owners,
            y_index_masses,
            y_index_owners
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


    def build_fragment_index(self, theoretical_ions_for_each_proteoform):
        # flatten every fragment into one mass-sorted array, keeping its proteoform
        masses = []
        owners = []
        for i, ions in enumerate(theoretical_ions_for_each_proteoform):
            ions = np.asarray(ions)
            ions = ions[ions > 0] # skip zero padding
            masses.append(ions)
            owners.append(np.full(len(ions), i))
        masses = np.concatenate(masses)
        owners = np.concatenate(owners)
        order = np.argsort(masses)
        return masses[order], owners[order]


    def count_matches_by_index(
        self, index_masses, index_owners, observed_ions, num_proteoforms, ppm_threshold
    ):
        # for each peak find the fragments within ppm tolerance (index position is a unique id)
        matched = []
        for peak in observed_ions:
            low = np.searchsorted(index_masses, peak * (1 - ppm_threshold / 1e6), side="left")
            high = np.searchsorted(index_masses, peak * (1 + ppm_threshold / 1e6), side="right")
            matched.append(np.arange(low, high))

        if matched:
            matched = np.concatenate(matched)
        else:
            matched = np.array([], dtype=int)

        # count each fragment once, then tally matches per proteoform
        matched = np.unique(matched)
        counts = np.bincount(index_owners[matched], minlength=num_proteoforms)
        return counts.tolist()


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

        # precompute each residue's b-ion mass once (fast_mass is slow to call per residue)
        unique_residues = set().union(*sequence_list)
        aa_b_mass = {aa: mass.fast_mass(aa, ion_type="b") for aa in unique_residues}

        for i, (sequence, n_terminal_modification_mass) in tqdm(
            enumerate(zip(sequence_list, n_terminal_modification_mass_list)),
            desc='Processing sequences',
            total=len(sequence_list)
        ):
            sequence_mass = mass.fast_mass(sequence, ion_type="M")
            aa_masses = np.array([aa_b_mass[aa] for aa in sequence])

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
        self, observed_ions, b_index_masses, b_index_owners,
        y_index_masses, y_index_owners, accuracy,
        proteoform_ids, proteoform_sequences, proteoform_n_term_mods, proteoform_mws,
        proteoform_disulfide_positions
    ):
        # count ion matches for each proteoform via the inverted fragment index
        num_proteoforms = len(proteoform_sequences)
        info('Searching b or c-type ions...')
        b_match_list = self.count_matches_by_index(
            b_index_masses, b_index_owners, observed_ions, num_proteoforms, accuracy
        )
        info('Searching y or z-type ions...')
        y_match_list = self.count_matches_by_index(
            y_index_masses, y_index_owners, observed_ions, num_proteoforms, accuracy
        )
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
