from itertools import product
import numpy as np
from Bio import SeqIO
from brainpy import isotopic_variants

# despite the ugly code, these function work better than networkx options

class DenovoFunctions():
    def find_aa_loss(self, peak_list, error):
        # amino acid mass lists
        aa_list = [
            "A", "R", "N", "D", "C",
            "E", "Q", "G", "H", "L",
            "K", "M", "F", "P", "S",
            "T", "W", "Y", "V"
        ]
        mass_list = [
            71.037114, 156.10111, 114.042928, 115.026944, 103.009185,
            129.042594, 128.058578, 57.021464, 137.058912, 113.084064,
            128.094963, 131.040485, 147.068414, 97.052764, 87.032029,
            101.047679, 186.079313, 163.063329, 99.068414
        ]

        # list of masses
        frag_masses = peak_list["monoisotopic_mw"].tolist()
        frag_mass_diff = np.subtract.outer(frag_masses, frag_masses)

        # for each amino acid, look for pairs of ions separated by that mass
        # store their indexes as a tuple
        aa_loss_pair_list = []
        for aa in range(len(aa_list)):
            mass = mass_list[aa]
            match_index = np.where((frag_mass_diff>=(mass-error))
                                   &(frag_mass_diff<=(mass+error)))
            num_match = len(match_index[0])
            ion_pairs = [] # list of all ion pairs for a particular aa
            for i in range(num_match):
                hi_mass_index = match_index[0][i]
                lo_mass_index = match_index[1][i]
                ion_pairs.append((hi_mass_index,
                                  lo_mass_index))
            aa_loss_pair_list.append(ion_pairs)
        return aa_loss_pair_list


    def find_peptides(self, pairs_list, length):

        def flatten(d):
            for i in d:
                yield from [i] if not isinstance(i, tuple) else flatten(i)

        pair_list = []
        for i in pairs_list:
            for j in i:
                pair_list.append(j)
        temp_peptide_list = []
        peptide_list = []
        peptide_list = pair_list
        for i in range(length-1):
            for pair1 in peptide_list:
                resid1 = pair1[-1]
                for pair2 in pair_list:
                    resid2 = pair2[0]
                    if resid1 == resid2:
                        joint = pair1, pair2[1]
                        temp_peptide_list.append(joint)
            peptide_list = temp_peptide_list
            temp_peptide_list = []
        peptides = []
        for data in peptide_list:
            temp = tuple(flatten(data))
            peptides.append(temp)
        # peptides contains a list of peak indices corresponding to sequence tags
        return peptides


    def find_aa_seq(self, peak_list, peptides, error):
        aa_dict = {"A" : [71.03711-error, 71.03711+error],
                   "R" : [156.10111-error, 156.10111+error],
                   "N" : [114.04293-error, 114.04293+error],
                   "D" : [115.02694-error, 115.02694+error],
                   "C" : [103.00919-error, 103.00919+error],
                   "E" : [129.04259-error, 129.04259+error],
                   "Q" : [128.05858-error, 128.05858+error],
                   "G" : [57.02146-error, 57.02146+error],
                   "H" : [137.05891-error, 137.05891+error],
                   "L" : [113.08406-error, 113.08406+error],
                   "K" : [128.09496-error, 128.09496+error],
                   "M" : [131.04049-error, 131.04049+error],
                   "F" : [147.06841-error, 147.06841+error],
                   "P" : [97.05276-error, 97.05276+error],
                   "S" : [87.03203-error, 87.03203+error],
                   "T" : [101.04768-error, 101.04768+error],
                   "W" : [186.07931-error, 186.07931+error],
                   "Y" : [163.06333-error, 163.06333+error],
                   "V" : [99.06841-error, 99.06841+error]}

        mz_list = []
        mass_list = []
        mass_diff_list = []
        aa_list = []
        seq_list = []
        ion_charge_list = []
        abundance_list = []

        masses = peak_list["monoisotopic_mw"].tolist()
        mzs = peak_list["monoisotopic_mz"].tolist()
        charges = peak_list["charge"].tolist()
        abundances = peak_list["abundance"].tolist()

        for index_list in peptides:
            charge_list = []
            # check ions have same charge
            for ion in index_list:
                charge_check = charges[ion]
                charge_list.append(charge_check)
            same_charge = charge_list.count(charge_list[0]) == len(charge_list)

            if same_charge:
            # extract mz, mass, and abundance for each ion
                mz_temp_list = []
                mass_temp_list = []
                mass_diff_temp_list = []
                aa_temp_list = []
                abundance_temp_list = []

                for m in index_list:
                    mz = mzs[m]
                    mz_temp_list.append(mz)
                    mass = masses[m]
                    mass_temp_list.append(mass)
                    abundance_temp_list.append(abundances[m])

                # look at each pair to assemble sequence from mass differences
                for j in range(len(index_list)-1):
                    m = index_list[j]
                    n = index_list[(j+1)]
                    mass_diff = masses[m]-masses[n]
                    mass_diff_temp_list.append(mass_diff)
                    for aa, bounds in aa_dict.items():
                        if mass_diff >= bounds[0] and mass_diff <= bounds[1]:
                            aa_temp_list.append(aa)

                seq = "".join(aa_temp_list)
                seq_list.append(seq)
                mz_list.append(mz_temp_list)
                mass_list.append(mass_temp_list)
                mass_diff_list.append(mass_diff_temp_list)
                aa_list.append(aa_temp_list)
                abundance_list.append(abundance_temp_list)
                ion_charge_list.append(charge_list[0])
        return(
            mz_list, mass_list, mass_diff_list, aa_list,
            seq_list, ion_charge_list, abundance_list
        )


    def averagine(self, mono_mass, charge):
        averagine_units = mono_mass/111.0543
        num_c = int(np.round(4.9384*averagine_units))
        num_h = int(np.round(7.7583*averagine_units))
        num_n = int(np.round(1.3577*averagine_units))
        num_o = int(np.round(1.4773*averagine_units))
        num_s = int(np.round(0.0417*averagine_units))
        molecule_mono_mass = 12*num_c+1.007825035*num_h+15.99491463*num_o\
                           +14.003074*num_n+31.9720707*num_s
        diff = np.round((mono_mass-molecule_mono_mass)/1.007825035)
        num_h = int(num_h+diff)
        molecule_mono_mass = 12*num_c+1.007825035*num_h+15.99491463*num_o\
                           +14.003074*num_n+31.9720707*num_s

        ion = {"H":num_h,
               "C":num_c,
               "O": num_o,
               "N": num_n,
               "S": num_s
              }

        theoretical_isotopic_cluster = isotopic_variants(
            ion,
            npeaks=10,
            charge=charge
        )

        x=[]
        y=[]

        for peak in theoretical_isotopic_cluster:
            x.append(peak.mz)
            y.append(peak.intensity)

        mono_theo = x[0]
        diff = mono_theo-(mono_mass+charge*1.007825035)/charge
        x2 = []
        for peak in x:
            x2.append(peak-diff)
        return(x2,y)


    def search_uniprot_fasta(self, xml_file, sequence):
        matches_prot = []
        matches_prot_name = []
        matches_start = []
        matches_end = []
        matches_sequence = []

        input_list = [sequence]
        # generate all possible replacements of "L" with "I" using product
        permutations = []
        for item in input_list:
            l_indices = [i for i, char in enumerate(item) if char == "L"]
            for replacement_combination in product("LI", repeat=len(l_indices)):
                replaced_item = list(item)
                for i, replacement in zip(l_indices, replacement_combination):
                    replaced_item[i] = replacement
                permutations.append("".join(replaced_item))
        # add the reverse complements of the permutations
        permutations += [seq[::-1] for seq in permutations]
        # create a set for faster membership checking
        sequence_set = set(permutations)

        xml_basename = xml_file.replace(".xml", "")
        fasta_file = f"{xml_basename}.fasta"

        with open(fasta_file, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                seqtest = str(record.seq)
                prot = str(record.id)
                prot_name = str(record.description)
                for seq in sequence_set:
                    if seq in seqtest:
                        start = seqtest.index(seq)
                        stop = start + len(seq)
                        matches_prot.append(prot)
                        matches_start.append(start)
                        matches_end.append(stop)
                        matches_prot_name.append(prot_name)
                        matches_sequence.append(seqtest)
                        break  # break after the first match for this sequence
        return matches_prot, matches_prot_name, matches_start, matches_end, matches_sequence
