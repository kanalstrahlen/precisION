import re
import bisect
import numpy as np
from pyteomics import mass, parser
from brainpy import isotopic_variants
from scipy.optimize import minimize

class ProteoformAnalysisRun():
    def generate_theoretical_fragments(self, sequence, max_charge, min_mz, max_mz, ion_type):
        unmod_seq, mod_arr = self.process_sequence_string(sequence)
        mod_arr = [self.mod_string_to_dict(mod) for mod in mod_arr]

        def formula_to_string(formula):
            formula_string=""
            for element, count in formula.items():
                if count >= 1:
                    formula_string += element
                    if count >= 2:
                        formula_string += str(count)
            return formula_string


        b_ion_formulae = []
        b_ion_masses = []
        b_ion_most_abundant_masses = []
        for i in range(1,len(unmod_seq)):
            ion_seq = unmod_seq[:i]
            formula = self.calculate_molecular_formula(ion_seq)
            if ion_type == "c/z•":
                formula["H"] += 3
                formula["N"] += 1
            mod_dicts = mod_arr[0:i]
            mod_sum = {key: 0 for d in mod_dicts for key in d}
            for d in mod_dicts:
                for key, value in d.items():
                    mod_sum[key] += value
            formula = {key: formula.get(key, 0) + mod_sum.get(key, 0) for key in set(formula) | set(mod_sum)}
            ion_mass = mass.calculate_mass(composition=formula)
            formula_string = formula_to_string(formula)
            most_abundant_formula = mass.most_probable_isotopic_composition(formula=formula_string)[0]
            most_abundant_mass = mass.calculate_mass(most_abundant_formula)
            b_ion_formulae.append(formula)
            b_ion_masses.append(ion_mass)
            b_ion_most_abundant_masses.append(most_abundant_mass)

        y_ion_formulae = []
        y_ion_masses = []
        y_ion_most_abundant_masses = []
        for i in range(len(unmod_seq)-1):
            ion_seq = unmod_seq[i+1:]
            formula = self.calculate_molecular_formula(ion_seq)
            formula["H"] += 2
            formula["O"] += 1
            if ion_type == "c/z•":
                formula["H"] += -2
                formula["N"] += -1
            mod_dicts = mod_arr[i+1:]
            mod_sum = {key: 0 for d in mod_dicts for key in d}
            for d in mod_dicts:
                for key, value in d.items():
                    mod_sum[key] += value
            formula = {key: formula.get(key, 0) + mod_sum.get(key, 0) for key in set(formula) | set(mod_sum)}
            ion_mass = mass.calculate_mass(composition=formula)
            formula_string = formula_to_string(formula)
            most_abundant_formula = mass.most_probable_isotopic_composition(formula=formula_string)[0]
            most_abundant_mass = mass.calculate_mass(most_abundant_formula)
            y_ion_formulae.append(formula)
            y_ion_masses.append(ion_mass)
            y_ion_most_abundant_masses.append(most_abundant_mass)

        y_ion_formulae = y_ion_formulae[::-1]
        y_ion_masses = y_ion_masses[::-1]
        y_ion_most_abundant_masses = y_ion_most_abundant_masses[::-1]

        ions = []

        for i, (ion_formula, ion_mass, most_abundant_mass) in enumerate(zip(b_ion_formulae, b_ion_masses, b_ion_most_abundant_masses)):
            for charge in range(1, max_charge):
                ion_name = f"b{i+1}"
                if ion_type == "c/z•":
                    ion_name = f"c{i+1}"
                mz = (ion_mass + charge * 1.00728) / charge
                most_abundant_mz = (most_abundant_mass + charge * 1.00728) / charge
                if min_mz <= mz <= max_mz:
                    ions.append((ion_name, ion_formula, ion_mass, charge, mz, most_abundant_mz))

        for i, (ion_formula, ion_mass, most_abundant_mass) in enumerate(zip(y_ion_formulae, y_ion_masses, y_ion_most_abundant_masses)):
            for charge in range(1, max_charge):
                ion_name = f"y{i+1}"
                if ion_type == "c/z•":
                    ion_name = f"z{i+1}"
                mz = (ion_mass + charge * 1.00728) / charge
                most_abundant_mz = (most_abundant_mass + charge * 1.00728) / charge              
                if min_mz <= mz <= max_mz:
                    ions.append((ion_name, ion_formula, ion_mass, charge, mz, most_abundant_mz))

        return ions


    def process_sequence_string(self, input_string):
        result_string = ""
        substrings_list = []

        i = 0
        while i < len(input_string):
            if input_string[i] == '[':
                closing_bracket_index = input_string.find(']', i)
                substring = input_string[i + 1:closing_bracket_index]
                substrings_list[-1] = substring
                i = closing_bracket_index + 1
            else:
                result_string += input_string[i]
                substrings_list.append(None)
                i += 1

        result_string = result_string.replace(" ", "")

        return result_string, substrings_list


    def mod_string_to_dict(self, mod):
        mod_formula = {
            'C': 0, 
            'H': 0,
            'N': 0,
            'O': 0,
            'S': 0,
            'P': 0,
            'Na': 0,
            'Ag': 0,
            'Al': 0,
            'As': 0,
            'Au': 0,
            'B': 0,
            'Ca': 0,
            'Cd': 0,
            'Cl': 0,
            'Co': 0,
            'Cr': 0,
            'F': 0,
            'Fe': 0,
            'Hg': 0,
            'I': 0,
            'K': 0,
            'Li': 0,
            'Mg': 0,
            'Mn': 0,
            'Mo': 0,
            'Ni': 0,
            'Pd': 0,
            'Pt': 0,
            'Ru': 0,
            'Se': 0,
            'Zn': 0
        }

        if mod is None:
            return mod_formula

        try:
            addition, loss = mod.split('-')
        except:
            if mod[0] == '-':
                addition = ''
                loss = mod[1:]
            else:
                addition = mod
                loss = ''

        addition_formulae = re.findall(r'([A-Z][a-z]*)(\d*)', addition)
        for formula in addition_formulae:
            element = formula[0]
            count = int(formula[1]) if formula[1] else 1
            mod_formula[element] += count

        subtraction_formulae = re.findall(r'([A-Z][a-z]*)(\d*)', loss)
        for formula in subtraction_formulae:
            element = formula[0]
            count = int(formula[1]) if formula[1] else 1
            mod_formula[element] -= count

        return mod_formula


    def calculate_molecular_formula(self, protein_sequence):
        amino_acid_formula = {
            'A': {'C': 3, 'H': 5, 'N': 1, 'O': 1},
            'R': {'C': 6, 'H': 12, 'N': 4, 'O': 1},
            'N': {'C': 4, 'H': 6, 'N': 2, 'O': 2},
            'D': {'C': 4, 'H': 5, 'N': 1, 'O': 3},
            'C': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1},
            'c': {'C': 3, 'H': 4, 'N': 1, 'O': 1, 'S': 1},
            'Q': {'C': 5, 'H': 8, 'N': 2, 'O': 2},
            'E': {'C': 5, 'H': 7, 'N': 1, 'O': 3},
            'G': {'C': 2, 'H': 3, 'N': 1, 'O': 1},
            'H': {'C': 6, 'H': 7, 'N': 3, 'O': 1},
            'I': {'C': 6, 'H': 11, 'N': 1, 'O': 1},
            'L': {'C': 6, 'H': 11, 'N': 1, 'O': 1},
            'K': {'C': 6, 'H': 12, 'N': 2, 'O': 1},
            'M': {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 1},
            'F': {'C': 9, 'H': 9, 'N': 1, 'O': 1},
            'P': {'C': 5, 'H': 7, 'N': 1, 'O': 1},
            'S': {'C': 3, 'H': 5, 'N': 1, 'O': 2},
            'T': {'C': 4, 'H': 7, 'N': 1, 'O': 2},
            'W': {'C': 11, 'H': 10, 'N': 2, 'O': 1},
            'Y': {'C': 9, 'H': 9, 'N': 1, 'O': 2},
            'V': {'C': 5, 'H': 9, 'N': 1, 'O': 1},
        }

        formula = {
            'C': 0, 
            'H': 0,
            'N': 0,
            'O': 0,
            'S': 0,
            'P': 0,
            'Na': 0,
            'Ag': 0,
            'Al': 0,
            'As': 0,
            'Au': 0,
            'B': 0,
            'Ca': 0,
            'Cd': 0,
            'Cl': 0,
            'Co': 0,
            'Cr': 0,
            'F': 0,
            'Fe': 0,
            'Hg': 0,
            'I': 0,
            'K': 0,
            'Li': 0,
            'Mg': 0,
            'Mn': 0,
            'Mo': 0,
            'Ni': 0,
            'Pd': 0,
            'Pt': 0,
            'Ru': 0,
            'Se': 0,
            'Zn': 0
        }

        for amino_acid in protein_sequence:
            if amino_acid in amino_acid_formula:
                for element, count in amino_acid_formula[amino_acid].items():
                    formula[element] += count

        return formula


    def load_centroid_spectrum(self, file_path, ppm_offset):
        spectrum = np.loadtxt(file_path)
        for i, x in enumerate(spectrum[:,0]):
            spectrum[:,0][i] = x + (ppm_offset / 1000000 * x)

        return spectrum


    def remove_noise_peaks(self, spectrum, signal_to_noise):
        def indices_less_than_cutoff(lst, cutoff):
            return [index for index, value in enumerate(lst) if value > cutoff]

        intensity = spectrum[:,1]
        hist, edges = np.histogram(intensity, bins=1000)
        mode_bin_index = np.argmax(hist)
        noise_level = (edges[mode_bin_index] + edges[mode_bin_index + 1]) / 2 # modal bin

        signal_cutoff = signal_to_noise * noise_level

        keep_index = indices_less_than_cutoff(intensity, signal_cutoff)

        filtered_spectrum_x = []
        filtered_spectrum_y = []
        for i in keep_index:
            filtered_spectrum_x.append(spectrum[:,0][i])
            filtered_spectrum_y.append(spectrum[:,1][i])

        return filtered_spectrum_x, filtered_spectrum_y, noise_level, signal_cutoff


    def match_ion_mz(self, peaks, theo_ions, ppm_threshold):
        def binary_search(arr, target, ppm_threshold):
            low = bisect.bisect_left(arr, target * (1 - ppm_threshold / 1e6))
            high = bisect.bisect_right(arr, target * (1 + ppm_threshold / 1e6))
            return high - low

        matched_ions = []
        for ion in theo_ions:
            mz = ion[5]
            num_matches = binary_search(peaks, mz, ppm_threshold)
            if num_matches >=1:
                matched_ions.append(ion)

        return matched_ions


    def fit_ions(self, peaks, intens, matched_ions, matching_tol, isotope_tol, min_score):
        def residual(scale_factor, list1, list2):
            scaled_list1 = np.array(list1) * scale_factor
            norm_list_1 = [(item / max(scaled_list1)) for item in scaled_list1]
            return np.sum(((scaled_list1 - list2) ** 2) * norm_list_1)

        def binary_search(arr, target, ppm_threshold):
            low = bisect.bisect_left(arr, target * (1 - ppm_threshold / 1e6))
            high = bisect.bisect_right(arr, target * (1 + ppm_threshold / 1e6))
            return range(low, high)

        fit_ions = []
        fit_ion_envelopes = []
        fit_ion_properties = []

        for ion in matched_ions:
            ion_formula = ion[1]
            charge = ion[3]

            theoretical_envelope = isotopic_variants(ion_formula, charge=charge)
            theo_x = [peak.mz for peak in theoretical_envelope]
            theo_y = [peak.intensity for peak in theoretical_envelope]

            index_to_delete = self.get_indices_below_n_percent(theo_y, 0.01)
            theo_x = [mz for i, mz in enumerate(theo_x) if i not in index_to_delete]
            theo_y = [intens for i, intens in enumerate(theo_y) if i not in index_to_delete]

            most_abundant_index = theo_y.index(max(theo_y))
            theo_mz = theo_x[most_abundant_index]
            match_range = binary_search(peaks, theo_mz, matching_tol)

            match_index = 0
            min_diff = 1
            for i in match_range:
                diff = abs(peaks[i] - theo_mz)
                if diff <= min_diff:
                    min_diff = diff
                    match_index = i
            if min_diff == 1:
                continue

            obs_mz = [0] * len(theo_x)
            obs_intens = [0] * len(theo_y)            
            obs_mz[most_abundant_index] = [peaks[match_index]]
            error = peaks[match_index] - theo_mz
            obs_intens[most_abundant_index] = [intens[match_index]]

            for index in range(len(theo_x)):
                theo_value = theo_x[index] + error
                match_range = binary_search(peaks, theo_value, isotope_tol)
                match_index = 0
                min_diff = 1
                for i in match_range:
                    diff = abs(peaks[i] - theo_value)
                    if diff <= min_diff:
                        min_diff = diff
                        match_index = i
                if min_diff == 1:
                    obs_mz[index] = theo_value
                    obs_intens[index] = 0
                else:
                    obs_mz[index] = peaks[match_index]
                    obs_intens[index] = intens[match_index]

            guess = max(obs_intens) / max(theo_y)

            result = minimize(residual, guess, args=(theo_y, obs_intens))
            ratio = result.x[0]
            theo_y = [intens * ratio for intens in theo_y]

            y_array = np.array(theo_y)
            exp_y_array = np.array(obs_intens)

            numerator = 0
            denominator = 0
            for i, theo_value in enumerate(y_array):
                exp_value = exp_y_array[i]
                if exp_value == 0:
                    ratio1 = 0
                else:
                    ratio1 = theo_value / exp_value
                ratio2 = exp_value / theo_value
                ratio = min([ratio1,ratio2])
                numerator += ratio * theo_value
                denominator += theo_value
            score = numerator / denominator

            if score >= min_score:
                fit_ions.append(ion)
                fit_ion_envelopes.append((theo_x, theo_y))
                ppm_error = error/obs_mz[0] * 1e6
                intensity = np.sum(theo_y)
                min_envelope = min(theo_x)
                max_envelope = max(theo_x)
                max_abu = max(theo_y)
                properties = (ppm_error, score, intensity, min_envelope, max_envelope, max_abu)
                fit_ion_properties.append(properties)

        return fit_ions, fit_ion_envelopes, fit_ion_properties


    def get_indices_below_n_percent(self, lst, n):
        arr = np.array(lst)
        max_value = np.max(arr)
        threshold = n / 100 * max_value
        
        indices_below_threshold = np.where(arr < threshold)[0]

        return indices_below_threshold.tolist()
