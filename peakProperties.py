import os
import pandas as pd
import numpy as np
from brainpy import isotopic_variants
from scipy.optimize import minimize
from scipy.stats import pearsonr, chisquare

class PeakProperties:
    def __init__(
        self,
        file_path,
        directory_path,
        s2n_threshold
    ):
        self.file_path = file_path
        self.directory_path = directory_path
        self.basename = os.path.basename(self.file_path).replace('.txt', '')
        self.input_cluster_csv = os.path.join(self.directory_path, f"{self.basename}.clusters.csv")
        self.input_centroid = os.path.join(self.directory_path, f"{self.basename}.centroid.rl.txt")
        if os.path.exists(self.input_centroid) is False:
            self.input_centroid = os.path.join(
                self.directory_path,
                f"{self.basename}.centroid.cwt.txt"
            )
        self.output = os.path.join(self.directory_path, f"{self.basename}.fullPeakList.csv")

        print("Calculating envelope properties...")
        self.calculate_prop(self.input_cluster_csv, self.input_centroid, float(s2n_threshold))


    def calculate_prop(
        self,
        input_cluster_csv,
        input_centroid,
        s2n_threshold
    ):
        peak_list = pd.read_csv(input_cluster_csv)
        spectrum = np.loadtxt(input_centroid)
        spec_mz = spectrum[:,0]
        spec_intens = spectrum[:,1]

        noise_coeff = self.fit_noise_level(spec_mz, spec_intens)

        peak_list = self.calculate_mz(peak_list)
        peak_list = self.calculate_num_charge_states(peak_list)

        for index, row in peak_list.iterrows():
            mw = row['monoisotopic_mw']
            charge = row['charge']
            mz = row['monoisotopic_mz']
            #generate the envelope
            theo_x, theo_y = self.averagine(mw, charge)

            matched_bool, matched_peaks, matched_theo = self.match_peaks(spectrum, theo_x, theo_y)

            theo_x = [peak[0] for peak in matched_theo]
            theo_y = [peak[1] for peak in matched_theo]
            exp_x = [peak[0] for peak in matched_peaks]
            exp_y = [peak[1] for peak in matched_peaks]

            diff_list = [(exp_x[i+1] - exp_x[i]) for i in range(len(theo_x) - 1)]

            # num matched peaks
            peak_list.at[index, 'matched_bool'] = str(matched_bool)
            peak_list.at[index, 'num_missing_peaks'] = matched_bool.count(0)
            peak_list.at[index, 'pc_missing_peaks'] = matched_bool.count(0) / len(matched_bool) * 100

            #abundance and log abundance and total abudnance and total log abundance
            peak_list.at[index, 'abundance'] = max(theo_y)
            peak_list.at[index, 'log_abundance'] = np.log10(max(theo_y))
            peak_list.at[index, 'total_abundance'] = np.sum(theo_y)
            peak_list.at[index, 'log_total_abundance'] = np.log10(np.sum(theo_y))

            #calculate fitting parameters
            correlation_coefficient, correlation_p_value = pearsonr(theo_y, exp_y)
            peak_list.at[index, 'correlation_coefficient'] = correlation_coefficient
            peak_list.at[index, 'correlation_p_value'] = correlation_p_value

            norm_exp_y = [value / np.sum(exp_y) for value in exp_y]
            norm_theo_y = [value / np.sum(theo_y) for value in theo_y]
            chisq_stat, chisq_p_value = chisquare(f_obs=norm_exp_y, f_exp=norm_theo_y)
            peak_list.at[index, 'chisq_stat'] = chisq_stat
            peak_list.at[index, 'chisq_p_value'] = chisq_p_value

            fit_score = self.fit_scoring(exp_y, theo_y)
            peak_list.at[index, 'fit_score'] = fit_score

            #calculate intergerence and s2n
            interference = self.calculate_interference(spec_mz, spec_intens, exp_x, exp_y)
            noise = noise_coeff[0]
            s2n = max(theo_y) / noise

            peak_list.at[index, 'interference'] = interference
            peak_list.at[index, 's2n'] = s2n
            peak_list.at[index, 'log_s2n'] = np.log10(s2n)

            # standard deviation of m/z values
            std_dev_mz = self.calculate_mz_std_dev(exp_x, theo_x)
            peak_list.at[index, 'mass_error_std'] = std_dev_mz

            peak_list.at[index, 'mass_diff'] = str(diff_list)
            peak_list.at[index, 'matched_x'] = str(exp_x)
            peak_list.at[index, 'matched_y'] = str(exp_y)
            peak_list.at[index, 'theo_x'] = str(theo_x)
            peak_list.at[index, 'theo_y'] = str(theo_y)

            peak_list.at[index, 'validated'] = False

        column_order = [
            'monoisotopic_mw',
            'charge',
            'monoisotopic_mz',
            'abundance',
            'log_abundance',
            'total_abundance',
            'log_total_abundance',
            's2n',
            'log_s2n',
            'fit_score',
            'interference',
            'correlation_coefficient',
            'correlation_p_value',
            'chisq_stat',
            'chisq_p_value',
            'mass_error_std',
            'num_charge_states',
            'num_missing_peaks',
            'pc_missing_peaks',
            'detected_thrash_0_1?',
            'detected_thrash_0_2?',
            'detected_thrash_0_3?',
            'detected_thrash_0_4?',
            'detected_msdeconv?',
            'detected_envcnn?',
            'matched_bool',
            'mass_diff',
            'matched_x',
            'matched_y',
            'theo_x',
            'theo_y',
            'validated'
        ]
        
        peak_list = peak_list.reindex(columns=column_order)
        peak_list = peak_list.sort_values(by='abundance', ascending=False)
        peak_list = peak_list[peak_list['s2n'] >= s2n_threshold] # filter by the s2n
        peak_list.to_csv(
            self.output,
            index=False
        )


    def fit_noise_level(self, x_array, y_array):
        peaks = list(zip(x_array, y_array))

        num_splits = 50
        split_size = len(peaks) // num_splits
        splits = [peaks[i:i + split_size] for i in range(0, len(peaks), split_size)]

        average_mz_list = []
        modal_intens_list = []

        for split in splits:
            average_mz = np.mean([t[0] for t in split])
            intensities = [t[1] for t in split]
            hist, edges = np.histogram(intensities, bins=500)
            mode_bin_index = np.argmax(hist)
            mode_value = (edges[mode_bin_index] + edges[mode_bin_index + 1]) / 2

            average_mz_list.append(average_mz)
            modal_intens_list.append(mode_value)

        coefficients = np.polyfit(average_mz_list, modal_intens_list, 0)
        a = coefficients

        return [a]


    def calculate_mz(self, peak_list):
        peak_list['monoisotopic_mz'] = (
            (peak_list['monoisotopic_mw'] +
             peak_list['charge'] * 1.007276466812) /
            peak_list['charge']
        )

        return peak_list


    def calculate_num_charge_states(self, peak_list):
        ppm_tolerance = 5
        ppm_diff_matrix = np.abs(
            (
                peak_list['monoisotopic_mw'].values[:, np.newaxis] -
                peak_list['monoisotopic_mw'].values
            )
            / peak_list['monoisotopic_mw'].values[:, np.newaxis]
            * 1e6
        )
        same_mw_mask = ppm_diff_matrix <= ppm_tolerance
        diff_charge_counts = np.sum(
            same_mw_mask &
            (peak_list['charge'].values[:, np.newaxis] != peak_list['charge'].values),
            axis=1
        ) + 1
        peak_list['num_charge_states'] = diff_charge_counts

        return peak_list


    def match_peaks(self, spectrum, theo_x, theo_y):
        matched_peaks = [[0, 0] for _ in range(len(theo_x))]
        matched_bool = []
        for t, peak in enumerate(theo_x):
            error_ppm = (peak - spectrum[:, 0]) / peak * 1e6
            matching_indices = np.where(np.abs(error_ppm) <= 5)[0]
            if len(matching_indices) > 0:
                matched_bool.append(1)
                matched_index = matching_indices[np.argmin(np.abs(error_ppm[matching_indices]))]
                matched_peaks[t] = [spectrum[matched_index, 0], spectrum[matched_index, 1]]
            else:
                matched_bool.append(0)
                matched_peaks[t] = [peak, 0]

        matched_intens = [peak[1] for peak in matched_peaks]

        theo_y = self.fit_intens(theo_y, matched_intens)

        matched_theo = list(zip(theo_x, theo_y))

        return matched_bool, matched_peaks, matched_theo


    def fit_intens(self, theo_y, exp_y):
        def residual(scale_factor, list1, list2):
            scaled_list1 = np.array(list1) * scale_factor
            norm_list_1 = [(item / max(scaled_list1)) for item in scaled_list1]
            return np.sum(((scaled_list1 - list2) ** 2) * norm_list_1)

        guess = max(exp_y) / max(theo_y)
        result = minimize(residual, guess, args=(theo_y, exp_y))

        ratio = result.x[0]
        theo_y = [intens * ratio for intens in theo_y]

        return theo_y


    def fit_scoring(self, exp_y, theo_y):
        theo_y_array = np.array(theo_y)
        exp_y_array = np.array(exp_y)

        numerator = 0
        denominator = 0

        for i, theo_value in enumerate(theo_y_array):
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

        return score


    def calculate_interference(self, spec_x, spec_y, exp_x, exp_y):
        total_target_signal = np.sum(exp_y)

        low_mz = min(exp_x)
        high_mz = max(exp_x)
        indices_within_range = [index for index, mz in enumerate(spec_x) if low_mz <= mz <= high_mz]

        total_signal = 0
        for i in indices_within_range:
            total_signal += spec_y[i]

        interference = 1 - total_target_signal / total_signal

        return interference


    def calculate_mz_std_dev(self, exp_x, theo_x):
        diff_list = []

        for i, exp in enumerate(exp_x):
            theo = theo_x[i]
            diff_list.append((exp-theo) / exp * 1e6)

        std_dev = np.std(diff_list)

        return std_dev


    def averagine(self, monoisotopic_mass, charge, cutoff=0.05):
        num_averagine = monoisotopic_mass / 111.0543
        num_c = int(round(4.9384 * num_averagine))
        num_h = int(round(7.7583 * num_averagine))
        num_n = int(round(1.3577 * num_averagine))
        num_o = int(round(1.4773 * num_averagine))
        num_s = int(round(0.0417 * num_averagine))

        polyaveragine_monoisotopic_mass = (
            num_c * 12 +
            num_h * 1.007825035 +
            num_n * 14.003074 +
            num_o * 15.99491463 +
            num_s * 31.9720707
        )
        h_mass_compensation = round((
            monoisotopic_mass - polyaveragine_monoisotopic_mass) / 1.007825035
        )
        num_h += h_mass_compensation

        ion = {'C': num_c, 'H': num_h, 'N': num_n, 'O': num_o, 'S': num_s}
        theoretical_isotopic_cluster = isotopic_variants(ion, charge=charge)

        x = [peak.mz for peak in theoretical_isotopic_cluster]
        y = [peak.intensity for peak in theoretical_isotopic_cluster]

        theoretical_monoisotopic_mz = x[0]
        mz_shift = theoretical_monoisotopic_mz - (monoisotopic_mass + charge * 1.0072765) / charge
        x2 = [peak - mz_shift for peak in x]
        y2 = list(y)

        filtered_x2 = []
        filtered_y2 = []
        for i, y_value in enumerate(y2):
            x_value = x2[i]
            if y_value >= cutoff * max(y2):
                filtered_x2.append(x_value)
                filtered_y2.append(y_value)

        return filtered_x2, filtered_y2
