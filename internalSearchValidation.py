import bisect
import re
import os
import math
import statistics
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib import gridspec
from pyteomics import mass
from brainpy import isotopic_variants
from scipy.optimize import minimize
from tqdm import tqdm
import wx
import numpy as np
import pandas as pd

# to do; we should probably factor in the number of theoretical fragments
# considered into the calculation of E-values. Could calculate distrubution for
# proportion of theoretical ions matched, instead of absolute number? 
# Will do later. Not that much of a problem as won't lead to wrong assignmnets,
# may just miss some. Also most of the ones missed should be covered by other terminal
# sequence ions. Decrease with increase neg offset also observable for continous scan
# trx analysis -- could also be used to correct maybe but further confirms it just decreases
# sensitivity rather than leading to false positives.


class InternalSearchValidationWindow(wx.Frame):
    def __init__(
        self, parent, title, file_path, directory_path, name, sequence,
        fine_accuracy, auto_accuracy, min_score, e_value#, ion_type
    ):
        super().__init__(parent, title=title, size=(800, 800))

        self.file_path = file_path
        self.directory_path = directory_path
        self.basename = os.path.basename(file_path).replace(".txt", "")
        self.calibrated_spectrum_file = os.path.join(
            self.directory_path,
            f"{self.basename}.recalibratedSpectrum.txt"
        )
        self.peak_assignment_file = os.path.join(
            self.directory_path,
            f"{self.basename}.assignedPeakList.csv"
        )
        self.accuracy = fine_accuracy
        self.auto_accuracy = auto_accuracy
        self.sequence = sequence
        self.name = name
        self.min_score = min_score
        self.e_value = e_value

        self.assignment_df = pd.DataFrame()
        self.spectrum = np.empty((1, 1))
        self.spectrum_x = []
        self.spectrum_y = []
        self.iteration = 0
        self.stage = ""
        self.length_match_list = 0
        self.observed_ion = 0
        self.ion_name = ""
        self.frag_site1 = ""
        self.frag_site2 = ""
        self.adduct = ""
        self.loss = ""
        self.total_adduct = ""
        self.total_loss = ""
        self.frag_site = ""
        self.formula_string = ""
        self.ppm_error = 0.0
        self.score = ""
        self.mw = 0.0
        self.scatter_x = []
        self.scatter_y = []
        self.total_intensity = 0.0
        self.charge_scaled_abundance = 0.0
        self.peak_list = []
        self.matched_theo_loss_ions = []
        self.matched_observed_loss_ions = []
        column_names = [
            "theo_mass",
            "theo_start_index",
            "theo_end_index",
            "obs_index",
            "ion_type",
            "charge",
            "obs_mass"
        ]
        self.matched_envelopes = pd.DataFrame(columns=column_names)

        #screen width variable to allow for smaller screens
        screen_width, screen_height = wx.DisplaySize()
        self.SetSize((screen_width // 2, screen_height // 2))

        self.panel = wx.Panel(self)
        self.plot_panel = ValidationPlotsPanel(self.panel)

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # instructions for use
        instruction_text = wx.StaticText(
            self.panel,
            label="Enter 't' for true match, 'f'' for false match:"
        )
        instruction_text.SetFont(
            wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD)
        )

        # hidden buttons for screening peaks and keyboard shortcuts using accelerator table
        self.true_button = wx.Button(self.panel)
        self.true_button.Bind(wx.EVT_BUTTON, self.on_true_button)
        self.true_button.Hide()
        self.false_button = wx.Button(self.panel)
        self.false_button.Bind(wx.EVT_BUTTON, self.on_false_button)
        self.false_button.Hide()
        accelerator_entries = [(wx.ACCEL_NORMAL, ord("t"), self.true_button.GetId()),
                               (wx.ACCEL_NORMAL, ord("f"), self.false_button.GetId())]
        accelerator_table = wx.AcceleratorTable(accelerator_entries)
        self.SetAcceleratorTable(accelerator_table)

        main_sizer.Add(instruction_text, 0, wx.EXPAND | wx.ALL, border=10)
        main_sizer.Add(self.plot_panel, 1, wx.EXPAND | wx.TOP, border=-10)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)
        self.tested_ions = []
        self.internal_ion_search()

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Show()


    def on_size(self, event):
        size = self.panel.GetClientSize()
        self.plot_panel.SetSize(size.width // 2, size.height)
        self.panel.Layout()
        event.Skip()


    def internal_ion_search(self):
        self.assignment_df = pd.read_csv(self.peak_assignment_file)
        self.spectrum = np.loadtxt(self.calibrated_spectrum_file)
        self.spectrum_x = self.spectrum[:, 0]
        self.spectrum_x_array = np.array(self.spectrum_x)
        self.spectrum_y = self.spectrum[:, 1]

        print("Calculating background matching...")
        expected = self.calculate_background()
        print(f"Mean background matching = {expected}")
        print("Identifying sets of internal fragments...")
        n_term_set, c_term_set = self.identify_sets(expected)
        if (len(self.sequence)) in c_term_set:
            c_term_set.remove(len(self.sequence))
        if 0 in n_term_set:
            n_term_set.remove(0)
        if len(n_term_set) == 0 and len(c_term_set) == 0:
            print("No significant internal fragments detected.")
            self.Close()
            return
        self.construct_matched_df(n_term_set, c_term_set)

        # intitiate matching scheme
        self.iteration = 0
        self.stage = "internal"
        self.length_match_list = len(self.matched_envelopes)
        # prevent double testing
        for tested_ion in self.matched_envelopes["obs_index"].to_list():
            self.tested_ions.append(tested_ion)
        if self.length_match_list >= 1:
            self.iterate_over_peaks(self.iteration, self.matched_envelopes, self.stage)
        else:
            self.neutral_loss_search()


    def calculate_background(self):
        unassigned_peaks = self.assignment_df[pd.isna(self.assignment_df['name'])]
        unassigned_peak_list = unassigned_peaks["monoisotopic_mw"].tolist()
        theo_b_ions, _ = self.theoretical_mass_generator(self.sequence)
        expected = self.background_matching(unassigned_peak_list, theo_b_ions)
        return expected


    def identify_sets(self, expected):
        unassigned_peaks = self.assignment_df[pd.isna(self.assignment_df['name'])]
        unassigned_peak_list = unassigned_peaks["monoisotopic_mw"].tolist()

        n_term = []
        for j in range(len(self.sequence)):
            test_sequence = self.sequence[j:]
            theo_b_ions, _ = self.theoretical_mass_generator(test_sequence)
            num_matches = self.count_matched_ions(unassigned_peak_list, theo_b_ions)
            n_term.append(num_matches)

        n_term_array = np.array(n_term)
        n_term_p_score = np.zeros_like(n_term_array, dtype=float)
        for i, match_count in enumerate(n_term_array):
            n_term_p_score[i] = (
                (expected ** match_count) *
                math.exp(-1 * expected) /
                math.factorial(match_count)
            )
        n_term_e_value = n_term_p_score * len(n_term_p_score)
        for i, value in enumerate(n_term_e_value):
            n_term_e_value[i] = self.neg_log_and_round(value, 1)
        n_term_sets = self.find_sig_indices(n_term_e_value, self.e_value)

        c_term = []
        for j in range(len(self.sequence)):
            test_sequence = self.sequence[:j+1]
            _, theo_y_ions = self.theoretical_mass_generator(test_sequence)
            theo_y_ions = [ion_mass - 18.0105646863 for ion_mass in theo_y_ions]
            num_matches = self.count_matched_ions(unassigned_peak_list, theo_y_ions)
            c_term.append(num_matches)

        c_term_array = np.array(c_term)
        c_term_p_score = np.zeros_like(c_term_array, dtype=float)
        for i, match_count in enumerate(c_term_array):
            c_term_p_score[i] = (
                (expected ** match_count) *
                math.exp(-1 * expected) /
                math.factorial(match_count)
            )
        c_term_e_value = c_term_p_score * len(c_term_p_score)
        for i, value in enumerate(c_term_e_value):
            c_term_e_value[i] = self.neg_log_and_round(value, 1)

        c_term_e_value = c_term_e_value[:-1]

        c_term_sets = self.find_sig_indices(c_term_e_value, self.e_value)

        return n_term_sets, c_term_sets


    def construct_matched_df(self, n_term_set, c_term_set):
        peak_list  = self.assignment_df["monoisotopic_mw"].tolist()
        unassigned_peaks = self.assignment_df[pd.isna(self.assignment_df['name'])]
        unassigned_peak_list = unassigned_peaks["monoisotopic_mw"].tolist()

        row_list = []
        for n_term_res in n_term_set:
            theo_ions = self.gen_internal_fragments(self.sequence, n_term_res, True)
            theo_masses = [tup[0] for tup in theo_ions]
            matched_theo, matched_obs = self.report_matches(theo_masses, peak_list, self.accuracy)
            for theo_index, obs_index in zip(matched_theo, matched_obs):
                theo_ion = theo_ions[theo_index]
                theo_mass = theo_ion[0]
                theo_start_index = theo_ion[1]
                theo_end_index = theo_ion[2]
                ion_type = "I"
                charge = self.assignment_df["charge"][obs_index]
                obs_mass = self.assignment_df["monoisotopic_mw"][obs_index]

                new_row = {
                    "theo_mass": theo_mass,
                    "theo_start_index": theo_start_index,
                    "theo_end_index": theo_end_index,
                    "obs_index": obs_index,
                    "ion_type": ion_type,
                    "charge": charge,
                    "obs_mass": obs_mass
                }

                row_list.append(new_row)

        for c_term_res in c_term_set:
            theo_ions = self.gen_internal_fragments(self.sequence, c_term_res, False)
            theo_masses = [tup[0] for tup in theo_ions]
            matched_theo, matched_obs = self.report_matches(theo_masses, peak_list, self.accuracy)
            for theo_index, obs_index in zip(matched_theo, matched_obs):
                theo_ion = theo_ions[theo_index]
                theo_mass = theo_ion[0]
                theo_start_index = theo_ion[1]
                theo_end_index = theo_ion[2]
                ion_type = "I"
                charge = self.assignment_df["charge"][obs_index]
                obs_mass = self.assignment_df["monoisotopic_mw"][obs_index]

                new_row = {
                    "theo_mass": theo_mass,
                    "theo_start_index": theo_start_index,
                    "theo_end_index": theo_end_index,
                    "obs_index": obs_index,
                    "ion_type": ion_type,
                    "charge": charge,
                    "obs_mass": obs_mass
                }

                row_list.append(new_row)

        self.matched_envelopes = pd.DataFrame(row_list)
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)

        assigned_peaks = self.assignment_df[self.assignment_df["ion"].notna()].index
        self.matched_envelopes = (
            self.matched_envelopes[~self.matched_envelopes['obs_index'].isin(assigned_peaks)]
        )
        self.matched_envelopes = self.matched_envelopes.reset_index(drop=True)


    def iterate_over_peaks(self, i, matched_peaks, stage):
        def residual(scale_factor, list1, list2):
            scaled_list1 = np.array(list1) * scale_factor
            norm_list_1 = [(item / max(scaled_list1)) for item in scaled_list1]
            return np.sum(((scaled_list1 - list2) ** 2) * norm_list_1)

        self.observed_ion = matched_peaks["obs_index"][i]
        charge = int(self.assignment_df["charge"][self.observed_ion])

        if stage == "internal":
            ion_type = matched_peaks["ion_type"][i]
            theo_start_index = int(matched_peaks["theo_start_index"][i])
            theo_end_index = int(matched_peaks["theo_end_index"][i])
            if ion_type == "I":
                ion_sequence = self.sequence[theo_start_index:theo_end_index]
                ion_formula = self.calculate_molecular_formula(ion_sequence)
                self.ion_name = f"{ion_type}{theo_start_index}-{theo_end_index} ({charge}+)"
                n_site1 = self.sequence[theo_start_index-1]
                c_site1 = self.sequence[theo_start_index]
                self.frag_site1 = n_site1 + c_site1
                n_site2 = self.sequence[theo_end_index-1]
                c_site2 = self.sequence[theo_end_index]
                self.frag_site2 = n_site2 + c_site2

            self.adduct = ""
            self.loss = ""
            self.total_adduct = ""
            self.total_loss = ""

        elif stage == "loss":
            self.ion_name = matched_peaks["original_ion_name"][i]
            self.frag_site = matched_peaks["frag_site"][i]
            ion_formula = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0, 'P': 0, 'Na': 0}
            original_formula = self.parse_molecular_formula(
                matched_peaks["original_ion_formula"][i]
            )
            adduct_formula = self.parse_molecular_formula(matched_peaks["adduct"][i])
            loss_formula = self.parse_molecular_formula(matched_peaks["loss"][i])
            for element, count in original_formula.items():
                ion_formula[element] += count
            for element, count in adduct_formula.items():
                ion_formula[element] += count
            for element, count in loss_formula.items():
                ion_formula[element] += -1 * count

            total_adduct_formula = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0, 'P': 0, 'Na': 0}
            if type(matched_peaks["original_ion_adduct"][i]) == str:
                original_adduct_formula = self.parse_molecular_formula(
                    matched_peaks["original_ion_adduct"][i]
                )
                for element, count in original_adduct_formula.items():
                    total_adduct_formula[element] += count
            for element, count in adduct_formula.items():
                total_adduct_formula[element] += count
            total_loss_formula = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0, 'P': 0, 'Na': 0}
            if type(matched_peaks["original_ion_loss"][i]) == str:
                original_loss_formula = self.parse_molecular_formula(
                    matched_peaks["original_ion_loss"][i]
                )
                for element, count in original_loss_formula.items():
                    total_loss_formula[element] += count
            for element, count in loss_formula.items():
                total_loss_formula[element] += count
            self.adduct = self.formula_to_string(adduct_formula)
            self.loss = self.formula_to_string(loss_formula)
            self.total_adduct = self.formula_to_string(total_adduct_formula)
            self.total_loss = self.formula_to_string(total_loss_formula)

        self.formula_string = self.formula_to_string(ion_formula)

        mz = self.assignment_df["monoisotopic_mz"][self.observed_ion]
        mw = self.assignment_df ["monoisotopic_mw"][self.observed_ion]

        theoretical_envelope = isotopic_variants(ion_formula, charge=charge)
        x = [peak.mz for peak in theoretical_envelope]
        y = [peak.intensity for peak in theoretical_envelope]

        error = mz - x[0]
        self.ppm_error = error / mz * 1000000

        self.mw = min(x)

        index_to_delete = self.get_indices_below_n_percent(y, 0.002)
        x = [mz for i, mz in enumerate(x) if i not in index_to_delete]
        y = [intens for i, intens in enumerate(y) if i not in index_to_delete]

        calibrated_x = [j + error for j in x]
        exp_x = []
        exp_y = []

        for peak_mz in calibrated_x:
            closest_index = self.find_closest_index(peak_mz, self.spectrum_x_array)
            exp_x.append(self.spectrum_x[closest_index])
            exp_y.append(self.spectrum_y[closest_index])

        guess = max(exp_y) / max(y)

        result = minimize(residual, guess, args=(y, exp_y))
        ratio = result.x[0]
        y = [intens * ratio for intens in y]

        theo_y_array = np.array(y)
        exp_y_array = np.array(exp_y)

        numerator = 0
        denominator = 0

        for j, theo_value in enumerate(theo_y_array):
            exp_value = exp_y_array[j]
            if exp_value == 0:
                ratio1 = 0
            else:
                ratio1 = theo_value / exp_value
            ratio2 = exp_value / theo_value
            ratio = min([ratio1,ratio2])
            numerator += ratio * theo_value
            denominator += theo_value
        self.score = numerator / denominator


        self.scatter_x = x
        self.scatter_y = y
        self.total_intensity = sum(y)
        self.charge_scaled_abundance = self.total_intensity / charge

        # mass and intensity ranges
        abu = max(y) * 1.2
        selected_min_mz = min(x) - 1.5 * (x[1] - x[0])
        selected_max_mz = max(x) + 1.5 * (x[1] - x[0])

        peak_params = [
            self.score,
            self.ppm_error
        ]
        peak_params = [round(value, 3) for value in peak_params]

        peak_params.append(self.ion_name)
        peak_params.append(self.total_adduct)
        peak_params.append(self.total_loss)

        prev_ion = self.assignment_df["ion"][self.observed_ion]

        if type(prev_ion) == str:
            old_score = self.assignment_df["fit_score"][self.observed_ion]
            new_score = self.score
            new_frag1 = self.frag_site1
            new_frag2 = self.frag_site2
            new_prop = (
                self.calculate_propensity(new_frag1[0], new_frag1[1]) *
                self.calculate_propensity(new_frag2[0], new_frag2[1])
            )
            new_prop = round(math.log10(new_prop), 1)
            old_frag1 = self.assignment_df["frag_site"][self.observed_ion].split(';')[0]
            old_frag2 = self.assignment_df["frag_site"][self.observed_ion].split(';')[1]
            old_prop = (
                self.calculate_propensity(old_frag1[0], old_frag1[1]) *
                self.calculate_propensity(old_frag2[0], old_frag2[1])
            )
            old_prop = round(math.log10(old_prop), 1)
            if new_score >= old_score and new_prop >= old_prop:
                self.on_true_button(None)
                return
            elif new_score <= old_score and new_prop < old_prop:#
                self.on_false_button(None)
                return
            else:
                print("!!!")
                print(f"ATTENTION: Fragment has already been assigned as {prev_ion}")
                print(f"Old score: {round(old_score,2)} New score: {round(new_score,2)}")
                print(f"Old propensity score: {old_prop} New propensity score: {new_prop}")
                print("Propensity scores closer to 0 are more likely.")
                print("!!!")
        else:
            pass

        if (self.score >= self.min_score and
            abs(self.ppm_error) <= self.auto_accuracy and
            type(prev_ion) != str):
            self.plot_panel.ax_main.clear()
            self.plot_panel.plot_main_validation_spectrum(
                self.spectrum,
                selected_min_mz,
                selected_max_mz,
                abu,
                mw,
                x,
                y,
                peak_params
            )
            self.plot_panel.canvas.draw()
            self.on_true_button(None)
            return
        else:
            self.plot_panel.ax_main.clear()
            self.plot_panel.plot_main_validation_spectrum(
                self.spectrum,
                selected_min_mz,
                selected_max_mz,
                abu,
                mw,
                x,
                y,
                peak_params
            )
            self.plot_panel.canvas.draw()


    def neutral_loss_search(self):
        self.assignment_df = pd.read_csv(self.peak_assignment_file)
        self.peak_list  = self.assignment_df["monoisotopic_mw"].tolist()
        self.spectrum = np.loadtxt(self.calibrated_spectrum_file)
        self.spectrum_x = self.spectrum[:, 0]
        self.spectrum_y = self.spectrum[:, 1]

        assigned_peaks = self.assignment_df[self.assignment_df["ion"].notna()].index

        obs_match_list = []
        original_ion_match_list = []
        original_ion_name_list = []
        original_formula_list = []
        original_adduct_list = []
        original_loss_list = []
        loss_match_list = []
        add_match_list = []
        charge_list = []
        original_frag_site_list = []

        for i in assigned_peaks:
            original_protein = self.assignment_df["name"][i]
            original_charge = self.assignment_df["charge"][i]
            original_name = self.assignment_df["ion"][i]
            original_formula = self.assignment_df["fitter_formula"][i]
            original_adduct = self.assignment_df["adduct"][i]
            original_loss = self.assignment_df["loss"][i]
            original_frag_site = self.assignment_df["frag_site"][i]
            matched_peak = self.peak_list[i]
            mod_name = ["H2O", "NH3", "Na"]
            mods = [-18.010565, -17.026549, 21.981945]
            theo_ions = [matched_peak + adduct for adduct in mods]

            if original_name[0] == "I" and original_protein == self.name:
                self.matched_theo_loss_ions, self.matched_observed_loss_ions = self.report_matches(
                    theo_ions, self.peak_list, self.accuracy*1.4142
                )
                for theo, obs in zip(self.matched_theo_loss_ions, self.matched_observed_loss_ions):
                    charge = self.assignment_df["charge"][obs]
                    if charge == original_charge:
                        obs_match_list.append(obs)
                        original_ion_match_list.append(i)
                        original_ion_name_list.append(original_name)
                        original_formula_list.append(original_formula)
                        original_adduct_list.append(original_adduct)
                        original_loss_list.append(original_loss)
                        original_frag_site_list.append(original_frag_site)
                        if mod_name[theo] in ["H2O", "NH3"]:
                            loss_match_list.append(mod_name[theo])
                            add_match_list.append("")
                        elif mod_name[theo] in ["Na"]:
                            loss_match_list.append("H")
                            add_match_list.append(mod_name[theo])
                        charge_list.append(charge)

        self.matched_envelopes = pd.DataFrame()
        self.matched_envelopes["obs_index"] = obs_match_list
        self.matched_envelopes["original_ion_index"] = original_ion_match_list
        self.matched_envelopes["original_ion_name"] = original_ion_name_list
        self.matched_envelopes["original_ion_formula"] = original_formula_list
        self.matched_envelopes["original_ion_adduct"] = original_adduct_list
        self.matched_envelopes["original_ion_loss"] =  original_loss_list
        self.matched_envelopes["adduct"] = add_match_list
        self.matched_envelopes["loss"] = loss_match_list
        self.matched_envelopes["charge"] = charge_list
        self.matched_envelopes["frag_site"] = original_frag_site_list

        self.matched_envelopes = (
            self.matched_envelopes[~self.matched_envelopes['obs_index'].isin(assigned_peaks)]
        )
        self.matched_envelopes = (
            self.matched_envelopes[~self.matched_envelopes['obs_index'].isin(self.tested_ions)]
        )
        self.matched_envelopes = self.matched_envelopes.reset_index(drop=True)

        for tested_ion in self.matched_envelopes["obs_index"].to_list():
            self.tested_ions.append(tested_ion)

        self.iteration = 0
        self.stage = "loss"
        self.length_match_list = len(self.matched_envelopes)
        if len(self.matched_envelopes) >= 1:
            self.iterate_over_peaks(self.iteration, self.matched_envelopes, self.stage)


    def on_true_button(self, event):
        self.assignment_df.loc[self.observed_ion, "name"] = self.name
        self.assignment_df.loc[self.observed_ion, "sequence"] = self.sequence
        self.assignment_df.loc[self.observed_ion, "n_mod"] = "internal"
        self.assignment_df.loc[self.observed_ion, "c_mod"] = "internal"
        self.assignment_df.loc[self.observed_ion, "ion"] = self.ion_name
        self.assignment_df.loc[self.observed_ion, "loss"] = self.total_loss
        self.assignment_df.loc[self.observed_ion, "adduct"] = self.total_adduct
        self.assignment_df.loc[self.observed_ion, "ppm_error"] = self.ppm_error
        self.assignment_df.loc[self.observed_ion, "fit_score"] = self.score
        self.assignment_df.loc[self.observed_ion, "theoretical_mz"] = self.mw
        self.assignment_df.loc[self.observed_ion, "total_intensity"] = self.total_intensity
        self.assignment_df.loc[self.observed_ion, "charge_scaled_abundance"] = self.charge_scaled_abundance
        self.assignment_df.loc[self.observed_ion, "fitter_formula"] = self.formula_string
        self.assignment_df.loc[self.observed_ion, "fitter_theo_x"] = str(self.scatter_x)
        self.assignment_df.loc[self.observed_ion, "fitter_theo_y"] = str(self.scatter_y)
        if self.stage == "internal":
            self.assignment_df.loc[self.observed_ion, "frag_site"] = (
                f"{str(self.frag_site1).upper()};{str(self.frag_site2).upper()}"
            )
        elif self.stage == "loss":
            self.assignment_df.loc[self.observed_ion, "frag_site"] = self.frag_site
        print(f"Assigned {self.name} {self.ion_name} (+){self.total_adduct} (-){self.total_loss}")

        if self.iteration < len(self.matched_envelopes) - 1:
            self.iteration += 1
            self.iterate_over_peaks(self.iteration, self.matched_envelopes, self.stage)

        elif self.iteration >= len(self.matched_envelopes) - 1:
            self.assignment_df.to_csv(self.peak_assignment_file, index=False)
            print(f'Saved updated assigned peak list to {self.peak_assignment_file}')
            if self.length_match_list == 0:
                self.filter_results()
                self.Close()
            print("Continuing with next round of neutral loss scans.")
            self.neutral_loss_search()


    def on_false_button(self, event):
        if self.iteration < len(self.matched_envelopes) - 1:
            self.iteration += 1
            self.iterate_over_peaks(self.iteration, self.matched_envelopes, self.stage)
        elif self.iteration >= len(self.matched_envelopes) - 1:
            self.assignment_df.to_csv(self.peak_assignment_file, index=False)
            print(f'Saved updated assigned peak list to {self.peak_assignment_file}')
            if self.length_match_list == 0:
                self.filter_results()
                self.Close()
            print("Continuing with next round of neutral loss scans.")
            self.neutral_loss_search()


    def filter_results(self):
        assignment_list = pd.read_csv(self.peak_assignment_file)
        ppm_threshold = self.accuracy
        mask = abs(assignment_list['ppm_error']) > ppm_threshold
        assignment_list.loc[mask, "name"] = None
        assignment_list.loc[mask, "sequence"] = None
        assignment_list.loc[mask, "n_mod"] = None
        assignment_list.loc[mask, "c_mod"] = None
        assignment_list.loc[mask, "ion"] = None
        assignment_list.loc[mask, "variable_mod"] = None
        assignment_list.loc[mask, "variable_mod_gain"] = None
        assignment_list.loc[mask, "variable_mod_loss"] = None
        assignment_list.loc[mask, "adduct"] = None
        assignment_list.loc[mask, "loss"] = None
        assignment_list.loc[mask, "theoretical_mz"] = None
        assignment_list.loc[mask, "ppm_error"] = None
        assignment_list.loc[mask, "fit_score"] = None
        assignment_list.loc[mask, "total_intensity"] = None
        assignment_list.loc[mask, "charge_scaled_abundance"] = None
        assignment_list.loc[mask, "fitter_formula"] = None
        assignment_list.loc[mask, "fitter_theo_x"] = None
        assignment_list.loc[mask, "fitter_theo_y"] = None
        assignment_list.loc[mask, "frag_site"] = None
        assignment_list.to_csv(self.peak_assignment_file, index=False)
        print(f"Removed assignments with mass errors greater than {ppm_threshold} ppm.")


    def gen_internal_fragments(self, sequence, start_index, n_term=True):
        ion_masses = []

        if n_term:
            sequence = sequence[start_index:]
        else:
            sequence = sequence[:start_index+1]

        # replace c with custom aa B to enable disulfide calc
        mass.std_aa_comp["B"] = mass.Composition({'C': 3, 'H': 4, 'N': 1, 'O': 1, 'S': 1})
        sequence = sequence.replace('c', 'B')

        sequence_mass = mass.calculate_mass(sequence, ion_type="M") # intact mass
        aa_masses = np.array([mass.calculate_mass(aa, ion_type="b") for aa in sequence])

        if n_term:
            ion_masses = np.cumsum(aa_masses)
            ion_masses = ion_masses[:-1]
        else:
            ion_masses = sequence_mass - np.cumsum(aa_masses) - 18.0105646863
            ion_masses = ion_masses[:-1]
            ion_masses.sort()

        ion_start_indices = []
        ion_end_indices = []
        if n_term:
            i = start_index + 1
            for _ in ion_masses:
                ion_start_indices.append(start_index)
                ion_end_indices.append(i)
                i += 1
        else:
            i = start_index - 1
            for _ in ion_masses:
                ion_start_indices.append(i+1)
                i += -1
                ion_end_indices.append(start_index+1)

        theo_ions = list(zip(ion_masses, ion_start_indices, ion_end_indices))

        return theo_ions


    def report_matches(self, theoretical_ions, observed_ions, ppm_threshold):
        # binary search w/ ppm tolerance helper function
        def binary_search(arr, target, ppm_threshold):
            low = bisect.bisect_left(arr, target * (1 - ppm_threshold / 1e6))
            high = bisect.bisect_right(arr, target * (1 + ppm_threshold / 1e6))
            return high - low, list(range(low, high))

        matched_theoretical_ion_indices = []
        matched_observed_ion_indices = []

        for i, theo_ion in enumerate(theoretical_ions):
            num_matches, match = binary_search(observed_ions, theo_ion, ppm_threshold)
            if num_matches >= 1:
                for item in match:
                    matched_theoretical_ion_indices.append(i)
                    matched_observed_ion_indices.append(item)

        return matched_theoretical_ion_indices, matched_observed_ion_indices


    def find_sig_indices(self, arr, set_point):
        indexes = []
        for i, element in enumerate(arr):
            if element > set_point and i >= 1:
                indexes.append(i)
        return indexes


    def neg_log_and_round(self, number, dp):
        if number == 0:
            return 0.0
        neg_log = -1 * math.log10(abs(number))
        rounded = round(neg_log, dp)
        return rounded


    def theoretical_mass_generator(self, sequence):
        # lists of b and y ions
        theo_b_ions_list = []
        theo_y_ions_list = []

        # replace c with custom aa B to enable disulfide calc
        mass.std_aa_comp["B"] = mass.Composition({'C': 3, 'H': 4, 'N': 1, 'O': 1, 'S': 1})
        sequence = sequence.replace('c', 'B')

        sequence_mass = mass.calculate_mass(sequence, ion_type="M") # intact mass
        aa_masses = np.array([mass.calculate_mass(aa, ion_type="b") for aa in sequence])

        b_ion_masses = (
            np.cumsum(aa_masses)
        )
        y_ion_masses = (
            sequence_mass -
            np.cumsum(aa_masses)
        )

        theo_b_ions_list.append(b_ion_masses[:-1]) # exclude the dehydrated full protein
        theo_y_ions_list.append(y_ion_masses[:-1])

        theo_y_ions_list.sort()

        theo_b_ions = theo_b_ions_list[0]
        theo_y_ions = theo_y_ions_list[0]

        return theo_b_ions, theo_y_ions


    def background_matching(self, peak_list, theo_ions):
        scan_spacing = self.accuracy / 1000
        # scan across offset of 900 to 1000 to define background matches
        # assumes no true peaks in this reason but robust as mean is not
        # greatly affected
        num_points = int((1000 - 900) // scan_spacing)

        offset_list = np.linspace(900, 1000, num=num_points)
        peak_list = np.array(peak_list)
        theo_ions = np.array(theo_ions)

        match_count_list = []
        for offset in tqdm(offset_list, total=len(offset_list)):
            temp_theo_ions = [ion_mass + offset for ion_mass in theo_ions]
            match_count = self.count_matched_ions(peak_list, temp_theo_ions)
            match_count_list.append(match_count)

        mean_matches = statistics.mean(match_count_list)

        return mean_matches


    def count_matched_ions(self, peak_list, theo_ions):
        match_count = 0

        def binary_search(arr, target, ppm_threshold):
            low = bisect.bisect_left(arr, target * (1 - ppm_threshold / 1e6))
            high = bisect.bisect_right(arr, target * (1 + ppm_threshold / 1e6))
            return high - low

        for ion in theo_ions:
            num_matches = binary_search(peak_list, ion, self.accuracy)
            if num_matches >=1:
                match_count += 1

        return match_count


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

        formula = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0, 'P': 0, 'Na': 0}

        for amino_acid in protein_sequence:
            if amino_acid in amino_acid_formula:
                for element, count in amino_acid_formula[amino_acid].items():
                    formula[element] += count
        return formula


    def find_closest_index(self, value, array):
        closest_index = np.abs(array - value).argmin()
        return closest_index


    def formula_to_string(self, formula):
        formula_string=""
        for element, count in formula.items():
            if count >= 1:
                formula_string += element
                if count >= 2:
                    formula_string += str(count)
        return formula_string


    def get_indices_below_n_percent(self, lst, n):
        arr = np.array(lst)
        max_value = np.max(arr)
        threshold = n / 100 * max_value
        indices_below_threshold = np.where(arr < threshold)[0]
        return indices_below_threshold.tolist()


    def parse_molecular_formula(self, formula_string):
        formula_dict = {}
        matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula_string)
        for element, count in matches:
            count = int(count) if count else 1
            formula_dict[element] = count
        return formula_dict


    def calculate_propensity(self, n_res, c_res):
        n_res = n_res.upper()
        c_res = c_res.upper()
        # fragmentation propensites from ives 2020
        n_term_props = {
            "A": 6.6,
            "C": 0.1,
            "D": 37,
            "E": 6.9,
            "F": 2.3,
            "G": 2.2,
            "H": 1.2,
            "I": 5.7,
            "K": 7.4,
            "L": 8.3,
            "M": 1.3,
            "N": 1.4,
            "P": 0.7,
            "Q": 2.1,
            "R": 0.6,
            "S": 2.0,
            "T": 2.4,
            "V": 9.8,
            "Y": 0.5,
            "W": 2.1,
        }

        c_term_props = {
            "A": 5.0,
            "C": 0.9,
            "D": 5.3,
            "E": 5.1,
            "F": 2.4,
            "G": 13,
            "H": 2.2,
            "I": 4.8,
            "K": 6.4,
            "L": 6.8,
            "M": 1.2,
            "N": 2.4,
            "P": 25,
            "Q": 2.5,
            "R": 1.5,
            "S": 3.6,
            "T": 3.8,
            "V": 4.7,
            "Y": 0.7,
            "W": 2.6,
        }

        propensity = n_term_props[n_res] / 100 * c_term_props[c_res] / 100
        return propensity



class ValidationPlotsPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        self.max_intens = 0

        sizer = wx.BoxSizer(wx.VERTICAL)

        self.figure = Figure(facecolor="#f0f0f0")
        self.gs = gridspec.GridSpec(1, 1)
        self.ax_main = self.figure.add_subplot(self.gs[0, :])
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.toolbar = NavigationToolbar2Wx(self.canvas)

        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.EXPAND)

        self.SetSizer(sizer)
        self.Bind(wx.EVT_SIZE, self.on_size)
        self.Layout()
        self.toolbar.Realize()


    def plot_main_validation_spectrum(
        self, spectrum, selected_min_mz, selected_max_mz, abu,
        ion_mass, scatter_x, scatter_y, peak_params
    ):
        self.ax_main.clear()
        self.ax_main.plot(spectrum[:, 0], spectrum[:, 1], color="black", linewidth=1)
        self.ax_main.set_xlabel("m/z")
        self.ax_main.set_ylabel("Intensity")
        self.ax_main.set_facecolor("#f0f0f0")
        self.ax_main.set_xlim(selected_min_mz, selected_max_mz)
        self.max_intens = abu * 1.1
        self.ax_main.set_ylim(0, self.max_intens)

        score = peak_params[0]
        ppm_error = peak_params[1]
        ion_name = peak_params[2]
        adduct = peak_params[3]
        loss = peak_params[4]

        variant_string = ""
        if adduct != "" or loss != "":
            variant_string += " "
        if adduct != "":
            variant_string += "+"
            variant_string += adduct
            variant_string += " "
        if adduct != "" and loss != "":
            variant_string += " "
        if loss != "":
            variant_string += "-"
            variant_string += loss

        if scatter_x is not None:
            label_text = (
                f"{ion_name}{variant_string}: {round(ion_mass, 2)} \n"
                f"Score: {score} \n"
                f"Mass error: {ppm_error} ppm \n"
            )
            self.ax_main.annotate(
                label_text,
                xy=(0.995,0.99), xycoords="axes fraction", ha="right", va="top"
            )
            self.ax_main.scatter(
                scatter_x, scatter_y,
                c="darkviolet", marker="o", alpha=0.6
            )

        self.canvas.draw()


    def on_size(self, event):
        self.fit_plot_to_panel()
        event.Skip()


    def fit_plot_to_panel(self):
        size = self.GetClientSize()
        if size[0] >= 50:
            dpi = self.GetContentScaleFactor() * 100
            width = size.width / dpi
            height = size.height / dpi - 0.3
            self.figure.set_size_inches(width, height)
            self.figure.tight_layout(rect=[0, 0, 1, 1])
            self.canvas.draw()
