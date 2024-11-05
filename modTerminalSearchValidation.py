import bisect
import re
import os
import math
import statistics
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
from matplotlib import gridspec
from pyteomics import mass
from brainpy import isotopic_variants
from scipy.optimize import minimize
import wx
import numpy as np
import pandas as pd



class ModTerminalSearchValidationWindow(wx.Frame):
    def __init__(
        self, parent, title, file_path, directory_path, name, sequence,
        n_mod_add, n_mod_sub, c_mod_add, c_mod_sub, mod_details,
        accuracy, auto_accuracy, min_score, ion_type
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
        self.accuracy = accuracy
        self.auto_accuracy = auto_accuracy
        self.sequence = sequence
        self.name = name
        self.n_mod_add = n_mod_add
        self.n_mod_sub = n_mod_sub
        self.c_mod_add = c_mod_add
        self.c_mod_sub = c_mod_sub
        self.min_score = min_score
        self.ion_type = ion_type
        self.mod_name = mod_details[0]
        self.mod_add_formula = mod_details[1]
        self.mod_loss_formula = mod_details[2]
        self.mod_mass_shift = mod_details[3]

        print(f"Considering {self.mod_name}...")

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
        self.terminal_ion_search()

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Show()


    def on_size(self, event):
        size = self.panel.GetClientSize()
        self.plot_panel.SetSize(size.width // 2, size.height)
        self.panel.Layout()
        event.Skip()


    def terminal_ion_search(self):
        self.assignment_df = pd.read_csv(self.peak_assignment_file)
        self.peak_list  = self.assignment_df["monoisotopic_mw"].tolist()
        self.spectrum = np.loadtxt(self.calibrated_spectrum_file)
        self.spectrum_x = self.spectrum[:, 0]
        self.spectrum_x_array = np.array(self.spectrum_x)
        self.spectrum_y = self.spectrum[:, 1]

        self.construct_matched_df()

        # intitiate matching scheme
        self.iteration = 0
        self.stage = "terminal"
        self.length_match_list = len(self.matched_peaks)
        # prevent double testing
        for tested_ion in self.matched_peaks["observed_ion_index"].to_list():
            self.tested_ions.append(tested_ion)
        if len(self.matched_peaks) >= 1:
            self.iterate_over_peaks(self.iteration, self.matched_peaks, self.stage)
        else:
            self.neutral_loss_search()


    def construct_matched_df(self):
        n_mod_mass = (
            mass.calculate_mass(formula = self.n_mod_add) -
            mass.calculate_mass(formula = self.n_mod_sub) +
            self.mod_mass_shift
        )
        c_mod_mass = (
            mass.calculate_mass(formula = self.c_mod_add) -
            mass.calculate_mass(formula = self.c_mod_sub) +
            self.mod_mass_shift
        )
        theo_b_ions, theo_y_ions = self.theoretical_mass_generator(
            self.sequence,
            n_mod_mass,
            c_mod_mass,
            self.ion_type
        )
        self.matched_theoretical_b_ions, self.matched_observed_b_ions = (
            self.report_matches_for_protein(theo_b_ions[0], self.peak_list, self.accuracy)
        )
        self.matched_theoretical_y_ions, self.matched_observed_y_ions = (
            self.report_matches_for_protein(theo_y_ions[0], self.peak_list, self.accuracy)
        )


        temp_theo_match_list = []
        temp_obs_match_list = []
        temp_charge_list = []
        for theo, obs in zip(self.matched_theoretical_b_ions, self.matched_observed_b_ions):
            for obs_individual in obs:
                temp_theo_match_list.append(theo)
                temp_obs_match_list.append(obs_individual)
                temp_charge_list.append(self.assignment_df["charge"][obs_individual])

        matched_b_peaks = pd.DataFrame()
        matched_b_peaks["theoretical_ion_index"] = temp_theo_match_list
        matched_b_peaks["observed_ion_index"] = temp_obs_match_list
        if self.ion_type == "b/y":
            matched_b_peaks["ion_type"] = "b"
        elif self.ion_type == "c/z•":
            matched_b_peaks["ion_type"] = "c"
        matched_b_peaks["charge"] = temp_charge_list

        temp_theo_match_list = []
        temp_obs_match_list = []
        temp_charge_list = []
        for theo, obs in zip(self.matched_theoretical_y_ions, self.matched_observed_y_ions):
            for obs_individual in obs:
                temp_theo_match_list.append(theo)
                temp_obs_match_list.append(obs_individual)
                temp_charge_list.append(self.assignment_df["charge"][obs_individual])

        matched_y_peaks = pd.DataFrame()
        matched_y_peaks["theoretical_ion_index"] = temp_theo_match_list
        matched_y_peaks["observed_ion_index"] = temp_obs_match_list
        if self.ion_type == "b/y":
            matched_y_peaks["ion_type"] = "y"
        elif self.ion_type == "c/z•":
            matched_y_peaks["ion_type"] = "z"
        matched_y_peaks["charge"] = temp_charge_list

        self.matched_peaks = pd.concat([matched_b_peaks, matched_y_peaks], ignore_index=True)

        assigned_peaks = self.assignment_df[self.assignment_df["ion"].notna()].index
        self.matched_peaks = (
            self.matched_peaks[~self.matched_peaks['observed_ion_index'].isin(assigned_peaks)]
        )
        self.matched_peaks= self.matched_peaks.reset_index(drop=True)


    def iterate_over_peaks(self, i, matched_peaks, stage):
        def residual(scale_factor, list1, list2):
            scaled_list1 = np.array(list1) * scale_factor
            norm_list_1 = [(item / max(scaled_list1)) for item in scaled_list1]
            return np.sum(((scaled_list1 - list2) ** 2) * norm_list_1)

        self.observed_ion = matched_peaks["observed_ion_index"][i]
        charge = int(self.assignment_df["charge"][self.observed_ion])

        if stage == "terminal":
            ion_type = matched_peaks["ion_type"][i]
            theo_ion = int(matched_peaks["theoretical_ion_index"][i])

            if ion_type == "b":
                term_mod_add_formula = self.parse_molecular_formula(self.n_mod_add)
                term_mod_sub_formula = self.parse_molecular_formula(self.n_mod_sub)
                ion_sequence = self.sequence[0:theo_ion+1]
                ion_formula = self.calculate_molecular_formula(ion_sequence)
                self.ion_name = f"{ion_type}{theo_ion+1} ({charge}+)"

                n_site = self.sequence[theo_ion]
                c_site = self.sequence[theo_ion+1]
                self.frag_site = n_site + c_site

            elif ion_type == "c":
                term_mod_add_formula = self.parse_molecular_formula(self.n_mod_add)
                term_mod_sub_formula = self.parse_molecular_formula(self.n_mod_sub)
                ion_sequence = self.sequence[0:theo_ion+1]
                ion_formula = self.calculate_molecular_formula(ion_sequence)
                adjust = {'H': 3, 'N': 1}
                for element, count in adjust.items():
                    ion_formula[element] += count
                self.ion_name = f"{ion_type}{theo_ion+1} ({charge}+)"

                n_site = self.sequence[theo_ion]
                c_site = self.sequence[theo_ion+1]
                self.frag_site = n_site + c_site

            elif ion_type == "y":
                term_mod_add_formula = self.parse_molecular_formula(self.c_mod_add)
                term_mod_sub_formula = self.parse_molecular_formula(self.c_mod_sub)
                ion_sequence = self.sequence[theo_ion+1:]
                ion_formula = self.calculate_molecular_formula(ion_sequence)
                adjust = {'H': 2, 'O': 1}
                for element, count in adjust.items():
                    ion_formula[element] += count
                self.ion_name = f"{ion_type}{len(self.sequence) - theo_ion - 1} ({charge}+)"

                n_site = self.sequence[theo_ion]
                c_site = self.sequence[theo_ion+1]
                self.frag_site = n_site + c_site

            elif ion_type == "z":
                term_mod_add_formula = self.parse_molecular_formula(self.c_mod_add)
                term_mod_sub_formula = self.parse_molecular_formula(self.c_mod_sub)
                ion_sequence = self.sequence[theo_ion+1:]
                ion_formula = self.calculate_molecular_formula(ion_sequence)
                adjust = {'O': 1, 'N': -1}
                for element, count in adjust.items():
                    ion_formula[element] += count
                self.ion_name = f"{ion_type}{len(self.sequence) - theo_ion - 1} ({charge}+)"

                n_site = self.sequence[theo_ion]
                c_site = self.sequence[theo_ion+1]
                self.frag_site = n_site + c_site

            for element, count in term_mod_add_formula.items():
                ion_formula[element] += count
            for element, count in term_mod_sub_formula.items():
                ion_formula[element] += -1 * count

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
        x = [peak.mz + (self.mod_mass_shift / charge) for peak in theoretical_envelope]
        y = [peak.intensity for peak in theoretical_envelope]
        error = mz - x[0]
        self.ppm_error = error / mz * 1000000 # -1 corrects for sign error

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
            new_frag = self.frag_site
            new_prop = self.calculate_propensity(new_frag[0], new_frag[1])
            new_prop = round(math.log10(new_prop), 1)
            old_frag = self.assignment_df["frag_site"][self.observed_ion]
            old_prop = self.calculate_propensity(old_frag[0], old_frag[1])
            old_prop = round(math.log10(old_prop), 1)

            if new_score >= old_score and new_prop >= old_prop:
                self.on_true_button(None)
                return
            elif new_score <= old_score and new_prop < old_prop:
                self.on_false_button(None)
                return
            else:
                print("!!!")
                print(f"ATTENTION: Fragment has already been assigned as {prev_ion}")
                print(f"Old score: {round(old_score,2)} New score: {round(new_score,2)}")
                print(f"Old propensity score: {old_prop} New propensity score: {new_prop}")
                print("Propensity scores closer to 0 are more likely.")
                print("!!!")

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
        original_variable_mod_list = [] 
        original_variable_mod_gain_list = []
        original_variable_mod_loss_list = []

        for i in assigned_peaks:
            original_protein = self.assignment_df["name"][i]
            original_charge = self.assignment_df["charge"][i]
            original_name = self.assignment_df["ion"][i]
            original_formula = self.assignment_df["fitter_formula"][i]
            original_adduct = self.assignment_df["adduct"][i]
            original_loss = self.assignment_df["loss"][i]
            original_frag_site = self.assignment_df["frag_site"][i]
            original_variable_mod = self.assignment_df["variable_mod"][i]
            original_variable_mod_gain = self.assignment_df["variable_mod_gain"][i]
            original_variable_mod_loss = self.assignment_df["variable_mod_loss"][i]
            matched_peak = self.peak_list[i]
            mod_name = ["H2O", "NH3", "Na"]
            mods = [-18.010565, -17.026549, 21.981945]
            theo_ions = [matched_peak + adduct for adduct in mods]
            self.matched_theo_loss_ions, self.matched_observed_loss_ions = (
                self.report_matches_for_protein(
                    theo_ions, self.peak_list, self.accuracy*1.4142
                )# * sqrt(2) as it is a 2 way comparison
            )

            if original_name[0] != "I" and original_protein == self.name and original_variable_mod == self.mod_name:
                for theo, obs in zip(self.matched_theo_loss_ions, self.matched_observed_loss_ions):
                    for obs_individual in obs:
                        charge = self.assignment_df["charge"][obs_individual]
                        if charge == original_charge:
                            obs_match_list.append(obs_individual)
                            original_ion_match_list.append(i)
                            original_ion_name_list.append(original_name)
                            original_formula_list.append(original_formula)
                            original_adduct_list.append(original_adduct)
                            original_loss_list.append(original_loss)
                            original_frag_site_list.append(original_frag_site)
                            original_variable_mod_list.append(original_variable_mod)
                            original_variable_mod_gain_list.append(original_variable_mod_gain)
                            original_variable_mod_loss_list.append(original_variable_mod_loss)
                            if mod_name[theo] in ["H2O", "NH3"]:
                                loss_match_list.append(mod_name[theo])
                                add_match_list.append("")
                            elif mod_name[theo] in ["Na"]:
                                loss_match_list.append("H")
                                add_match_list.append(mod_name[theo])
                            charge_list.append(charge)

        self.matched_peaks = pd.DataFrame()
        self.matched_peaks["observed_ion_index"] = obs_match_list
        self.matched_peaks["original_ion_index"] = original_ion_match_list
        self.matched_peaks["original_ion_name"] = original_ion_name_list
        self.matched_peaks["original_ion_formula"] = original_formula_list
        self.matched_peaks["original_ion_adduct"] = original_adduct_list
        self.matched_peaks["original_ion_loss"] =  original_loss_list
        self.matched_peaks["adduct"] = add_match_list
        self.matched_peaks["loss"] = loss_match_list
        self.matched_peaks["charge"] = charge_list
        self.matched_peaks["frag_site"] = original_frag_site_list
        self.matched_peaks["variable_mod"] = original_variable_mod_list
        self.matched_peaks["variable_mod_gain"] = original_variable_mod_gain_list
        self.matched_peaks["variable_mod_loss"] = original_variable_mod_loss_list

        self.matched_peaks = (
            self.matched_peaks[~self.matched_peaks['observed_ion_index'].isin(assigned_peaks)]
        )
        self.matched_peaks = (
            self.matched_peaks[~self.matched_peaks['observed_ion_index'].isin(self.tested_ions)]
        )
        self.matched_peaks = self.matched_peaks.reset_index(drop=True)

        for tested_ion in self.matched_peaks["observed_ion_index"].to_list():
            self.tested_ions.append(tested_ion)

        self.iteration = 0
        self.stage = "loss"
        self.length_match_list = len(self.matched_peaks)
        if len(self.matched_peaks) >= 1:
            self.iterate_over_peaks(self.iteration, self.matched_peaks, self.stage)



    def on_true_button(self, event):
        self.assignment_df.loc[self.observed_ion, "name"] = self.name
        self.assignment_df.loc[self.observed_ion, "sequence"] = self.sequence
        self.assignment_df.loc[self.observed_ion, "n_mod"] = (
            f"(+) {self.n_mod_add} (-) {self.n_mod_sub}"
        )
        self.assignment_df.loc[self.observed_ion, "c_mod"] = (
            f"(+) {self.c_mod_add} (-) {self.c_mod_sub}"
        )
        self.assignment_df.loc[self.observed_ion, "ion"] = self.ion_name
        self.assignment_df.loc[self.observed_ion, "variable_mod"] = self.mod_name
        self.assignment_df.loc[self.observed_ion, "variable_mod_gain"] = self.mod_add_formula
        self.assignment_df.loc[self.observed_ion, "variable_mod_loss"] = self.mod_loss_formula 
        self.assignment_df.loc[self.observed_ion, "loss"] = self.total_loss
        self.assignment_df.loc[self.observed_ion, "adduct"] = self.total_adduct
        self.assignment_df.loc[self.observed_ion, "ppm_error"] = self.ppm_error
        self.assignment_df.loc[self.observed_ion, "fit_score"] = self.score
        self.assignment_df.loc[self.observed_ion, "theoretical_mz"] = self.mw
        self.assignment_df.loc[self.observed_ion, "total_intensity"] = self.total_intensity
        self.assignment_df.loc[self.observed_ion, "charge_scaled_abundance"] = (
            self.charge_scaled_abundance
        )
        self.assignment_df.loc[self.observed_ion, "fitter_formula"] = self.formula_string
        self.assignment_df.loc[self.observed_ion, "fitter_theo_x"] = str(self.scatter_x)
        self.assignment_df.loc[self.observed_ion, "fitter_theo_y"] = str(self.scatter_y)
        self.assignment_df.loc[self.observed_ion, "frag_site"] = str(self.frag_site).upper()
        print(
            f"Assigned {self.name} {self.ion_name} (+){self.total_adduct} (-){self.total_loss} with {self.mod_name}"
        )

        if self.iteration < len(self.matched_peaks) - 1:
            self.iteration += 1
            self.iterate_over_peaks(self.iteration, self.matched_peaks, self.stage)

        elif self.iteration >= len(self.matched_peaks) - 1:
            self.assignment_df.to_csv(self.peak_assignment_file, index=False)
            if self.length_match_list == 0:
                self.filter_results()
                self.Close()
            print("Continuing with next round of neutral loss scans.")
            self.neutral_loss_search()


    def on_false_button(self, event):
        if self.iteration < len(self.matched_peaks) - 1:
            self.iteration += 1
            self.iterate_over_peaks(self.iteration, self.matched_peaks, self.stage)
        elif self.iteration >= len(self.matched_peaks) - 1:
            self.assignment_df.to_csv(self.peak_assignment_file, index=False)
            print(f'Saved updated assigned peak list to {self.peak_assignment_file}')
            if self.length_match_list == 0:
                self.filter_results()
                self.Close()
            print("Continuing with next round of neutral loss scans.")
            self.neutral_loss_search()


    def theoretical_mass_generator(
        self,
        sequence,
        n_terminal_modification_mass,
        c_terminal_modification_mass,
        ion_type
    ):

        # lists of b and y ions
        theo_b_ions_list = []
        theo_y_ions_list = []

        # replace c with custom aa B to enable disulfide calc
        mass.std_aa_comp["B"] = mass.Composition({'C': 3, 'H': 4, 'N': 1, 'O': 1, 'S': 1})
        sequence = sequence.replace('c', 'B')

        sequence_mass = mass.calculate_mass(sequence, ion_type="M")
        aa_masses = np.array([mass.calculate_mass(aa, ion_type="b") for aa in sequence])

        if ion_type == "b/y":
            b_ion_masses = (
                np.cumsum(aa_masses) +
                n_terminal_modification_mass
            )
            y_ion_masses = (
                sequence_mass +
                c_terminal_modification_mass -
                np.cumsum(aa_masses)
            )
        elif ion_type == "c/z•":
            b_ion_masses = (
                np.cumsum(aa_masses) +
                n_terminal_modification_mass +
                17.026549
            )
            y_ion_masses = (
                sequence_mass +
                c_terminal_modification_mass -
                np.cumsum(aa_masses) -
                16.018724
            )            

        theo_b_ions_list.append(b_ion_masses[:-1]) # exclude the dehydrated full protein
        theo_y_ions_list.append(y_ion_masses[:-1])

        theo_y_ions_list.sort()

        return theo_b_ions_list, theo_y_ions_list


    def report_matches_for_protein(self, theoretical_ions, observed_ions, ppm_threshold):
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
                matched_theoretical_ion_indices.append(i)
                matched_observed_ion_indices.append(match)

        return matched_theoretical_ion_indices, matched_observed_ion_indices


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


    def parse_molecular_formula(self, formula_string):
        formula_dict = {}
        matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula_string)
        for element, count in matches:
            count = int(count) if count else 1
            formula_dict[element] = count
        return formula_dict


    def formula_to_string(self, formula):
        formula_string=""
        for element, count in formula.items():
            if count >= 1:
                formula_string += element
                if count >= 2:
                    formula_string += str(count)
        return formula_string


    def find_closest_index(self, value, array):
        closest_index = np.abs(array - value).argmin()
        return closest_index


    def get_indices_below_n_percent(self, lst, n):
        arr = np.array(lst)
        max_value = np.max(arr)
        threshold = n / 100 * max_value
        indices_below_threshold = np.where(arr < threshold)[0]

        return indices_below_threshold.tolist()


    def calibrate_spectrum(self):
        self.assignment_df = pd.read_csv(self.peak_assignment_file)
        self.spectrum = np.loadtxt(self.calibrated_spectrum_file)
        self.spectrum_x = self.spectrum[:, 0]
        self.spectrum_y = self.spectrum[:, 1]

        assigned_peaks = self.assignment_df[self.assignment_df["ion"].notna()]
        assigned_peaks.loc[:, "error"] = (
            assigned_peaks["monoisotopic_mz"] * (assigned_peaks["ppm_error"]/1000000)
        )
        x = assigned_peaks["monoisotopic_mz"].tolist()
        y = assigned_peaks["error"].tolist()

        coefficients = np.polyfit(x, y, 1)
        m, b = coefficients
        equation_of_curve = f"y = {m}x + {b}"
        regression_line = np.polyval(coefficients, x)

        residuals = y - regression_line

        ppm_residuals = []
        max_deviation = 0
        for mz, residual in zip(x, residuals):
            ppm_residual = residual / mz * 1000000
            ppm_residuals.append(ppm_residual)
            if abs(ppm_residual) >= max_deviation:
                max_deviation = abs(ppm_residual)
        max_deviation = round(max_deviation) + 1

        self.spectrum_x = [mz - (m * mz + b) for mz in self.spectrum_x]
        self.spectrum = np.column_stack((self.spectrum_x, self.spectrum_y))
        np.savetxt(self.calibrated_spectrum_file, self.spectrum, delimiter=" ")

        self.assignment_df["monoisotopic_mz"] = (
            self.assignment_df["monoisotopic_mz"] -
            (m * self.assignment_df["monoisotopic_mz"] + b)
        )
        self.assignment_df["monoisotopic_mw"] = (
            self.assignment_df["monoisotopic_mz"] *
            self.assignment_df["charge"] - 1.00727647 * self.assignment_df["charge"]
        )
        self.assignment_df["ppm_error"] = (
            (self.assignment_df["monoisotopic_mz"] - self.assignment_df["theoretical_mz"]) /
            self.assignment_df["monoisotopic_mz"] * 1000000
        )

        self.assignment_df.to_csv(self.peak_assignment_file, index=False)

        _, (ax1, ax2) = plt.subplots(
            2, 1, sharex=True,
            figsize=(8, 6), gridspec_kw={'height_ratios': [3, 1]},
            num="Spectral Calibration Summary"
        )

        ax1.scatter(x, y, color='darkviolet', marker='o', alpha=0.5)
        ax1.plot(x, regression_line, color='mediumseagreen')
        ax1.set_ylabel('Error (Th)')

        ax2.scatter(x, ppm_residuals, color='darkviolet', marker='o', alpha = 0.5)
        ax2.set_ylim(-max_deviation, max_deviation)
        ax2.axhline(0, color='black', linestyle='--', linewidth=1)
        ax2.set_xlabel('m/z')
        ax2.set_ylabel('Residuals (ppm)')

        plt.tight_layout()
        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.SetIcon(wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO))

        plt.show()

        calibration_info = np.column_stack((x, y, regression_line))
        np.savetxt(
            self.calibration_info_file,
            calibration_info,
            delimiter=" ",
            header=equation_of_curve
        )


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


class ValidationPlotsPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)

        self.max_intens = 0.0

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
            # peak label at top right
            label_text = (
                f"{ion_name}{variant_string}: {round(ion_mass, 2)} \n"
                f"Score: {score} \n"
                f"Mass error: {ppm_error} ppm"
            )
            self.ax_main.annotate(
                label_text,
                xy=(0.995,0.99), xycoords="axes fraction", ha="right", va="top"
            )
            # scatter plot w/ theoretical distribution
            self.ax_main.scatter(scatter_x, scatter_y, c="darkviolet", marker="o", alpha=0.6)

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
