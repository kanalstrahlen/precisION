import bisect
import os
import wx
import wx.grid
import pandas as pd
from pyteomics import mass
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import re


class DefineTerminiWindow(wx.Frame):
    def __init__(self, parent, title, sequence, file_path, directory_path):
        super().__init__(
            parent,
            title=title,
            size=(415, 560),
            style=wx.DEFAULT_FRAME_STYLE & ~(wx.RESIZE_BORDER | wx.MAXIMIZE_BOX)
        )

        self.basename = os.path.basename(file_path).replace(".txt", "")
        self.directory_path = directory_path
        self.sequence = sequence

        panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        layout_sizer = wx.BoxSizer(wx.VERTICAL)

        title = wx.StaticText(panel, label="Define Termini")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext1 = wx.StaticText(
            panel,
            label="Scan sequences to identify the native protein termini."
        )
        subtext2 = wx.StaticText(
            panel,
            label="Define up to five modifications to consider."
        )

        # n_term_mod_sizer: options for N term mods to consider
        n_term_mod_staticbox = wx.StaticBox(panel)
        n_term_mod_sizer = wx.StaticBoxSizer(n_term_mod_staticbox, wx.VERTICAL)
        n_term_mod_header = wx.StaticText(panel, label="N-terminal Modifications")
        subtitle_font = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)
        n_term_mod_header.SetFont(subtitle_font)
        self.n_term_grid = self.create_modification_grid(panel, "N")

        n_term_mod_sizer.Add(n_term_mod_header, 0, wx.EXPAND | wx.BOTTOM, 5)
        n_term_mod_sizer.Add(self.n_term_grid, 0, wx.EXPAND | wx.BOTTOM, 0)
        n_term_mod_sizer.SetMinSize((400, -1))

        # c_term_mod_sizer: options for C term mods to consider
        c_term_mod_staticbox = wx.StaticBox(panel)
        c_term_mod_sizer = wx.StaticBoxSizer(c_term_mod_staticbox, wx.VERTICAL)
        c_term_mod_header = wx.StaticText(panel, label="C-terminal Modifications")
        c_term_mod_header.SetFont(subtitle_font)
        self.c_term_grid = self.create_modification_grid(panel, "C")

        c_term_mod_sizer.Add(c_term_mod_header, 0, wx.EXPAND | wx.BOTTOM, 5)
        c_term_mod_sizer.Add(self.c_term_grid, 0, wx.EXPAND | wx.BOTTOM, 0)
        c_term_mod_sizer.SetMinSize((400, -1))

        # options_sizer: options for search
        options_staticbox = wx.StaticBox(panel)
        options_sizer = wx.StaticBoxSizer(options_staticbox, wx.VERTICAL)
        options_header = wx.StaticText(panel, label="Search Options")
        options_header.SetFont(subtitle_font)

        options_subsizer = wx.BoxSizer(wx.HORIZONTAL)
        ppm_threshold_text = wx.StaticText(panel, label="Mass accuracy (ppm): ")
        self.ppm_threshold = wx.TextCtrl(panel, size=(30, 20), value="10")
        ion_type_text = wx.StaticText(panel, label = "Ion type: ")
        ion_type_choices = ["b/y", "c/z•"]
        self.ion_type = wx.ComboBox(
            panel,
            choices = ion_type_choices,
            style = wx.CB_READONLY
        )
        self.ion_type.SetValue("b/y")

        options_subsizer.Add(ppm_threshold_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        options_subsizer.Add(self.ppm_threshold, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)
        options_subsizer.Add(ion_type_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        options_subsizer.Add(self.ion_type, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)

        options_sizer.Add(options_header, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.Add(options_subsizer, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.SetMinSize((400, -1))

        define_termini_button = wx.Button(panel, label="Run Search",  size=(100, 40))
        define_termini_button.Bind(wx.EVT_BUTTON, self.on_define_termini_button)

        layout_sizer.Add(title, 0, wx.EXPAND | wx.BOTTOM, 5)
        layout_sizer.Add(subtext1, 0,  wx.EXPAND | wx.BOTTOM, 5)
        layout_sizer.Add(subtext2, 0,  wx.EXPAND | wx.BOTTOM, 0)
        layout_sizer.Add(n_term_mod_sizer, 0, wx.EXPAND | wx.BOTTOM, 0)
        layout_sizer.Add(c_term_mod_sizer, 0, wx.EXPAND | wx.BOTTOM, 0)
        layout_sizer.Add(options_sizer, 0, wx.EXPAND | wx.BOTTOM, 5)
        layout_sizer.Add(define_termini_button, 0, wx.BOTTOM |  wx.ALIGN_CENTER_HORIZONTAL, 0)

        main_sizer.Add(layout_sizer, 0, wx.EXPAND | wx.ALL, 5)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        panel.SetSizer(main_sizer)
        self.Show()


    def create_modification_grid(self, parent, terminus):
        grid = wx.grid.Grid(parent)
        grid.CreateGrid(5, 3)
        grid.SetDefaultCellBackgroundColour("#f0f0f0")

        grid.SetColLabelValue(0, "Modification Name")
        grid.SetColLabelValue(1, "Gain Formula")
        grid.SetColLabelValue(2, "Loss Formula")

        if terminus == "N":
            default_data = [
                ("No Modification", "", ""),
                ("Acetylation", "C2H2O", ""),
                ("Formylation", "CO", ""),
                ("", "", ""),
                ("", "", ""),
            ]
        elif terminus == "C":
            default_data = [
                ("No Modification", "", ""),
                ("Deamidation", "O", "NH"),
                ("Methylation", "CH2", ""),
                ("", "", ""),
                ("", "", ""),
            ]

        for i, (name, add_formula, subtract_formula) in enumerate(default_data):
            grid.SetCellValue(i, 0, name)
            grid.SetCellValue(i, 1, add_formula)
            grid.SetCellValue(i, 2, subtract_formula)

        for col in range(3):
            grid.AutoSizeColumn(col)

        return grid


    def load_peaks(self, df=True):
        # use the recalibrated values if they exist
        peak_file_path = os.path.join(self.directory_path, f"{self.basename}.assignedPeakList.csv")
        if os.path.exists(peak_file_path):
            peak_list = pd.read_csv(peak_file_path)
            peak_list = peak_list.sort_values(by="monoisotopic_mw")
        else:        
            peak_file_path = os.path.join(self.directory_path, f"{self.basename}.filteredPeakList.csv")  
            peak_list = pd.read_csv(peak_file_path)
            peak_list = peak_list.sort_values(by="monoisotopic_mw")
            peak_list = peak_list[peak_list["prediction"] == True]

        # vary if just the list of masses or full dataframe is required
        if df:
            pass
        else:
            peak_list = peak_list["monoisotopic_mw"].tolist()

        return peak_list


    def extract_modification_data(self, grid):
        modification_data = []
        for row in range(grid.GetNumberRows()):
            name = grid.GetCellValue(row, 0)
            add_formula = grid.GetCellValue(row, 1)
            subtract_formula = grid.GetCellValue(row, 2)
            modification_data.append((name, add_formula, subtract_formula))
        return modification_data


    def on_define_termini_button(self, event):
        n_term_mods = self.extract_modification_data(self.n_term_grid)
        c_term_mods = self.extract_modification_data(self.c_term_grid)

        observed_peaks = self.load_peaks(df=False)

        ppm_error = float(self.ppm_threshold.GetValue())
        ion_type = self.ion_type.GetValue()

        run = DefineTerminiFunctions()
        run.truncation_scan(
            self.sequence,
            n_term_mods,
            c_term_mods,
            observed_peaks,
            ppm_error,
            ion_type
        )

        self.Close()



class DefineTerminiFunctions():
    def truncation_scan(self, sequence, n_mods, c_mods, observed_peaks, ppm_error, ion_type):
        sequence = self.process_sequence_string(sequence)

        n_mod_count = len([mod for mod in n_mods if mod[0] != ""])
        c_mod_count = len([mod for mod in c_mods if mod[0] != ""])
        num_columns = max((n_mod_count, c_mod_count))

        _, ax = plt.subplots(2, num_columns, num="Define Termini")
        max_n_match_per_res = 0
        max_c_match_per_res = 0

        for i, mod in tqdm(
            enumerate(n_mods),
            desc='N-terminal modifications',
            total=n_mod_count
        ):
            mod_name = mod[0]
            if mod_name != "":
                mod_addition = mod[1]
                mod_subtraction = mod[2]
                mod_mass = (
                    mass.calculate_mass(formula=mod_addition) -
                    mass.calculate_mass(formula=mod_subtraction)
                )

                match_per_res = []
                for j in range(len(sequence)):
                    test_sequence = sequence[j:]
                    theo_b_ions, _ = self.theoretical_mass_generator(
                        test_sequence,
                        mod_mass,
                        0,
                        ion_type
                    )

                    matched_theoretical_b_ions, _ = self.report_matches_for_protein(
                        theo_b_ions[0],
                        observed_peaks,
                        ppm_error
                    )

                    match_per_res.append(len(matched_theoretical_b_ions))

                # set axis limit to be max of all plots and plot data
                res_index = [m+1 for m in list(range(len(sequence)))]
                max_value = max(match_per_res)
                max_n_match_per_res = max(max_n_match_per_res, max_value)
                max_index = res_index[match_per_res.index(max_value)]
                ax[0, i].plot(res_index, match_per_res, color="#205283")
                ax[0, i].set_title(mod_name)
                ax[0, i].annotate(f"Res:{max_index}", xy = (max_index, max_value+1))
                #print(sum(match_per_res))
                #print(sum(match_per_res))
                #print(sum(match_per_res)/len(match_per_res))

        # same for y ions
        for i, mod in tqdm(
            enumerate(c_mods),
            desc='C-terminal modifications',
            total=c_mod_count
        ):
            mod_name = mod[0]
            if mod_name != "":
                mod_addition = mod[1]
                mod_subtraction = mod[2]
                mod_mass = (
                    mass.calculate_mass(formula=mod_addition) -
                    mass.calculate_mass(formula=mod_subtraction)
                )

                match_per_res = []
                for j in range(len(sequence)):
                    test_sequence = sequence[:j+1]
                    _, theo_y_ions = self.theoretical_mass_generator(
                        test_sequence,
                        0,
                        mod_mass,
                        ion_type
                    )
                    matched_theoretical_y_ions, _ = self.report_matches_for_protein(
                        theo_y_ions[0],
                        observed_peaks,
                        ppm_error
                    )
                    match_per_res.append(len(matched_theoretical_y_ions))

                res_index = [m+1 for m in list(range(len(sequence)))]
                max_value = max(match_per_res)
                max_c_match_per_res = max(max_c_match_per_res, max_value)
                max_index = res_index[match_per_res.index(max_value)]

                ax[1, i].plot(res_index, match_per_res, color="#b42920")
                ax[1, i].set_title(mod_name)
                ax[1, i].annotate(f"Res:{max_index}", xy = (max_index, max_value+1), ha='right')
                #print(sum(match_per_res))
                #print(sum(match_per_res))
                #print(sum(match_per_res)/len(match_per_res))


        # general plotting parameters
        for i in range(num_columns):
            ax[0, i].set_ylim(0, max_n_match_per_res + 3)
            ax[1, i].set_ylim(0, max_c_match_per_res + 3)
            ax[0, i].set_xlim(-10, len(sequence)+12)
            ax[1, i].set_xlim(-10, len(sequence)+12)
            ax[0, 0].set_ylabel("Number of b-ion matches")
            ax[1, 0].set_ylabel("Number of y-ion matches")
            ax[1, i].set_xlabel("Residue index")

        # make plot window look like rest of applications
        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.SetIcon(wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO))

        plt.tight_layout()
        plt.show()


    def theoretical_mass_generator(
        self,
        sequence,
        n_terminal_modification_mass,
        c_terminal_modification_mass,
        ion_type = "b/y"
    ):

        if ion_type == "c/z•":
            cz_ions = True
        elif ion_type == "b/y":
            cz_ions = False

        # lists of b and y ions
        theo_b_ions_list = []
        theo_y_ions_list = []

        # replace c with B to enable disulfide calc
        mass.std_aa_comp["B"] = mass.Composition({'C': 3, 'H': 4, 'N': 1, 'O': 1, 'S': 1})
        sequence = sequence.replace('c', 'B')

        sequence_mass = mass.calculate_mass(sequence, ion_type="M") # intact mass
        aa_masses = np.array(
        [mass.calculate_mass(aa, ion_type="b") for aa in sequence]
        ) # individual residue masses

        if cz_ions:
            b_ion_masses = (
                np.cumsum(aa_masses) +
                n_terminal_modification_mass +
                17.026549 # + nh3
            )
            y_ion_masses = (
                sequence_mass +
                c_terminal_modification_mass -
                np.cumsum(aa_masses) -
                17.026549 + 1.00783 # - nh3 + extra proton to balance electron charge + electron
            )
        else:
            b_ion_masses = (
                np.cumsum(aa_masses) +
                n_terminal_modification_mass
            )
            y_ion_masses = (
                sequence_mass +
                c_terminal_modification_mass -
                np.cumsum(aa_masses)
            )

        theo_b_ions_list.append(b_ion_masses[:-1]) # exclude the dehydrated full protein
        theo_y_ions_list.append(y_ion_masses[:-1])

        # sort for binary search
        theo_y_ions_list.sort()

        return theo_b_ions_list, theo_y_ions_list


    def process_sequence_string(self, input_string):
        result_string = ""
        i = 0
        while i < len(input_string):
            if input_string[i] == '[':
                closing_bracket_index = input_string.find(']', i)
                i = closing_bracket_index + 1
            else:
                result_string += input_string[i]
                i += 1

        result_string = result_string.replace(" ", "")

        return result_string


    def report_matches_for_protein(self, theoretical_ions, observed_ions, ppm_threshold):
        # binary search w/ ppm tolerance helper function
        def binary_search(arr, target, ppm_threshold):
            low = bisect.bisect_left(arr, target * (1 - ppm_threshold / 1e6))
            high = bisect.bisect_right(arr, target * (1 + ppm_threshold / 1e6))
            return high - low, low

        matched_theoretical_ion_indices = []
        matched_observed_ion_indices = []

        # sort both lists for binary search; nearly always sorted but this is still there
        theoretical_ions.sort()
        observed_ions.sort()

        for i, theo_ion in enumerate(theoretical_ions):
            num_matches, match = binary_search(observed_ions, theo_ion, ppm_threshold)
            if num_matches >= 1:
                matched_theoretical_ion_indices.append(i)
                matched_observed_ion_indices.append(match)

        return matched_theoretical_ion_indices, matched_observed_ion_indices
