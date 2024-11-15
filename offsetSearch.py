import math
import bisect
import re
import wx
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from pyteomics import mass
from tqdm import tqdm
from scipy.stats import poisson
import matplotlib.pyplot as plt
from fastOffsetSearch import FastOffsetShiftCalculator


matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


# to do; internal offset search or other way to 
# identify  modified internal fragments other than massdiff
# could also add window to show which ions mass diff includes to help ID
# for

# to do; explore effect of mass defect on matching. if important
# could correct

class OffsetSearchWindow(wx.Frame):
    def __init__(self, parent, title, file_path, sequence, options):
        super().__init__(parent, title=title)
        self.file_path = file_path
        # screen width variable to allow for smaller screens
        screen_width, screen_height = wx.DisplaySize()
        self.SetSize((screen_width * 2 // 4, screen_height * 2 // 4))

        self.panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.plot_panel = OffsetSearchPlotPanel(self.panel)

        grid_sizer = wx.BoxSizer(wx.VERTICAL)

        subheader_font = wx.Font(
            10,
            wx.FONTFAMILY_DEFAULT,
            wx.FONTSTYLE_ITALIC,
            wx.FONTWEIGHT_BOLD
        )
        b_ion_grid_subheader = wx.StaticText(
            self.panel,
            label="Most common b-type ion matches"
        )
        y_ion_grid_subheader = wx.StaticText(
            self.panel,
            label="Most common y-type ion matches"
        )
        b_ion_grid_subheader.SetFont(subheader_font)
        y_ion_grid_subheader.SetFont(subheader_font)

        self.b_ion_grid = wx.grid.Grid(self.panel)
        self.b_ion_grid.SetDefaultCellBackgroundColour("#f0f0f0")

        self.y_ion_grid = wx.grid.Grid(self.panel)
        self.y_ion_grid.SetDefaultCellBackgroundColour("#f0f0f0")

        self.b_ion_num_rows = 1
        self.y_ion_num_rows = 1
        self.num_cols = 3

        self.b_ion_grid.CreateGrid(self.b_ion_num_rows, self.num_cols)
        self.y_ion_grid.CreateGrid(self.y_ion_num_rows, self.num_cols)

        self.b_ion_grid.Bind(wx.grid.EVT_GRID_SELECT_CELL, self.on_b_ion_cell_select)
        self.y_ion_grid.Bind(wx.grid.EVT_GRID_SELECT_CELL, self.on_y_ion_cell_select)

        column_names = [
            "Mass shift (Da)",
            "Number of matches",
            "-log(E-value)"
            ]
        for col, column_name in enumerate(column_names):
            self.b_ion_grid.SetColLabelValue(col, column_name)
            self.y_ion_grid.SetColLabelValue(col, column_name)

        button_sizer = wx.BoxSizer(wx.HORIZONTAL)

        reset_plots_button = wx.Button(
            self.panel,
            label="Reset plots",
            size=(170, 30)
        )
        reset_plots_button.Bind(wx.EVT_BUTTON, self.on_reset_plots_button)

        cz_ion_button = wx.Button(
            self.panel,
            label="Calculate c/z•-type ions",
            size=(170, 30)
        )
        cz_ion_button.Bind(wx.EVT_BUTTON, self.on_cz_ion_button)

        button_sizer.Add(reset_plots_button, 0, wx.RIGHT, 10)
        button_sizer.Add(cz_ion_button, 0, wx.RIGHT, 0)

        grid_sizer.Add(b_ion_grid_subheader, 0, wx.ALL, 0)
        grid_sizer.Add(self.b_ion_grid, 0, wx.ALL, 5)
        grid_sizer.Add(y_ion_grid_subheader, 0, wx.ALL, 0)
        grid_sizer.Add(self.y_ion_grid, 0, wx.ALL, 5)
        grid_sizer.Add(button_sizer, 0, wx.TOP | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.b_ion_grid.AutoSizeColumns()
        self.b_ion_grid.Refresh()
        self.y_ion_grid.AutoSizeColumns()
        self.y_ion_grid.Refresh()

        main_sizer.Add(self.plot_panel, 1, wx.EXPAND | wx.RIGHT, border=20)
        main_sizer.Add(grid_sizer, 0, wx.ALL, border=15)

        icon = wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        calculator = OffsetShiftCalculator(options)
        theo_b_ions, theo_y_ions = calculator.gen_theo_fragment_masses(sequence)
        peak_list = calculator.gen_peak_list(file_path)

        print("Searching b-type ions...")
        b_ion_offset, b_ion_num_matches = calculator.offset_scan(peak_list, theo_b_ions)
        print("Searching y-type ions...")
        y_ion_offset, y_ion_num_matches = calculator.offset_scan(peak_list, theo_y_ions)

        b_ion_df = calculator.find_max_peaks(b_ion_offset, b_ion_num_matches)
        y_ion_df = calculator.find_max_peaks(y_ion_offset, y_ion_num_matches)
        self.update_grids(b_ion_df, y_ion_df)

        self.plot_panel.plot_offset_scan(
            b_ion_offset,
            b_ion_num_matches,
            y_ion_offset,
            y_ion_num_matches
        )

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Show()


    def on_size(self, event):
        self.panel.Layout()
        event.Skip()


    def on_b_ion_cell_select(self, event):
        try:
            selected_row = event.GetRow()
            selected_col = event.GetCol()
            ion_mass = float(self.b_ion_grid.GetCellValue(selected_row, 0))
            max_value = float(self.b_ion_grid.GetCellValue(selected_row, 1))
            self.plot_panel.b_ion_zoom(ion_mass, max_value)
            event.Skip()
        except:
            event.Skip()


    def on_y_ion_cell_select(self, event):
        try:
            selected_row = event.GetRow()
            selected_col = event.GetCol()
            ion_mass = float(self.y_ion_grid.GetCellValue(selected_row, 0))
            max_value = float(self.y_ion_grid.GetCellValue(selected_row, 1))
            self.plot_panel.y_ion_zoom(ion_mass, max_value)
            event.Skip()
        except:
            event.Skip()


    def on_reset_plots_button(self, event):
        self.plot_panel.reset_zoom()
        event.Skip()


    def update_grids(self, b_ion_df, y_ion_df):
        self.b_ion_grid.ClearGrid()

        num_rows, num_cols = b_ion_df.shape

        diff = self.b_ion_num_rows - num_rows
        self.b_ion_grid.AppendRows(-1 * diff)

        for row in range(num_rows):
            for col in range(num_cols):
                value = str(b_ion_df.iloc[row, col])
                self.b_ion_grid.SetCellValue(row, col, value)
        self.b_ion_grid.AutoSizeColumns()
        self.b_ion_grid.Refresh()
        self.b_ion_grid.Show()

        self.y_ion_grid.ClearGrid()

        num_rows, num_cols = y_ion_df.shape

        diff = self.y_ion_num_rows - num_rows
        self.y_ion_grid.AppendRows(-1 * diff)

        for row in range(num_rows):
            for col in range(num_cols):
                value = str(y_ion_df.iloc[row, col])
                self.y_ion_grid.SetCellValue(row, col, value)
        self.y_ion_grid.AutoSizeColumns()
        self.y_ion_grid.Refresh()
        self.y_ion_grid.Show()

        self.panel.Layout()


    def on_cz_ion_button(self, event):
        num_rows = self.b_ion_grid.GetNumberRows()
        for row in range(num_rows):
            value = self.b_ion_grid.GetCellValue(row, 0)
            corrected_value = str(round(float(value) - 17.026549, 4))
            self.b_ion_grid.SetCellValue(row, 0, corrected_value)
        self.b_ion_grid.AutoSizeColumns()
        self.b_ion_grid.Refresh()
        self.b_ion_grid.Show()

        num_rows = self.y_ion_grid.GetNumberRows()
        for row in range(num_rows):
            value = self.y_ion_grid.GetCellValue(row, 0)
            corrected_value = str(round(float(value) + 17.026549 - 1.00783, 4))
            self.y_ion_grid.SetCellValue(row, 0, corrected_value)
        self.y_ion_grid.AutoSizeColumns()
        self.y_ion_grid.Refresh()
        self.y_ion_grid.Show()

        print("Transformed mass values in table to c/z•-type ions")



class OffsetShiftCalculator():
    def __init__(self, options):
        self.lower_range, self.upper_range, self.threshold, self.count_mode = options


    def find_max_peaks(self, offset_list, match_count_list):
        offset_array = np.array(offset_list)
        match_count_array = np.array(match_count_list)
        match_count_array = match_count_array.astype(int)

        expected = np.mean(match_count_array)

        p_score = np.zeros_like(match_count_array, dtype=float)

        for i, match_count in enumerate(match_count_array):
            p_score[i] = (
                (expected ** match_count) *
                math.exp(-1 * expected) /
                math.factorial(match_count)
            )

        e_value = p_score * len(p_score)

        unique, counts = np.unique(match_count_array, return_counts=True)
        relative_freq = counts / np.sum(counts)

        #plt.scatter(unique, relative_freq)
        #x = np.arange(
        #    poisson.ppf(0.0000001, expected),
        #    poisson.ppf(0.9999999, expected)
        #)
        #plt.plot(x, poisson.pmf(x, expected))
        ##print(expected)
        #plt.show()

        for i, value in enumerate(e_value):
            e_value[i] = self.neg_log_and_round(value, 1)

        sorted_indices = np.argsort(match_count_array)[::-1]
        top_indices = []

        for index in sorted_indices:
            # Check if the selected indices are at least min_distance apart
            if all(abs(index - selected_index) >= 30 for selected_index in top_indices):
                top_indices.append(index)
            # Stop when you have selected 10 indices
            if len(top_indices) == 10:
                break

        data = {
            'offset': [round(offset, 4) for offset in offset_array[top_indices]],
            'values': match_count_array[top_indices],
            'e_values': e_value[top_indices]
        }
        df = pd.DataFrame(data)

        return df


    def neg_log_and_round(self, number, dp):
        if number == 0:
            return 0.0

        neg_log = -1 * math.log10(abs(number))
        rounded = round(neg_log, dp)

        return rounded


    def offset_scan(self, peak_list, theo_ions):

        calculator = FastOffsetShiftCalculator(
            threshold=self.threshold,
            count_mode=self.count_mode,
            lower_range=self.lower_range,
            upper_range=self.upper_range
        )
        offset_list, match_count_list = calculator.offset_scan(peak_list, theo_ions)
        #scan_spacing = self.threshold / 1000 / 1.5 # 2/3 the threshold at m/z 1000
#
        #num_points = int((self.upper_range - self.lower_range) // scan_spacing)
#
        #offset_list = np.linspace(self.lower_range, self.upper_range, num=num_points)
        #peak_list = np.array(peak_list)
        #theo_ions = np.array(theo_ions)
#
        #match_count_list = []
#
        #for offset in tqdm(offset_list, total=len(offset_list)):
        #    temp_theo_ions = [ion_mass + offset for ion_mass in theo_ions]
        #    if len(theo_ions) > 800:
        #        match_count = self.count_matches(peak_list, temp_theo_ions)
        #    else:
        #        match_count = self.count_matched_ions(peak_list, temp_theo_ions)
        #    match_count_list.append(match_count)

        return offset_list, match_count_list


    def count_matches(self, exp_ions, theo_ions):
        exp_pointer = 0
        theo_pointer = 0
        matches = 0

        while exp_pointer < len(exp_ions) and theo_pointer < len(theo_ions):
            exp_mass = exp_ions[exp_pointer]
            theo_mass = theo_ions[theo_pointer]
            ppm_difference = abs(exp_mass - theo_mass) / theo_mass * 1e6

            if ppm_difference < self.threshold:
                matches += 1
                exp_pointer += 1
                if self.count_mode == "Theoretical ions":
                    theo_pointer += 1
                elif self.count_mode == "Observed ions":
                    theo_pointer = theo_pointer
            elif exp_mass < theo_mass:
                exp_pointer += 1
            else:
                theo_pointer += 1

        return matches


    def count_matched_ions(self, peak_list, theo_ions):
        match_count = 0

        def binary_search(arr, target, ppm_threshold):
            low = bisect.bisect_left(arr, target * (1 - ppm_threshold / 1e6))
            high = bisect.bisect_right(arr, target * (1 + ppm_threshold / 1e6))
            return high - low

        for ion in theo_ions:
            num_matches = binary_search(peak_list, ion, self.threshold)
            if num_matches >=1:
                if self.count_mode == "Theoretical ions":
                    match_count += 1
                elif self.count_mode == "Observed ions":
                    match_count += num_matches

        return match_count


    def gen_peak_list(self, assigned_peaks_file):
        all_peaks = pd.read_csv(assigned_peaks_file)
        unassigned_peaks = all_peaks[pd.isna(all_peaks['name'])]
        peak_list = unassigned_peaks["monoisotopic_mw"].to_list()
        return peak_list


    def gen_theo_fragment_masses(self, sequence):
        unmod_seq, mod_arr = self.process_sequence_string(sequence)
        mod_arr = [self.mod_string_to_dict(mod) for mod in mod_arr]

        b_ion_masses = []
        for i in range(1,len(unmod_seq)):
            ion_seq = unmod_seq[:i]
            formula = self.calculate_molecular_formula(ion_seq)
            mod_dicts = mod_arr[0:i]
            mod_sum = {key: 0 for d in mod_dicts for key in d}
            for d in mod_dicts:
                for key, value in d.items():
                    mod_sum[key] += value
            formula = {
                key: formula.get(key, 0) + mod_sum.get(key, 0)
                for key in set(formula) | set(mod_sum)
            }
            ion_mass = mass.calculate_mass(composition=formula)
            b_ion_masses.append(ion_mass)

        y_ion_masses = []
        for i in range(len(unmod_seq)-1):
            ion_seq = unmod_seq[i+1:]
            formula = self.calculate_molecular_formula(ion_seq)
            formula["H"] += 2
            formula["O"] += 1
            mod_dicts = mod_arr[i+1:]
            mod_sum = {key: 0 for d in mod_dicts for key in d}
            for d in mod_dicts:
                for key, value in d.items():
                    mod_sum[key] += value
            formula = {
                key: formula.get(key, 0) + mod_sum.get(key, 0)
                for key in set(formula) | set(mod_sum)
            }
            ion_mass = mass.calculate_mass(composition=formula)
            y_ion_masses.append(ion_mass)

        y_ion_masses = y_ion_masses[::-1]

        return b_ion_masses, y_ion_masses


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

        if mod == None:
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


class OffsetSearchPlotPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        self.figure = Figure(facecolor='#f0f0f0')

        self.ax1 = self.figure.add_subplot(211)
        self.ax1.set_facecolor('#f0f0f0')

        self.ax2 = self.figure.add_subplot(212)
        self.ax2.set_facecolor('#f0f0f0')

        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.toolbar = NavigationToolbar2Wx(self.canvas)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.SetSizer(sizer)
        self.Layout()

        self.toolbar.Show()
        self.Bind(wx.EVT_SIZE, self.on_size)

        self.b_ion_offset = []
        self.b_ion_num_matches = []
        self.y_ion_offset = []
        self.y_ion_num_matches = []


    def plot_offset_scan(
        self, b_ion_offset, b_ion_num_matches, y_ion_offset, y_ion_num_matches
    ):
        self.ax1.plot(b_ion_offset, b_ion_num_matches, color="#205283")
        self.ax2.plot(y_ion_offset, y_ion_num_matches, color="#b42920")

        self.b_ion_offset = b_ion_offset
        self.b_ion_num_matches = b_ion_num_matches
        self.y_ion_offset = y_ion_offset
        self.y_ion_num_matches = y_ion_num_matches

        self.ax1.set_ylim(0, max(b_ion_num_matches) * 1.2)
        self.ax2.set_ylim(0, max(y_ion_num_matches) * 1.2)
        self.ax1.set_xlim(min(b_ion_offset), max(b_ion_offset))
        self.ax2.set_xlim(min(y_ion_offset), max(y_ion_offset))

        self.ax1.set_ylabel("# matches")
        self.ax2.set_ylabel("# matches")
        self.ax2.set_xlabel("Mass difference (Da)")

        self.canvas.draw()


    def on_size(self, event):
        self.fit_plot_to_panel()
        event.Skip()


    def fit_plot_to_panel(self):
        try:
            size = self.GetClientSize()
            dpi = self.GetContentScaleFactor() * 100
            width = size.width / dpi
            height = size.height / dpi - 0.3
            self.figure.set_size_inches(width, height)
            self.figure.tight_layout(rect=[0, 0, 1, 1])
            self.canvas.draw()
        except:
            pass

    def b_ion_zoom(self, ion_mass, max_value):
        self.ax1.set_ylim(0, max_value * 1.2)
        self.ax1.set_xlim(ion_mass - 0.5, ion_mass + 0.5)
        self.canvas.draw()


    def y_ion_zoom(self, ion_mass, max_value):
        self.ax2.set_ylim(0, max_value * 1.2)
        self.ax2.set_xlim(ion_mass - 0.5, ion_mass + 0.5)
        self.canvas.draw()


    def reset_zoom(self):
        self.ax1.set_ylim(0, max(self.b_ion_num_matches) * 1.2)
        self.ax2.set_ylim(0, max(self.y_ion_num_matches) * 1.2)
        self.ax1.set_xlim(min(self.b_ion_offset), max(self.b_ion_offset))
        self.ax2.set_xlim(min(self.y_ion_offset), max(self.y_ion_offset))
        self.canvas.draw()
