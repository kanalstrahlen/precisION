import bisect
import re
import math
import xml.etree.ElementTree as ET
import wx
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from pyteomics import mass


matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

class OffsetIDWindow(wx.Frame):
    def __init__(self, parent, title, file_path, sequence, options):
        super().__init__(parent, title=title)
        self.file_path = file_path
        screen_width, screen_height = wx.DisplaySize()
        self.SetSize((screen_width * 1 // 2, screen_height * 1 // 2))

        calculator = OffsetIDCalculator(options)
        theo_b_ions, theo_y_ions = calculator.gen_theo_fragment_masses(sequence)
        peak_list = calculator.gen_peak_list(file_path)
        b_ion_matches = calculator.offset_search(peak_list, theo_b_ions)
        b_ion_matches = [match_index + 1 for match_index in b_ion_matches]
        y_ion_matches = calculator.offset_search(peak_list, theo_y_ions)
        unimod_df = calculator.find_unimod_matches()

        self.panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.plot_panel = OffsetSearchPlotPanel(
            self.panel,
            sequence,
            b_ion_matches,
            y_ion_matches
        )

        grid_sizer = wx.BoxSizer(wx.VERTICAL)

        unimod_grid_subheader = wx.StaticText(self.panel, label="Closest modifications")
        subheader_font = wx.Font(
            10,
            wx.FONTFAMILY_DEFAULT,
            wx.FONTSTYLE_ITALIC,
            wx.FONTWEIGHT_BOLD
        )
        unimod_grid_subheader.SetFont(subheader_font)

        self.unimod_grid = wx.grid.Grid(self.panel)
        self.unimod_grid.SetMinSize((560, 500))
        self.unimod_grid.SetDefaultCellBackgroundColour("#f0f0f0")
        self.num_rows = 1
        self.num_cols = 4
        self.unimod_grid.CreateGrid(self.num_rows, self.num_cols)
        unimod_column_names = [
            "Modification",
            "Formula",
            "Exact mass (Da)",
            "Error (Da)"
        ]
        for col, unimod_column_name in enumerate(unimod_column_names):
            self.unimod_grid.SetColLabelValue(col, unimod_column_name)

        grid_sizer.Add(unimod_grid_subheader, 0, wx.ALL, 0)
        grid_sizer.Add(self.unimod_grid, 0, wx.ALL, 5)

        self.unimod_grid.AutoSizeColumns()
        self.unimod_grid.Refresh()

        main_sizer.Add(self.plot_panel, 1, wx.EXPAND | wx.RIGHT, border=20)
        main_sizer.Add(grid_sizer, 0, wx.ALL, border=15)

        icon = wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.update_grid(unimod_df)

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Bind(wx.EVT_CLOSE, self.on_close)

        width = self.plot_panel.figure.get_figwidth()
        height = self.plot_panel.figure.get_figheight()
        dpi = self.plot_panel.figure.get_dpi()
        self.SetSize(((dpi * width + 560) * 1.05, (dpi * height) * 1.2))
        self.panel.Layout()

        self.Show()


    def on_close(self, event):
        plt.close("Offset ID")
        event.Skip()


    def on_size(self, event):
        self.panel.Layout()
        event.Skip()


    def update_grid(self, unimod_df):
        self.unimod_grid.ClearGrid()
        num_rows, num_cols = unimod_df.shape
        diff = self.num_rows - num_rows
        self.unimod_grid.AppendRows(-1 * diff)
        for row in range(num_rows):
            for col in range(num_cols):
                value = unimod_df.iloc[row, col]
                if type(value) != str:
                    value = round(value, 4)
                self.unimod_grid.SetCellValue(row, col, str(value))
        self.unimod_grid.SetColSize(0, 180)
        self.unimod_grid.SetColSize(1, 100)
        self.unimod_grid.SetColSize(2, 100)
        self.unimod_grid.SetColSize(3, 70)
        self.unimod_grid.Refresh()
        self.unimod_grid.Show()
        self.panel.Layout()



class OffsetIDCalculator():
    def __init__(self, options):
        self.offset, self.threshold = options


    def offset_search(self, peak_list, theo_ions):
        temp_theo_ions = [ion_mass + self.offset for ion_mass in theo_ions]
        matched_ion_indices = self.id_matched_ions(peak_list, temp_theo_ions)
        return matched_ion_indices


    def id_matched_ions(self, peak_list, theo_ions):
        matched_ion_indices = []

        def binary_search(arr, target, ppm_threshold):
            low = bisect.bisect_left(arr, target * (1 - ppm_threshold / 1e6))
            high = bisect.bisect_right(arr, target * (1 + ppm_threshold / 1e6))
            return high - low

        for i, ion in enumerate(theo_ions):
            num_matches = binary_search(peak_list, ion, self.threshold)
            if num_matches >=1:
                matched_ion_indices.append(i)

        return matched_ion_indices


    def gen_peak_list(self, assigned_peaks_file):
        all_peaks = pd.read_csv(assigned_peaks_file)
        unassigned_peaks = all_peaks[pd.isna(all_peaks['name'])]
        peak_list = all_peaks["monoisotopic_mw"].to_list()
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


    def find_unimod_matches(self):
        amino_acids = [
            ('Dehydrated Alanine', 71.03711, 'Internal'),
            ('Dehydrated Cysteine', 103.00919, 'Internal'),
            ('Dehydrated Aspartic Acid', 115.02694, 'Internal'),
            ('Dehydrated Glutamic Acid', 129.04259, 'Internal'),
            ('Dehydrated Phenylalanine', 147.06841, 'Internal'),
            ('Dehydrated Glycine', 57.02146, 'Internal'),
            ('Dehydrated Histidine', 137.05891, 'Internal'),
            ('Dehydrated Isoleucine', 113.08406, 'Internal'),
            ('Dehydrated Lysine', 128.09496, 'Internal'),
            ('Dehydrated Leucine', 113.08406, 'Internal'),
            ('Dehydrated Methionine', 131.04049, 'Internal'),
            ('Dehydrated Asparagine', 114.04293, 'Internal'),
            ('Dehydrated Proline', 97.05276, 'Internal'),
            ('Dehydrated Glutamine', 128.05858, 'Internal'),
            ('Dehydrated Arginine', 156.10111, 'Internal'),
            ('Dehydrated Serine', 87.03203, 'Internal'),
            ('Dehydrated Threonine', 101.04768, 'Internal'),
            ('Dehydrated Valine', 99.06841, 'Internal'),
            ('Dehydrated Tryptophan', 186.07931, 'Internal'),
            ('Dehydrated Tyrosine', 163.06333, 'Internal'),
            ('Dehydrated Alanine Loss', -71.03711, 'Internal'),
            ('Dehydrated Cysteine Loss', -103.00919, 'Internal'),
            ('Dehydrated Aspartic Acid Loss', -115.02694, 'Internal'),
            ('Dehydrated Glutamic Acid Loss', -129.04259, 'Internal'),
            ('Dehydrated Phenylalanine Loss', -147.06841, 'Internal'),
            ('Dehydrated Glycine Loss', -57.02146, 'Internal'),
            ('Dehydrated Histidine Loss', -137.05891, 'Internal'),
            ('Dehydrated Isoleucine Loss', -113.08406, 'Internal'),
            ('Dehydrated Lysine Loss', -128.09496, 'Internal'),
            ('Dehydrated Leucine Loss', -113.08406, 'Internal'),
            ('Dehydrated Methionine Loss', -131.04049, 'Internal'),
            ('Dehydrated Asparagine Loss', -114.04293, 'Internal'),
            ('Dehydrated Proline Loss', -97.05276, 'Internal'),
            ('Dehydrated Glutamine Loss', -128.05858, 'Internal'),
            ('Dehydrated Arginine Loss', -156.10111, 'Internal'),
            ('Dehydrated Serine Loss', -87.03203, 'Internal'),
            ('Dehydrated Threonine Loss', -101.04768, 'Internal'),
            ('Dehydrated Valine Loss', -99.06841, 'Internal'),
            ('Dehydrated Tryptophan Loss', -186.07931, 'Internal'),
            ('Dehydrated Tyrosine Loss', -163.06333, 'Internal')
        ]


        xml_file_path = './unimod.xml'
        mod_list = self.parse_xml(xml_file_path) + amino_acids

        set_mass = self.offset

        filtered_mod_list = [
            (index, mass, comp)
            for index, mass, comp in mod_list
            if (set_mass - 0.1) <= mass <= (set_mass + 0.1)
        ]

        differences = [abs(mass - set_mass) for _, mass, _ in filtered_mod_list]

        combined_data = list(zip(filtered_mod_list, differences))
        sorted_combined_data = sorted(combined_data, key=lambda x: x[1])

        sorted_mod_tuples = [item[0] for item in sorted_combined_data]
        differences = [item[1] for item in sorted_combined_data]

        mod_names = [item[0] for item in sorted_mod_tuples]
        mod_masses = [item[1] for item in sorted_mod_tuples]
        mod_comps = [item[2] for item in sorted_mod_tuples]

        if not mod_names:
            df = pd.DataFrame({
                'mod': ["No matching modifications"],
                'mono_mass': [0],
                'comp': [0],
                'difference': [0]
            })
        else:
            df = pd.DataFrame({
                'mod': mod_names,
                'comp': mod_comps,
                'mono_mass': mod_masses,
                'difference': differences
            })

        return df


    def concat_unique_masses(self, list1, list2):
        unique_masses = {}
        result = list1

        for name, ion_mass in list1:
            if ion_mass not in unique_masses:
                unique_masses[ion_mass] = True

        # Process the second list
        for name, ion_mass in list2:
            if ion_mass not in unique_masses:
                unique_masses[ion_mass] = True
                result.append((name, ion_mass))

        return result


    def parse_xml(self, xml_file):
        tree = ET.parse(xml_file)
        root = tree.getroot()

        mod_list = []

        for elem in root.findall(
            './/umod:mod',
            namespaces={'umod': 'http://www.unimod.org/xmlns/schema/unimod_2'}
        ):
            accession = elem.get('record_id')
            full_name = elem.get('full_name')

            delta_elem = elem.find(
                './/umod:delta',
                namespaces={'umod': 'http://www.unimod.org/xmlns/schema/unimod_2'}
            )
            mono_mass = (
                float(delta_elem.get('mono_mass'))
                if delta_elem.get('mono_mass') is not None else None
            )
            comp = (
                str(delta_elem.get('composition'))
                if delta_elem.get('mono_mass') is not None else None
            )
            mod_list.append((f"{full_name} [{accession}]", mono_mass, comp))

        return mod_list



class OffsetSearchPlotPanel(wx.Panel):
    def __init__(self, parent, sequence, b_fragment_positions, y_fragment_positions):
        super().__init__(parent)

        run = OffsetIDCalculator(["",""])
        protein_sequence, _ = run.process_sequence_string(sequence)

        internals_per_row = round(len(protein_sequence) / 100) * 6
        font_size = min(10, 6000 / len(protein_sequence))
        num_rows = math.ceil(len(protein_sequence) / internals_per_row)

        self.figure, ax = plt.subplots(
            num_rows,
            1,
            figsize=(font_size/20*10*internals_per_row/20, num_rows * (font_size / 20) * 0.5),
            num="Offset ID"
        )

        self.figure.set_facecolor("#f0f0f0")

        for i, ax_row in enumerate(ax):
            ax_row.set_facecolor("#f0f0f0")
            start = i * internals_per_row
            end = min((i + 1) * internals_per_row, len(protein_sequence))
            sequence_chunk = protein_sequence[start:end]

            if i == 0:
                ax_row.text(
                    -1, 0.5, "N",
                    ha="center", va="center",
                    fontsize=font_size, fontname="Courier New", color="red"
                )
            else:
                ax_row.text(
                    -0.5, 0.5, start+1,
                    ha="right", va="center",
                    fontsize=font_size, fontname="Courier New", color="grey"
                )

            if i == len(ax)-1:
                pass
            else:
                ax_row.text(
                    internals_per_row-0.5, 0.5, end,
                    ha="left", va="center",
                    fontsize=font_size, fontname="Courier New", color="grey"
                )

            pos = 0
            for j, res in enumerate(sequence_chunk):
                ax_row.text(
                    pos, 0.5, res,
                    ha="center", va="center",
                    fontsize=font_size, fontname="Courier New", weight = "bold"
                )
                pos += 1
                if i == len(ax)-1:
                    if j == len(sequence_chunk)-1:
                        ax_row.text(
                            pos, 0.5, "C",
                            ha="center", va="center",
                            fontsize=font_size, fontname="Courier New", color="red"
                        )

            ax_row.set_xlim(-1, internals_per_row+1)
            ax_row.axis("off")

            for position in b_fragment_positions:
                if start <= position-1 < end:
                    x = position - start - 0.5
                    ax_row.text(
                        x, 0.7, "⎤",
                        ha="right", va="center",
                        fontsize=font_size, fontname="monospace", color="#205283"
                    )

            for position in y_fragment_positions:
                if start <= position+1 < end:
                    x = position - start + 0.5
                    ax_row.text(
                        x, 0.3, "⎣",
                        ha="left", va="center",
                        fontsize=font_size, fontname="monospace", color="#b42920"
                    )


        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.toolbar = NavigationToolbar2Wx(self.canvas)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.SetSizer(sizer)
        self.Layout()

        self.figure.tight_layout(rect=[0, 0, 1, 1])

        self.toolbar.Show()
