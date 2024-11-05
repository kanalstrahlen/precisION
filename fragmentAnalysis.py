import os
import ast
import math
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import colors
import matplotlib.cm
from matplotlib.colors import LinearSegmentedColormap
import wx
import wx.grid
from Bio import PDB, pairwise2
import textalloc as ta


matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42



class FragmentAnalysisWindow(wx.Frame):
    def __init__(self, parent, title, file_path, directory_path):
        super().__init__(
            parent,
            title=title,
            size=(375, 380),
            style=wx.DEFAULT_FRAME_STYLE
        )

        self.file_path = file_path
        self.directory_path = directory_path

        self.panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # header_sizer: title and key information
        header_sizer = wx.BoxSizer(wx.VERTICAL)
        title = wx.StaticText(self.panel, label="Fragment Analysis")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext1 = wx.StaticText(
            self.panel,
            label="Analyse the effect of protein structure and sequence on "
        )
        subtext2 = wx.StaticText(
            self.panel,
            label="fragmentation using plots and summary statistics."
        )

        header_sizer.Add(title, 0, wx.BOTTOM, 5)
        header_sizer.Add(subtext1, 0, wx.BOTTOM, 3)
        header_sizer.Add(subtext2, 0, wx.BOTTOM, 5)

        # button_sizer: all buttons in nx2 grid
        button_sizer = wx.BoxSizer(wx.VERTICAL)

        button_subsizer1 = wx.BoxSizer(wx.HORIZONTAL)
        histogram_button = wx.Button(self.panel, label="Intensity Histogram",  size=(170, 30))
        histogram_button.Bind(wx.EVT_BUTTON, self.on_histogram_button)
        frag_position_button = wx.Button(self.panel, label="Fragment Position",  size=(170, 30))
        frag_position_button.Bind(wx.EVT_BUTTON, self.on_frag_position_button)

        button_subsizer1.Add(histogram_button, 0, wx.RIGHT, 5)
        button_subsizer1.Add(frag_position_button, 0, wx.RIGHT, 0)

        button_subsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        sequence_map_button = wx.Button(self.panel, label="Sequence Map",  size=(170, 30))
        sequence_map_button.Bind(wx.EVT_BUTTON, self.on_sequence_map_button)
        structure_button = wx.Button(self.panel, label="3D Structure",  size=(170, 30))
        structure_button .Bind(wx.EVT_BUTTON, self.on_structure_button )

        button_subsizer2.Add(sequence_map_button, 0, wx.RIGHT, 5)
        button_subsizer2.Add(structure_button, 0, wx.RIGHT, 0)

        button_subsizer3 = wx.BoxSizer(wx.HORIZONTAL)
        heatmap_button = wx.Button(self.panel, label="Residue Pair Heatmap",  size=(170, 30))
        heatmap_button.Bind(wx.EVT_BUTTON, self.on_heatmap_button)
        mod_finder_button = wx.Button(self.panel, label="Localise Modifications",  size=(170, 30))
        mod_finder_button.Bind(wx.EVT_BUTTON, self.on_mod_finder_button)

        button_subsizer3.Add(heatmap_button, 0, wx.RIGHT, 5)
        button_subsizer3.Add(mod_finder_button, 0, wx.RIGHT, 0)

        button_subsizer4 = wx.BoxSizer(wx.HORIZONTAL)
        spectrum_button = wx.Button(self.panel, label="Annotated Spectrum",  size=(170, 30))
        spectrum_button.Bind(wx.EVT_BUTTON, self.on_spectrum_button)
        statistics_button = wx.Button(self.panel, label="Summary Statistics",  size=(170, 30))
        statistics_button.Bind(wx.EVT_BUTTON, self.on_statistics_button)

        button_subsizer4.Add(spectrum_button, 0, wx.RIGHT, 5)
        button_subsizer4.Add(statistics_button, 0, wx.RIGHT, 0)

        button_sizer.Add(button_subsizer1, 0, wx.BOTTOM, 5)
        button_sizer.Add(button_subsizer2, 0, wx.BOTTOM, 5)
        button_sizer.Add(button_subsizer3, 0, wx.BOTTOM, 5)
        button_sizer.Add(button_subsizer4, 0, wx.BOTTOM, -10)

        # options_sizer: load pdb file and select options
        options_staticbox = wx.StaticBox(self.panel)
        options_sizer = wx.StaticBoxSizer(options_staticbox, wx.VERTICAL)
        options_sizer_header = wx.StaticText(self.panel, label="PDB File and Options")
        options_sizer_header.SetFont(
            wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)
        )

        # pdb_subsizer: select pdb file
        pdb_subsizer = wx.BoxSizer(wx.HORIZONTAL)
        pdb_label = wx.StaticText(self.panel, label='PDB File:')
        self.pdb_text = wx.TextCtrl(
            self.panel,
            size=(232, 20),
            style=wx.TE_READONLY | wx.TE_CENTER
        )
        pdb_button = wx.Button(self.panel, label="Browse",  size=(50, 25))
        pdb_button.Bind(wx.EVT_BUTTON, self.on_pdb_button)

        pdb_subsizer.Add(pdb_label, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        pdb_subsizer.Add(self.pdb_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        pdb_subsizer.Add(pdb_button, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)

        # name_subsizer: select protein in assigned peaks file for analysis
        name_subsizer = wx.BoxSizer(wx.HORIZONTAL)
        name_label = wx.StaticText(self.panel, label='Protein Name:')

        # read protein names from assigned peak list
        peak_list = pd.read_csv(file_path)
        peak_list = peak_list.dropna(subset=['name'])
        name_choices = peak_list['name'].unique()
        self.name_dropdown = wx.ComboBox(
            self.panel,
            choices=name_choices,
            style=wx.CB_READONLY
        )
        self.name_dropdown.SetValue(name_choices[0])
        self.name_dropdown.Bind(wx.EVT_COMBOBOX, self.on_name_select)

        color_label = wx.StaticText(self.panel, label='Color:')
        self.color_dropdown = wx.ComboBox(
            self.panel,
            choices=list(colors.CSS4_COLORS),
            style=wx.CB_READONLY
        )
        self.color_dropdown.SetValue("black")

        name_subsizer.Add(name_label, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        name_subsizer.Add(self.name_dropdown, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        name_subsizer.Add(color_label, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        name_subsizer.Add(self.color_dropdown, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)

        mod_subsizer = wx.BoxSizer(wx.HORIZONTAL)

        mod_label = wx.StaticText(self.panel, label="Modification:")
        self.current_name = self.name_dropdown.GetValue()
        name_data = peak_list[peak_list['name'] == self.current_name]
        mod_choices = name_data['variable_mod'].dropna().unique()
        if len(mod_choices) == 0:
            mod_choices = ["None"]

        self.mod_dropdown = wx.ComboBox(
            self.panel,
            choices=mod_choices,
            style=wx.CB_READONLY
        )
        self.mod_dropdown.SetValue(mod_choices[0])

        mod_label2 = wx.StaticText(self.panel, label="On residues:")
        self.mod_res = wx.TextCtrl(self.panel, size=(110, 20))
        self.mod_res.SetValue("S,T,Y")

        mod_subsizer.Add(mod_label, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        mod_subsizer.Add(self.mod_dropdown, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        mod_subsizer.Add(mod_label2, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        mod_subsizer.Add(self.mod_res, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)

        options_sizer.Add(options_sizer_header, 0, wx.BOTTOM, 5)
        options_sizer.Add(pdb_subsizer, 0, wx.BOTTOM, 5)
        options_sizer.Add(name_subsizer, 0, wx.BOTTOM, 5)
        options_sizer.Add(mod_subsizer, 0, wx.BOTTOM, 0)

        main_sizer.Add(header_sizer, 0, wx.ALL, 5)
        main_sizer.Add(button_sizer, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 5)
        main_sizer.Add(options_sizer, 0, wx.ALL,  5)

        icon = wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Show()


    def on_size(self, event):
        self.panel.Layout()
        event.Skip()


    def on_name_select(self, event):
        self.current_name = self.name_dropdown.GetValue()
        peak_list = pd.read_csv(self.file_path)

        name_data = peak_list[peak_list['name'] == self.current_name]
        mod_choices = name_data['variable_mod'].dropna().unique()
        if len(mod_choices) == 0:
            mod_choices = ["None"]
        self.mod_dropdown.Set(mod_choices)
        self.mod_dropdown.SetValue(mod_choices[0])


    def on_pdb_button(self, event):
        wildcard = "PDB files (*.pdb)|*.pdb|" "All files (*.*)|*.*"
        dlg = wx.FileDialog(self, "Choose a file", wildcard=wildcard, style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            file_path = dlg.GetPath()
            self.pdb_text.SetValue(file_path)
            print("Selected file:", file_path)
        dlg.Destroy()


    def on_histogram_button(self, event):
        name = self.name_dropdown.GetValue()
        color = self.color_dropdown.GetValue()
        pdb_file_path = self.pdb_text.GetValue()
        peak_list = pd.read_csv(self.file_path)
        histogram = IntensityHistogram(name, peak_list, color, pdb_file_path)
        histogram.gen_histogram()


    def on_frag_position_button(self, event):
        name = self.name_dropdown.GetValue()
        peak_list = pd.read_csv(self.file_path)
        color = self.color_dropdown.GetValue()
        frag_pos_plot = FragmentPosition(name, peak_list, color)
        frag_pos_plot.gen_frag_positions()


    def on_sequence_map_button(self, event):
        name = self.name_dropdown.GetValue()
        peak_list = pd.read_csv(self.file_path)
        seq_map = SequenceMap(name, peak_list)
        seq_map.gen_sequence_map()


    def on_structure_button(self, event):
        name = self.name_dropdown.GetValue()
        pdb_file_path = self.pdb_text.GetValue()
        peak_list = pd.read_csv(self.file_path)

        if os.path.exists(pdb_file_path):
            pass
        else:
            print("No pdb file provided.")
            return

        structure = Structure(name, peak_list, pdb_file_path)
        structure.gen_structure()


    def on_heatmap_button(self, event):
        name = self.name_dropdown.GetValue()
        color = self.color_dropdown.GetValue()
        peak_list = pd.read_csv(self.file_path)
        heatmap = Heatmap(name, peak_list, color)
        heatmap.gen_heatmap()


    def on_mod_finder_button(self, event):
        name = self.name_dropdown.GetValue()
        color = self.color_dropdown.GetValue()
        peak_list = pd.read_csv(self.file_path)
        mod = self.mod_dropdown.GetValue()
        res = self.mod_res.GetValue()

        if mod == "None":
            print("No modification selected.")
            return
        mod_finder_plot = ModFinder(name, peak_list, mod, res, color)
        mod_finder_plot.gen_plot()


    def on_spectrum_button(self, event):
        name = self.name_dropdown.GetValue()
        colour = self.color_dropdown.GetValue()
        plot = AnnotatedSpectrum()
        plot.create_plot(self.file_path, self.directory_path, name, colour)


    def on_statistics_button(self, event):
        name = self.name_dropdown.GetValue()
        basename = os.path.basename(self.file_path).replace(".assignedPeakList.csv", "")
        spectrum_path = os.path.join(
            self.directory_path,
            f"{basename}.centroid.cwt.txt"
        )
        run = FragmentStatistics(self.file_path, spectrum_path, name)
        print("\nSpectrum statistics")
        run.calc_pc_peaks_assigned()
        run.calc_pc_intens_assigned()
        run.calc_pc_tic_assigned()
        run.protein_summary()
        print(f"\n\nFor {name}...")
        run.pc_signal()
        run.sequence_coverage()
        print(f"Ion type breakdown for {name}:")
        run.ion_type_breakdown()



class IntensityHistogram():
    def __init__(self, name, peak_list, color, pdb_file_path):
        self.name = name
        self.peak_list = peak_list
        self.color = color
        self.pdb_file_path = pdb_file_path


    def gen_histogram(self):
        name_data = self.peak_list[self.peak_list['name'] == self.name]
        sequence = name_data['sequence'].tolist()[0]

        resid_array = np.zeros(len(sequence))
        for i, _ in enumerate(resid_array):
            resid_array[i] = i + 1

        b_intensity_array = np.zeros(len(sequence))
        y_intensity_array = np.zeros(len(sequence))

        b_ion_index = []
        b_ion_charge = []
        y_ion_index = []
        y_ion_charge = []

        for _, row in name_data.iterrows():
            ion_name = row["ion"]
            ion_serie = ion_name[0]
            ion_abundance = row["charge_scaled_abundance"]
            ion_charge = row["charge"]

            if ion_serie == "I":
                continue

            end_index = int(ion_name.find(' '))
            ion_index = int(ion_name[1:end_index])

            if ion_serie in ("b", "c"):
                b_intensity_array[ion_index - 1] += ion_abundance
                b_ion_index.append(ion_index)
                b_ion_charge.append(ion_charge)
            elif ion_serie in ("y", "z"):
                y_intensity_array[len(sequence) - ion_index - 1] += ion_abundance
                y_ion_index.append(len(sequence) - ion_index)
                y_ion_charge.append(ion_charge)

        self.plot_histogram(
            b_intensity_array, b_ion_index, b_ion_charge,
            y_intensity_array, y_ion_index, y_ion_charge, resid_array
        )


    def plot_histogram(
        self, b_intensity_array, b_ion_index, b_ion_charge,
        y_intensity_array, y_ion_index, y_ion_charge, resid_array
    ):
        #if os.path.exists(self.pdb_file_path):
        #    fig, ax = plt.subplots(
        #        3, 1,
        #        gridspec_kw={'height_ratios': [0.2, 1, 0.1], 'hspace': 0.05},
        #        num=f"precisION - Intensity Histogram [{self.name}]",
        #        sharex=True
        #    )
        #else:
        _, ax = plt.subplots(
            2, 1,
            gridspec_kw={'height_ratios': [0.2, 1], 'hspace': 0.05},
            num=f"precisION - Intensity Histogram [{self.name}]",
            sharex=True
        )

        ax[0].scatter(
            b_ion_index, b_ion_charge,
            facecolor='#205283', edgecolors='black',
            marker="s", alpha=0.3
        )
        ax[0].scatter(
            y_ion_index, y_ion_charge,
            facecolor='#b42920', edgecolors='black',
            marker="s", alpha=0.3
        )

        ax[0].set_ylabel("Charge")
        ax[0].set_yticks(
            [int(min(b_ion_charge + y_ion_charge)), int(max(b_ion_charge + y_ion_charge))]
        )
        ax[0].set_ylim(0, max(b_ion_charge + y_ion_charge) + 1)
        ax[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax[0].set_xlim(min(resid_array), max(resid_array))

        summed_intensity_array = b_intensity_array + y_intensity_array
        b_intensity_array = b_intensity_array / max(summed_intensity_array) * 100
        y_intensity_array = y_intensity_array / max(summed_intensity_array) * 100

        ax[1].bar(resid_array, b_intensity_array, width = 1, color="#205283")
        ax[1].bar(
            resid_array, y_intensity_array,
            bottom=b_intensity_array, width = 1, color='#b42920'
        )
        ax[1].set_xlim(min(resid_array), max(resid_array))
        ax[1].set_xlabel("Residue number")
        ax[1].set_ylabel("Relative intensity (%)")

# to do; fix errors and incorporate this sometime
#        if os.path.exists(self.pdb_file_path):
#            self.renumber_and_delete(pdb_file_path, sequence)
#            self.atomnum = {' N  ':0, ' CA ': 1, ' C  ': 2, ' O  ': 3}
#            coord = torch.tensor(self.read_pdbtext(open(pdb_file_path, 'r').read()))
#            numbers = pydssp.assign(coord, out_type='index').tolist()
#
#            axs[1].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
#            axs[1].set_xlabel("")
#
#            time_points = np.arange(0, len(numbers))
#            for i, (t, num) in enumerate(zip(time_points, numbers)):
#                len_plot = len(numbers)
#                prev_number = numbers[i-1]
#                try:
#                    next_num = numbers[i+1]
#                except:
#                    pass
#                if num == 0:
#                    if prev_number == 1 and next_num == 1:
#                        x = np.linspace(t, t + 1, 1000)
#                        y = 0.1 * np.sin(np.pi * (x - t_zero))  # Sine wave
#                        axs[2].plot(x, y, color='black', linewidth=1)
#                    elif prev_number == 2 and next_num == 2:
#                        axs[2].add_patch(plt.Rectangle((t, -0.05), 1, 0.1, color='black'))
#                    else:
#                        axs[2].plot([t, t+1], [0, 0], color='black', linewidth=1)  # Straight line
#                elif num == 1 and prev_number != 1:
#                    t_zero = t
#                    x = np.linspace(t, t + 1, 1000)
#                    y = 0.1 * np.sin(np.pi * (x - t_zero))  # Sine wave
#                    axs[2].plot(x, y, color='black', linewidth=1)
#                elif num == 1 and prev_number == 1:
#                    x = np.linspace(t, t + 1, 1000)
#                    y = 0.1 * np.sin(np.pi * (x - t_zero))  # Sine wave
#                    axs[2].plot(x, y, color='black', linewidth=1)
#                elif num == 2:
#                    axs[2].add_patch(plt.Rectangle((t, -0.05), 1, 0.1, color='black'))  # Rectangle
#
#            ax[2].set_xlim(0, len(numbers))
#            ax[2].set_ylim(-0.12, 0.12)  # Adjust the y-axis limits as needed
#            ax[2].set_xlabel("Residue number")
#
#            ax[2].spines['top'].set_visible(False)
#            ax[2].spines['right'].set_visible(False)
#            ax[2].spines['bottom'].set_visible(False)
#            ax[2].spines['left'].set_visible(False)
#            ax[2].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.SetIcon(wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO))

        plt.show()


class FragmentPosition():
    def __init__(self, name, peak_list, color):
        self.name = name
        self.peak_list = peak_list
        self.color = color


    def gen_frag_positions(self):
        name_data = self.peak_list[self.peak_list['name'] == self.name]
        sequence = name_data['sequence'].tolist()[0]

        resid_array = np.zeros(len(sequence))
        for i, _ in enumerate(resid_array):
            resid_array[i] = i + 1

        b_ion_list = []
        i_ion_list = []
        y_ion_list = []

        for _, row in name_data.iterrows():
            ion_name = row["ion"]
            ion_serie = ion_name[0]

            end_index = int(ion_name.find(' '))
            ion_index = ion_name[1:end_index]

            if ion_serie in ("b", "c"):
                pos_tuple = (1, int(ion_index))
                b_ion_list.append(pos_tuple)
            elif ion_serie == "I":
                indices = ion_index.split("-")
                start = int(indices[0]) + 1
                end = int(indices[1])
                pos_tuple = (start, end)
                i_ion_list.append(pos_tuple)
            elif ion_serie in ("y", "z"):
                pos_tuple = (len(sequence) - int(ion_index) + 1, len(sequence))
                y_ion_list.append(pos_tuple)

        i_ion_list = sorted(i_ion_list, key=lambda x: x[0])

        b_ion_counts = Counter(b_ion_list)
        i_ion_counts = Counter(i_ion_list)
        y_ion_counts = Counter(y_ion_list)

        b_unique_ions = list(b_ion_counts.keys())
        i_unique_ions = list(i_ion_counts.keys())
        y_unique_ions = list(y_ion_counts.keys())

        b_ion_occurrences = list(b_ion_counts.values())
        i_ion_occurrences = list(i_ion_counts.values())
        y_ion_occurrences = list(y_ion_counts.values())

        sorted_unique_ions = b_unique_ions + i_unique_ions + y_unique_ions[::-1]

        sorted_ion_occurences = b_ion_occurrences + i_ion_occurrences + y_ion_occurrences[::-1]

        self.plot_fragment_positions(sorted_unique_ions, sorted_ion_occurences, resid_array)


    def plot_fragment_positions(self, sorted_unique_ions, sorted_ion_occurences, resid_array):
        def interpolate_color(color, factor=0.5):
            rgb_color = colors.to_rgba(color)[:3]
            light_interpolated_color = [(1 - factor) + factor * c for c in rgb_color]
            dark_interpolated_color = [factor * c for c in rgb_color]

            return light_interpolated_color, dark_interpolated_color

        light_color, dark_color = interpolate_color(self.color, 0.25)
        custom_cmap = LinearSegmentedColormap.from_list(
            'custom_colormap', [light_color, dark_color], N=max(sorted_ion_occurences)
        )

        norm = plt.Normalize(min(sorted_ion_occurences)-0.5, max(sorted_ion_occurences)+0.5)

        # create a scalarmappable to map values to colors
        sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
        sm.set_array([])

        _, ax = plt.subplots(1, 1, num=f"precisION - Fragment Position [{self.name}]")
        index =  1
        for pos_tuple, occurrence in zip(sorted_unique_ions, sorted_ion_occurences):
            bar_color = sm.to_rgba(occurrence)
            ax.plot([pos_tuple[0], pos_tuple[1]], [index, index], color=bar_color)
            index += 1

        ax.set_xlim(1, max(resid_array))
        ax.set_ylim(0, len(sorted_unique_ions)+1)
        ax.set_xlabel("Residue number")
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('# of ions observed')
        cbar.set_ticks(range(min(sorted_ion_occurences), max(sorted_ion_occurences)+1))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.SetIcon(wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO))
        plt.tight_layout()

        plt.show()



class SequenceMap():
    def __init__(self, name, peak_list):
        self.name = name
        self.peak_list = peak_list


    def gen_sequence_map(self):
        name_data = self.peak_list[self.peak_list['name'] == self.name]

        sequence = name_data['sequence'].tolist()[0]
        n_mod = name_data['n_mod'].tolist()[0]
        c_mod = name_data['c_mod'].tolist()[0]

        b_fragment_positions = []
        y_fragment_positions = []
        for _, row in name_data.iterrows():
            ion_name = row["ion"]
            ion_serie = ion_name[0]
            if ion_serie != "I":
                end_index = int(ion_name.find(' '))
                ion_index = int(ion_name[1:end_index])

            if ion_serie in ("b", "c"):
                b_fragment_positions.append(ion_index)
            elif ion_serie in ("y", "z"):
                y_fragment_positions.append(len(sequence) - ion_index - 1)

        b_fragment_positions = list(set(b_fragment_positions))
        y_fragment_positions = list(set(y_fragment_positions))

        if n_mod != "(+)  (-) ":
            n_mod = True
        else:
            n_mod = False

        if c_mod != "(+)  (-) ":
            c_mod = True
        else:
            c_mod = False

        if len(sequence) >= 400:
            residues_per_row = 40
        else:
            residues_per_row = 20

        self.plot_sequence_map(
            self.name,
            sequence,
            b_fragment_positions,
            y_fragment_positions,
            n_mod,
            c_mod,
            residues_per_row
        )


    def plot_sequence_map(
        self, name, protein_sequence, b_fragment_positions, y_fragment_positions,
        n_term_mod=False, c_term_mod=False, residues_per_row=20
    ):
        num_rows = math.ceil(len(protein_sequence) / residues_per_row)

        font_size = 20
        _, ax = plt.subplots(
            num_rows, 1,
            figsize=(font_size/20*10*residues_per_row/20, num_rows * (font_size / 20) * 0.5),
            num=f"precisION - Sequence Map [{name}]"
        )

        for i, ax_row in enumerate(ax):
            start = i * residues_per_row
            end = min((i + 1) * residues_per_row, len(protein_sequence))
            sequence_chunk = protein_sequence[start:end]

            if c_term_mod:
                c_term_colour = "red"
            else:
                c_term_colour = "grey"

            if i == 0:
                if n_term_mod:
                    n_term_colour = "red"
                else:
                    n_term_colour = "grey"
                ax_row.text(
                    -1, 0.5, "N",
                    ha="center", va="center",
                    fontsize=font_size, fontname="Courier New",
                    color=n_term_colour
                )
            else:
                ax_row.text(
                    -0.5, 0.5, start+1,
                    ha="right", va="center",
                    fontsize=font_size, fontname="Courier New",
                    color="grey"
                )

            if i == len(ax)-1:
                pass
            else:
                ax_row.text(
                    residues_per_row-0.5, 0.5, end,
                    ha="left", va="center",
                    fontsize=font_size, fontname="Courier New",
                    color="grey"
                )

            pos = 0
            for j, res in enumerate(sequence_chunk):
                ax_row.text(
                    pos, 0.5, res,
                    ha="center", va="center",
                    fontsize=font_size, fontname="Courier New",
                    weight="bold"
                )
                pos += 1
                if i == len(ax)-1:
                    if j == len(sequence_chunk)-1:
                        ax_row.text(
                            pos, 0.5, "C",
                            ha="center", va="center",
                            fontsize=font_size, fontname="Courier New",
                            color=c_term_colour
                        )

            ax_row.set_xlim(-1, residues_per_row+1)
            ax_row.axis("off")

            for position in b_fragment_positions:
                if start <= position-1 < end:
                    x = position - start - 0.5
                    ax_row.text(
                        x, 0.7, "⎤",
                        ha="right", va="center",
                        fontsize=font_size, fontname="monospace",
                        color="#205283"
                    )

            for position in y_fragment_positions:
                if start <= position+1 < end:
                    x = position - start + 0.5
                    ax_row.text(
                        x, 0.3, "⎣",
                        ha="left", va="center",
                        fontsize=font_size, fontname="monospace",
                        color="#b42920"
                    )

        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.SetIcon(wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO))

        plt.tight_layout()
        plt.show()



class Structure():
    def __init__(self, name, peak_list, pdb_file_path):
        self.name = name
        self.peak_list = peak_list
        self.pdb_file_path = pdb_file_path


    def gen_structure(self):
        name_data = self.peak_list[self.peak_list['name'] == self.name]
        sequence = name_data['sequence'].tolist()[0]

        self.renumber_and_delete(self.pdb_file_path, sequence)

        b_intensity_array = np.zeros(len(sequence) - 1)
        y_intensity_array = np.zeros(len(sequence) - 1)

        for _, row in name_data.iterrows():
            ion_name = row["ion"]
            ion_serie = ion_name[0]
            ion_abundance = row["charge_scaled_abundance"]
            if ion_serie == "I":
                continue
            end_index = int(ion_name.find(' '))
            ion_index = int(ion_name[1:end_index])

            if ion_serie in ("b", "c"):
                b_intensity_array[ion_index - 1] += ion_abundance
            elif ion_serie in ("y", "z"):
                y_intensity_array[len(sequence) - ion_index - 1] += ion_abundance

        intensity_array = np.zeros(len(sequence) - 1)

        for i, b_ion_intens in enumerate(b_intensity_array):
            if b_ion_intens >= y_intensity_array[i]:
                y_intensity_array[i] = 0.0
            elif b_ion_intens <= y_intensity_array[i]:
                b_intensity_array[i] = 0.0

            intensity_array[i] = y_intensity_array[i] - b_intensity_array[i]

        min_value = abs(min(intensity_array))
        max_value = abs(max(intensity_array))
        norm_value = max([min_value, max_value])
        intensity_array = intensity_array / norm_value * 100

        self.edit_bfactor(self.pdb_file_path, intensity_array)


    def renumber_and_delete(self, pdb_file_path, sequence):
        # load PDB structure
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file_path)

        # extract sequence from PDB structure
        pdb_sequence = ''.join(
            [PDB.Polypeptide.index_to_one(
                PDB.Polypeptide.three_to_index(residue.get_resname())
            ) for model in structure for chain in model for residue in chain]
        )

        sequence = sequence.upper()

        if len(sequence) > len(pdb_sequence):
            print("ERROR. Sequence should be a subset of the structure.")
            return

        # perform sequence alignment
        alignments = pairwise2.align.globalxx(
            pdb_sequence,
            sequence,
            one_alignment_only=True,
            gap_char='-'
        )
        best_alignment = alignments[0]

        protein_seq = best_alignment.seqB

        if protein_seq.startswith('-'):
            leading_gap_length = 0
            for char in protein_seq:
                if char == '-':
                    leading_gap_length += 1
                else:
                    break
            alignment_start_sequence = leading_gap_length + 1
        else:
            alignment_start_sequence = 1

        alignment_end_sequence = len(best_alignment.seqB.rstrip('-'))

        new_residue_number = 1
        delete_residue_number = 999
        to_be_deleted = set()

        for model in structure:
            for chain in model:
                for residue in chain:
                    if alignment_start_sequence <= residue.id[1] <= alignment_end_sequence:
                        residue.id = (' ', int(new_residue_number), ' ')
                        new_residue_number += 1
                    else:
                        residue.id = (' ', int(delete_residue_number), ' ')
                        delete_residue_number += 1

        new_structure = PDB.Structure.Structure("filtered_structure")
        for model in structure:
            new_model = PDB.Model.Model(model.id)
            for chain in model:
                new_chain = PDB.Chain.Chain(chain.id)
                for residue in chain:
                    if residue.id[1] < 999:
                        new_chain.add(residue)
                new_model.add(new_chain)
            new_structure.add(new_model)

        io = PDB.PDBIO()
        io.set_structure(new_structure)
        output_path = pdb_file_path.replace(".pdb", ".renumbered.pdb")
        io.save(output_path)

        print(f"Renumbered residues in pdb file to match sequence. Saved as {output_path}.")


    def read_pdbtext(self, pdbstring: str):
        lines = pdbstring.split("\n")
        coords, atoms, resid_old, check = [], None, None, []
        for l in lines:
            if l.startswith('ATOM'):
                atomnum = {' N  ':0, ' CA ': 1, ' C  ': 2, ' O  ': 3}
                iatom = atomnum.get(l[12:16], None)
                resid = l[21:26]
                if resid != resid_old:
                    if atoms is not None:
                        coords.append(atoms)
                        check.append(atom_check)
                    atoms, resid_old, atom_check = [], resid, []
                if iatom is not None:
                    xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    atoms.append(xyz)
                    atom_check.append(iatom)
        if atoms is not None:
            coords.append(atoms)
            check.append(atom_check)
        coords = np.array(coords)
        # check
        assert len(coords.shape) == 3, "Some required atoms [N,CA,C,O] are missing in the input PDB file"
        check = np.array(check)
        assert np.all(check[:,0]==0), "Order of PDB line may be broken. It's required to be N->CA->C->O w/o any duplicate or lack"
        assert np.all(check[:,1]==1), "Order of PDB line may be broken. It's required to be N->CA->C->O w/o any duplicate or lack"
        assert np.all(check[:,2]==2), "Order of PDB line may be broken. It's required to be N->CA->C->O w/o any duplicate or lack"
        assert np.all(check[:,3]==3), "Order of PDB line may be broken. It's required to be N->CA->C->O w/o any duplicate or lack"

        return coords


    def edit_bfactor(self, pdb_file_path, b_factor_array):
        parser = PDB.PDBParser(QUIET=True)
        renumbered_path = pdb_file_path.replace(".pdb", ".renumbered.pdb")
        structure = parser.get_structure('protein', renumbered_path)
        for _, model in enumerate(structure):
            for _, chain in enumerate(model):
                for k, residue in enumerate(chain):
                    for _, atom in enumerate(residue):
                        if k < len(b_factor_array):
                            atom.bfactor = b_factor_array[k]
                        else:
                            atom.bfactor = 0
        output_path = pdb_file_path.replace(".pdb", ".bfactor.pdb")
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_path)
        print(
            f"Saved output as {output_path}. "
            "Please view in the molecular graphics package of your choosing and color by B-factor."
        )



class Heatmap():
    def __init__(self, name, peak_list, color):
        self.name = name
        self.peak_list = peak_list
        self.color = color


    def gen_heatmap(self):
        name_data = self.peak_list[self.peak_list['name'] == self.name]
        sequence = name_data['sequence'].tolist()[0]

        ions = name_data['ion'].tolist()
        frag_sites = name_data['frag_site'].tolist()
        abundances = name_data['charge_scaled_abundance'].tolist()

        amino_acid_pairs = []
        corrected_abundances = []

        for i, ion in enumerate(ions):
            ion_type = ion[0]
            if ion_type != 'i':
                amino_acid_pairs.append(frag_sites[i])
                corrected_abundances.append(abundances[i])
            elif ion_type == "I":
                frag_site_list = frag_sites[i].split(";")
                frag_site_1 = frag_site_list[0]
                frag_site_2 = frag_site_list[1]
                amino_acid_pairs.append(frag_site_1)
                corrected_abundances.append(abundances[i] / 2)
                amino_acid_pairs.append(frag_site_2)
                corrected_abundances.append(abundances[i] / 2)

        tuple_list = list(zip(amino_acid_pairs, corrected_abundances))

        all_amino_acids = [
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
        ]

        heatmap_data = np.zeros((len(all_amino_acids), len(all_amino_acids)))
        for pair, intensity in tuple_list:
            if type(pair) != float:
                i, j = all_amino_acids.index(pair[0]), all_amino_acids.index(pair[1])
                heatmap_data[i, j] += intensity

        # correct for the number of potential sites
        for i, _ in enumerate(all_amino_acids):
            for j, _ in enumerate(all_amino_acids):
                pair = all_amino_acids[i] + all_amino_acids[j]
                num_pair = sequence.count(pair)
                if num_pair >= 1:
                    heatmap_data[i, j] = heatmap_data[i, j] / sequence.count(pair)

        # normalise intensities
        max_intensity = 0
        for i, _ in enumerate(all_amino_acids):
            for j, _ in enumerate(all_amino_acids):
                if heatmap_data[i, j] >= max_intensity:
                    max_intensity = heatmap_data[i, j]
        for i, _ in enumerate(all_amino_acids):
            for j, _ in enumerate(all_amino_acids):
                heatmap_data[i, j] = heatmap_data[i, j] / max_intensity * 100

        self.plot_heatmap(heatmap_data, all_amino_acids)


    def plot_heatmap(self, heatmap_data, all_amino_acids):
        _, ax = plt.subplots(1, 1, num=f"precisION - Residue Pair Heatmap [{self.name}]")
        custom_cmap = LinearSegmentedColormap.from_list(
            'custom_colormap', ['white', self.color], N=256
        )
        heat_map = ax.imshow(heatmap_data, cmap=custom_cmap, vmin=0)
        ax.set_xticks(range(len(all_amino_acids)), all_amino_acids)
        ax.set_yticks(range(len(all_amino_acids)), all_amino_acids)
        ax.set_xlabel('C-terminal residue')
        ax.set_ylabel('N-terminal residue')

        cbar = plt.colorbar(heat_map, ax=ax, orientation='vertical', fraction = 0.05)
        cbar.set_label('Relative, normalised intensity (%)')
        plt.tight_layout(rect=[0, 0, 1, 0.97])

        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.SetIcon(wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO))
        plt.tight_layout()

        plt.show()



class ModFinder():
    def __init__(self, name, peak_list, mod, res, color):
        self.name = name
        self.peak_list = peak_list
        self.mod = mod
        self.res = res
        self.color = color


    def gen_plot(self):
        name_data = self.peak_list[self.peak_list['name'] == self.name]
        no_mod_data = name_data[name_data["variable_mod"].isna()]
        mod_data = name_data[name_data["variable_mod"] == self.mod]
        sequence = name_data['sequence'].tolist()[0]
        sequence = sequence.upper()

        resid_array = np.zeros(len(sequence))
        for i, _ in enumerate(resid_array):
            resid_array[i] = i + 1

        unmod_b_ion_list = []
        unmod_i_ion_list = []
        unmod_y_ion_list = []

        for _, row in no_mod_data.iterrows():
            ion_name = row["ion"]
            ion_serie = ion_name[0]
            end_index = int(ion_name.find(' '))
            ion_index = ion_name[1:end_index]
            if ion_serie in ("b", "c"):
                pos_tuple = (1, int(ion_index))
                unmod_b_ion_list.append(pos_tuple)
            elif ion_serie == "I":
                indices = ion_index.split("-")
                start = int(indices[0]) + 1
                end = int(indices[1])
                pos_tuple = (start, end)
                unmod_i_ion_list.append(pos_tuple)
            elif ion_serie in ("y", "z"):
                pos_tuple = (len(sequence) - int(ion_index) + 1, len(sequence))
                unmod_y_ion_list.append(pos_tuple)

        unmod_i_ion_list = sorted(unmod_i_ion_list, key=lambda x: x[0])
        unmod_b_ion_counts = Counter(unmod_b_ion_list)
        unmod_i_ion_counts = Counter(unmod_i_ion_list)
        unmod_y_ion_counts = Counter(unmod_y_ion_list)
        unmod_b_unique_ions = list(unmod_b_ion_counts.keys())
        unmod_i_unique_ions = list(unmod_i_ion_counts.keys())
        unmod_y_unique_ions = list(unmod_y_ion_counts.keys())
        sorted_unmod_ions = unmod_b_unique_ions + unmod_i_unique_ions + unmod_y_unique_ions[::-1]

        mod_b_ion_list = []
        mod_i_ion_list = []
        mod_y_ion_list = []

        for _, row in mod_data.iterrows():
            ion_name = row["ion"]
            ion_serie = ion_name[0]
            end_index = int(ion_name.find(' '))
            ion_index = ion_name[1:end_index]
            if ion_serie in ("b", "c"):
                pos_tuple = (1, int(ion_index))
                mod_b_ion_list.append(pos_tuple)
            elif ion_serie == "I":
                indices = ion_index.split("-")
                start = int(indices[0]) + 1
                end = int(indices[1])
                pos_tuple = (start, end)
                mod_i_ion_list.append(pos_tuple)
            elif ion_serie in ("y", "z"):
                pos_tuple = (len(sequence) - int(ion_index) + 1, len(sequence))
                mod_y_ion_list.append(pos_tuple)

        mod_b_ion_counts = Counter(mod_b_ion_list)
        mod_i_ion_counts = Counter(mod_i_ion_list)
        mod_y_ion_counts = Counter(mod_y_ion_list)
        mod_b_unique_ions = list(mod_b_ion_counts.keys())
        mod_i_unique_ions = list(mod_i_ion_counts.keys())
        mod_y_unique_ions = list(mod_y_ion_counts.keys())
        sorted_mod_ions = mod_b_unique_ions + mod_i_unique_ions + mod_y_unique_ions[::-1]

        if len(self.res) > 0:
            residues = self.res.split(",")
        else:
            residues = []

        def find_residues(string, letters):
            indices = []
            for letter in letters:
                index = -1
                while True:
                    index = string.find(letter, index + 1)
                    if index == -1:
                        break
                    indices.append(index + 1)
            return indices

        if len(residues) > 0:
            residue_pos = find_residues(sequence, residues)
        else:
            residue_pos = []

        self.plot_fragment_positions(sorted_unmod_ions, sorted_mod_ions, resid_array, residue_pos)


    def plot_fragment_positions(self, sorted_unmod_ions, sorted_mod_ions, resid_array, residue_pos):
        ratio = len(sorted_mod_ions) / len(sorted_unmod_ions)

        _, (ax1, ax2) = plt.subplots(
            2, 1,
            num=f"precisION - Localise Modifications [{self.name}]",
            sharex=True,
            gridspec_kw={'height_ratios': [ratio, 1]}
        )

        ax1.text(1.00, 0.02, self.mod, transform=ax1.transAxes, ha='right', va='bottom')
        ax2.text(1.00, 0.02, f"No Modification", transform=ax2.transAxes, ha='right', va='bottom')

        index =  1
        for pos_tuple in sorted_unmod_ions:
            ax2.plot([pos_tuple[0], pos_tuple[1]], [index, index], color='black')
            index += 1

        if len(residue_pos) >0:
            for pos in residue_pos:
                ax2.plot([pos, pos], [0, index+1], color="#b42920", linestyle="dotted", alpha=0.5)

        index = 1
        for pos_tuple in sorted_mod_ions:
            ax1.plot([pos_tuple[0], pos_tuple[1]], [index, index], color=self.color)
            index += 1

        if len(residue_pos) > 0:
            for pos in residue_pos:
                ax1.plot([pos, pos], [0, index+1], color="#b42920", linestyle="dotted", alpha=0.5)

        ax1.set_xlim(1, max(resid_array))
        ax1.set_ylim(0, len(sorted_mod_ions) + 1)
        
        ax2.set_xlim(1, max(resid_array))
        ax2.set_ylim(0, len(sorted_unmod_ions) + 1)
        ax2.set_xlabel("Residue number")

        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.SetIcon(wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO))
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.03)

        plt.show()



class FragmentStatistics():
    def __init__(self, assignment_file_path, spectrum_path, name):
        self.assignment_df = pd.read_csv(assignment_file_path)
        self.spectrum = np.loadtxt(spectrum_path)
        self.name = name


    def calc_pc_peaks_assigned(self):
        assigned_peaks = self.assignment_df.dropna(subset=['name'])
        pc_assigned = round(len(assigned_peaks) / len(self.assignment_df) * 100, 1)
        print(
            f"{pc_assigned}% ({len(assigned_peaks)}/{len(self.assignment_df)}) "
            "of the envelopes are assigned."
        )


    def calc_pc_intens_assigned(self):
        decon_intens = self.assignment_df['abundance'].sum()
        assigned_peaks = self.assignment_df.dropna(subset=['name'])
        assigned_intens = assigned_peaks['abundance'].sum()
        pc_assigned = round(assigned_intens / decon_intens * 100, 1)
        print(f"Approx. {pc_assigned}% of the deconvoluted signal is assigned.")


    def calc_pc_tic_assigned(self):
        spectrum_y = np.array(self.spectrum[:, 1])
        tic = np.sum(spectrum_y)

        assigned_peaks = self.assignment_df.dropna(subset=['name'])
        total_assigned_intens = assigned_peaks['total_intensity'].sum()

        pc_assigned = round(total_assigned_intens / tic * 100, 1)
        print(f"{pc_assigned}% of the total ion current is assigned.")


    def protein_summary(self):
        protein_df = self.assignment_df.dropna(subset=['name'])
        protein_names = protein_df["name"].unique()
        print(f"{len(protein_names)} protein/s identified.")


    def sequence_coverage(self):
        protein_peaks = self.assignment_df[self.assignment_df['name'] == self.name]
        sequence = protein_peaks["sequence"].iloc[0]
        sequence_length = len(sequence)
        num_potential_frag_sites = sequence_length - 1
        binary_list = [0] * num_potential_frag_sites
        ion_names = protein_peaks["ion"].unique()
        for ion in ion_names:
            ion_name = ion.split()[0]
            ion_type = ion_name[0]
            ion_index = ion_name[1:]
            if ion_type in ("b", "c"):
                binary_list[int(ion_index) - 1] = 1
            elif ion_type in ("y", "z"):
                index = sequence_length - int(ion_index) - 1
                binary_list[index] = 1
            elif ion_type == "I":
                indices = ion_index.split("-")
                start = int(indices[0])
                end = int(indices[1])
                binary_list[start - 1] = 1
                binary_list[end - 1] = 1

        num_hits = binary_list.count(1)
        potential_hits = len(binary_list)
        coverage = round(num_hits / potential_hits * 100, 1)

        print(f"Sequence coverage = {coverage}%")


    def pc_signal(self):
        assigned_peaks = self.assignment_df.dropna(subset=['name'])
        protein_peaks = assigned_peaks[assigned_peaks['name'] == self.name]

        total_assigned_intens = assigned_peaks['total_intensity'].sum()
        total_protein_intens = protein_peaks['total_intensity'].sum()
        pc_protein = round(total_protein_intens / total_assigned_intens * 100, 1)
        print(f"{pc_protein}% of the total assigned signal corresponds with {self.name}.")


    def ion_type_breakdown(self):
        protein_peaks = self.assignment_df[self.assignment_df['name'] == self.name]
        ion_types = list(set([ion[0] for ion in protein_peaks['ion'].dropna().unique()]))
        ions = [ion[0] for ion in protein_peaks['ion']]

        for ion_type in ion_types:
            num = ions.count(ion_type)
            print(f"{num} {ion_type}-type ions.")



class AnnotatedSpectrum():
    def create_plot(self, file_path, directory_path, name, colour):
        basename = os.path.basename(file_path).replace(".assignedPeakList.csv", "")
        spectrum_path = os.path.join(
            directory_path,
            f"{basename}.recalibratedSpectrum.txt"
        )

        screen_width, screen_height = wx.DisplaySize()

        peak_list = pd.read_csv(file_path)
        peak_list = peak_list.dropna(subset=['name'])

        name_data = peak_list[peak_list['name'] == name]
        name_data = name_data.sort_values(by="abundance", ascending=False)

        spectrum = np.loadtxt(spectrum_path)
        spectrum_x = np.array(spectrum[:, 0])
        spectrum_y = np.array(spectrum[:, 1])
        spectrum_y = spectrum_y / max(spectrum_y) * 100

        print("Creating plot...")
        fig, ax = plt.subplots(
            1, 1,
            num=f"precisION - Annotated Spectrum [{name}]",
            figsize=(screen_width/150, screen_height/150),
            dpi=100
        )
        ax.set_xlim(min(spectrum_x), max(spectrum_x))
        ax.set_ylim(0, 105)
        ax.set_yticks([0, 50, 100])
        ax.set_xlabel('m/z')
        ax.set_ylabel('Relative intensity (%)')
        def format_thousands(x, pos):
            return "{:,}".format(int(x))
        plt.gca().xaxis.set_major_formatter(FuncFormatter(format_thousands))

        ax.plot(spectrum_x, spectrum_y, color='black')

        labels = []
        label_x = []
        label_y = []
        num = 0
        for _, row in name_data.iterrows():
            ion_name = row['ion']
            adduct = row['adduct']
            loss = row['loss']
            num += 1

            x_plot = np.array([])
            y_plot = np.array([])

            theo_x = ast.literal_eval(row['fitter_theo_x'])
            for x in theo_x:
                peak_index = self.find_closest_index(x, spectrum_x)
                left_limit, right_limit = self.find_increasing_indexes(spectrum_y, peak_index)

                x_plot = np.concatenate((
                    x_plot,
                    spectrum_x[left_limit-1:left_limit],
                    spectrum_x[left_limit:right_limit],
                    spectrum_x[right_limit:right_limit+1]
                ))

                y_plot = np.concatenate((
                    y_plot,
                    np.zeros(1),
                    spectrum_y[left_limit:right_limit],
                    np.zeros(1)
                ))

            if colour != "black":
                ax.plot(x_plot, y_plot, color=colour)

            # label top 10 peaks
            if type(adduct) is float and type(loss) is float and num <= 10:
                labels.append(ion_name)
                label_x.append(min(x_plot))
                label_y.append(max(y_plot))

        print("Optimising annotation placement...")
        ta.allocate_text(fig, ax, label_x, label_y,
            labels,
            avoid_label_lines_overlap=True,
            max_distance=0.15,
            linecolor="grey",
            textcolor="#b42920",
            x_scatter=[spectrum_x], y_scatter=[spectrum_y],
            scatter_sizes = [5] * len(spectrum_x)
        )

        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.SetIcon(wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO))

        plt.tight_layout()
        plt.show()


    def find_closest_index(self, value, array):
        closest_index = np.abs(array - value).argmin()

        return closest_index


    def find_increasing_indexes(self, arr, index):
        left_index = index - 2
        while left_index > 0 and arr[left_index - 1] <= arr[left_index]:
            left_index -= 1
        if index - left_index >= 3:
            right_index = index + (index - left_index)
        else:
            right_index = index + 2
            while right_index <= len(arr) - 1 and arr[right_index] >= arr[right_index + 1]:
                right_index += 1
            left_index = index - (right_index - left_index)

        return left_index, right_index
