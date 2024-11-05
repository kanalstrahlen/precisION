import os
import math
import pickle
import warnings
import numpy as np
import pandas as pd
import wx
import wx.grid
import matplotlib
from matplotlib import colors
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
import matplotlib.cm
from matplotlib.ticker import FuncFormatter

from proteoformAnalysisFunctions import ProteoformAnalysisRun

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
warnings.simplefilter(action='ignore', category=FutureWarning)

class ProteoformAnalysisWindow(wx.Frame):
    def __init__(self, parent, title, file_path, directory_path):
        super().__init__(
            parent,
            title=title,
            size=(375, 360),
            style=wx.DEFAULT_FRAME_STYLE
        )
        self.file_path = file_path
        self.directory_path = directory_path
        self.basename = os.path.basename(file_path).replace(".txt", "")

        #screen width variable to allow for smaller screens
        screen_width, screen_height = wx.DisplaySize()
        self.SetSize((screen_width * 3 // 4, screen_height * 3 // 4))

        self.panel = wx.Panel(self)
        self.plot_panel = ProteoformAnalysisSpectrumPanel(self.panel)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # header_sizer: title and key information
        header_sizer = wx.BoxSizer(wx.VERTICAL)
        title = wx.StaticText(self.panel, label="Proteoform Analysis")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext = wx.StaticText(
            self.panel,
            label=(
                "Fit theoretical isotopic envelopes to fragment spectra without deconvolution. "
                "Useful for low intensity fragment assignment and targeted searching."
            )
        )

        header_sizer.Add(title, 0, wx.BOTTOM, 5)
        header_sizer.Add(subtext, 0, wx.BOTTOM, -5)

        # top_sizer: proteoform input, proteoform selector, data grid containing details
        top_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # proteoform_input_sizer: box for proteoform sequence
        protein_input_staticbox = wx.StaticBox(self.panel)
        proteoform_input_sizer = wx.StaticBoxSizer(protein_input_staticbox, wx.VERTICAL)

        proteoform_input_header = wx.StaticText(self.panel, label="Proteoform input")
        subheader_font = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)
        proteoform_input_header.SetFont(subheader_font)

        proteoform_name_sizer = wx.BoxSizer(wx.HORIZONTAL)
        proteoform_name_text = wx.StaticText(self.panel, label="Proteoform name:")
        self.proteoform_name = wx.TextCtrl(self.panel, size=(250, 20))
        self.proteoform_name.SetValue("Proteoform A")

        proteoform_name_sizer.Add(proteoform_name_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        proteoform_name_sizer.Add(self.proteoform_name, 0, wx.RIGHT| wx.ALIGN_CENTER_VERTICAL, 0)

        self.sequence_text = wx.TextCtrl(
            self.panel,
            style=wx.TE_MULTILINE | wx.TE_CHARWRAP | wx.VSCROLL,
            size=(350, 150)
        )
        self.sequence_text.SetValue('')

        proteoform_input_sizer.Add(proteoform_input_header, 0, wx.BOTTOM, 5)
        proteoform_input_sizer.Add(proteoform_name_sizer, 0, wx.BOTTOM, 5)
        proteoform_input_sizer.Add(self.sequence_text, 0, wx.BOTTOM, 0)

        # search_options_sizer: options for searching, and search button
        search_options_staticbox = wx.StaticBox(self.panel)
        search_options_sizer = wx.StaticBoxSizer(search_options_staticbox, wx.VERTICAL)

        search_options_header = wx.StaticText(self.panel, label="Search options")
        search_options_header.SetFont(subheader_font)

        options_text_input_sizer = wx.BoxSizer(wx.HORIZONTAL)

        options_text_sizer = wx.BoxSizer(wx.VERTICAL)
        signal_noise_text = wx.StaticText(self.panel, label="Min. signal to noise")
        matching_tol_text = wx.StaticText(self.panel, label="m/z tol. (ppm)")
        isotope_tol_text = wx.StaticText(self.panel, label="Isotopologue m/z tol. (ppm)")
        min_score_text = wx.StaticText(self.panel, label="Min. fitting score")
        max_charge_text = wx.StaticText(self.panel, label="Max. charge state")
        spectrum_offset_text = wx.StaticText(self.panel, label="Error correction (ppm)")
        ion_type_text = wx.StaticText(self.panel, label="Ion type")

        options_text_sizer.Add(signal_noise_text, 0, wx.BOTTOM, 9)
        options_text_sizer.Add(matching_tol_text, 0, wx.BOTTOM, 9)
        options_text_sizer.Add(isotope_tol_text, 0, wx.BOTTOM, 9)
        options_text_sizer.Add(min_score_text, 0, wx.BOTTOM, 9)
        options_text_sizer.Add(max_charge_text, 0, wx.BOTTOM, 9)
        options_text_sizer.Add(spectrum_offset_text, 0, wx.BOTTOM, 9)
        options_text_sizer.Add(ion_type_text, 0, wx.BOTTOM, 9)

        options_input_sizer = wx.BoxSizer(wx.VERTICAL)
        self.signal_noise = wx.TextCtrl(self.panel, size=(40, 20))
        self.matching_tol = wx.TextCtrl(self.panel, size=(40, 20))
        self.isotope_tol = wx.TextCtrl(self.panel, size=(40, 20))
        self.min_score = wx.TextCtrl(self.panel, size=(40, 20))
        self.max_charge = wx.TextCtrl(self.panel, size=(40, 20))
        self.spectrum_offset = wx.TextCtrl(self.panel, size=(40, 20))
        ion_types = ["b/y", "c/z•"]
        self.ion_type = wx.ComboBox(self.panel, choices=ion_types, style=wx.CB_READONLY)

        self.signal_noise.SetValue("3")
        self.matching_tol.SetValue("5")
        self.isotope_tol.SetValue("3")
        self.min_score.SetValue("0.5")
        self.max_charge.SetValue("10")
        self.spectrum_offset.SetValue("0.0")
        self.ion_type.SetValue("b/y")

        options_input_sizer.Add(self.signal_noise, 0, wx.BOTTOM, 5)
        options_input_sizer.Add(self.matching_tol, 0, wx.BOTTOM, 5)
        options_input_sizer.Add(self.isotope_tol, 0, wx.BOTTOM, 5)
        options_input_sizer.Add(self.min_score, 0, wx.BOTTOM, 5)
        options_input_sizer.Add(self.max_charge, 0, wx.BOTTOM, 5)
        options_input_sizer.Add(self.spectrum_offset, 0, wx.BOTTOM, 5)
        options_input_sizer.Add(self.ion_type, 0, wx.BOTTOM, 5)

        options_text_input_sizer.Add(options_text_sizer, 0, wx.RIGHT, 10)
        options_text_input_sizer.Add(options_input_sizer, 0, wx.RIGHT, 0)

        search_button = wx.Button(self.panel, label="Search spectrum",  size=(170, 30))
        search_button.Bind(wx.EVT_BUTTON, self.on_search_button)

        search_options_sizer.Add(search_options_header, 0, wx.BOTTOM, 5)
        search_options_sizer.Add(options_text_input_sizer, 0, wx.BOTTOM, 5)
        search_options_sizer.Add(search_button, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 0)

        # view_options_sizer: options for viewing spectra
        view_options_staticbox = wx.StaticBox(self.panel)
        view_options_sizer = wx.StaticBoxSizer(view_options_staticbox, wx.VERTICAL)

        view_options_header = wx.StaticText(self.panel, label="View options")
        view_options_header.SetFont(subheader_font)

        view_text = wx.StaticText(self.panel, label="Spectrum to view")
        self.view_choices = ["", "View all"]
        self.view_dropdown = wx.ComboBox(
            self.panel,
            choices=self.view_choices,
            style=wx.CB_READONLY
        )
        self.view_dropdown.SetValue("")
        self.Bind(wx.EVT_COMBOBOX, self.on_view_button, self.view_dropdown)

        color_text = wx.StaticText(self.panel, label="Color")
        self.color_dropdown = wx.ComboBox(
            self.panel,
            choices=list(colors.CSS4_COLORS),
            style=wx.CB_READONLY
        )
        self.color_dropdown.SetValue("dodgerblue")

        map_button = wx.Button(
            self.panel,
            label="View Sequence Map",
            size=(170, 30)
        )
        map_button.Bind(wx.EVT_BUTTON, self.on_map_button)

        save_button = wx.Button(
            self.panel,
            label="Export Output",
            size=(170, 30)
        )
        save_button.Bind(wx.EVT_BUTTON, self.on_save_button)

        view_options_sizer.Add(view_options_header, 0, wx.BOTTOM, 5)
        view_options_sizer.Add(view_text, 0, wx.BOTTOM, 3)
        view_options_sizer.Add(self.view_dropdown, 0, wx.BOTTOM, 5)
        view_options_sizer.Add(color_text, 0, wx.BOTTOM, 3)
        view_options_sizer.Add(self.color_dropdown, 0, wx.BOTTOM, 10)
        view_options_sizer.Add(map_button, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 5)
        view_options_sizer.Add(save_button, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 0)

        self.data_grid = wx.grid.Grid(self.panel)
        self.data_grid.SetDefaultCellBackgroundColour("#f0f0f0")
        self.data_grid.Bind(wx.grid.EVT_GRID_SELECT_CELL, self.on_cell_select)

        self.num_rows = 1
        self.num_cols = 8
        self.data_grid.CreateGrid(self.num_rows, self.num_cols)
        column_names = [
            "Proteoform",
            "Ion",
            "Mass (Da)", 
            "Charge",
            "m/z",
            "Intensity",
            "Error (ppm)",
            "Score"
            ]
        for col, column_name in enumerate(column_names):
            self.data_grid.SetColLabelValue(col, column_name)

        self.data_grid.AutoSizeColumns()
        self.data_grid.Refresh()

        top_sizer.Add(proteoform_input_sizer, 0, wx.RIGHT, 5)
        top_sizer.Add(search_options_sizer, 0, wx.RIGHT, 5)
        top_sizer.Add(view_options_sizer, 0, wx.RIGHT, 5)
        top_sizer.Add(self.data_grid, 1, wx.ALL, 5)

        top_sizer.SetItemMinSize(self.data_grid, -1, 220)

        main_sizer.Add(header_sizer, 0, wx.ALL, 5)
        main_sizer.Add(top_sizer, 0, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(self.plot_panel, 1, wx.EXPAND | wx.ALL, 5)

        icon = wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Bind(wx.EVT_CLOSE, self.on_close)
        try:
            self.load_empty_spectrum(file_path, directory_path)
            no_spec = False
        except:
            print("The spectrum must be analysed using the deconvolution module first.")
            no_spec = True
            return

        self.matched_proteoforms = []
        self.matched_envelopes = []
        self.noise = []
        self.cutoff = []
        self.sequences = []
        self.offsets = []

        # load session from pickle file if it exists
        pickle_path = os.path.join(self.directory_path, f"{self.basename}.proteoformAnalysis.pkl")
        if os.path.exists(pickle_path):
            with open(pickle_path, "rb") as f:
                state_dict = pickle.load(f)
                self.matched_proteoforms = state_dict.get("matched_proteoforms", [])
                self.matched_envelopes = state_dict.get("matched_envelopes", [])
                self.noise = state_dict.get("noise", [])
                self.cutoff = state_dict.get("cutoff", [])
                self.sequences = state_dict.get("sequences", [])
                self.offsets = state_dict.get("offsets", [])
                self.view_choices = state_dict.get("view_choices", ["", "View all"])
                self.view_dropdown.SetItems(self.view_choices)
                self.view_dropdown.SetValue("")
                self.view_dropdown.Refresh()
                self.panel.Layout()

        self.Show()
        if no_spec:
            self.Close()


    def load_empty_spectrum(self, file_path, directory_path):
        save_path = os.path.join(directory_path, f"{self.basename}.spectrum.txt")
        if not os.path.exists(save_path):
            rl_file_path = os.path.join(directory_path, f"{self.basename}.profile.rl.txt")
            if os.path.exists(rl_file_path):
                print('Loading spectrum...')
                output = np.loadtxt(rl_file_path)
                output = output[::2]
                np.savetxt(save_path, output, delimiter=" ")
                print(
                    f"Created {save_path} using RL spectrum. "
                    f"Please delete {self.basename}.profile.rl.txt and {self.basename}.spectrum.txt"
                    " and open the module again if you want to use the raw spectrum."
                )
            else:
                output = np.loadtxt(file_path)
                save_path = os.path.join(directory_path, f"{self.basename}.spectrum.txt")
                np.savetxt(save_path, output, delimiter=" ")

        self.spectrum = np.loadtxt(save_path)
        self.min_mz = min(self.spectrum[:, 0])
        self.max_mz = max(self.spectrum[:, 0])
        self.abu = max(self.spectrum[:,1])
        self.plot_panel.plot_spectrum(
            self.spectrum,
            self.min_mz,
            self.max_mz,
            self.abu
        )


    def on_map_button(self, event):
        selected_index = self.view_dropdown.GetSelection()
        if selected_index >= 2: # specific proteoform
            df = self.matched_proteoforms[selected_index-2]
            sequence = self.sequences[selected_index-2]
            name = self.view_dropdown.GetValue()
            ions = df["name"].unique()
            run = ProteoformAnalysisRun()
            sequence, _ = run.process_sequence_string(sequence)
            self.gen_sequence_map(name, sequence, ions)

        else:
            print("Please select a specific proteoform.")
            return


    def gen_sequence_map(self, name, sequence, ions):
        b_fragment_positions = []
        y_fragment_positions = []

        for ion in ions:
            ion_serie = ion[0]
            ion_index = int(ion[1:])
            if ion_serie in ("b", "c"):
                b_fragment_positions.append(ion_index)
            elif ion_serie in ("y", "z"):
                y_fragment_positions.append(len(sequence) - ion_index - 1)
        b_fragment_positions = list(set(b_fragment_positions))
        y_fragment_positions = list(set(y_fragment_positions))

        if len(sequence) >= 400:
            residues_per_row = 40
        else:
            residues_per_row = 20

        self.plot_sequence_map(
            name,
            sequence,
            b_fragment_positions,
            y_fragment_positions,
            residues_per_row
        )


    def plot_sequence_map(
        self, name, protein_sequence,
        b_fragment_positions, y_fragment_positions, residues_per_row=20
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

            c_term_colour = "grey"

            if i == 0:
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


    def on_save_button(self, event):
        save_path = os.path.join(
            self.directory_path,
            f"{self.basename}.proteoformAnalysis.output.csv"
        )

        df1 = pd.DataFrame(columns=[
            'proteoform',
            'name',
            'mass',
            'charge',
            'mz',
            'intensity',
            'error',
            'score',
            'min_envelope',
            'max_envelope',
            'max_abu'
        ])

        for i in range(2, len(self.view_choices)):
            df2 = self.matched_proteoforms[i-2]
            df1 = pd.concat([df1, df2], ignore_index=True)

        df1 = df1.drop(columns=['min_envelope', 'max_envelope', 'max_abu'])
        df1.to_csv(save_path, index=False)

        print(f'Saved combined results to {save_path}')


    def on_search_button(self, event):
        self.Freeze()
        name = self.proteoform_name.GetValue()
        sequence = self.sequence_text.GetValue()
        matching_tol = float(self.matching_tol.GetValue())
        isotope_tol = float(self.isotope_tol.GetValue())
        min_score = float(self.min_score.GetValue())
        max_charge = int(self.max_charge.GetValue())
        spectrum_offset = float(self.spectrum_offset.GetValue())
        signal_noise = float(self.signal_noise.GetValue())
        ion_type = self.ion_type.GetValue()

        run = ProteoformAnalysisRun()

        fragments = run.generate_theoretical_fragments(
            sequence,
            max_charge,
            self.min_mz,
            self.max_mz,
            ion_type
        )

        centroid_file_path = os.path.join(self.directory_path, f"{self.basename}_centroid.txt")

        spectrum = run.load_centroid_spectrum(centroid_file_path, spectrum_offset)
        spectrum_x, spectrum_y, noise, cutoff = run.remove_noise_peaks(spectrum, signal_noise)

        matched_ions = run.match_ion_mz(spectrum_x, fragments, matching_tol)

        fit_ions, fit_ion_envelopes, fit_ion_properties = run.fit_ions(
            spectrum[:,0],
            spectrum[:,1],
            matched_ions,
            matching_tol,
            isotope_tol,
            min_score
        )

        data = {
            'proteoform': [name] * len(fit_ions),
            'name': [ion[0] for ion in fit_ions],
            'mass': [round(ion[2], 2) for ion in fit_ions],
            'charge': [ion[3] for ion in fit_ions],
            'mz': [round(ion[4], 4) for ion in fit_ions],
            'intensity': [round(properties[2], 2) for properties in fit_ion_properties],
            'error': [round(properties[0], 2) for properties in fit_ion_properties],
            'score': [round(properties[1], 2) for properties in fit_ion_properties],
            'min_envelope': [round(properties[3], 2) for properties in fit_ion_properties],
            'max_envelope': [round(properties[4], 2) for properties in fit_ion_properties],
            'max_abu': [round(properties[5], 2) for properties in fit_ion_properties]
        }

        df = pd.DataFrame(data)

        if name in self.view_choices:
            index = self.view_choices.index(name)
            self.view_choices[index] = name
            self.matched_proteoforms[index-2] = df
            self.matched_envelopes[index-2] = fit_ion_envelopes
            self.noise[index-2] = noise
            self.cutoff[index-2] = cutoff
            self.sequences[index-2] = sequence
            self.offsets[index-2] = spectrum_offset
        else:
            self.view_choices.append(name)
            self.matched_proteoforms.append(df)
            self.matched_envelopes.append(fit_ion_envelopes)
            self.noise.append(noise)
            self.cutoff.append(cutoff)
            self.sequences.append(sequence)
            self.offsets.append(spectrum_offset)

        self.view_dropdown.SetItems(self.view_choices)
        self.view_dropdown.SetValue(name)
        self.view_dropdown.Refresh()
        self.panel.Layout()

        self.on_view_button(None)

        self.Thaw()


    def on_view_button(self, event):
        self.Freeze()
        selected_index = self.view_dropdown.GetSelection()

        if selected_index >= 2: # specific proteoform
            color = self.color_dropdown.GetValue()

            self.df = self.matched_proteoforms[selected_index-2]
            fit_ion_envelopes = self.matched_envelopes[selected_index-2]
            noise = self.noise[selected_index-2]
            cutoff = self.cutoff[selected_index-2]
            sequence = self.sequences[selected_index-2]
            spectrum_offset = self.offsets[selected_index-2]
            name = self.view_dropdown.GetValue()
            self.proteoform_name.SetValue(name)
            self.sequence_text.SetValue(sequence)

            self.plot_panel.plot_spectrum(
                self.spectrum,
                self.min_mz,
                self.max_mz,
                self.abu,
                spectrum_offset
            )

            self.plot_panel.plot_matches(fit_ion_envelopes, color)
            self.plot_panel.plot_noise(self.min_mz, self.max_mz, noise, cutoff)

            self.update_grid(self.df)


        elif selected_index == 0: #None
            self.plot_panel.plot_spectrum(
                self.spectrum,
                self.min_mz,
                self.max_mz,
                self.abu
            )

            self.clear_grid()

        elif selected_index == 1: #View all
            mean_offset = sum(self.offsets) / len(self.offsets)

            self.plot_panel.plot_spectrum(
                self.spectrum,
                self.min_mz,
                self.max_mz,
                self.abu,
                mean_offset
            )

            tab_colors = plt.cm.get_cmap('tab10')

            self.df = pd.DataFrame(columns=[
                'proteoform',
                'name',
                'mass',
                'charge',
                'mz',
                'intensity',
                'error',
                'score',
                'min_envelope',
                'max_envelope',
                'max_abu'
            ])

            for i in range(2, len(self.view_choices)):
                fit_ion_envelopes = self.matched_envelopes[i-2]
                df = self.matched_proteoforms[i-2]
                self.df = pd.concat([self.df, df], ignore_index=True)
                self.plot_panel.plot_matches(fit_ion_envelopes, tab_colors(i-2))

            self.update_grid(self.df)
            self.plot_panel.change_zoom(self.min_mz, self.max_mz, self.abu)

        self.Thaw()


    def update_grid(self, df):
        self.data_grid.ClearGrid()

        num_rows, num_cols = df.shape
        diff = self.num_rows - num_rows

        if diff > 0:
            self.data_grid.DeleteRows(self.num_rows - diff, diff)
            self.num_rows = num_rows
        elif diff < 0:
            self.data_grid.AppendRows(-1 * diff)
            self.num_rows = num_rows

        for row in range(self.num_rows):
            for col in range(8):
                value = str(df.iloc[row, col])
                self.data_grid.SetCellValue(row, col, value)

        self.data_grid.AutoSizeColumns()
        self.data_grid.Refresh()
        self.data_grid.Show()

        self.panel.Layout()


    def clear_grid(self):
        self.data_grid.ClearGrid()

        num_rows = 0
        diff = self.num_rows - num_rows

        if diff > 0:
            self.data_grid.DeleteRows(self.num_rows - diff, diff)
            self.num_rows = num_rows
        elif diff < 0:
            self.data_grid.AppendRows(-1 * diff)
            self.num_rows = num_rows

        self.data_grid.AutoSizeColumns()
        self.data_grid.Refresh()
        self.data_grid.Show()

        self.panel.Layout()


    def on_cell_select(self, event):
        self.selected_row = event.GetRow()
        charge = self.df.iloc[self.selected_row, self.df.columns.get_loc('charge')]
        border = 2 / charge
        min_mz = self.df.iloc[self.selected_row, self.df.columns.get_loc('min_envelope')] - border
        max_mz = self.df.iloc[self.selected_row, self.df.columns.get_loc('max_envelope')] + border
        max_abu = self.df.iloc[self.selected_row, self.df.columns.get_loc('max_abu')]
        self.plot_panel.change_zoom(min_mz, max_mz, max_abu)


    def on_size(self, event):
        size = self.panel.GetClientSize()
        self.panel.Layout()
        event.Skip()


    def on_close(self, event):
        self.save_state()
        event.Skip()


    def save_state(self):
        state_dict = {
            "matched_proteoforms": self.matched_proteoforms,
            "matched_envelopes": self.matched_envelopes,
            "noise": self.noise,
            "cutoff": self.cutoff,
            "sequences": self.sequences,
            "offsets": self.offsets,
            "view_choices": self.view_choices
            }

        output_path = os.path.join(self.directory_path, f"{self.basename}.proteoformAnalysis.pkl")
        with open(output_path, "wb") as f:
            pickle.dump(state_dict, f)

        print("Saved current progress.")



class ProteoformAnalysisSpectrumPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)

        sizer = wx.BoxSizer(wx.VERTICAL)

        self.figure = Figure(facecolor="#f0f0f0")
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.max_intens = 0

        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.EXPAND)

        self.SetSizer(sizer)
        self.Bind(wx.EVT_SIZE, self.on_size)
        self.Layout()
        self.toolbar.Realize()


    def plot_spectrum(
        self, spectrum, selected_min_mz, selected_max_mz, abu, spectrum_offset=0):
        self.ax.cla()

        x = [mz + mz * spectrum_offset / 1000000 for mz in spectrum[:, 0]]
        y = spectrum[:, 1]
        self.ax.plot(x, y, color="black", linewidth=1)

        self.ax.set_xlabel("m/z")
        self.ax.set_ylabel("Intensity")
        self.ax.set_facecolor("#f0f0f0")
        self.ax.set_xlim(selected_min_mz, selected_max_mz)
        self.ax.xaxis.set_major_formatter(FuncFormatter(self.format_with_commas))
        self.max_intens = abu * 1.1
        self.ax.set_ylim(0, self.max_intens)

        self.canvas.draw()


    def plot_points(
        self, spectrum_x, spectrum_y):
        self.ax.scatter(spectrum_x, spectrum_y)
        self.canvas.draw()


    def plot_noise(
        self, min_mz, max_mz, noise, cutoff):
        self.ax.plot([min_mz, max_mz], [noise, noise], color='red')
        self.ax.plot([min_mz, max_mz], [cutoff, cutoff], color='red', linestyle='dotted')
        self.canvas.draw()


    def plot_matches(
        self, fit_ions, color):
        for tup in fit_ions:
            scatter_x = tup[0]
            scatter_y = tup[1]
            self.ax.scatter(scatter_x, scatter_y, color=color, marker="o", alpha=0.5)
        self.canvas.draw()


    def change_zoom(self, min_mz, max_mz, abu):
        self.ax.set_xlim(min_mz, max_mz)
        self.ax.set_ylim(0, abu * 1.2)
        self.canvas.draw()


    def on_size(self, event):
        self.fit_plot_to_panel()
        event.Skip()


    def fit_plot_to_panel(self):
        size = self.GetClientSize()
        if size[0] >= 21:
            dpi = self.GetContentScaleFactor() * 100
            width = size.width / dpi
            height = size.height / dpi - 0.3
            self.figure.set_size_inches(width, height)
            self.figure.tight_layout(rect=[0, 0, 1, 1])
            self.canvas.draw()


    def format_with_commas(self, value, pos):
        formatted_value = "{:,.0f}".format(value)
        return formatted_value
