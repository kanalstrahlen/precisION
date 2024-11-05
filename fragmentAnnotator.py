import os
import ast
import requests
import wx
import wx.grid
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
import matplotlib.cm

from fragmentAnnotatorFunctions import AnnotatorFunctions
from defineTermini import DefineTerminiWindow
from unmodSearch import NoModSearchWindow
from modSearch import ModSearchWindow
from modDiscovery import ModDiscoveryWindow
from filterAssignments import FilterAssignmentsWindow

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


class FragmentAnnotatorWindow(wx.Frame):
    def __init__(self, parent, title, file_path, directory_path):
        super().__init__(
            parent,
            title=title,
            size=(375, 515),
            style=wx.DEFAULT_FRAME_STYLE & ~(wx.RESIZE_BORDER | wx.MAXIMIZE_BOX)
        )

        self.file_path = file_path
        self.directory_path = directory_path

        self.panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # header_sizer: title and key information
        header_sizer = wx.BoxSizer(wx.VERTICAL)
        title = wx.StaticText(self.panel, label="Fragment Assignment")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext1 = wx.StaticText(
            self.panel,
            label="Assign fragments using a hierachical approach that maximises"
        )
        subtext2 = wx.StaticText(
            self.panel,
            label="assignment likelihood."
        )

        header_sizer.Add(title, 0, wx.BOTTOM, 5)
        header_sizer.Add(subtext1, 0, wx.BOTTOM, 3)
        header_sizer.Add(subtext2, 0, wx.BOTTOM, -5)

        # protein_id_sizer: load protein from uniprot and vary sequence manually
        protein_id_staticbox = wx.StaticBox(self.panel)
        protein_id_sizer = wx.StaticBoxSizer(protein_id_staticbox, wx.VERTICAL)
        protein_id_sizer_header = wx.StaticText(self.panel, label="Protein ID and Sequence")
        protein_id_sizer_header.SetFont(
            wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)
        )

        # protein_id_subsizer1: load file using uniprot accession
        protein_id_subsizer1 = wx.BoxSizer(wx.HORIZONTAL)
        protein_id_label = wx.StaticText(self.panel, label='Protein ID:')
        self.id_text = wx.TextCtrl(self.panel, size=(220, 20))
        load_uniprot_button = wx.Button(self.panel, label="Load",  size=(50, 25))
        load_uniprot_button.Bind(wx.EVT_BUTTON, self.on_load_uniprot_button)

        protein_id_subsizer1.Add(protein_id_label, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        protein_id_subsizer1.Add(self.id_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        protein_id_subsizer1.Add(load_uniprot_button, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)

        # sequence_text: current sequence
        self.sequence_text = wx.TextCtrl(
            self.panel,
            style=wx.TE_MULTILINE | wx.VSCROLL,
            size=(350, 150)
        )

        # protein_id_subsizer2: define mods and update sequence plot
        protein_id_subsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        n_mod_label = wx.StaticText(self.panel, label='N-Mod:')
        self.n_mod_add = wx.TextCtrl(self.panel, size=(40, 20))
        add_label = wx.StaticText(self.panel, label='(+)')
        self.n_mod_sub = wx.TextCtrl(self.panel, size=(40, 20))
        sub_label = wx.StaticText(self.panel, label='(-)')
        c_mod_label = wx.StaticText(self.panel, label='C-Mod:')
        self.c_mod_add = wx.TextCtrl(self.panel, size=(40, 20))
        add_label2 = wx.StaticText(self.panel, label='(+)')
        self.c_mod_sub = wx.TextCtrl(self.panel, size=(40, 20))
        sub_label2 = wx.StaticText(self.panel, label='(-)')

        protein_id_subsizer2.Add(n_mod_label, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        protein_id_subsizer2.Add(self.n_mod_add, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 2)
        protein_id_subsizer2.Add(add_label, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 2)
        protein_id_subsizer2.Add(self.n_mod_sub, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 2)
        protein_id_subsizer2.Add(sub_label, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)
        protein_id_subsizer2.Add(c_mod_label, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        protein_id_subsizer2.Add(self.c_mod_add, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 2)
        protein_id_subsizer2.Add(add_label2, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 2)
        protein_id_subsizer2.Add(self.c_mod_sub, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 2)
        protein_id_subsizer2.Add(sub_label2, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)

        protein_id_sizer.Add(protein_id_sizer_header, 0, wx.BOTTOM, 5)
        protein_id_sizer.Add(protein_id_subsizer1, 0, wx.BOTTOM, 5)
        protein_id_sizer.Add(self.sequence_text, 0, wx.BOTTOM, 5)
        protein_id_sizer.Add(protein_id_subsizer2, 0, wx.BOTTOM, )

        # button_sizer: all buttons in 2xn grid
        button_sizer = wx.BoxSizer(wx.VERTICAL)

        button_subsizer1 = wx.BoxSizer(wx.HORIZONTAL)
        show_spectrum_button = wx.Button(self.panel, label="Show Spectrum",  size=(345, 30))
        show_spectrum_button.Bind(wx.EVT_BUTTON, self.on_show_spectrum_button)

        button_subsizer1.Add(show_spectrum_button, 0, wx.RIGHT, 0)

        line = wx.StaticLine(self.panel, wx.ID_ANY, style=wx.LI_HORIZONTAL)

        button_subsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        define_termini_button = wx.Button(self.panel, label="Define Termini",  size=(170, 30))
        define_termini_button.Bind(wx.EVT_BUTTON, self.on_define_termini_button)
        no_mod_search_button = wx.Button(self.panel, label="Unmodified Ion Search",  size=(170, 30))
        no_mod_search_button.Bind(wx.EVT_BUTTON, self.on_no_mod_search_button)

        button_subsizer2.Add(define_termini_button, 0, wx.RIGHT, 5)
        button_subsizer2.Add(no_mod_search_button, 0, wx.RIGHT, 0)

        button_subsizer3 = wx.BoxSizer(wx.HORIZONTAL)
        variable_search_button = wx.Button(
            self.panel,
            label="Modified Ion Search",
            size=(170, 30)
        )
        variable_search_button.Bind(wx.EVT_BUTTON, self.on_variable_search_button)
        mod_discovery_button = wx.Button(
            self.panel,
            label="Modification Discovery",
            size=(170, 30)
        )
        mod_discovery_button.Bind(wx.EVT_BUTTON, self.on_mod_discovery_button)

        button_subsizer3.Add(mod_discovery_button, 0, wx.RIGHT, 5)
        button_subsizer3.Add(variable_search_button, 0, wx.RIGHT, 0)

        button_subsizer4 = wx.BoxSizer(wx.HORIZONTAL)
        view_plots_button = wx.Button(self.panel, label="View Validation Plots", size=(170,30))
        view_plots_button.Bind(wx.EVT_BUTTON, self.on_view_plots_button)
        filter_assignments_button = wx.Button(self.panel, label="Filter Assignments", size=(170,30))
        filter_assignments_button.Bind(wx.EVT_BUTTON, self.on_filter_assignments_button)

        button_subsizer4.Add(view_plots_button, 0, wx.RIGHT, 5)
        button_subsizer4.Add(filter_assignments_button, 0, wx.RIGHT, 0)

        button_sizer.Add(button_subsizer1, 0, wx.BOTTOM, 5)
        button_sizer.Add(line, 0, wx.EXPAND | wx.BOTTOM, 5)
        button_sizer.Add(button_subsizer2, 0, wx.BOTTOM, 5)
        button_sizer.Add(button_subsizer3, 0, wx.BOTTOM, 5)
        button_sizer.Add(button_subsizer4, 0, wx.BOTTOM, 0)

        main_sizer.Add(header_sizer, 0, wx.ALL, 5)
        main_sizer.Add(protein_id_sizer, 0, wx.ALL,  5)
        main_sizer.Add(button_sizer, 0, wx.BOTTOM |  wx.ALIGN_CENTER_HORIZONTAL, 0)

        icon = wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        # set font size for sequence plot based on size of screen
        screen_width, _ = wx.DisplaySize()
        x_length = screen_width // 2
        self.sequence_font_size = x_length // 100

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Show()


    def on_show_spectrum_button(self, event):
        self.plot_window = FragmentAnnotatorSpectrumWindow(
            None,
            f"precisION - Fragment Matching [{os.path.basename(self.file_path)}]",
            self.file_path,
            self.directory_path
        )
        self.plot_window.Show()


    def on_load_uniprot_button(self, event):
        accession = self.id_text.GetValue()
        save_path = os.path.join(self.directory_path, f"{accession}.xml")

        url = "https://rest.uniprot.org/uniprotkb/stream"
        params = {
            "format": "xml",
            "query": f"accession:{accession}"
        }
        response = requests.get(url, params=params, timeout=10)

        if response.status_code == 200:
            xml_data = response.text
            with open(save_path, "w") as xml_file:
                xml_file.write(xml_data)
            print(f"Data downloaded and saved as {save_path}.")
        else:
            print(f"Failed to retrieve data. Status code: {response.status_code}")
            return

        sequence = AnnotatorFunctions.extract_protein_sequence_from_xml(save_path)
        self.sequence_text.SetValue(sequence)


    def on_define_termini_button(self, event):
        seq = self.sequence_text.GetValue()
        if len(seq) <= 1:
            print("No sequence input.")
            return

        define_termini_window = DefineTerminiWindow(
            self,
            f"precisION - Define Termini [{os.path.basename(self.file_path)}]",
            self.sequence_text.GetValue(),
            self.file_path,
            self.directory_path
        )
        define_termini_window.Show()


    def on_no_mod_search_button(self, event):
        seq = self.sequence_text.GetValue()
        if len(seq) == 0:
            print("No sequence input.")
            return
        name = self.id_text.GetValue()
        if len(name) == 0:
            print("No protein ID provided.")
            return    
        no_mod_search_window = NoModSearchWindow(
            self,
            f"precisION - Unmodified Ion Search [{os.path.basename(self.file_path)}]",
            self.id_text.GetValue(),
            self.sequence_text.GetValue(),
            self.file_path,
            self.directory_path,
            self.n_mod_add.GetValue(),
            self.n_mod_sub.GetValue(),
            self.c_mod_add.GetValue(),
            self.c_mod_sub.GetValue()
        )
        no_mod_search_window.Show()


    def on_mod_discovery_button(self, event):
        seq = self.sequence_text.GetValue()
        if len(seq) <= 1:
            print("No sequence input.")
            return

        assigned_peaks_file = (
            f"{os.path.splitext(os.path.basename(self.file_path))[0]}"
            ".assignedPeakList.csv"
        )
        assigned_peaks_path = os.path.join(self.directory_path, assigned_peaks_file)

        sequence = self.gen_sequence_with_modifications(
            self.sequence_text.GetValue(),
            self.n_mod_add.GetValue(),
            self.n_mod_sub.GetValue(),
            self.c_mod_add.GetValue(),
            self.c_mod_sub.GetValue()
        )

        if os.path.exists(assigned_peaks_path):
            mod_discovery_window = ModDiscoveryWindow(
                self,
                f"precisION - Modification Discovery [{os.path.basename(self.file_path)}]",
                self.file_path,
                self.directory_path,
                sequence
            )
            mod_discovery_window.Show()
        else:
            print(
                "The annotated peak list is not present in the working directory. "
                "Please complete the unmodified ion search."
            )


    def on_variable_search_button(self, event):
        seq = self.sequence_text.GetValue()
        if len(seq) == 0:
            print("No sequence input.")
            return
        name = self.id_text.GetValue()
        if len(name) == 0:
            print("No protein ID provided.")
            return    
        mod_search_window = ModSearchWindow(
            self,
            f"precisION - Modified Ion Search [{os.path.basename(self.file_path)}]",
            self.id_text.GetValue(),
            self.sequence_text.GetValue(),
            self.file_path,
            self.directory_path,
            self.n_mod_add.GetValue(),
            self.n_mod_sub.GetValue(),
            self.c_mod_add.GetValue(),
            self.c_mod_sub.GetValue()
        )
        mod_search_window.Show()        


    def on_view_plots_button(self, event):
        file_name = (
            f"{os.path.splitext(os.path.basename(self.file_path))[0]}"
            ".assignedPeakList.csv"
        )
        annotated_peaks_path = os.path.join(self.directory_path, file_name)
        cal_file_name = (
            f"{os.path.splitext(os.path.basename(self.file_path))[0]}"
            ".calibrationInformation.txt"
        )
        cal_file_path = os.path.join(self.directory_path, cal_file_name)
        if os.path.exists(annotated_peaks_path) and os.path.exists(cal_file_path):
            AnnotatorFunctions.validation_plots(
                annotated_peaks_path,
                cal_file_path
            )
        else:
            print(
                "The annotated peak list or calibration file is not present in the working directory. "
                "Please complete the unmodified ion search."
            )


    def on_filter_assignments_button(self, event):
        file_name = (
            f"{os.path.splitext(os.path.basename(self.file_path))[0]}"
            ".assignedPeakList.csv"
        )
        annotated_peaks_path = os.path.join(self.directory_path, file_name)
        if os.path.exists(annotated_peaks_path):
            FilterAssignmentsWindow(
                self,
                f"precisION - Filter Assignments [{os.path.basename(self.file_path)}]",
                annotated_peaks_path
            )
        else:
            print(
                "The annotated peak list is not present in the working directory. "
                "Please complete the unmodified ion search."
            )


    def on_size(self, event):
        self.panel.Layout()
        event.Skip()


    def gen_sequence_with_modifications(self, sequence, n_add, n_sub, c_add, c_sub):
        seq = sequence[0]
        if len(n_add)>=1 or len(n_sub)>=1:
            seq += "["
            if len(n_add)>=1:
                seq += n_add
            if len(n_sub)>=1:
                seq += f"-{n_sub}"
            seq += "]"
        seq += sequence[1:]
        if len(c_add)>=1 or len(c_sub)>=1:
            seq += "["
            if len(c_add)>=1:
                seq += c_add
            if len(c_sub)>=1:
                seq += f"-{c_sub}"
            seq += "]"

        seq = seq.replace("c", "C[-H]")

        return seq



class FragmentAnnotatorSpectrumWindow(wx.Frame):
    def __init__(self, parent, title, file_path, directory_path):
        super().__init__(parent, title=title)
        self.file_path = file_path
        self.directory_path = directory_path
        self.basename = os.path.basename(file_path).replace('.txt', '')

        # screen width variable to allow for smaller screens
        screen_width, screen_height = wx.DisplaySize()
        self.SetSize((screen_width * 3 // 4, screen_height * 2 // 4))

        self.panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.data_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.plot_panel = AnnotatorSpectrumPanel(self.panel)

        self.data_grid = wx.grid.Grid(self.panel)
        self.data_grid.SetDefaultCellBackgroundColour('#f0f0f0')
        self.data_grid.Bind(wx.grid.EVT_GRID_SELECT_CELL, self.on_cell_select)

        self.data_sizer.Add(self.data_grid, 1, wx.EXPAND)
        #self.data_sizer.SetMinSize(wx.Size(720, -1))

        main_sizer.Add(self.plot_panel, 0, wx.EXPAND | wx.RIGHT, border=-10)
        main_sizer.Add(self.data_sizer, 1, wx.EXPAND | wx.ALL, border=15)

        icon = wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.load_decon()
        self.load_spectrum(directory_path)

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Show()


    def on_size(self, event):
        self.panel.Layout()
        event.Skip()


    def load_spectrum(self, directory_path):
        save_path = os.path.join(directory_path, f"{self.basename}.recalibratedSpectrum.txt")
        self.spectrum = np.loadtxt(save_path)
        self.plot_panel.plot_spectrum(self.spectrum, self.peak_list, None)


    def load_decon(self):
        peaks_file = os.path.join(self.directory_path, f"{self.basename}.assignedPeakList.csv")

        self.peak_list = pd.read_csv(peaks_file)
        self.peak_list = self.peak_list.sort_values(by="monoisotopic_mz")
        column_filter = [
            'monoisotopic_mz',
            'charge',
            'abundance',
            'name',
            'ion',
            'variable_mod',
            'adduct',
            'loss',
            'ppm_error',
            'fit_score'
        ]
        peak_list_filtered = self.peak_list[column_filter].copy()

        # make color ids to help colour different proteins different colours
        unique_names = self.peak_list['name'].unique()
        value_to_index = {value: index for index, value in enumerate(unique_names)}
        self.peak_list['color_id'] = self.peak_list['name'].map(value_to_index)

        self.Freeze()
        self.data_grid.Destroy()
        self.data_grid = wx.grid.Grid(self.panel)
        self.data_grid.Hide()
        self.data_sizer.Add(self.data_grid, 1)
        self.data_grid.SetDefaultCellBackgroundColour('#f0f0f0')
        self.data_grid.Bind(wx.grid.EVT_GRID_SELECT_CELL, self.on_cell_select)

        num_rows, num_cols = peak_list_filtered.shape

        self.data_grid.CreateGrid(num_rows, num_cols)
        for row in range(num_rows):
            for col in range(num_cols):
                input_value = peak_list_filtered.iloc[row, col]
                try:
                    if np.isnan(input_value):
                        input_value = " "
                    elif isinstance(input_value, np.float64):
                        input_value = np.round(input_value, 2)
                except:
                    pass

                value = str(input_value)
                self.data_grid.SetCellValue(row, col, value)

        column_names = [
            'Monoisotopic m/z',
            'Charge', 
            'Abundance',
            'Protein',
            'Ion',
            'Modification',
            'Adduct',
            'Loss',
            'Error (ppm)',
            'Score'
            ]
        for col, column_name in enumerate(column_names):
            self.data_grid.SetColLabelValue(col, column_name)

        self.data_grid.AutoSizeColumns()
        self.data_grid.Refresh()
        self.data_grid.Show()
        self.panel.Layout()
        self.Thaw()


    def on_cell_select(self, event):
        selected_row = event.GetRow()
        self.plot_panel.plot_spectrum(self.spectrum, self.peak_list, selected_row)
        event.Skip()



class AnnotatorSpectrumPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        sizer = wx.BoxSizer(wx.VERTICAL)

        self.figure = Figure(facecolor='#f0f0f0')
        self.ax = self.figure.add_subplot(111)

        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.toolbar = NavigationToolbar2Wx(self.canvas)

        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.EXPAND)

        self.SetSizer(sizer)
        self.Layout()
        self.toolbar.Realize()
        self.Bind(wx.EVT_SIZE, self.on_size)

        self.max_intens = 0


    def plot_spectrum(self, spectrum, peak_list, selected_index):
        self.ax.cla()

        self.ax.plot(spectrum[:, 0], spectrum[:, 1], color='black', linewidth=1)

        self.ax.set_xlabel('m/z')
        self.ax.set_ylabel('Intensity')
        self.ax.set_facecolor('#f0f0f0')
        self.ax.set_xlim(min(spectrum[:, 0]), max(spectrum[:, 0]))
        self.max_intens = max(spectrum[:, 1]) * 1.25
        self.ax.set_ylim(0, self.max_intens)

        theo_x = [] # list of theo envelope mz values
        theo_y = [] # list of theo envelope intens values
        assignment = [] # list of assignment names
        color_id = [] # list of color id to color when multiple species

        for i in range(len(peak_list)):
            x_scatter = peak_list['fitter_theo_x'].iloc[i]
            y_scatter = peak_list['fitter_theo_y'].iloc[i]
            if type(x_scatter) != str:
                x_scatter = peak_list['theo_x'].iloc[i]
                y_scatter = peak_list['theo_y'].iloc[i]
            theo_x.append(ast.literal_eval(x_scatter))
            theo_y.append(ast.literal_eval(y_scatter))
            if type(peak_list['ion'].iloc[i]) == str:
                assignment.append(peak_list['ion'].iloc[i])
                if isinstance(peak_list['color_id'].iloc[i], np.int64):
                    color_id.append(peak_list['color_id'].iloc[i])
                else:
                    color_id.append(None)
            else:
                assignment.append(None)

        # breaking up envelopes into assigned and not assigned
        # extract max point for plotting
        matched_scatter_x = []
        matched_scatter_y = []
        unmatched_scatter_x = []
        unmatched_scatter_y = []

        for x, y, ass in zip(theo_x, theo_y, assignment):
            max_y = max(y)
            arg = y.index(max_y)
            label_x = x[arg]
            label_y = y[arg]
            if ass is None:
                unmatched_scatter_x.append(label_x)
                unmatched_scatter_y.append(label_y)
            elif ass == 'None':
                unmatched_scatter_x.append(label_x)
                unmatched_scatter_y.append(label_y)
            else:
                matched_scatter_x.append(label_x)
                matched_scatter_y.append(label_y)

        # set colors for each peak
        cmap = matplotlib.cm.tab10
        colors = [cmap(i) for i in color_id]

        # create zoomed out scatter plot (just major peak for each envelope)
        self.ax.scatter(unmatched_scatter_x, unmatched_scatter_y, c='grey', marker='o', alpha=0.3)
        self.ax.scatter(matched_scatter_x, matched_scatter_y, c=colors, marker='o', alpha=0.9)

        # zoom in when selecting a specific peak
        if selected_index is not None:
            selected_scatter_x = peak_list['fitter_theo_x'].iloc[selected_index]
            selected_scatter_y = peak_list['fitter_theo_y'].iloc[selected_index]
            if type(selected_scatter_x) != str:
                selected_scatter_x = peak_list['theo_x'].iloc[selected_index]
                selected_scatter_y = peak_list['theo_y'].iloc[selected_index]
            selected_scatter_x = ast.literal_eval(selected_scatter_x)
            selected_scatter_y = ast.literal_eval(selected_scatter_y)

            min_mz = min(selected_scatter_x) - 3 * (selected_scatter_x[1] - selected_scatter_x[0])
            max_mz = max(selected_scatter_x) + 4 * (selected_scatter_x[1] - selected_scatter_x[0])
            self.ax.set_xlim(min_mz, max_mz)
            self.ax.set_ylim(0, max(selected_scatter_y) * 1.25)

            selected_ass = peak_list['ion'].iloc[selected_index]
            if type(selected_ass) != str:
                self.ax.scatter(
                    selected_scatter_x,
                    selected_scatter_y,
                    c='grey',
                    marker='o',
                    alpha=1.0
                )

            else:
                if type(peak_list['adduct'].iloc[selected_index]) == str:
                    selected_ass += f" +{peak_list['adduct'].iloc[selected_index]}"

                if type(peak_list['loss'].iloc[selected_index]) == str:
                    selected_ass += f" -{peak_list['loss'].iloc[selected_index]}"

                color_id = peak_list['color_id'].iloc[selected_index]
                color_choice = cmap(color_id)

                self.ax.scatter(
                    selected_scatter_x,
                    selected_scatter_y,
                    c="red",
                    marker='o',
                    alpha=1.0
                )

                max_y = max(selected_scatter_y)
                arg = selected_scatter_y.index(max_y)
                label_x = selected_scatter_x[arg]
                label_y = selected_scatter_y[arg]

                self.ax.annotate(
                    selected_ass,
                    xy=(label_x, label_y * 1.02),
                    ha='center',
                    va='bottom',
                    c="red"
                )

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
