import statistics
import random
import os
import math
import threading
import wx
import wx.grid
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

from databaseSearchFunctions import DbSearchFunctions

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"


class DbSearchWindow(wx.Frame):
    def __init__(
        self,
        parent,
        title,
        file_path,
        directory_path,
        database_file
    ):
        super().__init__(parent, title=title)

        self.file_path = file_path
        self.directory_path = directory_path
        self.database_file = database_file

        #screen width variable to allow for smaller screens
        screen_width, screen_height = wx.DisplaySize()
        self.SetSize((screen_width * 3 // 4, screen_height * 3 // 4))

        self.panel = wx.Panel(self)
        self.plot_panel = DbSearchSpectrumPanel(self.panel)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # top_sizer: information, datagrid, and fragment map plotting options
        top_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # menu_sizer: information, search button, and search options
        menu_sizer = wx.BoxSizer(wx.VERTICAL)
        title= wx.StaticText(self.panel, label="Database Search")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext1 = wx.StaticText(
            self.panel,
            label="Match observed product ions to theoretical fragments using an open search."
        )
        subtext2 = wx.StaticText(
            self.panel,
            label="See documentation for full details."
        )
        search_button = wx.Button(self.panel, label="Run search",  size=(150, 30))
        search_button.Bind(wx.EVT_BUTTON, self.on_search_button)

        # options_sizer: options for db searching
        options_staticbox = wx.StaticBox(self.panel)
        options_sizer = wx.StaticBoxSizer(options_staticbox, wx.VERTICAL)
        options_text = wx.StaticText(self.panel, label="Search options")
        subtitle_font = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)
        options_text.SetFont(subtitle_font)

        # options_subsizer: option inputs
        options_subsizer = wx.BoxSizer(wx.HORIZONTAL)
        accuracy_text = wx.StaticText(self.panel, label="Accuracy (ppm): ")
        self.accuracy_text_ctrl = wx.TextCtrl(self.panel, size = (30, 20), value = "10")
        mw_text = wx.StaticText(self.panel, label="MW range (kDa): ")
        self.min_mw_text_ctrl = wx.TextCtrl(self.panel, size = (30, 20), value = "5")
        dash_text = wx.StaticText(self.panel, label = "-")
        self.max_mw_text_ctrl = wx.TextCtrl(self.panel, size = (30, 20), value = "100")
        ion_type_text = wx.StaticText(self.panel, label = "Ion type:")
        ion_type_choices = ["b/y", "c/z•"]
        self.ion_type_dropdown = wx.ComboBox(
            self.panel,
            choices=ion_type_choices,
            style=wx.CB_READONLY
        )
        self.ion_type_dropdown.SetValue("b/y")

        options_subsizer.Add(accuracy_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        options_subsizer.Add(self.accuracy_text_ctrl, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)
        options_subsizer.Add(mw_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        options_subsizer.Add(self.min_mw_text_ctrl, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        options_subsizer.Add(dash_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        options_subsizer.Add(self.max_mw_text_ctrl, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)

        options_subsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        disulfide_text = wx.StaticText(self.panel, label = "Disulfides:")
        disulfide_choices = ["Oxidised", "Reduced"]
        self.disulfide_dropdown = wx.ComboBox(
            self.panel,
            choices=disulfide_choices,
            style=wx.CB_READONLY
        )
        self.disulfide_dropdown.SetValue("Oxidised")

        options_subsizer2.Add(disulfide_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        options_subsizer2.Add(self.disulfide_dropdown, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)
        options_subsizer2.Add(ion_type_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 5)
        options_subsizer2.Add(self.ion_type_dropdown, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)

        options_sizer.Add(options_text, 0, wx.BOTTOM, 5)
        options_sizer.Add(options_subsizer, 0, wx.BOTTOM,  10)
        options_sizer.Add(options_subsizer2, 0, wx.BOTTOM,  5)

        menu_sizer.Add(title, 0, wx.BOTTOM, 5)
        menu_sizer.Add(subtext1, 0, wx.BOTTOM, 5)
        menu_sizer.Add(subtext2, 0, wx.BOTTOM, 10)
        menu_sizer.Add(search_button, 0, wx.BOTTOM |  wx.ALIGN_CENTER_HORIZONTAL, 5)
        menu_sizer.Add(options_sizer, 0, wx.EXPAND | wx.BOTTOM, 5)

        self.data_grid = wx.grid.Grid(self.panel)
        self.data_grid.SetDefaultCellBackgroundColour("#f0f0f0")
        self.data_grid.Bind(wx.grid.EVT_GRID_SELECT_CELL, self.on_protein_cell_select)

        # map_menu_sizer: fragment map options and plotting button
        map_menu_sizer = wx.BoxSizer(wx.VERTICAL)
        map_options_staticbox = wx.StaticBox(self.panel)
        map_options_sizer = wx.StaticBoxSizer(map_options_staticbox, wx.VERTICAL)
        map_options_text = wx.StaticText(self.panel, label="Plot options")
        map_options_text.SetFont(subtitle_font)

        map_rpr_sizer = wx.BoxSizer(wx.HORIZONTAL)
        map_rpr_text = wx.StaticText(self.panel, label="Residues per rows: ")
        self.map_rpr_text_ctrl = wx.TextCtrl(self.panel, size=(30, 20), value="20")

        map_rpr_sizer.Add(map_rpr_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        map_rpr_sizer.Add(self.map_rpr_text_ctrl, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)

        show_label_sizer = wx.BoxSizer(wx.HORIZONTAL)
        show_label_text = wx.StaticText(self.panel, label="Show ion labels: ")
        show_label_choices = ["True", "False"]
        self.show_label_dropdown = wx.ComboBox(
            self.panel,
            choices=show_label_choices,
            style=wx.CB_READONLY
        )
        self.show_label_dropdown.SetValue("True")

        show_label_sizer.Add(show_label_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        show_label_sizer.Add(self.show_label_dropdown, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)

        map_options_sizer.Add(map_options_text, 0, wx.BOTTOM, 5)
        map_options_sizer.Add(map_rpr_sizer, 0, wx.BOTTOM,  5)
        map_options_sizer.Add(show_label_sizer, 0, wx.BOTTOM,  5)

        map_button = wx.Button(self.panel, label="Open sequence map", size=(200, 30))
        map_button.Bind(wx.EVT_BUTTON, self.on_fragment_map_button)

        map_menu_sizer.Add(map_options_sizer, 0, wx.EXPAND | wx.BOTTOM, 5)
        map_menu_sizer.Add(map_button, 0, wx.BOTTOM |  wx.ALIGN_CENTER_HORIZONTAL, 5)

        top_sizer.Add(menu_sizer, 0, wx.RIGHT, 15)
        top_sizer.Add(self.data_grid, 1, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 15)
        top_sizer.Add(map_menu_sizer, 0, wx.RIGHT, 15)
        top_sizer.SetItemMinSize(self.data_grid, -1, 180)

        main_sizer.Add(top_sizer, 0, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(self.plot_panel, 1, wx.EXPAND | wx.TOP, -10)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.load_spectrum(file_path)
        self.run_tracker = False # used to check if search has been run- grid cant be easily updated

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Show()

        self.min_mw = 0
        self.max_mw = 0
        self.accuracy = 0
        self.ion_type = ""
        self.disulfide = ""

        self.matched_ids = []
        self.matched_sequences = []
        self.matched_n_term_mods =[]
        self.matched_mws = []
        self.matched_disulfide_positions = []
        self.num_matched_b_ions = []
        self.num_matched_y_ions = []


    def on_size(self, event):
        size = self.panel.GetClientSize()
        self.plot_panel.SetSize(size.width // 2, size.height)
        self.panel.Layout()
        event.Skip()


    def load_spectrum(self, file_path):
        self.basename = os.path.basename(file_path).replace(".txt", "")
        save_path = os.path.join(self.directory_path, f"{self.basename}.spectrum.txt")
        self.spectrum = np.loadtxt(save_path)
        self.min_mz = min(self.spectrum[:, 0])
        self.max_mz = max(self.spectrum[:, 0])
        self.abu = max(self.spectrum[:,1])
        self.plot_panel.plot_spectrum(
            self.spectrum,
            self.min_mz,
            self.max_mz,
            self.abu,
            None,
            None,
            None,
            None
        )


    def load_peaks(self, df = True):
        peak_file_path = os.path.join(self.directory_path, f"{self.basename}.filteredPeakList.csv")
        peak_list = pd.read_csv(peak_file_path)

        # sort here for binary search and subsequent indexing
        peak_list = peak_list.sort_values(by="monoisotopic_mw")

        peak_list = peak_list[peak_list["prediction"] == True]

        # vary if just the list of masses or full dataframe is required
        if df:
            pass
        else:
            peak_list = peak_list["monoisotopic_mw"].tolist()

        return peak_list


    def on_search_button(self, event):
        background_thread = threading.Thread(target=self.db_search)
        background_thread.start()


    def db_search(self):
        if self.run_tracker:
            print(
                "Already completed a database search. "
                "Please reopen the Protein ID window if you need to adjust settings."
            )
            return

        # load search parameters as variables
        self.min_mw = float(self.min_mw_text_ctrl.GetValue())*1000
        self.max_mw = float(self.max_mw_text_ctrl.GetValue())*1000
        self.accuracy = float(self.accuracy_text_ctrl.GetValue())
        self.ion_type = self.ion_type_dropdown.GetValue()
        self.disulfide = self.disulfide_dropdown.GetValue()

        if self.ion_type == "b/y":
            ion_type_file_name = "by"
        elif self.ion_type == "c/z•":
            ion_type_file_name = "cz"

        if self.disulfide == "Oxidised":
            disulfide = "ox"
        elif self.disulfide == "Reduced":
            disulfide = "red"

        h5_basename = os.path.basename(self.database_file).replace(".xml", "")
        h5_dir = os.path.dirname(self.database_file)
        h5_save_path = os.path.join(
            h5_dir,
            f"{h5_basename}.proteoforms.{int(self.min_mw//1000)}_{int(self.max_mw//1000)}kDa.cys_{disulfide}.{ion_type_file_name}_ions.h5"
        )

        db_search = DbSearchFunctions()

        # disulfide positions is a list of fixed mass shifts to apply to specific residues.
        # here it is just to factor in hydorgen loss when forming disulfide bonds
        # write_read_proteoform_file function could be edited to enable other
        # fixed modifications in search

        (
            proteoform_ids,
            proteoform_sequences,
            proteoform_n_term_mods,
            proteoform_mws,
            proteoform_disulfide_positions,
            theo_b_ions,
            theo_y_ions
        ) = db_search.write_read_proteoform_file(
            h5_save_path,
            self.database_file,
            self.min_mw,
            self.max_mw,
            self.ion_type,
            disulfide
        )

        # run database search
        print("Searching database...")
        observed_peaks = self.load_peaks(df = False) # load observed peaks from deconvolution


        (
            self.matched_ids,
            self.matched_sequences,
            self.matched_n_term_mods,
            self.matched_mws,
            self.matched_disulfide_positions,
            self.num_matched_b_ions,
            self.num_matched_y_ions
        ) = db_search.database_search(
            observed_peaks,
            theo_b_ions,
            theo_y_ions,
            self.accuracy,
            proteoform_ids,
            proteoform_sequences,
            proteoform_n_term_mods,
            proteoform_mws,
            proteoform_disulfide_positions
        )

        df = pd.DataFrame({
            "Proteoform ID": self.matched_ids,
            "Molecular weight (Da)": self.matched_mws,
            "# b ions": self.num_matched_b_ions,
            "# y ions": self.num_matched_y_ions,
        })
        df["Total # ions"] = df["# b ions"] + df["# y ions"]

        self.Freeze()
        num_rows, num_cols = df.shape
        self.data_grid.CreateGrid(num_rows, num_cols)
        for row in range(num_rows):
            for col in range(num_cols):
                value = str(df.iloc[row, col])
                self.data_grid.SetCellValue(row, col, value)
        column_names = [
            "ID",
            "Molecular weight", 
            "# b ions",
            "# y ions",
            "Total # ions",
            ]
        for col, column_name in enumerate(column_names):
            self.data_grid.SetColLabelValue(col, column_name)
        self.data_grid.AutoSizeColumns()
        self.data_grid.Refresh()
        self.data_grid.Show()

        self.panel.Layout()
        self.Thaw()

        self.run_tracker = True # db search completed


    def on_protein_cell_select(self, event):
        self.selected_row = event.GetRow()

        # theoretical ions
        sequence = self.matched_sequences[self.selected_row]
        n_terminal_modification = self.matched_n_term_mods[self.selected_row]
        disulfide_positions = self.matched_disulfide_positions[self.selected_row]

        # nFPS isn't applicable to this data
        if self.ion_type == "b/y":
            nfps = self.calculate_nfps(sequence, n_terminal_modification, disulfide_positions)
        else:
            nfps = None

        b_ions, y_ions = DbSearchFunctions.theoretical_frag_generator_single_search(
            sequence,
            n_terminal_modification,
            disulfide_positions,
            self.ion_type
        )

        # observed ions
        observed_peaks = self.load_peaks()
        mass_list = observed_peaks["monoisotopic_mw"].tolist()

        # matching
        (
            matched_theoretical_b_ions,
            matched_observed_b_ions
        ) = DbSearchFunctions.report_matches_for_protein(
            b_ions[0],
            mass_list,
            self.accuracy,
            all_charge=True
        )
        (
            matched_theoretical_y_ions,
            matched_observed_y_ions
        ) = DbSearchFunctions.report_matches_for_protein(
            y_ions[0],
            mass_list,
            self.accuracy,
            all_charge=True
        )

        # generate some lists to extract information from
        charge_list = observed_peaks["charge"].tolist()
        abundance_list = observed_peaks["abundance"].tolist()

        # lists to store matched ion information in
        b_ion_list = []
        b_ion_name_list = []
        y_ion_list = []
        y_ion_name_list = []

        # instead of editing plot, we can just change if we write ion names to list
        show_label = self.show_label_dropdown.GetValue()

        # obtain theoretical distribution and name of each matched ion
        for i, index in enumerate(matched_theoretical_b_ions):
            b_ion_mass = b_ions[0][index]
            observed_ion_index = matched_observed_b_ions[i]
            charge = charge_list[observed_ion_index]
            abundance = abundance_list[observed_ion_index]
            if self.ion_type == "b/y":
                name = f"b{index+1} {int(charge)}+"
            elif self.ion_type == "c/z•":
                name = f"c{index+1} {int(charge)}+"

            # calculate averagine distribtion using mass, charge, and abundance
            averagine_envelope = DbSearchFunctions.averagine(b_ion_mass, charge, abundance)
            b_ion_list.append(averagine_envelope)

            if show_label == "True":
                b_ion_name_list.append(name)
            else:
                b_ion_name_list.append("")

        for i, index in enumerate(matched_theoretical_y_ions):
            y_ion_mass = y_ions[0][index]
            observed_ion_index = matched_observed_y_ions[i]
            charge = charge_list[observed_ion_index]
            abundance = abundance_list[observed_ion_index]
            if self.ion_type == "b/y":
                name = f"y{index+1} {int(charge)}+"
            elif self.ion_type == "c/z•":
                name = f"z.{index+1} {int(charge)}+"

            averagine_envelope = DbSearchFunctions.averagine(y_ion_mass, charge, abundance)
            y_ion_list.append(averagine_envelope)

            if show_label == "True":
                y_ion_name_list.append(name)
            else:
                y_ion_name_list.append("")

        self.plot_panel.plot_spectrum(
            self.spectrum,
            self.min_mz,
            self.max_mz,
            self.abu,
            b_ion_list,
            b_ion_name_list,
            y_ion_list,
            y_ion_name_list,
            nfps
        )


    def calculate_nfps(self, forward_sequence, n_terminal_modification, disulfide_positions):
        # helper function to generate decoy sequences
        def generate_shuffled_string(input_string):
            characters = list(input_string)
            random.shuffle(characters)
            shuffled_string = ''.join(characters)
            return shuffled_string

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

        # load in observed data
        observed_peaks = self.load_peaks()
        mass_list = observed_peaks["monoisotopic_mw"].tolist()

        # determine matched ions and their positions in the sequence
        b_ions, y_ions = DbSearchFunctions.theoretical_frag_generator_single_search(
            forward_sequence,
            n_terminal_modification,
            disulfide_positions
        )
        matched_theoretical_b_ions, _ = DbSearchFunctions.report_matches_for_protein(
            b_ions[0],
            mass_list,
            self.accuracy
        )
        matched_theoretical_y_ions, _ = DbSearchFunctions.report_matches_for_protein(
            y_ions[0],
            mass_list,
            self.accuracy
        )
        matched_positions = [] # indices of locations that are c terminal to cleavage site
        for index in matched_theoretical_b_ions:
            matched_positions.append(index + 1)
        for index in matched_theoretical_y_ions:
            matched_positions.append(len(forward_sequence) - (index + 1))

        # calculate product of fragmentation propensities for the forward sequence
        propensities = []
        for position in matched_positions:
            n_term_res = forward_sequence[position - 1]
            c_term_res = forward_sequence[position]
            propensities.append(n_term_props[n_term_res] / 100 * c_term_props[c_term_res] / 100)
        product = 1
        for num in propensities:
            product *= num
        forward_product = product

        # calculate product of fragmentation propensities for 100 decoy sequences
        decoy_products = []
        for _ in range(100):
            decoy_sequence = generate_shuffled_string(forward_sequence)
            propensities = []
            for position in matched_positions:
                n_term_res = decoy_sequence[position - 1]
                c_term_res = decoy_sequence[position]
                propensities.append(n_term_props[n_term_res] / 100 * c_term_props[c_term_res] / 100)
            product = 1
            for num in propensities:
                product *= num
            decoy_products.append(product)

        # calculate nfps
        decoy_geo_mean = statistics.geometric_mean(decoy_products)
        nfps = math.log10(forward_product / decoy_geo_mean)

        return nfps


    def on_fragment_map_button(self, event):
        # theoretical ions
        sequence = self.matched_sequences[self.selected_row]
        n_terminal_modification = self.matched_n_term_mods[self.selected_row]
        disulfide_positions = self.matched_disulfide_positions[self.selected_row]

        b_ions, y_ions = DbSearchFunctions.theoretical_frag_generator_single_search(
            sequence,
            n_terminal_modification,
            disulfide_positions,
            self.ion_type
        )

        # observed ions
        observed_peaks = self.load_peaks()
        mass_list = observed_peaks["monoisotopic_mw"].tolist()

        # matching
        matched_theoretical_b_ions, _ = DbSearchFunctions.report_matches_for_protein(
            b_ions[0],
            mass_list,
            self.accuracy
        )
        matched_theoretical_y_ions, _ = DbSearchFunctions.report_matches_for_protein(
            y_ions[0],
            mass_list,
            self.accuracy
        )

        # convert ion indices to positions on the sequence
        b_fragment_positions = [index + 1 for index in matched_theoretical_b_ions]
        y_fragment_positions = [len(sequence) - (index + 2) for index in matched_theoretical_y_ions]

        # check to see if n terminus is modified for colouring
        n_term_mod = n_terminal_modification >= 1

        residues_per_row = int(self.map_rpr_text_ctrl.GetValue())

        # generate plot
        FragmentMapper.plot_fragment_map(
            sequence,
            b_fragment_positions,
            y_fragment_positions,
            n_term_mod,
            residues_per_row
        )



class DbSearchSpectrumPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)

        sizer = wx.BoxSizer(wx.VERTICAL)

        self.figure = Figure(facecolor="#f0f0f0")
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.toolbar = NavigationToolbar2Wx(self.canvas)

        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.EXPAND)

        self.SetSizer(sizer)
        self.Bind(wx.EVT_SIZE, self.on_size)
        self.Layout()
        self.toolbar.Realize()


    def plot_spectrum(
        self, spectrum, selected_min_mz, selected_max_mz,
        abu, b_ion_list, b_ion_name_list, y_ion_list, y_ion_name_list, nfps = None
    ):
        self.ax.cla()
        self.ax.plot(spectrum[:, 0], spectrum[:, 1], color="black", linewidth=1)
        self.ax.set_xlabel("m/z")
        self.ax.set_ylabel("Intensity")
        self.ax.set_facecolor("#f0f0f0")
        self.ax.set_xlim(selected_min_mz-100, selected_max_mz+100)
        self.max_intens = abu * 1.1
        self.ax.set_ylim(0, self.max_intens)

        if nfps == None:
            pass
        elif nfps > 5:
            label_text = f"nFPS: {round(nfps, 2)}"
            self.ax.annotate(
                label_text,
                xy=(0.995,0.99),
                xycoords="axes fraction",
                ha="right",
                va="top",
                color="green"
            )
        elif 2 <= nfps <= 5:
            label_text = f"nFPS: {round(nfps, 2)}"
            self.ax.annotate(
                label_text,
                xy=(0.995,0.99),
                xycoords="axes fraction",
                ha="right",
                va="top",
                color='orange'
            )
        elif nfps < 2:
            label_text = f"nFPS: {round(nfps, 2)}"
            self.ax.annotate(
                label_text,
                xy=(0.995,0.99),
                xycoords="axes fraction",
                ha="right",
                va="top",
                color='red'
            )

        # plotting b ion averagine envelopes
        if b_ion_list == None:
            pass
        else:
            for i, envelope in enumerate(b_ion_list):
                x = envelope[0]
                y = envelope[1]
                x_label = x[0]
                y_label = max(y) * 1.02
                label = b_ion_name_list[i]
                self.ax.scatter(x, y, c="#205283", marker="o", alpha=0.5)
                self.ax.text(x_label, y_label, label, c="#205283")

        # plotting y ion averagine envelopes
        if y_ion_list == None:
            pass
        else:
            for i, envelope in enumerate(y_ion_list):
                x = envelope[0]
                y = envelope[1]
                x_label = x[0]
                y_label = y[0] * 1.02
                label = y_ion_name_list[i]
                self.ax.scatter(x, y, c="#b42920", marker="o", alpha=0.5)
                self.ax.text(x_label, y_label, label, c="#b42920")

        self.canvas.draw()


    def on_size(self, event):
        self.fit_plot_to_panel()
        event.Skip()


    def fit_plot_to_panel(self):
        size = self.GetClientSize()
        if size[0] >= 51:
            dpi = self.GetContentScaleFactor() * 100
            width = size.width / dpi
            height = size.height / dpi - 0.3
            self.figure.set_size_inches(width, height)
            self.figure.tight_layout(rect=[0, 0, 1, 1])
            self.canvas.draw()



class FragmentMapper():
    def plot_fragment_map(
        protein_sequence, b_fragment_positions, y_fragment_positions,
        n_term_mod=False, residues_per_row=20):
        num_rows = math.ceil(len(protein_sequence) / residues_per_row)

        font_size = 20
        _, ax = plt.subplots(
            num_rows,
            1,
            figsize=(font_size/20*10*residues_per_row/20, num_rows * (font_size / 20) * 0.5),
            num="precisION - Fragment Map"
        )

        for i, ax_row in enumerate(ax):
            start = i * residues_per_row
            end = min((i + 1) * residues_per_row, len(protein_sequence))
            sequence_chunk = protein_sequence[start:end]

            if i == 0:
                if n_term_mod == True:
                    n_term_colour = "red"
                else:
                    n_term_colour = "grey"
                ax_row.text(
                    -1, 0.5, "N", ha="center", va="center",
                    fontsize=font_size, fontname="Courier New", color=n_term_colour
                )
            else:
                ax_row.text(-0.5, 0.5, start+1, ha="right", va="center",
                    fontsize=font_size, fontname="Courier New", color="grey"
                )

            if i == len(ax)-1:
                pass
            else:
                ax_row.text(
                    residues_per_row-0.5, 0.5, end, ha="left", va="center",
                    fontsize=font_size, fontname="Courier New", color="grey"
                )

            pos = 0
            for j, res in enumerate(sequence_chunk):
                ax_row.text(pos, 0.5, res, ha="center", va="center",
                    fontsize=font_size, fontname="Courier New", weight = "bold"
                )
                pos += 1
                if i == len(ax)-1:
                    if j == len(sequence_chunk)-1:
                        ax_row.text(pos, 0.5, "C", ha="center", va="center",
                            fontsize=font_size, fontname="Courier New", color="grey"
                        )

            ax_row.set_xlim(-1, residues_per_row+1)
            ax_row.axis("off")

            for position in b_fragment_positions:
                if start <= position-1 < end:
                    x = position - start - 0.5
                    ax_row.text(x, 0.7, "⎤", ha="right", va="center",
                        fontsize=font_size, fontname="monospace", color="#205283"
                    )

            for position in y_fragment_positions:
                if start <= position+1 < end:
                    x = position - start + 0.5
                    ax_row.text(x, 0.3, "⎣", ha="left", va="center",
                        fontsize=font_size, fontname="monospace", color="#b42920"
                    )

        fig_manager = plt.get_current_fig_manager()
        fig_manager.window.SetIcon(wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO))

        plt.tight_layout()
        plt.show()
