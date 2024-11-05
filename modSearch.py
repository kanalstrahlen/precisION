import os
import random
import wx
import pandas as pd
from pyteomics import mass
import time

from modTerminalSearchValidation import ModTerminalSearchValidationWindow
from modInternalSearchValidation import ModInternalSearchValidationWindow



class ModSearchWindow(wx.Frame):
    def __init__(
        self, parent, title, name, sequence,
        file_path, directory_path,
        n_mod_add="", n_mod_sub="", c_mod_add="", c_mod_sub=""
    ):
        super().__init__(
            parent,
            title=title,
            size=(490, 425),
            style=wx.DEFAULT_FRAME_STYLE
        )

        self.basename = os.path.basename(file_path).replace(".txt", "")
        self.directory_path = directory_path
        self.file_path = file_path
        self.sequence = sequence
        self.n_mod_add = n_mod_add
        self.n_mod_sub = n_mod_sub
        self.c_mod_add = c_mod_add
        self.c_mod_sub = c_mod_sub
        self.name = name

        panel = wx.Panel(self)
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # window title
        title = wx.StaticText(panel, label="Modified Ion Assignment")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext1 = wx.StaticText(
            panel,
            label="Assign sequence ions that have been biologically or chemically modified."
        )
        subtext2 = wx.StaticText(
            panel,
            label="Modifications can be input as chemical formulas or mass shifts."
        )


        # mod_sizer: mods to consider
        mod_staticbox = wx.StaticBox(panel)
        mod_sizer = wx.StaticBoxSizer(mod_staticbox, wx.VERTICAL)
        mod_header = wx.StaticText(panel, label="Modification")
        subtitle_font = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)
        mod_header.SetFont(subtitle_font)
        self.mod_grid = self.create_modification_grid(panel)

        mod_sizer.Add(mod_header, 0, wx.EXPAND | wx.BOTTOM, 5)
        mod_sizer.Add(self.mod_grid, 0, wx.EXPAND | wx.BOTTOM, 0)
        mod_sizer.SetMinSize((400, -1))

        options_staticbox = wx.StaticBox(panel)
        options_sizer = wx.StaticBoxSizer(options_staticbox, wx.VERTICAL)
        options_header = wx.StaticText(panel, label="Search Options")
        options_header.SetFont(
            wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)
        )

        # options_sizer: first row of matching options
        options_sizer_1 = wx.BoxSizer(wx.HORIZONTAL)

        ppm_label = wx.StaticText(panel, label="Maximum mass error (ppm):")
        self.ppm_input = wx.TextCtrl(panel, value="3", size=(30, -1))

        options_sizer_1.Add(ppm_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        options_sizer_1.Add(self.ppm_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)

        options_sizer_2 = wx.BoxSizer(wx.HORIZONTAL)

        automatic_matching_score_label = wx.StaticText(
            panel, label="Minimum envelope fit score for automatic assignment:"
        )
        self.automatic_matching_score_input = wx.TextCtrl(panel, value="0.75", size=(30, -1))

        options_sizer_2.Add(
            automatic_matching_score_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5
        )
        options_sizer_2.Add(
            self.automatic_matching_score_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 0
        )

        options_sizer_3 = wx.BoxSizer(wx.HORIZONTAL)

        auto_ppm_label = wx.StaticText(
            panel, label="Maximum mass error for automatic assignment (ppm):"
        )
        self.auto_ppm_input = wx.TextCtrl(panel, value="3", size=(30, -1))

        options_sizer_3.Add(auto_ppm_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        options_sizer_3.Add(self.auto_ppm_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 0)

        options_sizer_4 = wx.BoxSizer(wx.HORIZONTAL)

        shared_internal_label = wx.StaticText(
            panel,
            label="Minimum -log(E-value) for internal fragment assignment:"
        )
        self.internal_shared_termini_input = wx.TextCtrl(panel, value="3", size=(30, -1))

        options_sizer_4.Add(
            shared_internal_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5
        )
        options_sizer_4.Add(
            self.internal_shared_termini_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 0
        )

        options_sizer_5 = wx.BoxSizer(wx.HORIZONTAL)

        ion_type_label = wx.StaticText(panel, label="Ion types:")
        ion_type_choices = ["b/y", "c/z•"]
        self.ion_type_dropdown = wx.ComboBox(panel, choices=ion_type_choices, style=wx.CB_READONLY)
        self.ion_type_dropdown.SetValue("b/y")

        options_sizer_5.Add(ion_type_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        options_sizer_5.Add(self.ion_type_dropdown, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 0)

        options_sizer.Add(options_header, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.Add(options_sizer_1, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.Add(options_sizer_2, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.Add(options_sizer_3, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.Add(options_sizer_4, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.Add(options_sizer_5, 0, wx.EXPAND | wx.BOTTOM, 5)

        button_sizer = wx.BoxSizer(wx.HORIZONTAL)

        terminal_match_button = wx.Button(
            panel, label="Assisted Match - Terminals",  size=(170, 30)
        )
        terminal_match_button.Bind(wx.EVT_BUTTON, self.on_terminal_match_button)
        internal_match_button = wx.Button(
            panel, label="Assisted Match - Internals",  size=(170, 30)
        )
        internal_match_button.Bind(wx.EVT_BUTTON, self.on_internal_match_button)

        button_sizer.Add(terminal_match_button, 0, wx.RIGHT, 5)
        button_sizer.Add(internal_match_button, 0, wx.RIGHT, 0)

        main_sizer.Add(title, 0, wx.EXPAND | wx.BOTTOM, 5)
        main_sizer.Add(subtext1, 0, wx.EXPAND | wx.BOTTOM, 3)
        main_sizer.Add(subtext2, 0, wx.EXPAND | wx.BOTTOM, 0)
        main_sizer.Add(mod_sizer, 0, wx.EXPAND | wx.BOTTOM, 5)        
        main_sizer.Add(options_sizer, 0, wx.EXPAND | wx.BOTTOM, 5)
        main_sizer.Add(button_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM, 0)

        panel_sizer.Add(main_sizer, 0, wx.ALL, 5)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        panel.SetSizer(panel_sizer)
        self.Show()

        self.directory_path = directory_path


    def create_modification_grid(self, parent):
        grid = wx.grid.Grid(parent)
        grid.CreateGrid(1, 4)
        grid.SetDefaultCellBackgroundColour("#f0f0f0")

        grid.SetColLabelValue(0, "Modification Name")
        grid.SetColLabelValue(1, "Gain Formula")
        grid.SetColLabelValue(2, "Loss Formula")
        grid.SetColLabelValue(3, "Mass shift (Da)")

        default_data = [
            ("Monoisotope error", "", "", "-1.003"),
        ]

        for i, (name, add_formula, subtract_formula, mass_shift) in enumerate(default_data):
            grid.SetCellValue(i, 0, name)
            grid.SetCellValue(i, 1, add_formula)
            grid.SetCellValue(i, 2, subtract_formula)
            grid.SetCellValue(i, 3, mass_shift)

        for col in range(4):
            grid.AutoSizeColumn(col)

        return grid


    def on_terminal_match_button(self, event):
        peak_assignment_file = os.path.join(
            self.directory_path, f"{self.basename}.assignedPeakList.csv"
        )
        if not os.path.exists(peak_assignment_file):
            print("Unmodified terminal fragments should be assigned first.")
            return
        self.update_protein_name()
        mod_list = self.gen_mod_list()
        for mod in mod_list:
            self.terminal_matching(mod)


    def on_internal_match_button(self, event):
        peak_assignment_file = os.path.join(
            self.directory_path, f"{self.basename}.assignedPeakList.csv"
        )
        if not os.path.exists(peak_assignment_file):
            print("Unmodified terminal fragments should be assigned first.")
            return
        self.update_protein_name()
        mod_list = self.gen_mod_list()
        for mod in mod_list:
            self.internal_matching(mod)


    def update_protein_name(self):
        peak_assignment_file = os.path.join(
            self.directory_path, f"{self.basename}.assignedPeakList.csv"
        )
        peak_list = pd.read_csv(peak_assignment_file)
        names = peak_list['name'].unique()
        if self.name in names:
            index = peak_list.index[peak_list['name'] == self.name].min()
            prev_sequence = peak_list.loc[index, 'sequence']
            if self.sequence != prev_sequence:
                self.name = f'{self.name}_newIsoform{random.randint(10,99)}'
                print("Protein ID has already been used for another sequence!")
                print(f"Updating ID to {self.name}")


    def terminal_matching(self, mod_details):
        val_window = ModTerminalSearchValidationWindow(
            self,
            f"precisION - Modified Ions [{mod_details[0]}]",
            self.file_path,
            self.directory_path,
            self.name,
            self.sequence,
            self.n_mod_add,
            self.n_mod_sub,
            self.c_mod_add,
            self.c_mod_sub,
            mod_details,
            float(self.ppm_input.GetValue()),
            float(self.auto_ppm_input.GetValue()),
            float(self.automatic_matching_score_input.GetValue()),
            self.ion_type_dropdown.GetValue()
        )


    def internal_matching(self, mod_details):

        if self.ion_type_dropdown.GetValue() == "c/z•":
            print("Internal fragments produced by ExD are not currently supported in precisION.")
            print("Support will be added once mechanistic details are further explored.")
            return

        val_window = ModInternalSearchValidationWindow(
            self,
            f"precisION - Modified Internal Ions [{mod_details[0]}]",
            self.file_path,
            self.directory_path,
            self.name,
            self.sequence,
            float(self.ppm_input.GetValue()),
            float(self.auto_ppm_input.GetValue()),
            float(self.automatic_matching_score_input.GetValue()),
            float(self.internal_shared_termini_input.GetValue()),
            mod_details
        )

        val_window.Show()



    def gen_mod_list(self):
        modification_data = []
        for row in range(self.mod_grid.GetNumberRows()):
            name = self.mod_grid.GetCellValue(row, 0)
            add_formula = self.mod_grid.GetCellValue(row, 1)
            subtract_formula = self.mod_grid.GetCellValue(row, 2)
            mass_shift = self.mod_grid.GetCellValue(row, 3)
            if name != "":
                if add_formula == "" and subtract_formula == "":
                    mass_shift = float(mass_shift)
                    name = f"{name} ({round(mass_shift, 4)})"
                    add_formula = "undefined"
                    loss_formula = "undefined"
                elif add_formula != "" and subtract_formula == "":
                    mass_shift = mass.calculate_mass(formula=add_formula)
                    name = f"{name} ({round(mass_shift, 4)})"
                    add_formula = add_formula
                    loss_formula = ""                
                elif add_formula == "" and subtract_formula != "":
                    mass_shift = -1 * mass.calculate_mass(formula=subtract_formula)
                    name = f"{name} ({mass_shift})"
                    add_formula = ""
                    loss_formula = subtract_formula      
                elif add_formula != "" and subtract_formula != "":
                    mass_shift = mass.calculate_mass(formula=add_formula) - mass.calculate_mass(formula=subtract_formula)
                    name = f"{name} ({round(mass_shift, 4)})"
                    add_formula = add_formula
                    loss_formula = subtract_formula 

                modification_data.append([name, add_formula, loss_formula, mass_shift])

        return modification_data