import pandas as pd
import wx


class FilterAssignmentsWindow(wx.Frame):
    def __init__(self, parent, title, assignment_file):
        super().__init__(
            parent,
            title=title,
            size=(340, 365),
            style=wx.DEFAULT_FRAME_STYLE
        )

        self.assignment_file = assignment_file
        assignment_list = pd.read_csv(assignment_file)

        names = ["Any protein ID"] + list(assignment_list['name'].dropna().unique())

        fixed_mods = []
        try:
            fixed_mods += list(["N "] + assignment_list['n_mod'].dropna().unique())
        except:
            pass
        try:
            fixed_mods += list(["C "] + assignment_list['c_mod'].dropna().unique())
        except:
            pass
        fixed_mods = ["Any fixed modification"] + fixed_mods
        variable_mods = ["Any variable modification"]
        try:
            variable_mods += list(assignment_list['variable_mod'].dropna().unique())
        except:
            pass
        ions = list(set([ion[0] for ion in assignment_list['ion'].dropna().unique()]))
        ion_types = ["Any ion type"] + ions

        panel = wx.Panel(self)
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # window title
        title = wx.StaticText(panel, label="Filter Assignments")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext = wx.StaticText(
            panel,
            label="Filter assignments to exclude outliers and correct mistakes."
        )

        options_staticbox = wx.StaticBox(panel)
        options_sizer = wx.StaticBoxSizer(options_staticbox, wx.VERTICAL)
        options_header = wx.StaticText(panel, label="Filtering Options")
        options_header.SetFont(
            wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)
        )

        ppm_sizer = wx.BoxSizer(wx.HORIZONTAL)

        ppm_label = wx.StaticText(panel, label="Maximum mass error (ppm):")
        self.ppm_input = wx.TextCtrl(panel, value="3", size=(30, -1))
        ppm_button = wx.Button(panel, label="Filter",  size=(100, 25))
        ppm_button.Bind(wx.EVT_BUTTON, self.on_ppm_button)

        ppm_sizer.Add(ppm_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        ppm_sizer.Add(self.ppm_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        ppm_sizer.Add(ppm_button, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)

        fit_sizer = wx.BoxSizer(wx.HORIZONTAL)

        fit_label = wx.StaticText(panel, label="Minimum fit score:")
        self.fit_input = wx.TextCtrl(panel, value="0.0", size=(30, -1))
        fit_button = wx.Button(panel, label="Filter",  size=(100, 25))
        fit_button.Bind(wx.EVT_BUTTON, self.on_fit_button)

        fit_sizer.Add(fit_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        fit_sizer.Add(self.fit_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        fit_sizer.Add(fit_button, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)

        set_sizer = wx.BoxSizer(wx.VERTICAL)

        set_label = wx.StaticText(panel, label="Set of ions")
        set_label.SetFont(
            wx.Font(9, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD)
        )

        self.name_dropdown = wx.ComboBox(
            panel,
            choices=names,
            style=wx.CB_READONLY
        )
        self.name_dropdown.SetValue(names[0])
        self.fixed_mod_dropdown = wx.ComboBox(
            panel,
            choices=fixed_mods,
            style=wx.CB_READONLY
        )
        self.fixed_mod_dropdown.SetValue(fixed_mods[0])
        self.variable_mod_dropdown = wx.ComboBox(
            panel,
            choices=variable_mods,
            style=wx.CB_READONLY
        )
        self.variable_mod_dropdown.SetValue(variable_mods[0])
        self.ion_type_dropdown = wx.ComboBox(
            panel,
            choices=ion_types,
            style=wx.CB_READONLY
        )
        self.ion_type_dropdown.SetValue(ion_types[0])

        set_sizer.Add(set_label, 0, wx.BOTTOM, 5)
        set_sizer.Add(self.name_dropdown, 0, wx.BOTTOM, 5)
        set_sizer.Add(self.fixed_mod_dropdown, 0, wx.BOTTOM, 5)
        set_sizer.Add(self.variable_mod_dropdown, 0, wx.BOTTOM, 5)
        set_sizer.Add(self.ion_type_dropdown, 0, wx.BOTTOM, 5)

        set_button = wx.Button(panel, label="Filter",  size=(100, 25))
        set_button.Bind(wx.EVT_BUTTON, self.on_set_button)

        options_sizer.Add(options_header, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.Add(ppm_sizer, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.Add(fit_sizer, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.Add(set_sizer, 0, wx.EXPAND | wx.BOTTOM, 5)
        options_sizer.Add(set_button, 0, wx.BOTTOM, 5)

        main_sizer.Add(title, 0, wx.EXPAND | wx.BOTTOM, 5)
        main_sizer.Add(subtext, 0, wx.EXPAND | wx.BOTTOM, 0)
        main_sizer.Add(options_sizer, 0, wx.EXPAND | wx.BOTTOM, 0)

        panel_sizer.Add(main_sizer, 0, wx.ALL, 5)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        panel.SetSizer(panel_sizer)
        self.Show()


    def on_ppm_button(self, event):
        assignment_list = pd.read_csv(self.assignment_file)
        ppm_threshold = float(self.ppm_input.GetValue())
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
        assignment_list.to_csv(self.assignment_file, index=False)
        print(f"Removed assignments with mass errors greater than {ppm_threshold} ppm.")


    def on_fit_button(self, event):
        assignment_list = pd.read_csv(self.assignment_file)
        fit_threshold = float(self.fit_input.GetValue())
        mask = assignment_list['fit_score'] < fit_threshold
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
        assignment_list.to_csv(self.assignment_file, index=False)
        print(f"Removed assignments with fitting scores less than {fit_threshold}.")


    def on_set_button(self, event):
        assignment_list = pd.read_csv(self.assignment_file)
        name = str(self.name_dropdown.GetValue())
        fixed_mod = str(self.fixed_mod_dropdown.GetValue())
        variable_mod = str(self.variable_mod_dropdown.GetValue())
        ion_type = str(self.ion_type_dropdown.GetValue())

        if name == "Any protein ID":
            df_name = assignment_list
        else:
            df_name = assignment_list[assignment_list["name"] == name]

        if fixed_mod == "Any fixed modification":
            df_fixed = df_name
        else:
            fixed_mod_terminus = fixed_mod[0]
            fixed_mod = fixed_mod[2:]
            if fixed_mod_terminus == "N":
                df_fixed = df_name[df_name["n_mod"] == fixed_mod]
            elif fixed_mod_terminus == "C":
                df_fixed = df_name[df_name["c_mod"] == fixed_mod]

        if variable_mod == "Any variable modification":
            df_variable = df_fixed
        else:
            df_variable = df_fixed[df_fixed["variable_mod"] == variable_mod]

        if ion_type == "Any ion type":
            df_type = df_variable
        else:
            df_type = df_variable[df_variable["ion"].str[0] == ion_type]

        assignment_list.loc[df_type.index, "name"] = None
        assignment_list.loc[df_type.index, "sequence"] = None
        assignment_list.loc[df_type.index, "n_mod"] = None
        assignment_list.loc[df_type.index, "c_mod"] = None
        assignment_list.loc[df_type.index, "ion"] = None
        assignment_list.loc[df_type.index, "variable_mod"] = None
        assignment_list.loc[df_type.index, "variable_mod_gain"] = None
        assignment_list.loc[df_type.index, "variable_mod_loss"] = None
        assignment_list.loc[df_type.index, "adduct"] = None
        assignment_list.loc[df_type.index, "loss"] = None
        assignment_list.loc[df_type.index, "theoretical_mz"] = None
        assignment_list.loc[df_type.index, "ppm_error"] = None
        assignment_list.loc[df_type.index, "fit_score"] = None
        assignment_list.loc[df_type.index, "total_intensity"] = None
        assignment_list.loc[df_type.index, "charge_scaled_abundance"] = None
        assignment_list.loc[df_type.index, "fitter_formula"] = None
        assignment_list.loc[df_type.index, "fitter_theo_x"] = None
        assignment_list.loc[df_type.index, "fitter_theo_y"] = None
        assignment_list.loc[df_type.index, "frag_site"] = None

        assignment_list.to_csv(self.assignment_file, index=False)
        print("Removed selected set of ions.")
