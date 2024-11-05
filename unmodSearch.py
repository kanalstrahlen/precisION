import random
import os
import shutil
import wx
import pandas as pd

from unmodSearchValidation import UnmodSearchValidationWindow
from internalSearchValidation import InternalSearchValidationWindow

# to do; add c/z to everything

class NoModSearchWindow(wx.Frame):
    def __init__(
        self, parent, title, name, sequence,
        file_path, directory_path,
        n_mod_add="", n_mod_sub="", c_mod_add="", c_mod_sub=""
    ):
        super().__init__(
            parent,
            title=title,
            size=(450, 330),
            style=wx.DEFAULT_FRAME_STYLE & ~(wx.RESIZE_BORDER | wx.MAXIMIZE_BOX)
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
        title = wx.StaticText(panel, label="Unmodified Ion Assignment")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext1 = wx.StaticText(
            panel,
            label="Sequentially match unmodified ions, neutral losses, and then internal fragments."
        )
        subtext2 = wx.StaticText(
            panel,
            label="Spectra will be recalibrated using unmodified terminal ions and neutral losses."
        )

        options_staticbox = wx.StaticBox(panel)
        options_sizer = wx.StaticBoxSizer(options_staticbox, wx.VERTICAL)
        options_header = wx.StaticText(panel, label="Search Options")
        options_header.SetFont(
            wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)
        )

        # options_sizer: first row of matching options
        options_sizer_1 = wx.BoxSizer(wx.HORIZONTAL)

        ppm_label = wx.StaticText(panel, label="Maximum mass error (ppm):")
        self.rough_ppm_input = wx.TextCtrl(panel, value="10", size=(30, -1))
        pre_cal_label = wx.StaticText(panel, label="(rough)")
        self.fine_ppm_input = wx.TextCtrl(panel, value="3", size=(30, -1))
        post_cal_label = wx.StaticText(panel, label="(fine)")


        options_sizer_1.Add(ppm_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        options_sizer_1.Add(self.rough_ppm_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        options_sizer_1.Add(pre_cal_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        options_sizer_1.Add(self.fine_ppm_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        options_sizer_1.Add(post_cal_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)

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
        self.auto_ppm_input = wx.TextCtrl(panel, value="1.5", size=(30, -1))

        options_sizer_3.Add(auto_ppm_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        options_sizer_3.Add(self.auto_ppm_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 0)

        options_sizer_4 = wx.BoxSizer(wx.HORIZONTAL)

        shared_internal_label = wx.StaticText(
            panel,
            label="Minimum -log(E-value) for internal fragment assignment:"
        )
        self.internal_shared_termini_input = wx.TextCtrl(panel, value="3", size=(30, -1))

        options_sizer_4.Add(shared_internal_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
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
        button_sizer.Add(internal_match_button, 0, wx.RIGHT, 5)

        main_sizer.Add(title, 0, wx.EXPAND | wx.BOTTOM, 5)
        main_sizer.Add(subtext1, 0, wx.EXPAND | wx.BOTTOM, 3)
        main_sizer.Add(subtext2, 0, wx.EXPAND | wx.BOTTOM, 0)
        main_sizer.Add(options_sizer, 0, wx.EXPAND | wx.BOTTOM, 5)
        main_sizer.Add(button_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM, 0)

        panel_sizer.Add(main_sizer, 0, wx.ALL, 5)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        panel.SetSizer(panel_sizer)
        self.Show()

        self.directory_path = directory_path


    def on_terminal_match_button(self, event):
        self.generate_assignment_files()
        self.update_protein_name()
        self.terminal_matching()


    def on_internal_match_button(self, event):
        peak_assignment_file = os.path.join(
            self.directory_path, f"{self.basename}.assignedPeakList.csv"
        )
        if not os.path.exists(peak_assignment_file):
            print("Terminal fragments should be assigned first.")
            return
        self.update_protein_name()
        self.internal_matching()


    def generate_assignment_files(self):
        uncalibrated_spectrum_file = os.path.join(
            self.directory_path, f"{self.basename}.spectrum.txt"
        )
        calibrated_spectrum_file = os.path.join(
            self.directory_path, f"{self.basename}.recalibratedSpectrum.txt"
        )

        if os.path.exists(calibrated_spectrum_file):
            pass
        else:
            shutil.copy(uncalibrated_spectrum_file, calibrated_spectrum_file)

        filtered_peak_file = os.path.join(
            self.directory_path, f"{self.basename}.filteredPeakList.csv"
        )
        peak_assignment_file = os.path.join(
            self.directory_path, f"{self.basename}.assignedPeakList.csv"
        )

        if os.path.exists(peak_assignment_file):
            pass
        else:
            peak_list = pd.read_csv(filtered_peak_file)
            peak_list = peak_list[peak_list["prediction"] == True]

            columns_to_keep = [
                "monoisotopic_mw",
                "charge",
                "abundance",
                "monoisotopic_mz",
                "theo_x",
                "theo_y"
            ]

            peak_list = peak_list[columns_to_keep]
            peak_list["name"] = None
            peak_list["sequence"] = None
            peak_list["n_mod"] = None
            peak_list["c_mod"] = None
            peak_list["ion"] = None
            peak_list["variable_mod"] = None
            peak_list["variable_mod_gain"] = None
            peak_list["variable_mod_loss"] = None
            peak_list["adduct"] = None
            peak_list["loss"] = None
            peak_list["theoretical_mz"] = None
            peak_list["ppm_error"] = None
            peak_list["fit_score"] = None
            peak_list["total_intensity"] = None
            peak_list["charge_scaled_abundance"] = None
            peak_list["fitter_formula"] = None
            peak_list["fitter_theo_x"] = None
            peak_list["fitter_theo_y"] = None
            peak_list["frag_site"] = None

            # sort for binary search
            peak_list = peak_list.sort_values(by="monoisotopic_mw")
            peak_list.to_csv(peak_assignment_file, index=False)


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


    def terminal_matching(self):
        val_window = UnmodSearchValidationWindow(
            self,
            f"precisION - Assisted Terminal Matching [{os.path.basename(self.file_path)}]",
            self.file_path,
            self.directory_path,
            self.name,
            self.sequence,
            self.n_mod_add,
            self.n_mod_sub,
            self.c_mod_add,
            self.c_mod_sub,
            float(self.rough_ppm_input.GetValue()),
            float(self.fine_ppm_input.GetValue()),
            float(self.auto_ppm_input.GetValue()),
            float(self.automatic_matching_score_input.GetValue()),
            self.ion_type_dropdown.GetValue()
        )

        val_window.Show()


    def internal_matching(self):

        if self.ion_type_dropdown.GetValue() == "c/z•":
            print("Internal fragments produced by ExD are not currently supported in precisION.")
            print("Support will be added once mechanistic details are further explored.")
            return

        val_window = InternalSearchValidationWindow(
            self,
            f"precisION - Assisted Internal Matching [{os.path.basename(self.file_path)}]",
            self.file_path,
            self.directory_path,
            self.name,
            self.sequence,
            float(self.fine_ppm_input.GetValue()),
            float(self.auto_ppm_input.GetValue()),
            float(self.automatic_matching_score_input.GetValue()),
            float(self.internal_shared_termini_input.GetValue())
        )

        val_window.Show()
