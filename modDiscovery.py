import os
import wx

from peakDifferences import PeakDifferencesWindow
from offsetSearch import OffsetSearchWindow
from offsetID import OffsetIDWindow



class ModDiscoveryWindow(wx.Frame):
    def __init__(self, parent, title, file_path, directory_path, sequence):
        super().__init__(
            parent,
            title=title,
            size=(460, 750),
            style=wx.DEFAULT_FRAME_STYLE & ~(wx.RESIZE_BORDER | wx.MAXIMIZE_BOX)
        )

        column_width = 350
        subheader_font = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_BOLD)

        assigned_peaks_file = os.path.basename(file_path).replace(".txt", ".assignedPeakList.csv")
        self.directory_path = directory_path
        self.assigned_peaks_file_path = os.path.join(directory_path, assigned_peaks_file)

        self.panel = wx.Panel(self)
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # window title
        title = wx.StaticText(self.panel, label="Modification Discovery")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext1 = wx.StaticText(
            self.panel,
            label="Search for modifications in top-down mass spectra using untargeted approaches."
        )
        subtext2 = wx.StaticText(
            self.panel,
            label=(
                "Identified modifications can be incorporated "
                "into a variable search for assignment."
            )
        )

        # search_options_sizer: options for all three search modes
        search_options_sizer = wx.BoxSizer(wx.VERTICAL)

        # peak_diff_sizer: which peaks to compare?
        peak_diff_staticbox = wx.StaticBox(self.panel)
        peak_diff_sizer = wx.StaticBoxSizer(peak_diff_staticbox, wx.VERTICAL)
        peak_diff_header = wx.StaticText(self.panel, label="Offset finder")
        peak_diff_header.SetFont(subheader_font)
        peak_diff_info = wx.StaticText(
            self.panel,
            label="Identify common mass differences between fragment ions."
        )

        peak_diff_selection_choices= ["All peaks", "Assigned peaks", "Unassigned peaks"]

        peak_diff_selection_a_sizer = wx.BoxSizer(wx.HORIZONTAL)
        peak_diff_selection_a_text = wx.StaticText(self.panel, label="Peak list A:")
        self.peak_diff_selection_a_dropdown = wx.ComboBox(
            self.panel,
            choices=peak_diff_selection_choices,
            style=wx.CB_READONLY
        )
        self.peak_diff_selection_a_dropdown.SetValue("Assigned peaks")

        peak_diff_selection_a_sizer.Add(
            peak_diff_selection_a_text,
            0,
            wx.ALIGN_CENTER_VERTICAL | wx.RIGHT,
            2
        )
        peak_diff_selection_a_sizer.Add(
            self.peak_diff_selection_a_dropdown,
            0,
            wx.ALIGN_CENTER_VERTICAL | wx.RIGHT,
            5
        )

        peak_diff_selection_b_sizer = wx.BoxSizer(wx.HORIZONTAL)
        peak_diff_selection_b_text = wx.StaticText(self.panel, label="Peak list B:")
        self.peak_diff_selection_b_dropdown = wx.ComboBox(
            self.panel,
            choices=peak_diff_selection_choices,
            style=wx.CB_READONLY
        )
        self.peak_diff_selection_b_dropdown.SetValue("Unassigned peaks")

        peak_diff_selection_b_sizer.Add(
            peak_diff_selection_b_text,
            0,
            wx.ALIGN_CENTER_VERTICAL | wx.RIGHT,
            2
        )
        peak_diff_selection_b_sizer.Add(
            self.peak_diff_selection_b_dropdown,
            0,
            wx.ALIGN_CENTER_VERTICAL | wx.RIGHT,
            5
        )

        peak_diff_button = wx.Button(self.panel, label="Compute Differences (B-A)",  size=(180, 25))
        peak_diff_button.Bind(wx.EVT_BUTTON, self.on_peak_diff_button)

        peak_diff_sizer.Add(peak_diff_header, 0, wx.BOTTOM, 5)
        peak_diff_sizer.Add(peak_diff_info, 0, wx.BOTTOM, 5)
        peak_diff_sizer.Add(peak_diff_selection_a_sizer, 0, wx.BOTTOM, 5)
        peak_diff_sizer.Add(peak_diff_selection_b_sizer, 0, wx.BOTTOM, 10)
        peak_diff_sizer.Add(peak_diff_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM, 5)

        peak_diff_sizer.SetMinSize((column_width, -1))

        # offset_search_sizer: offset range, b_ions, y_ions or both
        offset_search_staticbox = wx.StaticBox(self.panel)
        offset_search_sizer = wx.StaticBoxSizer(offset_search_staticbox, wx.VERTICAL)
        offset_search_header = wx.StaticText(self.panel, label="Fragment-level open search")
        offset_search_header.SetFont(subheader_font)
        offset_search_info = wx.StaticText(
            self.panel,
            label=(
                "Identify statistically significant mass offsets that "
                "correspond with modifications."
            )
        )

        offset_range_sizer = wx.BoxSizer(wx.HORIZONTAL)
        offset_range_text = wx.StaticText(self.panel, label="Offset range (Da):")
        self.offset_range_lower = wx.TextCtrl(self.panel, size=(40, 20))
        self.offset_range_lower.SetValue("-100")
        offset_range_text2 = wx.StaticText(self.panel, label="-")
        self.offset_range_upper = wx.TextCtrl(self.panel, size=(40, 20))
        self.offset_range_upper.SetValue("500")

        offset_range_sizer.Add(offset_range_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        offset_range_sizer.Add(self.offset_range_lower, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 2)
        offset_range_sizer.Add(offset_range_text2, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 2)
        offset_range_sizer.Add(self.offset_range_upper, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)

        offset_tol_sizer = wx.BoxSizer(wx.HORIZONTAL)
        offset_tol_text = wx.StaticText(self.panel, label="Tolerance (ppm):")
        self.offset_tol_input = wx.TextCtrl(self.panel, size=(40, 20))
        self.offset_tol_input.SetValue("3")
        count_mode_text = wx.StaticText(self.panel, label="Counting mode:")
        count_mode_choices = ["Theoretical ions", "Observed ions"]
        self.count_mode_dropdown = wx.ComboBox(
            self.panel,
            choices=count_mode_choices,
            style=wx.CB_READONLY
        )
        self.count_mode_dropdown.SetValue(count_mode_choices[0])

        offset_tol_sizer.Add(offset_tol_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 2)
        offset_tol_sizer.Add(self.offset_tol_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        offset_tol_sizer.Add(count_mode_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 2)
        offset_tol_sizer.Add(self.count_mode_dropdown, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)

        offset_search_button = wx.Button(self.panel, label="Run Offset Search",  size=(180, 25))
        offset_search_button.Bind(wx.EVT_BUTTON, self.on_offset_search_button)

        offset_search_sizer.Add(offset_search_header, 0, wx.BOTTOM, 5)
        offset_search_sizer.Add(offset_search_info, 0, wx.BOTTOM, 5)
        offset_search_sizer.Add(offset_range_sizer, 0, wx.BOTTOM, 5)
        offset_search_sizer.Add(offset_tol_sizer, 0, wx.BOTTOM, 10)
        offset_search_sizer.Add(offset_search_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM, 5)

        offset_search_sizer.SetMinSize((column_width, -1))

        # offset_id: offset, view apply to both termini, search unimod
        offset_id_staticbox = wx.StaticBox(self.panel)
        offset_id_sizer = wx.StaticBoxSizer(offset_id_staticbox, wx.VERTICAL)
        offset_id_header = wx.StaticText(self.panel, label="Offset identification")
        offset_id_header.SetFont(subheader_font)
        offset_id_info = wx.StaticText(
            self.panel,
            label="Localise mass offsets and search UniMod to identify composition."
        )

        offset_id_offset_sizer = wx.BoxSizer(wx.HORIZONTAL)
        offset_id_offset_text = wx.StaticText(self.panel, label="Offset (Da):")
        self.offset_id_offset_input = wx.TextCtrl(self.panel, size=(40, 20))
        self.offset_id_offset_input.SetValue("0")

        offset_id_offset_sizer.Add(offset_id_offset_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 2)
        offset_id_offset_sizer.Add(
            self.offset_id_offset_input,
            0,
            wx.ALIGN_CENTER_VERTICAL | wx.RIGHT,
            5
        )

        offset_id_tol_sizer = wx.BoxSizer(wx.HORIZONTAL)
        offset_id_tol_text = wx.StaticText(self.panel, label="Tolerance (ppm):")
        self.offset_id_tol_input = wx.TextCtrl(self.panel, size=(40, 20))
        self.offset_id_tol_input.SetValue("3")

        offset_id_tol_sizer.Add(offset_id_tol_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 2)
        offset_id_tol_sizer.Add(self.offset_id_tol_input, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)

        offset_id_button = wx.Button(self.panel, label="Identify Offset",  size=(180, 25))
        offset_id_button.Bind(wx.EVT_BUTTON, self.on_offset_id_button)

        offset_id_sizer.Add(offset_id_header, 0, wx.BOTTOM, 5)
        offset_id_sizer.Add(offset_id_info, 0, wx.BOTTOM, 5)
        offset_id_sizer.Add(offset_id_offset_sizer, 0, wx.BOTTOM, 5)
        offset_id_sizer.Add(offset_id_tol_sizer, 0, wx.BOTTOM, 10)
        offset_id_sizer.Add(offset_id_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM, 5)

        offset_id_sizer.SetMinSize((column_width, -1))

        search_options_sizer.Add(peak_diff_sizer, 0, wx.EXPAND | wx.BOTTOM, 0)
        search_options_sizer.Add(offset_search_sizer, 0, wx.EXPAND | wx.BOTTOM, 0)
        search_options_sizer.Add(offset_id_sizer, 0, wx.EXPAND | wx.BOTTOM, 0)

        # sequence_sizer: input proteoform sequence
        sequence_staticbox = wx.StaticBox(self.panel)
        sequence_sizer = wx.StaticBoxSizer(sequence_staticbox, wx.VERTICAL)
        sequence_header = wx.StaticText(self.panel, label="Proteoform sequence")
        sequence_header.SetFont(subheader_font)
        self.sequence_text = wx.TextCtrl(
            self.panel,
            style=wx.TE_MULTILINE | wx.TE_CHARWRAP | wx.VSCROLL,
            size=(420, 150)
        )

        self.sequence_text.SetValue(sequence)

        sequence_sizer.Add(sequence_header, 0, wx.BOTTOM, 5)
        sequence_sizer.Add(self.sequence_text, 0, wx.BOTTOM, 5)

        main_sizer.Add(title, 0, wx.EXPAND | wx.BOTTOM, 5)
        main_sizer.Add(subtext1, 0, wx.EXPAND | wx.BOTTOM, 3)
        main_sizer.Add(subtext2, 0, wx.EXPAND | wx.BOTTOM, 0)
        main_sizer.Add(search_options_sizer, 0, wx.EXPAND | wx.BOTTOM, 0)
        main_sizer.Add(sequence_sizer, 0, wx.EXPAND | wx.BOTTOM, 0)

        panel_sizer.Add(main_sizer, 0, wx.ALL, 5)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.panel.SetSizer(panel_sizer)
        self.Show()


    def on_peak_diff_button(self, event):
        options = [
            self.peak_diff_selection_a_dropdown.GetValue(),
            self.peak_diff_selection_b_dropdown.GetValue(),
        ]

        peak_diff_window = PeakDifferencesWindow(
            self,
            "precisION - Offset Finder",
            self.assigned_peaks_file_path,
            options
            )

        peak_diff_window.Show()


    def on_offset_search_button(self, event):
        options = [
            float(self.offset_range_lower.GetValue()),
            float(self.offset_range_upper.GetValue()),
            float(self.offset_tol_input.GetValue()),
            self.count_mode_dropdown.GetValue()
        ]

        offset_search_window = OffsetSearchWindow(
            self,
            "precisION - Fragment-Level Open Search",
            self.assigned_peaks_file_path,
            self.sequence_text.GetValue(),
            options
            )

        offset_search_window.Show()


    def on_offset_id_button(self, event):
        options = [
            float(self.offset_id_offset_input.GetValue()),
            float(self.offset_id_tol_input.GetValue()),
        ]

        offset_id_window = OffsetIDWindow(
            self,
            "precisION - Offset Identification",
            self.assigned_peaks_file_path,
            self.sequence_text.GetValue(),
            options
            )

        offset_id_window.Show()
