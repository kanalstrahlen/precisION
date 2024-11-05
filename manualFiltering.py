import os
import wx
import pandas as pd


class FilterWindow(wx.Frame):
    def __init__(self, parent, title, file_path, directory_path):
        super().__init__(parent, title=title, size=(300, 510))

        panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # window title
        title = wx.StaticText(panel, label="Manual Filtering")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))

        # options_sizer: filtering option
        options_sizer = wx.BoxSizer(wx.VERTICAL)

        min_charge_label = wx.StaticText(panel, label="Minimum charge:")
        self.min_charge_input = wx.TextCtrl(panel, value="1")
        max_charge_label = wx.StaticText(panel, label="Maximum charge:")
        self.max_charge_input = wx.TextCtrl(panel, value="30")
        abu_label = wx.StaticText(panel, label="Minimum abundance (%):")
        self.abu_input = wx.TextCtrl(panel, value="0.001")
        fit_label = wx.StaticText(panel, label="Minimum fit score:")
        self.fit_input = wx.TextCtrl(panel, value="0.4")
        interference_label = wx.StaticText(panel, label="Maximum interference score:")
        self.interference_input = wx.TextCtrl(panel, value="0.8")
        s2n_label = wx.StaticText(panel, label="Minimum signal-to-noise ratio:")
        self.s2n_input = wx.TextCtrl(panel, value="4")
        pc_missing_label = wx.StaticText(panel, label="Maximum fraction of missing peaks (%):")
        self.pc_missing_input = wx.TextCtrl(panel, value="90")
        std_mass_label = wx.StaticText(panel, label="Maximum std(m/z error) (ppm):")
        self.std_mass_input = wx.TextCtrl(panel, value="3")

        options_sizer.Add(min_charge_label, 0, wx.ALIGN_LEFT)
        options_sizer.Add(self.min_charge_input, 0, wx.EXPAND | wx.TOP, 5)
        options_sizer.Add(max_charge_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        options_sizer.Add(self.max_charge_input, 0, wx.EXPAND | wx.TOP, 5)
        options_sizer.Add(abu_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        options_sizer.Add(self.abu_input, 0, wx.EXPAND | wx.TOP, 5)
        options_sizer.Add(fit_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        options_sizer.Add(self.fit_input, 0, wx.EXPAND | wx.TOP, 5)
        options_sizer.Add(interference_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        options_sizer.Add(self.interference_input, 0, wx.EXPAND | wx.TOP, 5)
        options_sizer.Add(s2n_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        options_sizer.Add(self.s2n_input, 0, wx.EXPAND | wx.TOP, 5)
        options_sizer.Add(pc_missing_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        options_sizer.Add(self.pc_missing_input, 0, wx.EXPAND | wx.TOP, 5)
        options_sizer.Add(std_mass_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        options_sizer.Add(self.std_mass_input, 0, wx.EXPAND | wx.TOP, 5)

        # filter_button: apply selected filters
        filter_button = wx.Button(panel, label="Apply filters")
        filter_button.Bind(wx.EVT_BUTTON, self.on_filter_button)

        main_sizer.Add(title, 0, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(options_sizer, 0, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(filter_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 7)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        panel.SetSizer(main_sizer)
        self.Show()

        self.directory_path = directory_path
        self.file_path = file_path
        self.basename = os.path.basename(self.file_path).replace(".txt", "")

    def on_filter_button(self, event):
        # read in filter parameters
        min_charge_value = float(self.min_charge_input.GetValue())
        max_charge_value = float(self.max_charge_input.GetValue())
        abu_value = float(self.abu_input.GetValue())
        fit_value = float(self.fit_input.GetValue())
        interference_value = float(self.interference_input.GetValue())
        s2n_value = float(self.s2n_input.GetValue())
        pc_missing_value = float(self.pc_missing_input.GetValue())
        std_mass_value = float(self.std_mass_input.GetValue())

        peaks_file = os.path.join(self.directory_path, f"{self.basename}.fullPeakList.csv")
        peak_list = pd.read_csv(peaks_file)
        peak_list['rel_abundance'] = peak_list["abundance"]/peak_list["abundance"][0] * 100

        # apply filter
        peak_list['prediction'] = (
            (peak_list['charge'] >= min_charge_value) &
            (peak_list['charge'] <= max_charge_value) &
            (peak_list['rel_abundance'] >= abu_value) &
            (peak_list['fit_score'] >= fit_value) &
            (peak_list['interference'] <= interference_value) &
            (peak_list['s2n'] >= s2n_value) &
            (peak_list['pc_missing_peaks'] <= pc_missing_value) &
            (peak_list['mass_error_std'] <= std_mass_value)
        )

        output_file = os.path.join(self.directory_path, f"{self.basename}.filteredPeakList.csv")
        peak_list.to_csv(output_file, index=False)
        print(f"Saved filtered peak list to: {output_file}")
        self.Close()
