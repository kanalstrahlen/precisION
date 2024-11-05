import wx
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
from scipy.signal import savgol_filter

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


class PeakDifferencesWindow(wx.Frame):
    def __init__(self, parent, title, file_path, options):
        super().__init__(parent, title=title)
        self.file_path = file_path
        screen_width, screen_height = wx.DisplaySize()
        self.SetSize((screen_width * 2 // 4, screen_height * 2 // 4))

        self.panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.plot_panel = HistogramPanel(self.panel)

        options_sizer = wx.BoxSizer(wx.VERTICAL)

        bin_size_text = wx.StaticText(self.panel, label="Bin size (Da):")
        options_header_font = wx.Font(
            10,
            wx.FONTFAMILY_DEFAULT,
            wx.FONTSTYLE_NORMAL,
            wx.FONTWEIGHT_BOLD
        )
        bin_size_text.SetFont(options_header_font)

        bin_size_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.bin_size_slider = wx.Slider(
            self.panel,
            value=-200,
            minValue=-300,
            maxValue=0.000,
            style=wx.SL_HORIZONTAL
        )
        self.bin_size_input = wx.TextCtrl(
            self.panel,
            value="0.01",
            style=wx.TE_READONLY | wx.NO_BORDER
        )
        self.Bind(wx.EVT_SCROLL, self.on_bin_size_slider_scroll, self.bin_size_slider)

        bin_size_sizer.Add(self.bin_size_slider, 0, wx.RIGHT, 2)
        bin_size_sizer.Add(self.bin_size_input, 0, wx.RIGHT, 5)

        score_mode_text = wx.StaticText(self.panel, label="Scoring mode:")
        score_mode_text.SetFont(options_header_font)
        score_mode_choices = ["Summed intensity", "Intensity × counts", "Counts"]
        self.score_mode_dropdown = wx.ComboBox(
            self.panel,
            choices=score_mode_choices,
            style=wx.CB_READONLY
        )
        self.score_mode_dropdown.SetValue("Counts")

        smoothing_text = wx.StaticText(self.panel, label="Smoothing width:")
        smoothing_text.SetFont(options_header_font)
        smoothing_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.smoothing_slider = wx.Slider(
            self.panel,
            value=1,
            minValue=1,
            maxValue=100,
            style=wx.SL_HORIZONTAL
        )
        self.smoothing_input = wx.TextCtrl(
            self.panel,
            value="1",
            style=wx.TE_READONLY | wx.NO_BORDER
        )
        self.Bind(wx.EVT_SCROLL, self.on_smoothing_slider_scroll, self.smoothing_slider)

        smoothing_sizer.Add(self.smoothing_slider, 0, wx.RIGHT, 2)
        smoothing_sizer.Add(self.smoothing_input, 0, wx.RIGHT, 5)

        options_sizer.Add(bin_size_text, 0, wx.BOTTOM, 4)
        options_sizer.Add(bin_size_sizer, 0, wx.BOTTOM, 10)
        options_sizer.Add(score_mode_text, 0, wx.BOTTOM, 4)
        options_sizer.Add(self.score_mode_dropdown, 0, wx.BOTTOM, 20)
        options_sizer.Add(smoothing_text, 0, wx.BOTTOM, 4)
        options_sizer.Add(smoothing_sizer, 0, wx.BOTTOM, 20)

        main_sizer.Add(self.plot_panel, 1, wx.EXPAND | wx.RIGHT, border=20)
        main_sizer.Add(options_sizer, 0, wx.ALL, border=15)

        icon = wx.Icon('./icons/icon.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        diff_run = DifferenceCalculator(file_path, options)
        self.difference_list = diff_run.calc_diff_list()
        self.plot_panel.bin_values(
            self.difference_list,
            float(self.bin_size_input.GetValue()),
            self.score_mode_dropdown.GetValue(),
            int(self.smoothing_input.GetValue())
        )

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Show()


    def on_size(self, event):
        self.panel.Layout()
        event.Skip()


    def on_bin_size_slider_scroll(self, event):
        self.update_bin_size()
        self.plot_panel.bin_values(
            self.difference_list,
            float(self.bin_size_input.GetValue()),
            self.score_mode_dropdown.GetValue(),
            int(self.smoothing_input.GetValue())
        )


    def update_bin_size(self):
        slider_value = self.bin_size_slider.GetValue()
        value = 10 ** (slider_value / 100)
        self.bin_size_input.SetValue(str(round(value,4)))


    def on_smoothing_slider_scroll(self, event):
        self.update_smoothing()
        self.plot_panel.bin_values(
            self.difference_list,
            float(self.bin_size_input.GetValue()),
            self.score_mode_dropdown.GetValue(),
            int(self.smoothing_input.GetValue())
        )


    def update_smoothing(self):
        slider_value = self.smoothing_slider.GetValue()
        self.smoothing_input.SetValue(str(slider_value))



class DifferenceCalculator():
    def __init__(self, assigned_peaks_file, options):
        peak_list = pd.read_csv(assigned_peaks_file)
        selection_a, selection_b = options

        self.peak_list_a = self.df_to_peak_tuples(peak_list, selection_a)
        self.peak_list_b = self.df_to_peak_tuples(peak_list, selection_b)


    def calc_diff_list(self):
        diff_list = self.calc_difference(self.peak_list_a, self.peak_list_b)
        return diff_list


    def df_to_peak_tuples(self, peak_list, selection):
        if selection == "All peaks":
            peak_tuples = list(
                peak_list[['monoisotopic_mw', 'abundance']].to_records(index=False)
            )
        elif selection == "Assigned peaks":
            temp_peak_list = peak_list.dropna(subset=['name'])
            peak_tuples = list(
                temp_peak_list[['monoisotopic_mw', 'abundance']].to_records(index=False)
            )
        elif selection == "Unassigned peaks":
            temp_peak_list = peak_list[pd.isna(peak_list['name'])]
            peak_tuples = list(
                temp_peak_list[['monoisotopic_mw', 'abundance']].to_records(index=False)
            )

        return peak_tuples


    def calc_difference(self, peak_list_a, peak_list_b):
        diff_tuple_list = []

        for b_tuple in peak_list_b:
            for a_tuple in peak_list_a:
                mass_diff = b_tuple[0] - a_tuple[0]
                intens_product = b_tuple[1] * a_tuple[1]
                if -200 <= mass_diff <= 1500:
                    diff_tuple_list.append((mass_diff, intens_product, 1))

        return diff_tuple_list



class HistogramPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        self.figure = Figure(facecolor='#f0f0f0')
        self.ax = self.figure.add_subplot(111)
        self.ax.set_facecolor('#f0f0f0')
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.SetSizer(sizer)
        self.Layout()
        self.toolbar.Show()
        self.Bind(wx.EVT_SIZE, self.on_size)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.cursor_position_text = self.ax.text(
            0.99,
            0.99,
            '',
            transform=self.ax.transAxes,
            ha='right',
            va='top'
        )


    def bin_values(self, data, bin_size, scoring_method, smooth_window_size):
        mass_differences, intensities, counts = zip(*data)

        mass_differences = np.array(mass_differences)
        intensities = np.array(intensities)
        counts = np.array(counts)

        bin_edges = np.arange(
            min(mass_differences),
            max(mass_differences) + bin_size,
            bin_size
        )
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        if scoring_method == "Intensity × counts":
            hist, _ = np.histogram(mass_differences, bins=bin_edges)
            bin_indices = np.searchsorted(bin_edges, mass_differences, side='right') - 1
            weighted_hist = np.bincount(bin_indices, weights=intensities) * hist

        elif scoring_method == "Summed intensity":
            weighted_hist, _ = np.histogram(
                mass_differences,
                bins=bin_edges,
                weights=intensities
            )

        elif scoring_method == "Counts":
            hist, _ = np.histogram(mass_differences, bins=bin_edges)
            bin_indices = np.searchsorted(bin_edges, mass_differences, side='right') - 1
            weighted_hist = np.bincount(bin_indices, weights=counts) * hist

        try:
            smoothed_hist = savgol_filter(
                weighted_hist,
                window_length=smooth_window_size,
                polyorder=2
            )
        except:
            smoothed_hist = weighted_hist
        smoothed_hist = smoothed_hist / max(smoothed_hist) * 100

        self.plot_spectrum(bin_centers, smoothed_hist)


    def plot_spectrum(self, bin_centers, smoothed_hist):
        self.ax.clear()
        self.ax.plot(bin_centers, smoothed_hist, color='black')
        self.ax.set_ylim(0, max(smoothed_hist) * 1.1)
        self.ax.set_xlim(min(bin_centers), max(bin_centers))
        self.ax.set_ylabel("Relative score (%)")
        self.ax.set_xlabel("Mass difference (Da)")
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


    def on_click(self, event):
        if event.inaxes:
            x, _ = event.xdata, event.ydata
            self.cursor_position_text.remove()
            self.cursor_position_text = self.ax.text(
                0.99,
                0.99,
                f'Peak shift = {x:.5f} Da',
                transform=self.ax.transAxes,
                ha='right',
                va='top'
            )
            self.canvas.draw()
