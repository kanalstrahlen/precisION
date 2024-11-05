import os
import ast
import numpy as np
import pandas as pd
import wx
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec

from classifier import PeakListClassifier

class ValidationWindow(wx.Frame):
    def __init__(self, parent, title, file_path, directory_path, num_shuffle):
        super().__init__(parent, title=title, size=(1500, 800))
        self.num_shuffle = num_shuffle

        self.min_mz = 0
        self.max_mz = 1
        self.abu = 0
        self.spectrum = []

        self.file_path = file_path
        self.directory_path = directory_path
        self.basename = os.path.basename(file_path).replace(".txt", "")
        self.output_file = os.path.join(
            self.directory_path,
            f"{self.basename}.filteredPeakList.csv"
        )

        #screen width variable to allow for smaller screens
        screen_width, screen_height = wx.DisplaySize()
        desired_width = screen_width * 2 // 4
        self.SetSize((desired_width, screen_height * 3 // 4))

        self.panel = wx.Panel(self)
        self.plot_panel = ValidationPlotsPanel(self.panel)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # instructions for use
        instruction_text = wx.StaticText(
            self.panel,
            label="Enter 't' for true positive, 'f'' for false positive, or 'x' to attempt classification:"
        )
        instruction_text.SetFont(
            wx.Font(
                12,
                wx.FONTFAMILY_DEFAULT,
                wx.FONTSTYLE_NORMAL,
                wx.FONTWEIGHT_BOLD
            )
        )

        # hidden buttons for screening peaks and keyboard shortcuts using accelerator table
        self.true_button = wx.Button(self.panel)
        self.true_button.Bind(wx.EVT_BUTTON, self.on_true_button)
        self.true_button.Hide()
        self.false_button = wx.Button(self.panel)
        self.false_button.Bind(wx.EVT_BUTTON, self.on_false_button)
        self.false_button.Hide()
        self.finished_button = wx.Button(self.panel)
        self.finished_button.Bind(wx.EVT_BUTTON, self.on_finished_button)
        self.finished_button.Hide()
        accelerator_entries = [(wx.ACCEL_NORMAL, ord("t"), self.true_button.GetId()),
                               (wx.ACCEL_NORMAL, ord("f"), self.false_button.GetId()),
                               (wx.ACCEL_NORMAL, ord("x"), self.finished_button.GetId())]
        accelerator_table = wx.AcceleratorTable(accelerator_entries)
        self.SetAcceleratorTable(accelerator_table)

        main_sizer.Add(instruction_text, 0, wx.EXPAND | wx.ALL, border=10)
        main_sizer.Add(self.plot_panel, 1, wx.EXPAND | wx.TOP, border=-10)

        self.setup_plots()

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)
        self.plot_peaks(self.i)

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Show()



    def setup_plots(self):
        #load peak list and spectrum
        self.load_spectrum(self.directory_path)
        peaks_file = os.path.join(self.directory_path, f"{self.basename}.fullPeakList.csv")
        peak_list = pd.read_csv(peaks_file)
        self.peak_list = peak_list
        self.peak_list["checked"] = False
        self.peak_list["IndexColumn"] = self.peak_list.index

        #shuffle top 100 entries and rejoin to ensure a mix of correct and not
        num_rows_to_shuffle = self.num_shuffle
        top_rows = self.peak_list.head(num_rows_to_shuffle)
        shuffled_top_rows = top_rows.sample(frac=1)
        shuffled_top_rows = shuffled_top_rows.set_index(top_rows.index)
        self.peak_list = pd.concat([shuffled_top_rows, self.peak_list.iloc[num_rows_to_shuffle:]])

        # spectrum counter; incremented by one after each peak screened
        self.i = 0

        # initialise plots
        # subplot_1: number of true and false entries
        self.true_peaks = 0
        self.false_peaks = 0
        self.plot_panel.plot_subplot_1(self.true_peaks, self.false_peaks)

        # subplot_2: abundance vs s2n scatter plot with true/false colouring
        self.plot_panel.ax_subplots[1].scatter(
            self.peak_list["interference"],
            self.peak_list["log_abundance"],
            color="grey",
            alpha=0.5
        )
        self.plot_panel.ax_subplots[1].set_ylabel("log10(Abundance)")
        self.plot_panel.ax_subplots[1].set_xlabel("Interference score")
        self.plot_panel.ax_subplots[1].xaxis.set_major_locator(plt.AutoLocator())

        # subplot_3: true/false for each deconvolution result
        self.thrash_1_true = 0
        self.thrash_1_false = 0
        self.thrash_2_true = 0
        self.thrash_2_false = 0
        self.thrash_3_true = 0
        self.thrash_3_false = 0
        self.thrash_4_true = 0
        self.thrash_4_false = 0
        self.msdeconv_true = 0
        self.msdeconv_false = 0
        self.envcnn_true = 0
        self.envcnn_false = 0
        categories = [
            "Thrash 0.1",
            "Thrash 0.2",
            "Thrash 0.3",
            "Thrash 0.4",
            "MS-Deconv",
            "EnvCNN"
        ]
        self.plot_panel.ax_subplots[2].bar(categories, [0,0,0,0,0,0])
        self.plot_panel.ax_subplots[2].set_xticks(range(len(categories)))
        self.plot_panel.ax_subplots[2].set_xticklabels(categories, rotation=90)
        self.plot_panel.ax_subplots[2].set_ylabel("Number of peaks")

        # subplot_4: comparative frequency histogram (fit vs count)
        self.plot_panel.ax_subplots[3].set_xlabel("Fit Score")
        self.plot_panel.ax_subplots[3].set_ylabel("Frequency")
        self.plot_panel.ax_subplots[3].set_xlim((-0.015,0.415))
        self.true_fits = []
        self.false_fits = []


    def on_size(self, event):
        size = self.panel.GetClientSize()
        self.plot_panel.SetSize(size.width // 2, size.height)
        self.panel.Layout()
        event.Skip()


    def load_spectrum(self, directory_path):
        #load spectrum from the .spectrum file
        save_path = os.path.join(directory_path, f"{self.basename}.spectrum.txt")
        self.spectrum = np.loadtxt(save_path)
        self.min_mz = min(self.spectrum[:, 0])
        self.max_mz = max(self.spectrum[:, 0])
        self.abu = max(self.spectrum[:,1])


    def plot_peaks(self, i):
        # read in parameters to send to the plot
        # mass and charge
        mass = self.peak_list.loc[i, "monoisotopic_mw"]
        charge = self.peak_list.loc[i, "charge"]
        mass_diff = self.peak_list.loc[i, "mass_diff"]
        mass_diff = ast.literal_eval(mass_diff)

        # theoretical isotopic envelope
        scatter_x = self.peak_list.loc[i, "theo_x"]
        scatter_x = ast.literal_eval(scatter_x)
        scatter_y = self.peak_list.loc[i, "theo_y"]
        scatter_y = ast.literal_eval(scatter_y)

        # mass and intensity ranges
        abu = self.peak_list.loc[i, "abundance"]
        selected_min_mz = min(scatter_x) - 3 * (scatter_x[1] - scatter_x[0])
        selected_max_mz = max(scatter_x) + 4 * (scatter_x[1] - scatter_x[0])

        # color options
        color_list = self.peak_list.loc[i, "matched_bool"]
        color_list = ast.literal_eval(color_list)
        matched_mz = self.peak_list.loc[i, "matched_x"]
        matched_mz = ast.literal_eval(matched_mz)
        matched_intens = self.peak_list.loc[i, "matched_y"]
        matched_intens = ast.literal_eval(matched_intens)
        color_mapping = {0: "#cd4120", 1: "#6ab4e6"}
        color_list = [color_mapping[item] for item in color_list]
        label_color = "#b44a52"

        self.plot_panel.ax_main.clear()
        self.plot_panel.plot_main_validation_spectrum(
            self.spectrum,
            selected_min_mz,
            selected_max_mz,
            abu,
            color_list,
            mass,
            charge,
            mass_diff,
            matched_mz,
            matched_intens,
            label_color,
            scatter_x,
            scatter_y
        )
        self.plot_panel.canvas.draw()


    def on_true_button(self, event):
        if self.i < len(self.peak_list):

            # update peak_list
            self.peak_list.loc[self.i, "checked"] = True
            self.peak_list.loc[self.i, "validated"] = True

            # update subplot_1 data
            self.true_peaks += 1

            # update subplot_3 data
            keys = [
                "detected_thrash_0_1?",
                "detected_thrash_0_2?",
                "detected_thrash_0_3?",
                "detected_thrash_0_4?"
            ]
            for thrash_key in keys:
                if self.peak_list.loc[self.i, thrash_key]:
                    thrash_number = int(thrash_key.split("_")[-1][0])
                    setattr(
                        self,
                        f"thrash_{thrash_number}_true",
                        getattr(self, f"thrash_{thrash_number}_true") + 1
                    )
            if self.peak_list.loc[self.i, "detected_msdeconv?"]:
                self.msdeconv_true += 1
            if self.peak_list.loc[self.i, "detected_envcnn?"]:
                self.envcnn_true += 1

            # update subplot_4 data
            self.true_fits.append(self.peak_list.loc[self.i, "fit_score"])

            # update all plots
            self.plot_panel.plot_subplot_1(self.true_peaks, self.false_peaks)
            self.plot_panel.plot_subplot_2(
                self.peak_list.loc[self.i, "interference"],
                self.peak_list.loc[self.i, "log_abundance"],
                True
            )
            self.plot_panel.plot_subplot_3(
                self.thrash_1_true,
                self.thrash_1_false,
                self.thrash_2_true,
                self.thrash_2_false,
                self.thrash_3_true,
                self.thrash_3_false,
                self.thrash_4_true,
                self.thrash_4_false,
                self.msdeconv_true,
                self.msdeconv_false,
                self.envcnn_true,
                self.envcnn_false
            )
            self.plot_panel.plot_subplot_4(self.true_fits, self.false_fits)

            # move to the next peak and update plot
            self.i += 1
            self.plot_peaks(self.i)

        else:
            self.on_finishing_classification()


    def on_false_button(self, event):
        if self.i < len(self.peak_list):

            # update peak_list
            self.peak_list.loc[self.i, "checked"] = True
            self.peak_list.loc[self.i, "validated"] = False

            # update subplot_1 data
            self.false_peaks += 1

            # update subplot_3 data
            keys = [
                "detected_thrash_0_1?",
                "detected_thrash_0_2?",
                "detected_thrash_0_3?",
                "detected_thrash_0_4?"
            ]
            for thrash_key in keys:
                if self.peak_list.loc[self.i, thrash_key]:
                    thrash_number = int(thrash_key.split("_")[-1][0])
                    setattr(
                        self,
                        f"thrash_{thrash_number}_false",
                        getattr(self, f"thrash_{thrash_number}_false") + 1
                    )
            if self.peak_list.loc[self.i, "detected_msdeconv?"]:
                self.msdeconv_false += 1
            if self.peak_list.loc[self.i, "detected_envcnn?"]:
                self.envcnn_false += 1

            # update subplot_4 data
            self.false_fits.append(self.peak_list.loc[self.i, "fit_score"])

            # update all plots
            self.plot_panel.plot_subplot_1(self.true_peaks, self.false_peaks)
            self.plot_panel.plot_subplot_2(
                self.peak_list.loc[self.i, "interference"],
                self.peak_list.loc[self.i, "log_abundance"],
                False
            )
            self.plot_panel.plot_subplot_3(
                self.thrash_1_true,
                self.thrash_1_false,
                self.thrash_2_true,
                self.thrash_2_false,
                self.thrash_3_true,
                self.thrash_3_false,
                self.thrash_4_true,
                self.thrash_4_false,
                self.msdeconv_true,
                self.msdeconv_false,
                self.envcnn_true,
                self.envcnn_false
            )
            self.plot_panel.plot_subplot_4(self.true_fits, self.false_fits)

            # move to the next peak and update plot
            self.i += 1
            self.plot_peaks(self.i)

        else:
            print("Classification complete.")
            print("Please press 't'.")


    def on_finished_button(self, event):
        print("Starting classification...")
        self.peak_list.sort_values(by="IndexColumn").to_csv(self.output_file, index=False)
        PeakListClassifier(self.output_file)
        self.Close()


    def on_finishing_classification(self):
        print("Finished classification")
        self.peak_list["prediction"] = self.peak_list["validated"]
        self.peak_list.sort_values(by="IndexColumn").to_csv(self.output_file, index=False)
        print(f"Saved filtered peak list as: {self.output_file}")
        self.Close()



class ValidationPlotsPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        self.max_intens = 0
        sizer = wx.BoxSizer(wx.VERTICAL)

        self.figure = Figure(facecolor="#f0f0f0")
        self.gs = gridspec.GridSpec(2, 4, height_ratios=[2, 1])
        self.ax_main = self.figure.add_subplot(self.gs[0, :])
        self.ax_subplots = [
            self.figure.add_subplot(self.gs[1, i], facecolor="#f0f0f0") for i in range(4)
        ]
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.toolbar = NavigationToolbar2Wx(self.canvas)

        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.EXPAND)

        self.SetSizer(sizer)
        self.Bind(wx.EVT_SIZE, self.on_size)
        self.Layout()
        self.toolbar.Realize()




    def plot_main_validation_spectrum(
        self, spectrum, selected_min_mz, selected_max_mz, abu,
        color_list, mass, charge, mass_diff, matched_mz,
        matched_intens, label_color, scatter_x=None, scatter_y=None
    ):
        self.ax_main.clear()
        self.ax_main.plot(spectrum[:, 0], spectrum[:, 1], color="black", linewidth=1)
        self.ax_main.set_xlabel("m/z")
        self.ax_main.set_ylabel("Intensity")
        self.ax_main.set_facecolor("#f0f0f0")
        self.ax_main.set_xlim(selected_min_mz, selected_max_mz)
        self.max_intens = abu * 1.1
        self.ax_main.set_ylim(0, self.max_intens)

        if scatter_x is not None:
            # peak label at top right
            label_text = f"{mass} ({charge}+)"
            self.ax_main.annotate(
                label_text,
                xy=(0.995,0.99),
                xycoords="axes fraction",
                ha="right",
                va="top"
            )

            # scatter plot w/ theoretical distribution
            self.ax_main.scatter(
                scatter_x,
                scatter_y,
                c=color_list,
                marker="o",
                alpha=0.8
            )

            # vertical line at each peak
            for peak in zip(matched_mz, matched_intens):
                x_coords = [peak[0], peak[0]]
                y_coords = [0, peak[1]]
                self.ax_main.plot(
                    x_coords,
                    y_coords,
                    color=label_color,
                    linestyle="--",
                    linewidth=1
                )

            # mass differences between adjacent peaks
            for i, diff in enumerate(mass_diff):
                diff = round(float(diff), 3)
                diff = f"{diff:.3f}"
                x_pos = (matched_mz[i] + matched_mz[i+1])/2
                y_pos = min([matched_intens[i], matched_intens[i+1]])
                if y_pos > (self.max_intens*0.01):
                    x_coords = [matched_mz[i], matched_mz[i+1]]
                    y_coords = [y_pos, y_pos]
                    self.ax_main.plot(
                        x_coords,
                        y_coords,
                        color=label_color,
                        linestyle=":",
                        linewidth=1
                    )
                    self.ax_main.text(
                        x_pos,
                        y_pos,
                        diff,
                        ha="center",
                        va="bottom",
                        color=label_color
                    )

        self.canvas.draw()


    def plot_subplot_1(self, true_peaks, false_peaks):
        ax = self.ax_subplots[0]
        ax.clear()
        ax.set_ylabel("Number of peaks")
        x_labels = ["True", "False"]
        y_values = [true_peaks, false_peaks]
        ax.bar(x_labels, y_values, color=["#6ab4e6", "#cd4120"], alpha=0.8)
        self.canvas.draw()


    def plot_subplot_2(self, x, y, selection):
        ax = self.ax_subplots[1]
        if selection == True:
            ax.scatter(x, y, color="#84c2ea")
        elif selection == False:
            ax.scatter(x, y, color="#d4644c")
        self.canvas.draw()


    def plot_subplot_3(
        self, thrash_1_true, thrash_1_false, thrash_2_true,
        thrash_2_false, thrash_3_true, thrash_3_false,
        thrash_4_true, thrash_4_false, msdeconv_true, msdeconv_false,
        envcnn_true, envcnn_false
    ):
        ax = self.ax_subplots[2]
        ax.clear()
        ax.set_ylabel("Number of peaks")
        categories = [
            "Thrash 0.1",
            "Thrash 0.2",
            "Thrash 0.3",
            "Thrash 0.4",
            "MS-Deconv",
            "EnvCNN"
        ]
        true = [
            thrash_1_true,
            thrash_2_true,
            thrash_3_true,
            thrash_4_true,
            msdeconv_true,
            envcnn_true
        ]
        false = [
            thrash_1_false,
            thrash_2_false,
            thrash_3_false,
            thrash_4_false,
            msdeconv_false,
            envcnn_false
        ]
        ax.bar(categories, true, color="#6ab4e6", alpha=0.8)
        ax.bar(categories, false, bottom=true, color="#cd4120", alpha=0.8)
        ax.set_xticks(range(len(categories)))
        ax.set_xticklabels(categories, rotation=90)
        self.canvas.draw()


    def plot_subplot_4(self, true_fits, false_fits):
        ax = self.ax_subplots[3]
        ax.clear()
        bins = [
            0,
            0.05,
            0.1,
            0.15,
            0.2,
            0.25,
            0.3,
            0.35,
            0.4,
            0.45,
            0.5,
            0.55,
            0.6,
            0.65,
            0.7,
            0.75,
            0.8,
            0.85,
            0.9,
            0.95,
            1.0
        ]
        ax.set_ylabel("Frequency")
        ax.set_xlabel("Fit Score")
        weight_list = [-1 for _ in false_fits] # make the false fit histogram negative
        ax.hist(
            true_fits,
            bins=bins,
            color="#6ab4e6",
            alpha=0.8,
            label="True"
        )
        ax.hist(
            false_fits,
            bins=bins,
            color="#cd4120",
            weights=weight_list,
            alpha=0.8,
            label="False"
        )
        lim = max(ax.get_ylim(), key=abs)
        ax.set_ylim((-abs(lim), abs(lim)))
        ax.axhline(y=0, color="black", linestyle="dashed", linewidth=0.5)
        ax.legend(facecolor="#f0f0f0")
        self.canvas.draw()


    def on_size(self, event):
        self.fit_plot_to_panel()
        event.Skip()


    def fit_plot_to_panel(self):
        size = self.GetClientSize()
        if size[0] >= 50:
            dpi = self.GetContentScaleFactor() * 100
            width = size.width / dpi
            height = size.height / dpi - 0.3
            self.figure.set_size_inches(width, height)
            self.figure.tight_layout(rect=[0, 0, 1, 1])
            self.canvas.draw()
