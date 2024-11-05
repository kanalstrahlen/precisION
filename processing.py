import os
import ast
import wx
import wx.grid
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
import joblib
from sklearn.preprocessing import StandardScaler

# precisION modules
from validation import ValidationWindow
from manualFiltering import FilterWindow

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

class ProcessingWindow(wx.Frame):
    def __init__(self, parent, title, file_path, directory_path):
        super().__init__(parent, title=title)

        self.file_path = file_path
        self.directory_path = directory_path
        self.basename = os.path.basename(file_path).replace(".txt", "")

        #screen width variable to allow for smaller screens
        screen_width, _ = wx.DisplaySize()
        desired_width = screen_width * 3 // 4
        self.SetSize((desired_width, 800))

        self.panel = wx.Panel(self)
        self.plot_panel = ProcessingSpectrumPanel(self.panel)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # data_sizer: title, buttons and grid
        self.data_sizer = wx.BoxSizer(wx.VERTICAL)

        # menu_sizer: title and buttons
        menu_sizer = wx.BoxSizer(wx.VERTICAL)
        title = wx.StaticText(self.panel, label="Envelope Filtering")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext1 = wx.StaticText(
            self.panel,
            label="View and filter deconvolution results using "
            "supervised classification or manual filtering."
        )

        # button_sizer: buttons for user input
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        load_decon_button = wx.Button(self.panel, label="Load full deconvolution")
        load_decon_button.Bind(wx.EVT_BUTTON, self.on_load_decon_button)
        val_button = wx.Button(self.panel, label="Supervised classification")
        val_button.Bind(wx.EVT_BUTTON, self.on_val_button)
        shuffle_peak_text = wx.StaticText(self.panel, label="Num. shuffled envelopes for training:")
        self.shuffle_peak_input = wx.TextCtrl(self.panel, value="1000",  size=(40, -1))
        existing_classifier_button = wx.Button(self.panel, label="Classify using exisiting model")
        existing_classifier_button.Bind(wx.EVT_BUTTON, self.on_existing_classifier_button)
        filter_button = wx.Button(self.panel, label="Apply manual filters")
        filter_button.Bind(wx.EVT_BUTTON, self.on_filter_button)
        load_validated_decon_button = wx.Button(self.panel, label="Load filtered deconvolution")
        load_validated_decon_button.Bind(wx.EVT_BUTTON, self.on_load_validated_decon_button)

        button_sizer.Add(load_decon_button, 0, wx.ALL, 2)
        button_sizer.Add(val_button, 0, wx.ALL, 2)
        button_sizer.Add(shuffle_peak_text, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 2)
        button_sizer.Add(self.shuffle_peak_input, 0, wx.ALL, 2)
        button_sizer.Add(existing_classifier_button, 0, wx.ALL, 2)
        button_sizer.Add(filter_button, 0, wx.ALL, 2)
        button_sizer.Add(load_validated_decon_button, 0, wx.ALL, 2)

        menu_sizer.Add(title, 0, wx.BOTTOM, 5)
        menu_sizer.Add(subtext1, 0, wx.BOTTOM, 5)
        menu_sizer.Add(button_sizer, 0, wx.BOTTOM, 5)

        # data_grid: grid for loading in deconvolution results, interactive
        self.data_grid = wx.grid.Grid(self.panel)
        self.data_grid.SetDefaultCellBackgroundColour("#f0f0f0")
        self.data_grid.Bind(wx.grid.EVT_GRID_SELECT_CELL, self.on_cell_select)

        self.data_sizer.Add(menu_sizer, 0)
        self.data_sizer.Add(self.data_grid, 1)
        self.data_sizer.SetItemMinSize(self.data_grid, -1, 220)

        main_sizer.Add(self.data_sizer, 0, wx.EXPAND | wx.ALL, border=10)
        main_sizer.Add(self.plot_panel, 1, wx.EXPAND | wx.TOP, border=-10)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.load_empty_spectrum(file_path, directory_path)
        self.Show()

        self.peak_list = pd.DataFrame()


    def on_size(self, event):
        size = self.panel.GetClientSize()
        self.plot_panel.SetSize(size.width // 2, size.height)
        self.panel.Layout()
        event.Skip()


    def load_empty_spectrum(self, file_path, directory_path):
        save_path = os.path.join(directory_path, f"{self.basename}.spectrum.txt")

        # make a spectrum file inside the directory
        if not os.path.exists(save_path):
            rl_file_path = os.path.join(directory_path, f"{self.basename}.profile.rl.txt")
            if os.path.exists(rl_file_path):
                print('Loading spectrum...')
                output = np.loadtxt(rl_file_path)
                output = output[::2]
                save_path = os.path.join(directory_path, f"{self.basename}.spectrum.txt")
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
            self.abu,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None
        )



    def on_filter_button(self, event):
        try:
            filter_window = FilterWindow(
                self,
                f"precisION - Manual Filtering [{os.path.basename(self.file_path)}]",
                self.file_path,
                self.directory_path
            )
            filter_window.Show()
        except FileNotFoundError:
            print("No deconvolution file is present in the working directory!")


    def on_existing_classifier_button(self, event):
        dlg = wx.FileDialog(self, "Choose a .pk1 file for logistic regression")
        if dlg.ShowModal() == wx.ID_OK:
            lr_path = dlg.GetPath()
            print("Selected file:", lr_path)
        dlg.Destroy()
        dlg = wx.FileDialog(self, "Choose a .pk1 file for gradient boosting classification")
        if dlg.ShowModal() == wx.ID_OK:
            gb_path = dlg.GetPath()
            print("Selected file:", gb_path)
        dlg.Destroy()
        self.classify_from_model(lr_path, gb_path)

    def classify_from_model(self, lr_file, gb_file):
        peaks_file = os.path.join(self.directory_path, f"{self.basename}.fullPeakList.csv")
        peak_list = pd.read_csv(peaks_file)
        self.peak_list = peak_list
        self.peak_list['checked'] = False
        output_file = os.path.join(self.directory_path, f"{self.basename}.filteredPeakList.csv")
        self.peak_list.to_csv(output_file, index=False)
        df = pd.read_csv(output_file)

        # define features and target
        features = [
            "charge",
            "fit_score",
            "interference",
            "log_s2n",
            "pc_missing_peaks",
            "mass_error_std",
            "chisq_p_value",
            "correlation_p_value",
            "correlation_coefficient",
            "chisq_stat"
        ]

        checked_data = df[df['checked'] == True]
        other_data = df[df['checked'] != True]
        combined_data = pd.concat([checked_data, other_data], ignore_index=True)
        x_combined = combined_data[features]
        scaler = StandardScaler()
        x_combined_scaled = scaler.fit_transform(x_combined)
        lr_model = joblib.load(lr_file)
        gb_model = joblib.load(gb_file)

        # run models to predict unchecked data
        logistic_predictions = lr_model.predict_proba(x_combined_scaled)
        log_true_prob = [sublist[1] for sublist in logistic_predictions]
        combined_data["logistic_prediction_prob"] = log_true_prob

        gradient_boosting_predictions = gb_model.predict_proba(x_combined_scaled)
        gradient_true_prob = [sublist[1] for sublist in gradient_boosting_predictions]
        combined_data["gradient_boosting_prediction_prob"] = gradient_true_prob
        combined_data["prediction"] = self.voting(lr_model, gb_model, x_combined_scaled)

        combined_data.to_csv(output_file, index=False)
        print(f"Saved filtered peak list to {output_file}")


    def voting(self, logistic_model, gradient_boosting_model, x):
        logistic_proba = logistic_model.predict_proba(x)
        gradient_boosting_proba = gradient_boosting_model.predict_proba(x)

        vote_list = []
        for i in range(len(logistic_proba)):
            log_prob = logistic_proba[i][1]
            grad_prob = gradient_boosting_proba[i][1]
            if (log_prob > 0.5) or (grad_prob > 0.5):
                vote_list.append(True)
            else:
                vote_list.append(False)

        return vote_list


    def load_decon_button_handler(self, event, validated=False):
        try:
            file_suffix = "filtered" if validated else "full"
            peaks_file = os.path.join(
                self.directory_path,
                f"{self.basename}.{file_suffix}PeakList.csv"
            )
            peak_list = pd.read_csv(peaks_file)
            self.peak_list = peak_list
            if not validated:
                peak_list["prediction"] = False # just for plotting purposes
            list_filter = [
                "monoisotopic_mw",
                "charge",
                "monoisotopic_mz",
                "abundance",
                "fit_score",
                "interference",
                "s2n",
                "num_missing_peaks",
                "mass_error_std",
                "prediction"
            ]
            peak_list_filtered = peak_list[list_filter].copy()

            self.Freeze()
            self.data_grid.Destroy()
            self.data_grid = wx.grid.Grid(self.panel)
            self.data_grid.Hide()
            self.data_sizer.Add(self.data_grid, 1)
            self.data_sizer.SetItemMinSize(self.data_grid, 2000, 220)
            self.data_grid.SetDefaultCellBackgroundColour("#f0f0f0")
            self.data_grid.Bind(wx.grid.EVT_GRID_SELECT_CELL, self.on_cell_select)

            num_rows, num_cols = peak_list_filtered.shape
            self.data_grid.CreateGrid(num_rows, num_cols)
            for row in range(num_rows):
                for col in range(num_cols):
                    value = peak_list_filtered.iloc[row, col]
                    if col in [0, 2, 4, 5, 6, 8]:
                        value = round(value, 3)
                    if col in [3]:
                        value = round(value, 1)
                    value = str(value)
                    self.data_grid.SetCellValue(row, col, value)
                    self.data_grid.SetReadOnly(row, col)

            column_names = [
                "Monoisotopic MW (Da)",
                "Charge",
                "Monoisotopic m/z",
                "Intensity",
                "Fit",
                "Interference",
                "s/n",
                "# Missing Peaks",
                "std(Mass error)",
                "Validated?"
            ]
            for col, column_name in enumerate(column_names):
                self.data_grid.SetColLabelValue(col, column_name)

            self.data_grid.AutoSizeColumns()
            self.data_grid.Refresh()
            self.data_grid.Show()
            self.panel.Layout()
            self.Thaw()
            event.Skip()

        except FileNotFoundError:
            print(
                f"No {'filtered ' if validated else ''}"
                "deconvolution file is present in the working directory!"
            )


    def on_load_decon_button(self, event):
        self.load_decon_button_handler(event, validated=False)


    def on_load_validated_decon_button(self, event):
        self.load_decon_button_handler(event, validated=True)


    def on_cell_select(self, event):
        selected_row = event.GetRow()

        # read in parameters to send to the plot
        # mass and charge
        mass = self.peak_list.loc[selected_row, "monoisotopic_mw"]
        charge = self.peak_list.loc[selected_row, "charge"]
        mass_diff = self.peak_list.loc[selected_row, "mass_diff"]
        mass_diff = ast.literal_eval(mass_diff)

        # theoretical isotopic envelope
        scatter_x = self.peak_list.loc[selected_row, "theo_x"]
        scatter_x = ast.literal_eval(scatter_x)
        scatter_y = self.peak_list.loc[selected_row, "theo_y"]
        scatter_y = ast.literal_eval(scatter_y)

        # mass and intensity ranges
        abu = float(self.data_grid.GetCellValue(selected_row, 3))*1.2
        selected_min_mz = min(scatter_x) - 3 * (scatter_x[1] - scatter_x[0])
        selected_max_mz = max(scatter_x) + 4 * (scatter_x[1] - scatter_x[0])

        # color options
        color_list = self.peak_list.loc[selected_row, "matched_bool"]
        color_list = ast.literal_eval(color_list)
        matched_mz = self.peak_list.loc[selected_row, "matched_x"]
        matched_mz = ast.literal_eval(matched_mz)
        matched_intens = self.peak_list.loc[selected_row, "matched_y"]
        matched_intens = ast.literal_eval(matched_intens)
        validated = self.peak_list.loc[selected_row, "prediction"]

        if validated == True:
            color_mapping = {0: "#cd4120", 1: "#6ab4e6"}
            color_list = [color_mapping[item] for item in color_list]
            label_color = "#b44a52"
        else:
            color_mapping = {0: "#5a3962", 1: "#009462"}
            color_list = [color_mapping[item] for item in color_list]
            label_color = "#b44a52"

        self.plot_panel.plot_spectrum(
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

        event.Skip()


    def on_val_button(self, event):
        try:
            num_shuffle = int(self.shuffle_peak_input.GetValue())
            val_window = ValidationWindow(
                self,
                f"precisION - Supervised Classification [{os.path.basename(self.file_path)}]",
                self.file_path,
                self.directory_path,
                num_shuffle
            )
            val_window.Show()
        except FileNotFoundError:
            print("No deconvolution file is present in the working directory!")



class ProcessingSpectrumPanel(wx.Panel):
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
        self.max_intens = 0


    def plot_spectrum(
        self, spectrum, selected_min_mz, selected_max_mz, abu,
        color_list, mass, charge, mass_diff, matched_mz,
        matched_intens, label_color, scatter_x=None, scatter_y=None
    ):
        self.ax.cla()
        self.ax.plot(spectrum[:, 0], spectrum[:, 1], color="black", linewidth=1)
        self.ax.set_xlabel("m/z")
        self.ax.set_ylabel("Intensity")
        self.ax.set_facecolor("#f0f0f0")
        self.ax.set_xlim(selected_min_mz, selected_max_mz)
        self.max_intens = abu * 1.1
        self.ax.set_ylim(0, self.max_intens)

        if scatter_x is not None:
            # peak label at top right
            label_text = f"{mass} ({charge}+)"
            self.ax.annotate(
                label_text,
                xy=(0.995,0.99),
                xycoords="axes fraction",
                ha="right",
                va="top"
            )

            # scatter plot w/ theoretical distribution
            self.ax.scatter(scatter_x, scatter_y, c=color_list, marker="o", alpha=0.8)

            # vertical line at each peak
            for peak in zip(matched_mz, matched_intens):
                x_coords = [peak[0], peak[0]]
                y_coords = [0, peak[1]]
                self.ax.plot(x_coords, y_coords, color=label_color, linestyle="--", linewidth=1)

            # mass differences between adjacent peaks
            for i, diff in enumerate(mass_diff):
                diff = round(float(diff), 3)
                diff = f"{diff:.3f}"
                x_pos = (matched_mz[i] + matched_mz[i+1])/2
                y_pos = min([matched_intens[i], matched_intens[i+1]])
                if y_pos > (self.max_intens*0.01):
                    x_coords = [matched_mz[i], matched_mz[i+1]]
                    y_coords = [y_pos, y_pos]
                    self.ax.plot(
                        x_coords,
                        y_coords,
                        color=label_color,
                        linestyle=":",
                        linewidth=1
                    )
                    self.ax.text(
                        x_pos,
                        y_pos,
                        diff,
                        ha="center",
                        va="bottom",
                        color=label_color
                    )

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
