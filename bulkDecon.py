import os
import threading
import shutil
import subprocess
import wx
import pandas as pd

from peakProperties import PeakProperties
from peakPicking import RLPeakPicking, CWTPeakPicking
from decon import Thrash, TopFD, Cluster

class DeconvolutionWindow(wx.Frame):
    def __init__(self, parent, title, directory_path):
        super().__init__(parent, title=title, size=(430, 610))

        panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # window title
        title = wx.StaticText(panel, label="Spectral Deconvolution Settings")
        title.SetFont(
            wx.Font(
                16,
                wx.FONTFAMILY_DEFAULT,
                wx.FONTSTYLE_NORMAL,
                wx.FONTWEIGHT_BOLD
            )
        )
        subtext = wx.StaticText(
            panel,
            label="Inputs should be tab-delimited .txt files containing the m/z and intensity."
        )

        # options_sizer: all decon options
        options_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # left_sizer: UniThrash and decon_button
        left_sizer = wx.BoxSizer(wx.VERTICAL)

        # mode_sizer: general settings
        mode_staticbox = wx.StaticBox(panel)
        mode_sizer = wx.StaticBoxSizer(mode_staticbox, wx.VERTICAL)
        subtitle_font = wx.Font(
            12,
            wx.FONTFAMILY_DEFAULT,
            wx.FONTSTYLE_NORMAL,
            wx.FONTWEIGHT_BOLD
        )
        mode_header = wx.StaticText(panel, label="Deconvolution mode")
        mode_header.SetFont(subtitle_font)
        mode_label = wx.StaticText(panel, label="Preset:")
        mode_choices = ["Extensive", "Rapid"]
        self.mode_dropdown = wx.ComboBox(panel, choices=mode_choices, style=wx.CB_READONLY)
        self.mode_dropdown.SetValue("Rapid")

        mode_sizer.Add(mode_header, 0, wx.EXPAND | wx.BOTTOM, 5)
        mode_sizer.Add(mode_label, 0, wx.ALIGN_LEFT)
        mode_sizer.Add(self.mode_dropdown, 0, wx.EXPAND | wx.TOP, 5)
        mode_sizer.SetMinSize((200, -1))

        # thrash_sizer: all UniThrash options
        thrash_staticbox = wx.StaticBox(panel)
        thrash_sizer = wx.StaticBoxSizer(thrash_staticbox, wx.VERTICAL)
        thrash_header = wx.StaticText(panel, label="THRASH")

        thrash_header.SetFont(subtitle_font)
        s2n_label = wx.StaticText(panel, label="Signal to noise ratio:")
        self.s2n_input = wx.TextCtrl(panel, value="3")
        pbr_label = wx.StaticText(panel, label="Peak background ratio:")
        self.pbr_input = wx.TextCtrl(panel, value="3")
        threshold_label = wx.StaticText(panel, label="Data is thresholded:")
        threshold_choices = ["True", "False"]
        self.threshold_dropdown = wx.ComboBox(
            panel,
            choices=threshold_choices,
            style=wx.CB_READONLY
        )
        self.threshold_dropdown.SetValue("True")
        pbrt_label = wx.StaticText(panel, label="Peptide background ratio:")
        self.pbrt_input = wx.TextCtrl(panel, value="3")
        max_charge_label = wx.StaticText(panel, label="Maximum fragment ion \ncharge:")
        self.max_charge_input = wx.TextCtrl(panel, value="20")
        max_mw_label = wx.StaticText(panel, label="Maximum fragment ion \nmass (Da):")
        self.max_mw_input = wx.TextCtrl(panel, value="50000")

        thrash_sizer.Add(thrash_header, 0, wx.EXPAND | wx.BOTTOM, 5)
        thrash_sizer.Add(s2n_label, 0, wx.ALIGN_LEFT)
        thrash_sizer.Add(self.s2n_input, 0, wx.EXPAND | wx.TOP, 5)
        thrash_sizer.Add(pbr_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        thrash_sizer.Add(self.pbr_input, 0, wx.EXPAND | wx.TOP, 5)
        thrash_sizer.Add(threshold_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        thrash_sizer.Add(self.threshold_dropdown, 0, wx.EXPAND | wx.TOP, 5)
        thrash_sizer.Add(pbrt_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        thrash_sizer.Add(self.pbrt_input, 0, wx.EXPAND | wx.TOP, 5)
        thrash_sizer.Add(max_charge_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        thrash_sizer.Add(self.max_charge_input, 0, wx.EXPAND | wx.TOP, 5)
        thrash_sizer.Add(max_mw_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        thrash_sizer.Add(self.max_mw_input, 0, wx.EXPAND | wx.TOP, 5)
        thrash_sizer.SetMinSize((200, -1))

        left_sizer.Add(mode_sizer, 0)
        left_sizer.Add(thrash_sizer, 0)

        # right_sizer: uniTopFD and clustering
        right_sizer = wx.BoxSizer(wx.VERTICAL)

        # msdeconv_sizer: all uniTopFD options
        topfd_staticbox = wx.StaticBox(panel)
        topfd_sizer = wx.StaticBoxSizer(topfd_staticbox, wx.VERTICAL)
        topfd_sizer.SetMinSize((200, -1))
        topfd_header = wx.StaticText(panel, label="TopFD")
        topfd_header.SetFont(subtitle_font)
        topfd_activation_label = wx.StaticText(panel, label="Activation type:")
        topfd_activation_choices = ["CID", "ETD", "UVPD"]
        self.topfd_activation_dropdown = wx.ComboBox(
            panel,
            choices=topfd_activation_choices,
            style=wx.CB_READONLY
        )
        self.topfd_activation_dropdown.SetValue("CID")
        topfd_max_charge_label = wx.StaticText(panel, label="Maximum fragment ion \ncharge:")
        self.topfd_max_charge_input = wx.TextCtrl(panel, value="20")
        topfd_max_mass_label = wx.StaticText(panel, label="Maximum fragment ion \nmass (Da):")
        self.topfd_max_mass_input = wx.TextCtrl(panel, value="50000")
        topfd_mz_error_label = wx.StaticText(panel, label="m/z error (Da):")
        self.topfd_mz_error_input = wx.TextCtrl(panel, value="0.02")
        topfd_s2n_label = wx.StaticText(panel, label="Signal to noise ratio:")
        self.topfd_s2n_input = wx.TextCtrl(panel, value="2")

        topfd_sizer.Add(topfd_header, 0)
        topfd_sizer.Add(topfd_activation_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        topfd_sizer.Add(self.topfd_activation_dropdown, 0, wx.EXPAND | wx.TOP, 5)
        topfd_sizer.Add(topfd_max_charge_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        topfd_sizer.Add(self.topfd_max_charge_input, 0, wx.EXPAND | wx.TOP, 5)
        topfd_sizer.Add(topfd_max_mass_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        topfd_sizer.Add(self.topfd_max_mass_input, 0, wx.EXPAND | wx.TOP, 5)
        topfd_sizer.Add(topfd_mz_error_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        topfd_sizer.Add(self.topfd_mz_error_input, 0, wx.EXPAND | wx.TOP, 5)
        topfd_sizer.Add(topfd_s2n_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        topfd_sizer.Add(self.topfd_s2n_input, 0, wx.EXPAND | wx.TOP, 5)

        # cluster_sizer: clustering options
        cluster_staticbox = wx.StaticBox(panel)
        cluster_sizer = wx.StaticBoxSizer(cluster_staticbox, wx.VERTICAL)
        cluster_sizer.SetMinSize((200, -1))
        cluster_header = wx.StaticText(panel, label="Clustering")
        cluster_header.SetFont(subtitle_font)
        cluster_cutoff_label = wx.StaticText(panel, label="Clustering cutoff (ppm):")
        self.cluster_cutoff_input = wx.TextCtrl(panel, value="5")
        cluster_s2n_label = wx.StaticText(panel, label="Threshold signal-to-noise ratio:")
        self.cluster_s2n_input = wx.TextCtrl(panel, value="3")

        cluster_sizer.Add(cluster_header, 0)
        cluster_sizer.Add(cluster_cutoff_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        cluster_sizer.Add(self.cluster_cutoff_input, 0, wx.EXPAND | wx.TOP, 5)
        cluster_sizer.Add(cluster_s2n_label, 0, wx.ALIGN_LEFT | wx.TOP, 5)
        cluster_sizer.Add(self.cluster_s2n_input, 0, wx.EXPAND | wx.TOP, 5)

        right_sizer.Add(topfd_sizer, 0)
        right_sizer.Add(cluster_sizer)

        options_sizer.Add(left_sizer, 0, wx.LEFT, 5 | wx.CENTER, 0)
        options_sizer.Add(right_sizer, 0, wx.LEFT, 5 | wx.CENTER, 0)

        # decon_button: run deconvolution
        decon_button = wx.Button(panel, label="Deconvolve spectra", size=(400, 30))
        decon_button.Bind(wx.EVT_BUTTON, self.on_decon_button)

        #html_sizer: save html file option
        html_sizer = wx.BoxSizer(wx.HORIZONTAL)

        html_label = wx.StaticText(panel, label="Save deconvolved spectrum (large file size):")
        self.html_checkbox = wx.CheckBox(panel)
        self.html_checkbox.SetValue(True)

        html_sizer.Add(html_label, 0)
        html_sizer.Add(self.html_checkbox, 0, wx.LEFT, 5)

        main_sizer.Add(title, 0, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(subtext, 0, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(options_sizer, 0)
        main_sizer.Add(decon_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, 5)
        main_sizer.Add(html_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, 5)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        panel.SetSizer(main_sizer)
        self.Show()

        self.directory_path = directory_path


    def on_decon_button(self, event):
        background_thread = threading.Thread(target=self.decon)
        background_thread.start()


    def decon(self):
        mode = self.mode_dropdown.GetValue()
        thrash_decon_params = self.read_thrash_params()
        topfd_decon_params = self.read_topfd_params()
        cluster_ppm = self.cluster_cutoff_input.GetValue()
        cluster_s2n = self.cluster_s2n_input.GetValue()
        html_file = self.html_checkbox.GetValue()

        print("Deconvolution can take some time. Please be patient.")

        # loop over each file in the directory and run deconvolutions and clustering
        all_files = os.listdir(self.directory_path)
        txt_files = [file for file in all_files if file.endswith(".txt")]

        for i, file in enumerate(txt_files):
            file_path = os.path.join(self.directory_path, file)
            file_name = os.path.basename(file)
            folder_path = os.path.join(
                self.directory_path,
                f"{os.path.splitext(file_name)[0]}_precisION_files"
            )
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
                print("Created folder:", folder_path)

            if mode == 'Extensive':
                print(f"Deconvolving {file} using extensive settings")
                file_base = os.path.basename(file_path).replace(".txt", "")
                cwt_file = os.path.join(folder_path, f"{file_base}.centroid.cwt.txt")
                rl_file = os.path.join(folder_path, f"{file_base}.centroid.rl.txt")
                convolved_file = os.path.join(folder_path, f"{file_base}.profile.rl.txt")
                centroid_mzml = os.path.join(folder_path, f"{file_base}_centroid.mzML")
                centoid_txt = os.path.join(folder_path, f"{file_base}_centroid.txt")

                if all(os.path.exists(file) for file in (cwt_file, rl_file, convolved_file)):
                    print('Peak picking already complete. Rerunning deconvolution.')
                else:
                    print('Picking peaks using RL deconvolution.')
                    RLPeakPicking(file_path, folder_path, html_file)

                Thrash(convolved_file, folder_path, thrash_decon_params, file_base, True)
                TopFD(centroid_mzml, folder_path, topfd_decon_params)
                Cluster(file_path, folder_path, cluster_ppm)
                PeakProperties(file_path, folder_path, cluster_s2n)
                shutil.copy(rl_file, centoid_txt)
                print(f"\n{i+1}/{len(txt_files)} spectra completed.\n")

            elif mode == "Rapid":
                print(f"Deconvolving {file} using rapid settings")
                file_base = os.path.basename(file_path).replace(".txt", "")
                cwt_file = os.path.join(folder_path, f"{file_base}.centroid.cwt.txt")
                centroid_mzml = os.path.join(folder_path, f"{file_base}_centroid.mzML")
                centoid_txt = os.path.join(folder_path, f"{file_base}_centroid.txt")

                if os.path.exists(cwt_file):
                    print('Peak picking already complete. Rerunning deconvolution.')
                else:
                    print('Picking peaks using CWT.')
                    CWTPeakPicking(file_path, folder_path, html_file)

                Thrash(file_path, folder_path, thrash_decon_params, file_base, False)
                TopFD(centroid_mzml, folder_path, topfd_decon_params)
                Cluster(file_path, folder_path, cluster_ppm)
                PeakProperties(file_path, folder_path, cluster_s2n)
                shutil.copy(cwt_file, centoid_txt)
                print(f"\n{i+1}/{len(txt_files)} spectra completed.\n")

        print("Deconvolution complete. You can now load the clustered results.")
        wx.CallAfter(self.Close)


    def read_topfd_params(self):
        topfd_activation_value = self.topfd_activation_dropdown.GetValue()
        topfd_max_charge_value = self.topfd_max_charge_input.GetValue()
        topfd_max_mass_value = self.topfd_max_mass_input.GetValue()
        topfd_mz_error_value = self.topfd_mz_error_input.GetValue()
        topfd_s2n_value = self.topfd_s2n_input.GetValue()

        topfd_decon_params = [
            topfd_activation_value,
            topfd_max_charge_value,
            topfd_max_mass_value,
            topfd_mz_error_value,
            topfd_s2n_value
        ]

        return topfd_decon_params


    def read_thrash_params(self):
        s2n_value = self.s2n_input.GetValue()
        pbr_value = self.pbr_input.GetValue()
        threshold_value = self.threshold_dropdown.GetValue()
        pbrt_value = self.pbrt_input.GetValue()
        max_charge_value = self.max_charge_input.GetValue()
        max_mw_value = self.max_mw_input.GetValue()

        thrash_decon_params = [
            s2n_value,
            pbr_value,
            threshold_value,
            pbrt_value,
            max_charge_value,
            max_mw_value
        ]

        return thrash_decon_params

