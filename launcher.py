import os
import matplotlib
import wx
import wx.adv
from matplotlib.backends import backend_pdf, backend_svg

from bulkDecon import DeconvolutionWindow
from processing import ProcessingWindow
from identification import ProteinIDWindow
from fragmentAnnotator import FragmentAnnotatorWindow
from fragmentAnalysis import FragmentAnalysisWindow
from proteoformAnalysis import ProteoformAnalysisWindow

matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

# to do; alter peaklist to envelopelist and update all other names
# to do evaluate cosine score for socring peak fits

class LauncherWindow(wx.Frame):
    def __init__(self, parent):
        super().__init__(
            parent,
            title="precisION",
            size=(905, 300),
            style=wx.DEFAULT_FRAME_STYLE & ~(wx.RESIZE_BORDER | wx.MAXIMIZE_BOX)
        )

        print("precisION")
        print("If you use this software, please cite:")

        panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        #info_sizer: left panel w/ software information
        info_sizer = wx.BoxSizer(wx.VERTICAL)
        logo = wx.Image("./icons/icon.png", wx.BITMAP_TYPE_ANY)
        logo = logo.Rescale(100, 100)
        bitmap = wx.StaticBitmap(panel, bitmap=wx.Bitmap(logo))
        title_text = wx.StaticText(panel, label="precisION")
        title_text.SetFont(wx.Font(
            16,
            wx.FONTFAMILY_DEFAULT,
            wx.FONTSTYLE_NORMAL,
            wx.FONTWEIGHT_BOLD
            )
        )
        sub_text = wx.StaticText(panel, label="Precise and accurate identification of native ")
        sub_text2 = wx.StaticText(panel, label="proteoforms from top-down mass spectra.")
        cite_text = wx.StaticText(panel, label="Manual and References ")
        cite_text.SetFont(wx.Font(
            10,
            wx.FONTFAMILY_DEFAULT,
            wx.FONTSTYLE_NORMAL,
            wx.FONTWEIGHT_BOLD
            )
        )
        pub_text = wx.adv.HyperlinkCtrl(
            panel,
            label="GitHub",
            url="https://github.com/kanalstrahlen/precisION"
        )

        info_sizer.Add(bitmap, 0, wx.ALIGN_CENTER_HORIZONTAL)
        info_sizer.Add(title_text, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 10)
        info_sizer.Add(sub_text, 0, wx.ALIGN_CENTER_HORIZONTAL)
        info_sizer.Add(sub_text2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        info_sizer.Add(cite_text, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, 10)
        info_sizer.Add(pub_text, 0, wx.ALIGN_CENTER_HORIZONTAL)

        main_sizer.Add(info_sizer, 0, wx.ALL | wx.CENTER, 10)

        #input_sizer: right panel w/ module and input selection
        input_sizer = wx.BoxSizer(wx.VERTICAL)
        top_button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        bottom_button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        file_input_sizer = wx.BoxSizer(wx.HORIZONTAL)
        dir_input_sizer = wx.BoxSizer(wx.HORIZONTAL)

        #decon module button
        decon_button_label = (
            "Spectral Deconvolution\n"
            "Deconvolute fragmentation spectra using multiple algorithms and cluster results"
        )
        decon_button = wx.Button(panel, label=decon_button_label, size=(200, 85))
        decon_button.Bind(wx.EVT_BUTTON, self.on_decon_button)

        #decon processing module button
        processing_button_label = (
            "Envelope Filtering\n"
            "Filter deconvolution results using supervised classification or manual filters"
        )
        processing_button = wx.Button(panel, label=processing_button_label, size=(200, 85))
        processing_button.Bind(wx.EVT_BUTTON, self.on_processing_button)

        #protein id module
        protein_id_button_label = (
            "Protein Identification\n"
            "Identify unknown proteins from fragmentation spectra "
            "using sequence tag or database searches"
        )
        protein_id_button = wx.Button(panel, label=protein_id_button_label, size=(200, 85))
        protein_id_button.Bind(wx.EVT_BUTTON, self.on_protein_id_button)

        #fragment matching module
        assignment_button_label = (
            "Fragment Assignment\n"
            "Assign ions in fragmentation\n"
            "spectra using a hierarchical,\n"
            "semi-supervised approach"
        )
        assignment_button = wx.Button(panel, label=assignment_button_label, size=(200, 85))
        assignment_button.Bind(wx.EVT_BUTTON, self.on_assignment_button)

        #fragment analysis module
        analysis_button_label = (
            "Fragment Analysis\n"
            "Analyse the relationship between protein sequence, structure, and fragmentation"
        )
        analysis_button = wx.Button(panel, label=analysis_button_label, size=(200, 85))
        analysis_button.Bind(wx.EVT_BUTTON, self.on_analysis_button)

        #proteoform analysis module
        proteoform_analysis_button_label = (
            "Proteoform Analysis\n"
            "Assign fragments to specific proteoforms without "
            "deconvolution by fitting theoretical envelopes"
        )
        proteoform_analysis_button = wx.Button(
            panel,
            label=proteoform_analysis_button_label,
            size=(200, 85)
        )
        proteoform_analysis_button.Bind(wx.EVT_BUTTON, self.on_proteoform_analysis_button)

        top_button_sizer.Add(decon_button, 0, wx.ALL, 5)
        top_button_sizer.Add(processing_button, 0, wx.ALL, 5)
        top_button_sizer.Add(protein_id_button, 0, wx.ALL, 5)

        bottom_button_sizer.Add(assignment_button, 0, wx.ALL, 5)
        bottom_button_sizer.Add(analysis_button, 0, wx.ALL, 5)
        bottom_button_sizer.Add(proteoform_analysis_button, 0, wx.ALL, 5)

        #file selection
        default_file = "Spectra should be deconvoluted prior to selection."
        file_button = wx.Button(panel, label="Select File")
        file_button.Bind(wx.EVT_BUTTON, self.on_file_button)
        self.selected_file_text = wx.TextCtrl(
            panel,
            value=default_file,
            style=wx.TE_READONLY | wx.TE_CENTER
        )
        self.selected_file_text.SetMinSize((540, -1))

        file_input_sizer.Add(file_button, 0, wx.ALIGN_LEFT | wx.RIGHT, 5)
        file_input_sizer.Add(self.selected_file_text, 0)

        #directory selection
        default_dir = "The working directory will be set automatically after selecting file."
        dir_button = wx.Button(panel, label="Select Dir.")
        dir_button.Bind(wx.EVT_BUTTON, self.on_dir_button)
        self.selected_dir_text = wx.TextCtrl(
            panel,
            value=default_dir,
            style=wx.TE_READONLY | wx.TE_CENTER
        )
        self.selected_dir_text.SetMinSize((540, -1))

        dir_input_sizer.Add(dir_button, 0, wx.ALIGN_LEFT | wx.RIGHT, 5)
        dir_input_sizer.Add(self.selected_dir_text, 0)

        input_sizer.Add(top_button_sizer, 0, wx.ALL | wx.CENTER, 0)
        input_sizer.Add(bottom_button_sizer, 0, wx.ALL | wx.CENTER, 0)
        input_sizer.Add(file_input_sizer, 0, wx.ALL | wx.CENTER, 2)
        input_sizer.Add(dir_input_sizer, 0, wx.ALL | wx.CENTER, 5)

        main_sizer.Add(input_sizer, 0, wx.ALL, 5 | wx.CENTER, 0)

        panel.SetSizer(main_sizer)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.Show()


    def on_file_button(self, event):
        wildcard = "Text files (*.txt)|*.txt|" \
           "All files (*.*)|*.*"

        dlg = wx.FileDialog(
            self,
            "Choose a file",
            style=wx.FD_OPEN,
            wildcard=wildcard
        )

        if dlg.ShowModal() == wx.ID_OK:
            file_path = dlg.GetPath()
            self.selected_file_text.SetValue(file_path)
            print("Selected file:", file_path)

            # Create a folder for the processed files
            directory_path = os.path.dirname(file_path)
            file_name = os.path.basename(file_path)
            folder_path = os.path.join(
                directory_path,
                "{}_precisION_files".format(os.path.splitext(file_name)[0])
            )

            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
                print("Created folder:", folder_path)
            self.selected_dir_text.SetValue(folder_path)
            print("Selected directory:", folder_path)
        dlg.Destroy()


    def on_dir_button(self, event):
        dlg = wx.DirDialog(self, "Choose a directory")
        if dlg.ShowModal() == wx.ID_OK:
            directory_path = dlg.GetPath()
            self.selected_dir_text.SetValue(directory_path)
            print("Selected directory:", directory_path)
        dlg.Destroy()


    def on_decon_button(self, event):
        dlg = wx.DirDialog(self, "Choose a directory")
        if dlg.ShowModal() == wx.ID_OK:
            directory_path = dlg.GetPath()
            print("Selected directory:", directory_path)
        dlg.Destroy()
        decon_settings_window = DeconvolutionWindow(
            self,
            "precisION - Spectral Deconvolution",
            directory_path
        )
        decon_settings_window.Show()


    def on_processing_button(self, event):
        try:
            file_path = self.selected_file_text.GetValue()
            directory_path = self.selected_dir_text.GetValue()
            _, file_extension = os.path.splitext(file_path)
            if file_extension.lower() != ".txt":
                raise ValueError()
            processing_button_window = ProcessingWindow(
                self,
                f"precisION - Envelope Filtering [{os.path.basename(file_path)}]",
                file_path,
                directory_path
            )
            processing_button_window.Show()
        except:
            print((
                "Could not open file. "
                "Inputs should be tab-delimited .txt files "
                "containing the m/z and intensity. "
                "Deconvolution must be completed first."
                )
            )


    def on_protein_id_button(self, event):
        file_path = self.selected_file_text.GetValue()
        directory_path = self.selected_dir_text.GetValue()

        _, file_extension = os.path.splitext(file_path)
        if file_extension.lower() != ".txt":
            print((
                "Could not open file. "
                "Inputs should be tab-delimited .txt files "
                "containing the m/z and intensity. "
                "Deconvolution and filtering must be completed first."
                )
            )
            return

        filtered_peaklist_name = os.path.basename(file_path).replace(".txt", ".filteredPeakList.csv")
        filtered_peaklist_path = os.path.join(directory_path, filtered_peaklist_name)
        
        if os.path.exists(filtered_peaklist_path):
            print(f"Using {filtered_peaklist_path} as peak list for searching.")
        else:
            print("Peak lists must be filtered using the envelope filtering module.")
            return

        identification_button_window = ProteinIDWindow(
            self,
            f"precisION - Protein Identification [{os.path.basename(file_path)}]",
            file_path,
            directory_path
        )
        identification_button_window.Show()


    def on_assignment_button(self, event):
        file_path = self.selected_file_text.GetValue()
        directory_path = self.selected_dir_text.GetValue()
        _, file_extension = os.path.splitext(file_path)
        if file_extension.lower() != ".txt":
            print((
                "Could not open file. "
                "Inputs should be tab-delimited .txt files "
                "containing the m/z and intensity. "
                "Deconvolution and filtering must be completed first."
                )
            )
            return

        filtered_peaklist_name = os.path.basename(file_path).replace(".txt", ".filteredPeakList.csv")
        filtered_peaklist_path = os.path.join(directory_path, filtered_peaklist_name)

        if os.path.exists(filtered_peaklist_path):
            print(f"Using {filtered_peaklist_path} as peak list for assignment.")
        else:
            print("Peak lists must be filtered using the envelope filtering module.")
            return

        assignment_button_window = FragmentAnnotatorWindow(
            self,
            f"precisION - Fragment Assignment [{os.path.basename(file_path)}]",
            file_path,
            directory_path
        )
        assignment_button_window.Show()


    def on_analysis_button(self, event):
        file_path = self.selected_file_text.GetValue()
        directory_path = self.selected_dir_text.GetValue()

        file_name = os.path.basename(file_path).replace(".txt", ".assignedPeakList.csv")
        file_path = os.path.join(directory_path, file_name)

        if os.path.exists(file_path):
            analysis_button_window = FragmentAnalysisWindow(
                self,
                f"precisION - Fragment Analysis [{os.path.basename(file_path)}]",
                file_path,
                directory_path
            )
            analysis_button_window.Show()
        else:
            print("The assigned peak list file is not present.")


    def on_proteoform_analysis_button(self, event):
        file_path = self.selected_file_text.GetValue()
        directory_path = self.selected_dir_text.GetValue()

        if os.path.exists(file_path):
            proteoform_analysis_button_window = ProteoformAnalysisWindow(
                self,
                f"precisION - Proteoform Analysis [{os.path.basename(file_path)}]",
                file_path,
                directory_path
            )
        else:
            print((
                "Could not open file. "
                "The accepted input is a .txt file containing the m/z (x) and intensity (y)."
                )
            )


if __name__ == "__main__":
    app = wx.App()
    LauncherWindow(None)
    app.MainLoop()
