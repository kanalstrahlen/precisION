import os
import wx
import wx.richtext as rt
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
from Bio import SeqIO
from denovoFunctions import DenovoFunctions

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


class DenovoWindow(wx.Frame):
    def __init__(self, parent, title, file_path, directory_path, database_file):
        super().__init__(parent, title=title)

        self.file_path = file_path
        self.directory_path = directory_path
        self.database_file = database_file

        #screen width variable to allow for smaller screens
        screen_width, screen_height = wx.DisplaySize()
        self.SetSize((screen_width * 3 // 4, screen_height * 3 // 4))

        self.panel = wx.Panel(self)
        self.plot_panel = DeNovoSpectrumPanel(self.panel)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # top_sizer: all information and user inputs
        top_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # menu_sizer: all information, run button and options
        menu_sizer = wx.BoxSizer(wx.VERTICAL)
        title= wx.StaticText(self.panel, label="De Novo Sequencing")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext1 = wx.StaticText(
            self.panel,
            label="Identify sequence tags and search sequence databases for those tags."
        )
        subtext2 = wx.StaticText(
            self.panel,
            label="High mass accuracy is required to distinguish K/Q (<18 mDa)."
        )
        search_button = wx.Button(self.panel, label="Run search",  size=(150, 30))
        search_button.Bind(wx.EVT_BUTTON, self.on_search_button)

        # options_sizer: all options for searching
        options_staticbox = wx.StaticBox(self.panel)
        options_sizer = wx.StaticBoxSizer(options_staticbox, wx.VERTICAL)
        options_text = wx.StaticText(self.panel, label="Search options")
        options_text.SetFont(
            wx.Font(
                10,
                wx.FONTFAMILY_DEFAULT,
                wx.FONTSTYLE_ITALIC,
                wx.FONTWEIGHT_BOLD
            )
        )

        # options_subsizer: horizontal options input
        options_subsizer = wx.BoxSizer(wx.HORIZONTAL)
        accuracy_text = wx.StaticText(self.panel, label="Mass accuracy (mDa): ")
        self.accuracy_text_ctrl = wx.TextCtrl(self.panel, size=(30, 20), value="5")
        length_text = wx.StaticText(self.panel, label="Peptide length: ")
        self.length_text_ctrl = wx.TextCtrl(self.panel, size=(30, 20), value="5")

        options_subsizer.Add(accuracy_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        options_subsizer.Add(self.accuracy_text_ctrl, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)
        options_subsizer.Add(length_text, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        options_subsizer.Add(self.length_text_ctrl, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)

        options_sizer.Add(options_text, 0, wx.BOTTOM, 5)
        options_sizer.Add(options_subsizer, 0, wx.BOTTOM,  5)

        menu_sizer.Add(title, 0, wx.BOTTOM, 5)
        menu_sizer.Add(subtext1, 0, wx.BOTTOM, 5)
        menu_sizer.Add(subtext2, 0, wx.BOTTOM, 10)
        menu_sizer.Add(search_button, 0, wx.BOTTOM |  wx.ALIGN_CENTER_HORIZONTAL, 5)
        menu_sizer.Add(options_sizer, 0, wx.EXPAND | wx.BOTTOM, 5)

        # tag_sizer: list of sequence tags identified from denovo search
        tags_sizer = wx.BoxSizer(wx.VERTICAL)
        tag_title= wx.StaticText(self.panel, label="Sequence tags")
        subtitle_font = wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD)
        tag_title.SetFont(subtitle_font)
        self.sequence_tags = [] #list of tags currently displayed in the box
        self.tag_list_box = wx.ListBox(
            self.panel,
            choices=self.sequence_tags,
            size=(200, 150),
            style=wx.VSCROLL | wx.HSCROLL
        )
        bg_colour = "#f0f0f0"
        self.tag_list_box.SetBackgroundColour(bg_colour)
        self.tag_list_box.Bind(wx.EVT_LISTBOX, self.on_tag_select)

        tags_sizer.Add(tag_title, 0, wx.BOTTOM, 5)
        tags_sizer.Add(self.tag_list_box, 0, wx.BOTTOM, 0)

        # proteins_sizer: list of proteins identified by database search of selected tag
        proteins_sizer = wx.BoxSizer(wx.VERTICAL)
        protein_title= wx.StaticText(self.panel, label="Proteins")
        protein_title.SetFont(subtitle_font)
        self.proteins = [] #list of proteins currently displayed in the box
        self.proteins_list_box = wx.ListBox(
            self.panel,
            choices=self.proteins,
            size=(200, 150),
            style=wx.VSCROLL | wx.HSCROLL
        )
        self.proteins_list_box.SetBackgroundColour(bg_colour)
        self.proteins_list_box.Bind(wx.EVT_LISTBOX, self.on_protein_select)

        proteins_sizer.Add(protein_title, 0, wx.BOTTOM, 5)
        proteins_sizer.Add(self.proteins_list_box, 0, wx.BOTTOM, 0)

        # fullseq_sizer: display the full sequence of the protein that has returned a match
        fullseq_sizer = wx.BoxSizer(wx.VERTICAL)
        fullseq_title= wx.StaticText(self.panel, label="Full sequence")
        fullseq_title.SetFont(subtitle_font)
        self.fullseq = "" #sequence currently displayed in the box
        self.full_seq_text_box = rt.RichTextCtrl(
            self.panel,
            style=wx.VSCROLL | wx.HSCROLL,
            size=(200, 150)
        )
        self.full_seq_text_box.SetBackgroundColour(bg_colour)
        self.full_seq_text_box.SetValue(self.fullseq)
        
        start = 0
        end = 0
        range_style = rt.RichTextAttr()
        range_style.SetFontWeight(wx.FONTWEIGHT_BOLD)
        self.full_seq_text_box.SetStyle(start, end, range_style)

        fullseq_sizer.Add(fullseq_title, 0, wx.BOTTOM, 5)
        fullseq_sizer.Add(self.full_seq_text_box, 0, wx.EXPAND | wx.BOTTOM, 0 )

        # button_sizer: buttons for analysing reulsts
        button_sizer = wx.BoxSizer(wx.VERTICAL)
        db_search_button =  wx.Button(
            self.panel,
            label="Search database",
            size=(160, 30)
        )
        db_search_button.Bind(wx.EVT_BUTTON, self.on_search_db_button)
        add_match_button = wx.Button(
            self.panel,
            label="Add protein to .xml file",
            size=(160, 30)
        )
        add_match_button.Bind(wx.EVT_BUTTON, self.on_add_match_button)
        add_tag_matches_button =  wx.Button(
            self.panel,
            label="Add tag matches to .xml file",
            size=(160, 30)
        )
        add_tag_matches_button.Bind(wx.EVT_BUTTON, self.on_add_tag_matches_button)
        add_all_matches_button =  wx.Button(
            self.panel,
            label="Add all matches to .xml file",
            size=(160, 30)
        )
        add_all_matches_button.Bind(wx.EVT_BUTTON, self.on_add_all_matches_button)
        save_xml_button =  wx.Button(
            self.panel,
            label="Save .xml file",
            size=(160, 30)
        )
        save_xml_button.Bind(wx.EVT_BUTTON, self.on_save_xml_button)

        button_sizer.Add(db_search_button, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 5)
        button_sizer.Add(add_match_button, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 5)
        button_sizer.Add(add_tag_matches_button, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 5)
        button_sizer.Add(add_all_matches_button, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 5)
        button_sizer.Add(save_xml_button, 0, wx.BOTTOM | wx.ALIGN_CENTER_HORIZONTAL, 5)

        top_sizer.Add(menu_sizer, 0, wx.RIGHT, 15)
        top_sizer.Add(tags_sizer, 0, wx.RIGHT, 15)
        top_sizer.Add(proteins_sizer, 0, wx.RIGHT, 15)
        top_sizer.Add(fullseq_sizer, 0, wx.RIGHT, 15)
        top_sizer.Add(button_sizer, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)

        main_sizer.Add(top_sizer, 0, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(self.plot_panel, 1, wx.EXPAND | wx.TOP, -10)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.load_spectrum(self.file_path)
        self.protein_selection_id_list = [] #list of proteins selected to save in xml file

        # generate fasta file. have to do otherwise database search is incredibly slow
        xml_basename = self.database_file.replace(".xml", "")
        fasta_file = f"{xml_basename}.fasta"

        # if fasta file already exists, use already created file
        if os.path.exists(fasta_file):
            print(f"Using {fasta_file}...")
        else:
            print("Generating FASTA file...")
            with open(self.database_file, "r") as file, open(fasta_file, "w") as output_handle:
                records = SeqIO.parse(file, "uniprot-xml")
                SeqIO.write(records, output_handle, "fasta")

        self.panel.SetSizer(main_sizer)
        self.panel.Bind(wx.EVT_SIZE, self.on_size)
        self.Show()

        self.mz_list = []
        self.mass_list = []
        self.mass_diff_list = []
        self.aa_list =[]
        self.peptides = []
        self.charge_list = []
        self.abundance_list = []
        self.prot_id_list = []
        self.prot_name_list = []
        self.start_res_list = []
        self.end_res_list = []
        self.prot_sequence_list = []


    def on_size(self, event):
        size = self.panel.GetClientSize()
        self.plot_panel.SetSize(size.width // 2, size.height)
        self.panel.Layout()
        event.Skip()


    def load_spectrum(self, file_path):
        self.basename = os.path.basename(file_path).replace(".txt", "")
        save_path = os.path.join(self.directory_path, f"{self.basename}.spectrum.txt")
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
            None
        )


    def on_search_button(self, event):
        self.tag_list_box.Clear()

        # load in peak list and the relevant files
        peak_file_path = os.path.join(self.directory_path, f"{self.basename}.filteredPeakList.csv")
        peak_list = pd.read_csv(peak_file_path)
        accuracy = float(self.accuracy_text_ctrl.GetValue())
        length = int(self.length_text_ctrl.GetValue())

        # find pairs of losses, then assemble these into tags
        run = DenovoFunctions()
        aa_loss_pair_list = run.find_aa_loss(peak_list, accuracy/1000)
        tag_matches = run.find_peptides(aa_loss_pair_list, length)
        # compile all relevant information for plots and filtering
        (
            self.mz_list,
            self.mass_list,
            self.mass_diff_list,
            self.aa_list,
            self.peptides,
            self.charge_list,
            self.abundance_list
        ) = run.find_aa_seq(peak_list, tag_matches, accuracy/1000)

        try:
            self.tag_list_box.InsertItems(self.peptides, 0)
        except:
            print("No sequence tags found with the current search parameters.")


    def on_tag_select(self, event):
        self.proteins_list_box.Clear()
        self.full_seq_text_box.Clear()

        sequence_index = self.tag_list_box.GetSelection()

        # get all information about current sequence tag
        mz = self.mz_list[sequence_index]
        mass = self.mass_list[sequence_index]
        aa = self.aa_list[sequence_index]
        charge = [self.charge_list[sequence_index]]*len(mass)
        abundance = self.abundance_list[sequence_index]

        # calculate averagine distributions for all peaks and scale to data
        run = DenovoFunctions()
        averagine = [run.averagine(m, z) for m, z in zip(mass, charge)]
        averagine_x = []
        averagine_y=[]
        for i, peak_abundance in enumerate(abundance):
            peaks_mz = averagine[i][0]
            peaks_abu = averagine[i][1]
            max_peak = max(peaks_abu)
            ratio = peak_abundance/max_peak
            peaks_abu = [peak_abu * ratio for peak_abu in peaks_abu]
            averagine_x += peaks_mz
            averagine_y += peaks_abu

        # plot range
        min_mz = min(mz)
        max_mz = max(mz)

        # update plot
        self.plot_panel.plot_spectrum(
            self.spectrum,
            min_mz,
            max_mz,
            max(abundance),
            mz,
            aa,
            averagine_x,
            averagine_y
        )


    def on_search_db_button(self, event):
        # extract sequence tag from selection
        sequence_index = self.tag_list_box.GetSelection()
        sequence_tag = self.peptides[sequence_index]

        run = DenovoFunctions()
        # search uniprot database
        (
            self.prot_id_list,
            self.prot_name_list,
            self.start_res_list,
            self.end_res_list,
            self.prot_sequence_list
        ) = run.search_uniprot_fasta(self.database_file, sequence_tag)


        # update protein list
        if len(self.prot_name_list) >= 1:
            self.proteins_list_box.Clear()
            self.proteins_list_box.InsertItems(self.prot_name_list, 0)
        else:
            print("No matching proteins found!")


    def on_protein_select(self, event):
        # extract protein name from selection
        protein_index = self.proteins_list_box.GetSelection()
        sequence = self.prot_sequence_list[protein_index]

        self.full_seq_text_box.SetValue(sequence)
        #bold the matched region
        start = self.start_res_list[protein_index]
        end = self.end_res_list[protein_index]
        range_style = rt.RichTextAttr()
        range_style.SetFontWeight(wx.FONTWEIGHT_BOLD)
        self.full_seq_text_box.SetStyle(start, end, range_style)


    def on_add_match_button(self, event):
        protein_index = self.proteins_list_box.GetSelection()
        up_id = self.prot_id_list[protein_index]
        self.protein_selection_id_list.append(up_id)
        print("Added protein. Remember to save the file.")


    def on_add_tag_matches_button(self, event):
        for up_id in self.prot_id_list:
            self.protein_selection_id_list.append(up_id)
        print("Added proteins that match the selected tag. Remember to save the file.")


    def on_add_all_matches_button(self, event):
        run = DenovoFunctions()
        print("Searching for all matches:")
        for tag in self.peptides:
            prot_id_list, prot_name_list, _, _, _ = run.search_uniprot_fasta(
                self.database_file,
                tag
            )
            if len(prot_name_list) > 0:
                print(f"Matched {tag} with:")
                for name in prot_name_list:
                    print(name)
            for up_id in prot_id_list:
                self.protein_selection_id_list.append(up_id)
        print("Added proteins that match the list of tags. Remember to save the file.")


    def on_save_xml_button(self, event):
        selected_ids = set(self.protein_selection_id_list)
        print("Generating .xml file for the following records:")
        for identifier in selected_ids:
            print(identifier)
        self.write_xml(selected_ids)


    def write_xml(self, selected_ids):
        xml_basename = os.path.basename(self.database_file).replace(".xml", "")
        xml_dir = os.path.dirname(self.directory_path)
        save_path = os.path.join(xml_dir, f"{xml_basename}.denovoFiltered.{self.basename}.xml")

        with open(self.database_file, "r", encoding="utf-8") as input_file, \
             open(save_path, "w", encoding="utf-8") as output_file:
            for _ in range(2):
                line = input_file.readline()
                output_file.write(line)
            in_entry = False
            entry_lines = []
            last_line = []
            for line in input_file:
                if line.strip().startswith("<entry "):
                    in_entry = True
                    entry_lines = [line]
                elif in_entry:
                    entry_lines.append(line)
                    if line.strip() == "</entry>":
                        accession = None
                        for entry_line in entry_lines:
                            if "<accession>" in entry_line:
                                accession = entry_line.strip() \
                                .replace("<accession>", "") \
                                .replace("</accession>", "")
                                break
                        if accession in selected_ids:
                            output_file.writelines(entry_lines)
                        in_entry = False
                last_line.append(line)
                if len(last_line) > 1:
                    last_line.pop(0)
            output_file.writelines(last_line)
            print(f"Saved selected entries to {save_path}")



class DeNovoSpectrumPanel(wx.Panel):
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
        self.max_intens = 1


    def plot_spectrum(
        self, spectrum, selected_min_mz, selected_max_mz,
        abu, mz_lines, aa, scatter_x, scatter_y
    ):
        self.ax.cla()
        self.ax.plot(spectrum[:, 0], spectrum[:, 1], color="black", linewidth=1)
        self.ax.set_xlabel("m/z")
        self.ax.set_ylabel("Intensity")
        self.ax.set_facecolor("#f0f0f0")
        self.ax.set_xlim(selected_min_mz-100, selected_max_mz+100)
        self.max_intens = abu * 1.1
        self.ax.set_ylim(0, self.max_intens)

        # add in vertical lines at peak mz and annotate halfway between with aa id
        if mz_lines is not None:
            for i, peak in enumerate(mz_lines):
                x_coords = [peak, peak]
                y_coords = [0, abu*1.1]
                self.ax.plot(x_coords, y_coords, color= "#b42920", linestyle="--", linewidth=1)
                if i >= 1:
                    res_name = aa[i-1]
                    pos = (peak + prev_peak) / 2
                    self.ax.text(
                        pos,
                        abu,
                        res_name,
                        horizontalalignment="center",
                        color = "#b42920"
                    )
                prev_peak = peak

        # plot averagine distributions
        self.ax.scatter(scatter_x, scatter_y, c="#205283", marker="o", alpha=0.5)

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
