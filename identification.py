import os
import wx
import wx.grid

from denovo import DenovoWindow
from databaseSearch import DbSearchWindow


class ProteinIDWindow(wx.Frame):
    def __init__(self, parent, title, file_path, directory_path):
        super().__init__(
            parent,
            title=title,
            size=(406, 245),
            style=wx.DEFAULT_FRAME_STYLE & ~(wx.RESIZE_BORDER | wx.MAXIMIZE_BOX)
        )

        self.file_path = file_path
        self.directory_path = directory_path

        self.panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        #header_sizer: title and information
        header_sizer = wx.BoxSizer(wx.VERTICAL)
        title = wx.StaticText(self.panel, label="Protein Identification")
        title.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        subtext1 = wx.StaticText(
            self.panel,
            label="Identify unknown proteins from fragment data alone."
        )
        subtext2 = wx.StaticText(
            self.panel,
            label="Use either sequence tag (de novo) or open searches."
        )

        header_sizer.Add(title, 0, wx.BOTTOM, 5)
        header_sizer.Add(subtext1, 0, wx.BOTTOM, 5)
        header_sizer.Add(subtext2, 0, wx.BOTTOM, 5)

        #button_sizer: launchers for both modules
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        denovo_button_label = "De novo sequencing"
        denovo_button = wx.Button(self.panel, label=denovo_button_label,  size=(150, 30))
        denovo_button.Bind(wx.EVT_BUTTON, self.on_denovo_button)
        db_button_label = "Open search"
        db_button = wx.Button(self.panel, label=db_button_label,  size=(150, 30))
        db_button.Bind(wx.EVT_BUTTON, self.on_db_button)

        button_sizer.Add(db_button, 0, wx.RIGHT, 5)
        button_sizer.Add(denovo_button, 0, wx.RIGHT, 0)


        # options_sizer: options for filtering and database input
        options_staticbox = wx.StaticBox(self.panel)
        options_sizer = wx.StaticBoxSizer(options_staticbox, wx.VERTICAL)
        options_text = wx.StaticText(self.panel, label="General search options")
        options_text.SetFont(
            wx.Font(
                10,
                wx.FONTFAMILY_DEFAULT,
                wx.FONTSTYLE_ITALIC,
                wx.FONTWEIGHT_BOLD
            )
        )

        database_file_label = wx.StaticText(self.panel, label="Database file (.xml):")

        db_file_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.selected_file_text = wx.TextCtrl(self.panel, style=wx.TE_READONLY | wx.TE_CENTER)
        self.selected_file_text.SetMinSize((290, -1))
        file_button = wx.Button(self.panel, label="Browse")
        file_button.Bind(wx.EVT_BUTTON, self.on_file_button)

        db_file_sizer.Add(self.selected_file_text, 0, wx.RIGHT, 5)
        db_file_sizer.Add(file_button, 0, wx.RIGHT, 0)

        options_sizer.Add(options_text, 0, wx.BOTTOM, 5)
        options_sizer.Add(database_file_label, 0, wx.BOTTOM, 0)
        options_sizer.Add(db_file_sizer, 0, wx.BOTTOM, 5)

        main_sizer.Add(header_sizer, 0, wx.ALL, 5)
        main_sizer.Add(button_sizer, 0, wx.BOTTOM |  wx.ALIGN_CENTER_HORIZONTAL, 0)
        main_sizer.Add(options_sizer, 0, wx.EXPAND | wx.ALL, 5)

        self.panel.SetSizer(main_sizer)

        icon = wx.Icon("./icons/icon.ico", wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        self.Show()


    def on_file_button(self, event):
        wildcard = "XML files (*.xml)|*.xml|" \
                   "All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self,
            "Choose a database file (.xml)",
            style=wx.FD_OPEN,
            defaultDir=os.getcwd(),
            wildcard=wildcard
        )

        if dlg.ShowModal() == wx.ID_OK:
            file_path = dlg.GetPath()
            self.selected_file_text.SetValue(file_path)
            print("Selected file:", file_path)
        dlg.Destroy()


    def on_denovo_button(self, event):
        db_file = self.selected_file_text.GetValue()
        if db_file == "":
            print("No database file selected!")
            return
        denovo_button_window = DenovoWindow(
            self,
            f"precisION - De Novo Sequencing [{os.path.basename(self.file_path)}]",
            self.file_path,
            self.directory_path,
            db_file
        )
        denovo_button_window.Show()


    def on_db_button(self, event):
        db_file = self.selected_file_text.GetValue()
        if db_file == "":
            print("No database file selected!")
            return
        denovo_button_window = DbSearchWindow(
            self,
            f"precisION - Open Search [{os.path.basename(self.file_path)}]",
            self.file_path,
            self.directory_path,
            db_file
        )
        denovo_button_window.Show()
