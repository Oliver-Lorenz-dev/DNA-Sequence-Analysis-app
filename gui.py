import tkinter as tk
from tkinter import StringVar
import tkinter.ttk as ttk
from tkinter import simpledialog
from tkinter import messagebox
from tkinter.filedialog import askopenfile
from Bio import SeqIO
from Bio.SeqUtils import GC


class GUI(tk.Frame):
    def __init__(self, master=tk.Tk()):
        # create a window
        super().__init__(master)
        self.window = master
        self.text_out = tk.StringVar()
        # set size of GUI
        self.window.geometry("400x200")
        # create a title
        title = tk.Label(self.window, text="Bioinformatics app")
        # add title to window
        title.pack()
        # create a quit button and add to window
        tk.Button(self.window, text="Quit", command=self.window.destroy).pack()
        # create open file button
        tk.Button(self.window, text="Open file", command=lambda: self.open_file()).pack()
        # create a GC button
        tk.Button(self.window, text="GC content", command=lambda: self.gc()).pack()
        # create output label
        out_label = tk.Label(self.window, textvariable=self.text_out)
        # add out_label to window
        out_label.pack()
        # wait for user input
        self.mainloop()

    def open_file(self):
        """function which opens a file"""
        self.file = askopenfile(mode='r+', filetypes=[('FASTA Files', '*.fasta'), ('FASTA files', '*.fa')])
        # check user has selected a file
        try:
            self.content = SeqIO.read(self.file, 'fasta')
        # tell user no file was opened if they press cancel
        except AttributeError:
            self.text_out.set('You pressed cancel, no new file was opened')
            # note if a file has already been opened this will stay open

    def gc(self):
        """function which determines the GC content of the sequence in a FASTA file"""
        # check user has opened a file
        try:
            self.GC_content = GC(self.content.seq)
            self.text_out.set('GC content: ' + "{:.2f}".format(self.GC_content) + '%')
        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')


g = GUI()
g