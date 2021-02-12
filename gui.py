import tkinter as tk
from tkinter import StringVar
import tkinter.ttk as ttk
from tkinter import simpledialog
from tkinter import messagebox
from tkinter.filedialog import askopenfile
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import nt_search


class GUI(tk.Frame):
    def __init__(self, master=tk.Tk()):
        # create a window
        super().__init__(master)
        self.window = master
        # set string variable
        self.text_out = tk.StringVar()
        # create string variable for sequence search
        self.search = tk.StringVar()
        # set size of GUI
        self.window.geometry("500x250")
        # create a title
        title = tk.Label(self.window, text="Bioinformatics App")
        # add title to window
        title.pack()
        # create a quit button and add to window
        tk.Button(self.window, text="Quit", command=self.window.destroy).pack()
        # create open file button
        tk.Button(self.window, text="Open file", command=lambda: self.open_file()).pack()
        # create a GC button
        tk.Button(self.window, text="Sequence length", command=lambda: self.get_length()).pack()
        # create a GC button
        tk.Button(self.window, text="GC content", command=lambda: self.gc()).pack()
        # create a base frequency button
        tk.Button(self.window, text="Base frequency", command=lambda: self.base_frequency()).pack()
        # create entry
        tk.Entry(self.window, textvariable=self.search).place(x=20, y=70)
        # create a base frequency button
        tk.Button(self.window, text="Sequence search", command=lambda: self.seq_search()).place(x=30, y=92)
        # create a transcribe button
        tk.Button(self.window, text="Transcribe sequence", command=lambda: self.transcribe()).pack()
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
            # reset text_out
            self.text_out.set('')
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

    def base_frequency(self):
        """function which determines the percentage of each base in a DNA sequence"""
        # check user has opened a file
        try:
            L = len(self.content.seq)
            A = (self.content.seq.count('A') / L) * 100
            G = (self.content.seq.count('G') / L) * 100
            C = (self.content.seq.count('C') / L) * 100
            T = (self.content.seq.count('T') / L) * 100
            N = (self.content.seq.count('N') / L) * 100
            self.text_out.set('A: ' + "{:.2f}".format(A) + '%  ' + 'G: ' + "{:.2f}".format(G) +
                  '%  ' + 'C: ' + "{:.2f}".format(C) + '%  ' +
                  'T: ' + "{:.2f}".format(T) + '%  ' + 'N: ' + "{:.2f}".format(N) + '%  ')
        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')

    def seq_search(self):
        """function which allows user to search for a nucleotide sequence of their choice"""
        # check user has opened a file
        try:
            # check for search input
            if len(self.search.get()) > 0:
                # reset text_out
                self.text_out.set('')
                search_input = self.search.get()
                # avoid case sensitivity
                search_input = search_input.upper()
            else:
                self.text_out.set('No search input detected')
            # print to stdout because output can be very large for this function
            # check if search input is correct
            try:
                if len(self.search.get()) > 0:
                    # check for matches
                    if len(nt_search(str(self.content.seq), str(search_input))) > 1:
                        match_length = len(nt_search(str(self.content.seq), str(search_input)))
                        matches = match_length - 1
                        self.text_out.set('Sequence: ' + search_input + ': ' + str(matches) + ' matches')
                    else:
                        self.text_out.set('Sequence: ' + search_input + ' No matches found')
                else:
                    pass
            # tell user search input is incorrect
            except KeyError:
                self.text_out.set('Please search for the bases A,C,G and T only')

        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')

    def get_length(self):
        """function which returns the length of a sequence"""
        # check user has opened a file
        try:
            seq_len = len(self.content.seq)
            self.text_out.set('Sequence length: ' + str(seq_len))
        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')

    def transcribe(self):
        """function which transcribes DNA to mRNA"""
        # check user has opened a file
        try:
            mRNA = self.content.seq.transcribe()
            if len(mRNA) < 51:
                self.text_out.set('mRNA sequence: ' + mRNA)
            else:
                self.text_out.set('Your mRNA sequence is longer than 50 bases, the output is in the CLI.')
                print(mRNA)
        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')


g = GUI()
g