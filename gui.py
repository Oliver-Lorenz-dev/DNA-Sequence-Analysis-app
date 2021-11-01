import tkinter as tk
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
        self.window.geometry("700x350")
        # create a title
        title = tk.Label(self.window, text="FASTA File Analyser")
        # add title to window
        title.pack()
        # create open file button
        tk.Button(self.window, text="Open file", command=lambda: self.open_file()).pack()
        # create a sequence length button
        tk.Button(self.window, text="Sequence length", command=lambda: self.get_length()).pack()
        # create a GC button
        tk.Button(self.window, text="GC content", command=lambda: self.gc()).pack()
        # create a base frequency button
        tk.Button(self.window, text="Base frequency", command=lambda: self.base_frequency()).pack()
        # create entry
        tk.Entry(self.window, textvariable=self.search).place(x=50, y=140)
        # create a sequence search button
        tk.Button(self.window, text="Sequence search", command=lambda: self.seq_search()).place(x=60, y=162)
        # create a transcribe button
        tk.Button(self.window, text="Transcribe DNA sequence", command=lambda: self.transcribe()).pack()
        # create a translate button
        tk.Button(self.window, text="Translate DNA sequence", command=lambda: self.translate()).pack()
        # create a protein analysis button
        tk.Button(self.window, text="Protein sequence analysis", command=lambda: self.protein_analysis()).pack()
        # create a complement button
        tk.Button(self.window, text="Complement DNA sequence", command=lambda: self.complement()).pack()
        # create a reverse complement button
        tk.Button(self.window, text="Reverse Complement DNA sequence", command=lambda: self.rev_complement()).pack()
        # create a get record button
        tk.Button(self.window, text="Get sequence record", command=lambda: self.get_record()).pack()
        # create a quit button and add to window
        tk.Button(self.window, text="Quit", command=self.window.destroy).pack()
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
            # check if sequence is too long to be printed in the GUI
            if len(mRNA) < 101:
                self.text_out.set('mRNA sequence: ' + mRNA)
            else:
                self.text_out.set('Your mRNA sequence is longer than 100 bases, output directed to stdout.')
                print('mRNA: ' + mRNA)
        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')

    def translate(self):
        """function which transcribes DNA to amino acid sequence"""
        # check user has opened a file
        try:
            # sequence length must be a multiple of 3 for translation
            if len(self.content.seq) % 3 == 0:
                protein = self.content.seq.translate()
            elif (len(self.content.seq)+1) % 3 == 0:
                edited_seq_1 = self.content.seq + 'N'
                protein = edited_seq_1.translate()
            else:
                edited_seq_2 = self.content.seq + 'NN'
                protein = edited_seq_2.translate()
            # check if sequence is too long to be printed in the GUI
            if len(protein) < 101:
                self.text_out.set('Protein sequence: ' + protein)
            else:
                self.text_out.set('Your protein sequence is longer than 100 amino acids, output directed to stdout.')
                print('Protein sequence: ' + protein)
        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')

    def complement(self):
        """function which complements DNA sequence"""
        # check user has opened a file
        try:
            comp = self.content.seq.complement()
            # check if sequence is too long to be printed in the GUI
            if len(comp) < 101:
                self.text_out.set('Complement sequence: ' + comp)
            else:
                self.text_out.set('Your complement sequence is longer than 100 bases, output directed to stdout.')
                print('Complement: ' + comp)
        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')

    def rev_complement(self):
        """function which complements DNA sequence"""
        # check user has opened a file
        try:
            rev_comp = self.content.seq.reverse_complement()
            # check if sequence is too long to be printed in the GUI
            if len(rev_comp) < 101:
                self.text_out.set('Reverse Complement sequence: ' + rev_comp)
            else:
                self.text_out.set('Your reverse complement sequence is longer than 100 bases, output directed to stdout.')
                print('Reverse complement: ' + rev_comp)
        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')

    def get_record(self):
        """function which gets the record of the FASTA file"""
        # check user has opened a file
        try:
            # check file has a record
            if self.content.description:
                record = str(self.content.description)
                # check if record is too long to be printed in the GUI
                if len(record) < 101:
                    self.text_out.set(record)
                else:
                    self.text_out.set('Record is longer than 100 characters, output directed to stdout')
                    print(record)
            else:
                self.text_out.set('No record avaliable')
        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')

    def protein_analysis(self):
        """function which generates a simple analysis report of the translated DNA sequence"""
        # check user has opened a file
        try:
            # translate DNA sequence
            if len(self.content.seq) % 3 == 0:
                protein = self.content.seq.translate()
            elif (len(self.content.seq) + 1) % 3 == 0:
                edited_seq_1 = self.content.seq + 'N'
                protein = edited_seq_1.translate()
            else:
                edited_seq_2 = self.content.seq + 'NN'
                protein = edited_seq_2.translate()
            # get length of sequence
            protein_length = len(protein)
            # count types of amino acids
            hydrophobic = (protein.count('A') + protein.count('C') + protein.count('I')
            + protein.count('L') + protein.count('M') + protein.count('F') + protein.count('W')
            + protein.count('V'))
            neutral = (protein.count('G') + protein.count('H') + protein.count('P')
            + protein.count('S') + protein.count('T') + protein.count('Y'))
            hydrophilic = (protein.count('R') + protein.count('N') + protein.count('D')
            + protein.count('Q') + protein.count('E') + protein.count('K'))
            self.text_out.set('Sequence length: ' + str(protein_length) + ' Hydrophobic AAs: ' + str(hydrophobic)
                              + ' Neutral AAs: ' + str(neutral) + ' Hydrophilic AAs: ' + str(hydrophilic))
        # tell user to open a file
        except AttributeError:
            self.text_out.set('Please open a FASTA file before using other functions of this application')


