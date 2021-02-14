Advanced Programming for Digital Biology
========================================
This package is built as a part of the CSC8330: Advanced Programming
for Digital Biology.

Bioinformatics GUI App
----------------------
Usage
-----
Type `python main.py` to run this application.

This application is designed to work on DNA sequences in the form of a FASTA file.

Options
-------
Open file: Click the `Open file` button to open the file browser and then click on the FASTA file you want to open. Note: you must open a file
before using other functions of this application.

Sequence length: Click on the `Sequence length` button to find out the length of your DNA sequence

GC content: Click on the `GC content` button to find out the GC content of your DNA sequence.

Base frequency: Click on the `Base frequency` button to find out the number of each type of base in your DNA sequence:
(A,C,T,G and N): N represents an unknown base. 

Transcribed DNA sequence: Click on the `Transcribe DNA sequence` button to transcribe your DNA sequence into mRNA. 
Note: if your sequence is longer than 100 bases it will be printed to stdout (the command line), if it is less than 100
bases long it will be printed in the GUI.

Translate DNA sequence: Click on the `Translate DNA sequence` button to translate your DNA sequence into a protein
sequence. Note: if your sequence is longer than 100 amino acids it will be printed to stdout (the command line), 
if it is less than 100 amino acids long it will be printed in the GUI.

Protein sequence analysis: Click on the `Protein sequence analysis` button to generate an analysis of the
translated DNA sequence. This analysis includes: the length of the protein sequence, and the number of 
hydrophobic, neutral and hydrophilic amino acids in the sequence.
The definitions of what amino acids are hydrophobic, neutral and hydrophilic are from:
http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html

Complement DNA sequence: Click on the `Complement DNA sequence` button to generate the complement of your DNA sequence.
Note: if your sequence is longer than 100 bases it will be printed to stdout (the command line), 
if it is less than 100 bases long it will be printed in the GUI.

Reverse complement DNA sequence: Click on the `Reverse Complement DNA sequence` button 
to generate the reverse complement of your DNA sequence.
Note: if your sequence is longer than 100 bases it will be printed to stdout (the command line), 
if it is less than 100 bases long it will be printed in the GUI.

Get sequence record: Click on the `Get sequence record` button to get the record of your sequence. 
Note: if your record is longer than 100 characters it will be printed to stdout (the command line), 
if it is less than 100 characters long it will be printed in the GUI.

Sequence search: Type your search input into the input box above the `Sequence search` button and then
press the sequence search button to perform the search. This will return the number of times the sequence 
you searched for is in your DNA sequence (from your FASTA file). Your search must contain only the following characters:
A, C, G and T. 

Quit: Press the `Quit` button to exit the application.
