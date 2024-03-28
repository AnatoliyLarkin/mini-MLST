# Repository1

Offline extention for MiniMLST program allows free and effective opportunity to study and interpret multilocus sequence typing data for Neisseria gonorrhoeae.

Current version of offline-MiniMLST supports 3 modes:

1) Obtain MLST from nucleotide profile
This mode allows to determine which sequence-type will be characterized by entered nucleotide profile. In addition, in case of unknown nucleotide profile, the most similar MLSTs will be shown. Output of this mode is accompanied with sequences alignment. Notice, that nucleotide profiles, that differ from the entered one by more than two mismatches will not be shown. Mode input supports sequences, containing 18 nucleotides in capital notation (ATGC). In case of doubtful nucleotide position or gaps letter N is required.

2) Obtain nucleotide profile from MLST
This mode allows to extract the reference nucleotide profile for specific MLST. Input supports numbers [0-9]. MLST of interest will be graphically charactized by nucleotide profile with included coordinates and SNPs.

3) Determine which genogroup the sequence type belongs to
This mode allows to obtain the genogroup, that will be characteristic for MLST of interest. Input supports numbers [0-9]. In addition to genogroup determination, output also provides the list of other      MLSTs inside the specific genogroup.

Before program execution ensure, that you have following packages installed and updated: Pandas, regex, biopython. You can download them manually or type 'pip install -r requirements.txt'.

To execute the program type 'python3 MiniMLST.py'.

Information about program navigation is available by writing command 'help' or '4' in main menu. To ensure the stable and effective program usage we recommend you to follow console guidelines and input appropriate data only.

If any problems, errors or questions arose, please contact us - larkinanatanat@yandex.ru 
 


