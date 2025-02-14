# Mini-MLST typing tool for N. gonorrhoeae

The program is designed to genotype N. gonorrhoeae isolates by analyzing 18 informative single nucleotide polymorphisms in seven housekeeping genes abcZ, adk, aroE, fumC, gdh, pdhC, pgm.

## The current version of Offline Mini MLST supports 3 modes:

1. Obtain MLST from nucleotide profile. This mode allows you to determine which sequence type is characterized by the entered nucleotide profile. Additionally, if the nucleotide profile is unknown, the most similar MLSTs are displayed. The output of this mode is accompanied by a sequence alignment. Note that nucleotide profiles differing by more than two mismatches will not be displayed. The input mode supports sequences containing 18 nucleotides in uppercase (ATGC). In case of doubtful nucleotide positions or gaps, the letter N is required.

2. Obtain Nucleotide Profile from MLST. This mode allows you to extract the reference nucleotide profile for a specific MLST. Input supports numbers [0-9]. The MLST of interest will be graphically characterized by the nucleotide profile with included coordinates and SNPs.

3. This mode allows you to determine the genogroup that will be characteristic of the MLST of interest. The input supports numbers [0-9]. In addition to genogroup determination, the output also provides a list of other MLSTs within the specified genogroup.

Before running the program, make sure you have installed and updated the following packages: Pandas, regex, Biopython. You can download them manually or type 'pip install -r requirements.txt'.

To run the program, type 'python3 MiniMLST.py'.

Information about program navigation is available by typing "help" or "4" in the main menu. To ensure stable and effective use of the program, we recommend that you follow the console guidelines and enter only appropriate data.

## Copyright:
© This is an open software distributed under the terms of the Creative Commons Attribution License , which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.

## Funding:

This work was supported by the Russian Science Foundation Grant No. 24-25-20084

## Citation: 
Kandinov, I.; Shaskolskiy, B.; Kravtsov, D.; Filippova, M.; Larkin, A.; Gryadunov, D. Mini-Multilocus Sequence Typing Scheme for the Global Population of Neisseria gonorrhoeae. Int. J. Mol. Sci. 2024, 25, 5781. https://doi.org/10.3390/ijms25115781 


