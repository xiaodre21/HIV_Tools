# HIV_Tools
This repository serves to store tools for various objectives related with HIV sequencing data.

The tools present in this repository will grow and suffer changes along time, as new problems arise and new solutions get discovered.
As such, this README file will be always updated in order to reflect descriptions, instructions and documentation for all the tools present.

## Typing with COMET: (COntext-based Modeling for Expeditious Typing)

COMET is a platform provided by the Luxemburg Institute of Health, for HIV-1 and HIV-2 typing.
It is currently only available in their web platform, which can be practical for some but impractical for others.

**typing_with_comet.py** is a python script created to be ran in the terminal (commandline) to type all sequences present in a fasta file and produce an excel or csv report with the typing results. By default, after typing, the script will automatically create a new fasta file where the headers of each sequences will have the identified type, separated by a comma, and then the previous header.
<br />
<br />
E.g.
(original fasta file)
> IMCJ_KR020_1 
tggcgcccgaacagggacttgaggaagagtgagagtcttcggagcacggctgagtgagggcagtaagggcggcaggaatc
aaccacgacggagagctcctgtaaaagcgcaggccggtaccaggcagcgtgaggagcgggaggagaagaggcctccggga

(new fasta file)
> B.IMCJ_KR020_1 
tggcgcccgaacagggacttgaggaagagtgagagtcttcggagcacggctgagtgagggcagtaagggcggcaggaatc
aaccacgacggagagctcctgtaaaagcgcaggccggtaccaggcagcgtgaggagcgggaggagaagaggcctccggga
<br />

It has the following arguments:
- **--hiv_type**: specify the type of sequences to type (1 for HIV-1, 2 for HIV-2).
- **--output_type**: specify the output type, xlsx or csv.
- **--fasta_file**: the directory for the fasta (or fasta.gz).
- **--output_directory**: the directory for the output. A folder called "Comit_Output" will be created in this directory, storing the results.
- **--disable_auto_rename**: Optional argument, "Y" to disable the creation of the new fasta file.
<br />

### Dependencies
- needs to be filled

