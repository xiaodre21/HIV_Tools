# HIV_Tools
This repository serves to store tools for various objectives related with HIV sequencing data.

The tools present in this repository will grow and suffer changes along time, as new problems arise and new solutions get discovered.
As such, this README file will be always updated in order to reflect descriptions, instructions and documentation for all the tools present.

## Typing with COMET: (COntext-based Modeling for Expeditious Typing)

COMET is a platform provided by the Luxemburg Institute of Health (https://comet.lih.lu/), for HIV-1 and HIV-2 typing.
It is currently only available on their web platform, which can be practical for some but impractical for others.

**typing_with_comet.py** is a python script created to be ran in the terminal (command-line) to type all sequences present in a fasta file and produce an excel or csv report with the typing results. By default, after typing, the script will automatically create a new fasta file where the headers of each sequences will have the identified type, separated by a comma, and then the previous header. This can be turned off using an optional argument described below.<br />
The script contains a -h or --Help with instructions and descriptions for the arguments.
<br />
<br />
E.g.
(original fasta file)
>IMCJ_KR020_1 
tggcgcccgaacagggacttgaggaagagtgagagtcttcggagcacggctgagtgagggcagtaagggcggcaggaatc
aaccacgacggagagctcctgtaaaagcgcaggccggtaccaggcagcgtgaggagcgggaggagaagaggcctccggga

(new fasta file)
>B.IMCJ_KR020_1 
tggcgcccgaacagggacttgaggaagagtgagagtcttcggagcacggctgagtgagggcagtaagggcggcaggaatc
aaccacgacggagagctcctgtaaaagcgcaggccggtaccaggcagcgtgaggagcgggaggagaagaggcctccggga
<br />

It has the following arguments:
- **--hiv_type**: specify the type of sequences to type (**1** for HIV-1, **2** for HIV-2).
- **--output_type**: specify the output type, **xlsx** or **csv**.
- **--fasta_file**: the directory for the fasta (or fasta.gz).
- **--output_directory**: the directory for the output. A folder called "Comit_Output" will be created in this directory, storing the results.
- **--disable_auto_rename**: Optional argument, **Y** to disable the creation of the new fasta file.
<br />

***Example of usage:***
<br />
**[Windows]**
<br />
python typing_using_comet.py --hiv_type 1 --fasta_file unknown_seqs.fasta --output_directory C:\Users\results --output_type xlsx
<br />

**[MAC]**
<br />
python3 typing_using_comet.py --hiv_type 2 --fasta_file unknown_seqs.fasta.gz --output_directory C:\Users\results --output_type csv --disable_auto_rename Y

### Dependencies

Its always advisable to use a Virtual Environment.<br />
- Bio=1.6.2<br />
- biopython=1.83<br />
- pandas=2.0.3<br />
- Requests=2.28.1<br />
- selenium=4.4.0<br />
- webdriver_manager=4.0.2<br />

**[Windows]**
<br />
python -m venv venv<br />
venv\Scripts\activate<br />
pip install -r requirements.txt

**[MAC]**
<br />
python3 -m venv venv<br />
source venv/bin/activate<br />
pip3 install -r requirements.txt<br />

