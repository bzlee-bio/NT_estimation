# Prediction of spider novel neurotoxic peptides & peptide data augmentation
There are two tools that are proposed in the research.
1. Prediction of spider neurotoxic peptide from amino acid sequences only.<br>
2. Data augmentation tool for peptide data

<br><br>
# Prediction tools for neurotoxic peptide 
## Dependencies
To install python dependencies for neurotoxic peptides, run: `pip install -r requirements_neurotox_pred.txt`

## Input file 
Input file type is fasta format in which amino acids are represented using single-letter codes.

Detailed information of fasta links: https://en.wikipedia.org/wiki/FASTA_format

## Output file
Output file only contains peptide list which are predicted as neurotoxin.

## Running prediction of spider novel neurotoxic peptides

`python ntest.py --fasta <input_fasta_file.fasta> --output <output_file_name.csv>`
<br><br>
# Data augmentation for proteins
## Dependencies
To install python dependencies, run: `pip install -r requirements.txt`<br>
Also, NCBI BLAST should be installed.<br><br>
After install NCBI blast, BLAST bin path in aug_config.ini file should be modified.

<br><br>
## Installation of NCBI BLAST
### Ubuntu / Debian
Run
`$ sudo apt-get install ncbi-blast`
### CenotOS
1. Download latest BLAST version (<a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download">Link</a>)
2. Decompress
`$ tar -zxvf ncbi-blast-version_info-linux.tar.gz`

## Input file 
Input file type is fasta format in which amino acids are represented using single-letter codes.

Detailed information of fasta links: https://en.wikipedia.org/wiki/FASTA_format

## Output file
Output file contains information about probabilities of four respective ion-channel modulability.

Probability with >=0.5 predicts as modulator peptides for respective ion channels.

## Running a data augmentation tool
`python peptide_augmentation --fasta <input_fasta_file.fasta> --cpu <max_cpu_usage> --eval <E-value cutoff>`



## Citation
Lee, Byungjo, Min K. Shin, In-Wook Hwang, Junghyun Jung, Yu J. Shim, Go W. Kim, Seung T. Kim, Wonhee Jang, and Jung-Suk Sung. 2021. "A Deep Learning Approach with Data Augmentation to Predict Novel Spider Neurotoxic Peptides" International Journal of Molecular Sciences 22, no. 22: 12291. https://doi.org/10.3390/ijms222212291
