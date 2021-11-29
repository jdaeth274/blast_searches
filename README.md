# BLAST searching isolates #

## Installation ##

This requires conda, please install conda first [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
Once installed clone the repo:

    git clone https://github.com/jdaeth274/blast_searches.git

Then use the environment.yml file to install the dependencies with conda

    cd blast_searches 
    conda env create --file=environment.yml

Activate this environment using 

    conda activate blast_env

## Usage ##

To run you need a `.txt` file with the full path to each of the fasta files on each line and a fasta file of the 
query sequence you want to blast for:

    python blast_searches/blast-runner.py --seqs <fasta txt file> --query <fasta file of query> --output <output name for csv results>

This will produce a csv file of the BLAST results. 






