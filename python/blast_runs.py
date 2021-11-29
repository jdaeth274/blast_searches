import subprocess
import os
import sys


def blast_runs(tmp_dna_dir, data_dir, output):
    """ Function to blast the mfa file against the conserved ends
        of the AbaRs in the data directory """

    ## Make the blast db for the input .mfa files

    mfa_file = tmp_dna_dir + "/output.mfa"
    mfa_blast = mfa_file + ".nin"
    mfa_blast_2 = mfa_file + ".nal"
    if not os.path.exists(mfa_blast) and not os.path.exists(mfa_blast_2):
        print("No BLASTDB file exists")
        make_db_cmd = "makeblastdb -dbtype nucl -out " + tmp_dna_dir + "/output.mfa -max_file_sz 2GB -in " + mfa_file
        try:
            subprocess.check_output(make_db_cmd, shell=True)
        except subprocess.SubprocessError:
            sys.exit("Failed making the BLASTDB for the output mfa file, ensure BLAST is installed")



    ## Now we'll blast the left end
        ## Now we'll blast for intact comM
    comM_file = data_dir
    blast_cmd = "blastn -db " + mfa_file + " -query " + comM_file + " -outfmt 10 -out " + output + " -evalue 0.001 -num_alignments 1000000"

    try:
        subprocess.check_output(blast_cmd, shell=True)
    except subprocess.SubprocessError:
        sys.exit("Failed running comM blast search")

