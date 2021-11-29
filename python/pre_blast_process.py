import os
import subprocess
import sys
import re
import pandas
import numpy
import tqdm
import shutil

def make_dna_db(seq_loc):
    """ Function to take in a list of fasta sequences
        Check if each sequence is already a chromosome
        if so copy to new dna_dir else make dna and move to
        dna_dir, then combine dna files into one mfa file
    """

    python_dir_name = os.path.dirname(os.path.realpath(__file__))
    perl_dir = re.sub("python", "perl/", python_dir_name)

    if os.path.isdir('tmp_dna_lib'):
        shutil.rmtree('contig_bounds')
        print("Deleting current dna directory")
        os.mkdir("contig_bounds")

    else:
        os.mkdir("tmp_dna_lib")

    for seq in range(len(seq_loc)):
        current_seq = seq_loc[seq]
        grep_cmd = "grep -c '^>' " + current_seq
        grep_out = subprocess.check_output(grep_cmd, shell=True)
        grep_num = int(grep_out.strip().decode(encoding="utf-8"))

        if grep_num > 1:
            dna_cmd = "perl " + perl_dir + "converting_velvet_contigs_to_dna.pl " + current_seq
            try:
                dna_file = subprocess.check_output(dna_cmd, shell=True)
            except subprocess.SubprocessError:
                sys.exit("Failed creating the DNA file for %s" % os.path.basename(current_seq))
            dna_file = re.sub("^.*is\s","",dna_file.strip().decode(encoding="utf-8"))
            mv_cmd = "mv " + dna_file + " ./tmp_dna_lib"
            try:
                subprocess.check_output(mv_cmd, shell=True)
            except subprocess.SubprocessError:
                sys.exit("Failed moving the newly created DNA file into the tmp dna dir for iso %s" % os.path.basename(current_seq))
        elif grep_num == 1:
            new_dna_name = os.path.splitext(os.path.basename(current_seq))[0] + ".dna"
            cp_cmd = "cp " + current_seq + " ./tmp_dna_lib/" + new_dna_name
            new_fast_header_cmd = "sed -i '1s/.*/\>" + os.path.basename(current_seq)+ "/' ./tmp_dna_lib/" + new_dna_name
            print(new_fast_header_cmd)
            try:
                subprocess.check_output(cp_cmd, shell=True)
            except subprocess.SubprocessError:
                sys.exit("Failed copying the chromosome file into the tmp dna dir for iso %s" % os.path.basename(current_seq))
            try:
                subprocess.check_output(new_fast_header_cmd, shell=True)
            except subprocess.SubprocessError:
                sys.exit("Failed renaming the fasta header for the chromosome file iso %s" % os.path.basename(current_seq))

        else:
            sys.exit("Non-positive number of contigs in fasta file: %s" % os.path.basename(current_seq))


    ## Now creating the mfa file

    mfa_cmd = "cd ./tmp_dna_lib && cat *.dna > output.mfa"
    print()
    print("Running MFA concatenation step now")
    print()
    try:
        subprocess.check_output(mfa_cmd, shell=True)
    except subprocess.SubprocessError:
        sys.exit("Failed creating the mfa DNA file")

    return "./tmp_dna_lib"

def contig_bounds(seqs):
    """ Function to check the bounds of contigs in a file"""
    if os.path.isdir('contig_bounds'):
        shutil.rmtree('contig_bounds')
        print("Deleting current contig bounds directory")
        os.mkdir("contig_bounds")
    else:
        os.mkdir("contig_bounds")
    print("This many seqs to get through for contig bounds: " + str(len(seqs)))
    for seq in tqdm.tqdm(range(len(seqs))):

        grp_cmd = "grep -n '^>' " + seqs[seq] + " > grep_output_temp"

        try:
            subprocess.check_output(grp_cmd, shell=True)
        except subprocess.SubprocessError:
            sys.exit("Failed trying to search for contigs in %s " % os.path.basename(seqs[seq]))

        grep_file = open("grep_output_temp", "r")

        grep_file_lines = grep_file.readlines()

        contig_starts = []

        for k in range(len(grep_file_lines)):
            current_line = grep_file_lines[k]
            current_line_start = re.split(":", current_line)[0]
            current_line_start = re.sub("\'", "", current_line_start)
            current_line_start = int(current_line_start)
            contig_starts.append(current_line_start)

        ###############################################################################
        ## Now we'll get the fasta file in ############################################
        ###############################################################################

        fasta_file = open(seqs[seq], "r")

        temp = fasta_file.read().splitlines()

        contig_bounds = pandas.DataFrame(data=numpy.zeros(shape=(len(contig_starts), 2)))

        for k in range(len(contig_starts)):
            if k == 0:
                contig_bounds.iloc[0, 0] = 1
            else:
                contig_bounds.iloc[k, 0] = contig_bounds.iloc[k - 1, 1] + 1

            current_start = contig_starts[k]
            if (k + 1) != len(contig_starts):
                current_end = contig_starts[k + 1] - 1
            else:
                current_end = len(temp)

            if (k + 1) != len(contig_starts):
                dist = sum(len(i) for i in temp[current_start:current_end])
                contig_bounds.iloc[k, 1] = contig_bounds.iloc[k, 0] + dist - 1
            else:
                dist = sum(len(i) for i in temp[current_start:current_end])
                contig_bounds.iloc[k, 1] = contig_bounds.iloc[k, 0] + dist - 1


        contig_name = os.path.basename(seqs[seq])
        file_path = "./contig_bounds/" + contig_name + "#contig_bounds.csv"

        contig_bounds.to_csv(path_or_buf=file_path, index=False)

    return "./contig_bounds"













