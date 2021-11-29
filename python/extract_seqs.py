import subprocess
import sys
import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import re
import time

def parser():
    purpose = ''' This is a scipt to extract specific positions from fasta sequences Usage:
            python acba_mlst_runner.py --seq <fasta_file> --pos <csv of pos to extract> --output <output_prefixes> '''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='acba_mlst_runner.py')

    parser.add_argument('--seq', required=True, help='List of seqeuence files (FASTA) (required)', type=str)
    parser.add_argument('--output', required=True, help='Prefix of output files  (required)', type=str)
    parser.add_argument('--pos', required=True, help='csv of positions to extract', type=str)
    parser.add_argument('--length-cutoff', default=0, help='minimum length of fragment to extract [default = 0]', type=int)

    args = parser.parse_args()

    return args

def whole_fasta_creator(seq_file):
    """ Take in a fasta file with multiple contigs and then extract a DNA sequence out """
    seqs_list = []
    with open(seq_file) as input_handle:
        seqs = SeqIO.parse(input_handle, "fasta")
        for seq in seqs:
            seqs_list.append(str(seq.seq))

    seqsy = "".join(seqs_list)

    return seqsy
def get_pos_list(pos_file):
    """ take in csv file and get list of lists out"""
    liszt = []
    with open(pos_file, "r") as input_handle:
        listo = input_handle.read().splitlines()
        for item in listo:
            items = re.split(",", item)
            current_list = []
            counter = 1
            for num in items:
                if counter == 1:
                    num_c = int(num) - 1
                else:
                    num_c = int(num)
                counter += 1
                current_list.append(num_c)
            liszt.append(current_list)

    return liszt

def write_out_locs(sequence, locs, output, length):
    """ Take in the sequence and the positions, write out each sequence as a separate item"""
    counter = 1
    tot_length = len(locs)
    number = tot_length
    count = 0
    while (number > 0):
        number = number//10
        count = count + 1

    num_digits = "{0:0=" + str(count) + "d}"
    print("Beginning write out")

    for loc in tqdm.tqdm(locs):
        iter_val = num_digits.format((counter))
        # print("Writing out %s of %s rows." % (iter_val, tot_length),
        #       end="\r", flush=True)

        if (loc[1] - loc[0]) < length:
            continue
        current_out = output + "_" + str(counter)
        current_seq = Seq(sequence[loc[0]:loc[1]])
        current_record = SeqRecord(current_seq, id = current_out,
                                   annotations={"molecule_type":"DNA"})
        SeqIO.write(current_record,(current_out + ".fasta"), "fasta")
        time.sleep(0.01)

        counter += 1
    print("")

    concat_cmd = "cat " + output + "_*.fasta > " + output + "_total.fasta"
    try:
        subprocess.check_call(concat_cmd, shell=True)
    except subprocess.SubprocessError:
        rm_cmd = "rm " + output + "_total.fasta"
        try:
            subprocess.check_call(rm_cmd, shell=True)
            subprocess.check_call(concat_cmd, shell=True)
        except subprocess.SubprocessError:
            sys.exit("Failed concating the output fragments")

    rm_cmd = "rm " + output + "_[0-9]*.fasta"
    try:
        subprocess.check_call(rm_cmd, shell=True)
    except subprocess.SubprocessError:
        sys.exit("Failed removing the individual fragment files")


if __name__ == '__main__':
    # Get input args

    input_args = parser()

    # Get sequence

    seqcess = whole_fasta_creator(input_args.seq)

    # Get the locs

    locs = get_pos_list(input_args.pos)

    #Write out the sequences

    write_out_locs(seqcess, locs, input_args.output, input_args.length_cutoff)

    print("Done")




