import os
from python import pre_blast_process
from python import blast_runs
from python import post_blast_process
import pandas
import re

def main(input_args):
    """ Main function to run through steps of ABar and comM identification
    """
    print(input_args)
    with open(input_args.seqs, "r") as input_handle:
        seq_lines = input_handle.read().splitlines()

    if input_args.dna_dir is None:
        dna_dir = pre_blast_process.make_dna_db(seq_lines)
    else:
        dna_dir = input_args.dna_dir

    if input_args.contig_bounds is None:
        contig_bounds = pre_blast_process.contig_bounds(seq_lines)
    else:
        contig_bounds = input_args.contig_bounds

    ## BLAST the fragments now
    python_dir_name = os.path.dirname(os.path.realpath(__file__))
    data_dir = re.sub("python", "data", python_dir_name)
    #tmp_dna_dir = "./tmp_dna_lib"

    blast_runs.blast_runs(dna_dir, input_args.query, input_args.output)
    blast_csv = pandas.read_csv(input_args.output,
                                names = ['query','subject','pid',
                                         'align','gap','mismatch',
                                         'qstart','qend','sstart',
                                         'send','eval','bitscore'])
    blast_csv.to_csv(input_args.output, index = False)

    ## Post-BLAST merging and hit idents
    # R_dir = re.sub("python", "R", python_dir_name)
    # post_blast_process.merge_blast_hits(contig_dir=contig_bounds,R_dir=R_dir)
    #
    # out_name = input_args.output + "_hits.csv"
    # post_blast_process.extract_hits(out_name, input_args.no_contigs)
    #
    # ## Check if want to search for comM
    # if input_args.comM:
    #     out_name = input_args.output + "_complete_comM.csv"
    #     post_blast_process.extract_comM(out_name)


    return True

