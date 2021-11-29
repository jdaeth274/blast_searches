# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import argparse
import sys
import time

from python.common import main

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.

def parse_args():
    purpose = ''' This is a scipt to blast a sequence against a collection of fasta sequences Usage:
        python blast-runner.py --seqs <list_of_fastas> --query <fasta_file_of_query> --output <output_prefixes> --threads <num_cores_to_use>'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='blast-runner.py')

    parser.add_argument('--seqs', required=True, help='List of seqeuence files (FASTA) (required)', type=str)
    parser.add_argument('--query', required=True, help='Sequence of query for BLAST (FASTA) (required)', type=str)
    parser.add_argument('--output', required=True, help='Prefix of output files  (required)', type=str)
    parser.add_argument('--dna-dir', default=None, help= 'Location of directory of DNA files with output.mfa file present',
                        type=str)
    parser.add_argument('--contig-bounds', default=None, help="Location of the contig bounds of files", type=str)
    parser.add_argument('--threads', default=1, help='Number of threads to use for ORF finder', type=int)
    parser.add_argument('--no-contigs', default=True, action='store_false', help='Do not use contig bounds when defining a hit [default = False]')

    args = parser.parse_args()

    return args
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    start = time.perf_counter()
    if main(parse_args()):
        end = time.perf_counter()
        print("Took this long to complete: %s (s)" % (end - start))
    sys.exit("Done")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
