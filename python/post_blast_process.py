import subprocess
import os
import pandas
pandas.set_option('display.max_columns', None)

def merge_blast_hits(contig_dir, R_dir):
    """ Function to run the R scripts for merging the conserved region hits
        directly after the BLAST search"""

    merge_cmd = "Rscript --vanilla " + R_dir + "/testing_out_functions.R left_end_blast.csv right_end_blast.csv " + contig_dir

    try:
        subprocess.check_output(merge_cmd, shell=True)
    except:
        subprocess.SubprocessError

def extract_hits(out_name, contig_use):
    """ Function to take the merged blast dbs and extracate the hits
        using the criteria that a left and right end must be on the
        same contig and witin 150k of each other on that contig"""

    left_end_csv = pandas.read_csv(filepath_or_buffer="./left_end_merged.csv",header=0)
    right_end_csv = pandas.read_csv(filepath_or_buffer="./right_end_merged.csv",header=0)

    ## run through unique values in left end (if they're not in left end they're not
    ## a hit)

    left_end_vals = left_end_csv.subject.unique()

    out_df = pandas.DataFrame(columns=['id','hit_start','hit_end','ori','contig', 'index'])
    for iso_index,subject in enumerate(left_end_vals):

        current_left_row = left_end_csv[left_end_csv['subject'] == subject]
        current_left_row = current_left_row.sort_values(by = 'sstart', ascending=True)
        current_right_row = right_end_csv[right_end_csv['subject'] == subject]
        if current_right_row.empty:
            continue
        for index, row in current_left_row.iterrows():
            if contig_use:
                right_match = current_right_row[(current_right_row['ori'] == current_left_row.loc[index,'ori']) & \
                                                (current_right_row['contig'] == row.loc['contig'])]
                if right_match.empty:
                    continue
                if row.loc['ori'] == "forward":
                    right_hits = right_match[(right_match['sstart'] <= (row.loc['send'] + 150000)) &
                                             right_match['sstart'] >= (row.loc['send'] - 50)]
                    if right_hits.empty:
                        continue
                    right_hits = right_hits.sort_values(by='sstart', ascending=False)
                    right_hits = right_hits.reset_index(drop = True)
                elif row.loc['ori'] == "reverse":

                    right_hits = right_match[(right_match['sstart'] >= (row.loc['send'] - 150000)) &
                                             (right_match['sstart'] <= (row.loc['send'] + 50))]
                    if right_hits.empty:
                        continue
                    right_hits =  right_hits.sort_values(by='sstart', ascending=True)
                    right_hits = right_hits.reset_index(drop=True)

                new_row = {'id': row.loc['subject'],
                           'hit_start':row.loc['sstart'],
                           'hit_end':right_hits.loc[0, "send"],
                           'ori':row.loc['ori'],
                           'contig':row.loc['contig'],
                           'index':(str(iso_index) + "_" + str(index))}

                if row.loc['ori'] == "forward":
                    other_longer_hits = out_df[(out_df['id'] == row.loc['subject']) &
                                               (out_df['contig'] == row.loc['contig']) &
                                               (out_df['ori'] == "forward") &
                                               (out_df['hit_start'] < new_row["hit_start"]) &
                                               (out_df['hit_end'] >= new_row["hit_end"])]
                    if other_longer_hits.empty:
                        ## Check for any hits that will be nested within this new one
                        other_shorter_hits = out_df[(out_df['id'] == row.loc['subject']) &
                                                    (out_df['contig'] == row.loc['contig']) &
                                                    (out_df['ori'] == "forward") &
                                                    (out_df['hit_start'] > new_row['hit_start'] )&
                                                    (out_df['hit_end'] <= new_row['hit_end'])]
                        if not other_shorter_hits.empty:
                            ## get their index and remove them
                            shorter_indies = other_shorter_hits['index'].to_list()
                            ## remove these from the out_df
                            indy_rows = out_df[out_df['index'].isin(shorter_indies)].index
                            out_df = out_df.drop(indy_rows)
                            out_df = out_df.reset_index(drop=True)
                        out_df = out_df.append(new_row, ignore_index=True)
                elif row.loc['ori'] == "reverse":
                    other_longer_hits = out_df[(out_df['id'] == row.loc['subject']) &
                                               (out_df['contig'] == row.loc['contig']) &
                                               (out_df['ori'] == "reverse") &
                                               (out_df['hit_start'] >= new_row["hit_start"]) &
                                               (out_df['hit_end'] < new_row["hit_end"])]
                    if other_longer_hits.empty:
                        ## Check for any hits that will be nested within this new one
                        other_shorter_hits = out_df[(out_df['id'] == row.loc['subject']) &
                                                    (out_df['contig'] == row.loc['contig']) &
                                                    (out_df['ori'] == "reverse") &
                                                    (out_df['hit_start'] <= new_row['hit_start']) &
                                                    (out_df['hit_end'] >= new_row['hit_end'])]
                        if not other_shorter_hits.empty:
                            ## get their index and remove them
                            shorter_indies = other_shorter_hits['index'].to_list()
                            ## remove these from the out_df
                            indy_rows = out_df[out_df['index'].isin(shorter_indies)].index
                            out_df = out_df.drop(indy_rows)
                            out_df = out_df.reset_index(drop=True)
                        out_df = out_df.append(new_row, ignore_index=True)



            else:
                right_match = current_right_row[(current_right_row['ori'] == current_left_row.loc[index, 'ori'])]
                if right_match.empty:
                    continue
                if row.loc['ori'] == "forward":
                    right_hits = right_match[right_match['sstart'] <= (row.loc['send'] + 150000)]
                    if right_hits.empty:
                        continue
                    right_hits = right_hits.sort_values(by='sstart', ascending=True)
                    right_hits = right_hits.reset_index(drop=True)
                elif row.loc['ori'] == "reverse":
                    right_hits = right_match[right_match['sstart'] >= (row.loc['send'] - 150000)]
                    if right_hits.empty:
                        continue
                    right_hits = right_hits.sort_values(by='sstart', ascending=False)
                    right_hits = right_hits.reset_index(drop=True)

                new_row = {'id': row.loc['subject'],
                           'hit_start': row.loc['sstart'],
                           'hit_end': right_hits.loc[0, "send"],
                           'ori': row.loc['ori'],
                           'contig': row.loc['contig'],
                           'index': (str(iso_index) + "_" + str(index))}



    out_df.to_csv(out_name, index=False)

def extract_comM(out_name):
    '''Function to check if the comM hits in an isolate are intact '''

    blast_cols = ['query', 'subject', 'pid', 'align', 'gap', 'mismatch',
                  'qstart', 'qend', 'sstart', 'send', 'eval', 'bitscore']

    comM_hits = pandas.read_csv("./comM_hits.csv", names=blast_cols, header=None)

    comM_hitters = comM_hits[comM_hits['align'] >= 1450]

    comM_hitters.to_csv(out_name, index=False)






