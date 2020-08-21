import numpy as np
import pandas as pd
import os
import csv
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def consiliate(mydict, fna_record):
    output = {}
    addit = False
    for key in mydict:
        for output_key in output:
            if (mydict[key]['Seqence name'].split(' ')[0] == output[output_key]['Seqence name']
                    and not addit):
                if (abs(mydict[key]['CRISPR start'] - output[output_key]['CRISPR end']) < 10000
                        or abs(mydict[key]['CRISPR end'] - output[output_key]['CRISPR start']) < 10000):
                    output[output_key]['CRISPR start'] = min(output[output_key]['CRISPR start'],
                                                             mydict[key]['CRISPR start'])
                    output[output_key]['CRISPR end'] = max(output[output_key]['CRISPR end'], mydict[key]['CRISPR end'])
                    if mydict[key]['Program'] not in output[output_key]['Program']:
                        output[output_key]['Program'].append(mydict[key]['Program'])
                    if mydict[key]['Repeat sequence'] not in output[output_key]['Repeat sequence']:
                        output[output_key]['Repeat sequence'].append(mydict[key]['Repeat sequence'])
                    output[output_key]['Number of repeats'] = max(output[output_key]['Number of repeats'],
                                                                  mydict[key]['Number of repeats'])
                    #get the sequence
                    seq = fna_record[output[output_key]['Seqence name']].seq
                    if (output[output_key]['CRISPR start'] - 20000) < 0:
                        start = 0
                    else:
                        start = (output[output_key]['CRISPR start'] - 20000)
                    if (output[output_key]['CRISPR end'] + 20000) > len(seq):
                        end = len(seq)
                    else:
                        end = (output[output_key]['CRISPR end'] + 20000)
                    output[output_key]['20K DNA sequence'] = seq[start:end]

                    addit = True

        if not addit:
            CRISPR_num = 'CRISPR_' + str(len(output) + 1)
            output[CRISPR_num] = {}
            output[CRISPR_num]['Seqence name'] = mydict[key]['Seqence name'].split(' ')[0]
            output[CRISPR_num]['CRISPR start'] = mydict[key]['CRISPR start']
            output[CRISPR_num]['CRISPR end'] = mydict[key]['CRISPR end']
            output[CRISPR_num]['Program'] = [mydict[key]['Program']]
            output[CRISPR_num]['Repeat sequence'] = [mydict[key]['Repeat sequence']]
            output[CRISPR_num]['Number of repeats'] = mydict[key]['Number of repeats']
            # get the sequence
            seq = fna_record[output[CRISPR_num]['Seqence name']].seq
            if (output[CRISPR_num]['CRISPR start'] - 20000) < 0:
                start = 0
            else:
                start = (output[CRISPR_num]['CRISPR start'] - 20000)
            if (output[CRISPR_num]['CRISPR end'] + 20000) > len(seq):
                end = len(seq)
            else:
                end = (output[CRISPR_num]['CRISPR end'] + 20000)
            output[CRISPR_num]['20K DNA sequence'] = seq[start:end]
            #print(output[CRISPR_num]['20K DNA sequence'])
        addit = False
    return output

parser = argparse.ArgumentParser(description = "Consolidate CRISPR array and extract sequence information")
parser.add_argument("-i", "-input", type = str, help = "input file relative directory")
parser.add_argument("-o", "-output", type = str, help = "output relative directory")
parser.add_argument("-f", "-fna", type = str, help = "genomic sequence in fasta format")

args = parser.parse_args()


dirname = os.getcwd()
input_relative = args.i
fna_relative = args.f
out_relative = args.o
input_file = os.path.join(dirname, input_relative)
fna_file = os.path.join(dirname, fna_relative)
out_csv_file = os.path.join(dirname, out_relative, "CRISPR_consolidate.csv")
out_fasta_path = os.path.join(dirname, out_relative)
fna_record = SeqIO.index(fna_file, "fasta")

df = pd.read_csv(input_file)
pre_consolidate_dict = df.transpose().to_dict()
summary = {}
summary = consiliate(pre_consolidate_dict, fna_record)

#The out put file should be a panda data frame
dataframe = pd.DataFrame(columns = ["Seqence name", "CRISPR start", "CRISPR end", "Program", "Repeat sequence", "Number of repeats", "20K DNA sequence"])
for key in summary:
    dataframe = dataframe.append(summary[key], ignore_index=True)
outfile_exists = os.path.isfile(out_csv_file)
with open(out_csv_file, "a") as csvfile:
    header = ["Seqence name", "CRISPR start", "CRISPR end", "Program", "Repeat sequence", "Number of repeats", "20K DNA sequence"]
    if not outfile_exists:
        headwriter = csv.DictWriter(csvfile, fieldnames=header)
        headwriter.writeheader()
    dataframe.to_csv(csvfile, mode = "a", header = False, index = None)


#Output sequence as fasta
for key in summary:
    if len(summary[key]["20K DNA sequence"]) > 5000:
        record= SeqRecord(seq = summary[key]["20K DNA sequence"], id = key, description = summary[key]["Seqence name"])
        with open(out_fasta_path + "/" + key + "_sequence.fasta", "a") as output_handle:
            SeqIO.write(record, output_handle, "fasta")


