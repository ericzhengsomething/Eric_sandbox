import os
import argparse
from Bio import SeqIO
from protein_annotation_funcions import *
import pandas as pd
import csv

#pd.set_option('display.max_columns', None)

#def output_script(cdd_file, crispr_file, bit_score, index_file, orf_file, current_genome_summary_file, protein_summary_dic, subtype_summary_location):

parser = argparse.ArgumentParser(description = "Summary CRISPR array and the pritein samples.")
parser.add_argument("--cdd_file", type = str, help = "input cdd annotation")
parser.add_argument("--crispr_file", type = str, help = "input crispr annotation")
parser.add_argument("--index_file", type = str, help = "input crispr index file")
parser.add_argument("--pfam_tab_file", type = str, help = "input pfam tab file")
parser.add_argument("--fasta_file", type = str, help = "input orf fasta file")
parser.add_argument("--output_folder", type = str, help = "output directory")
parser.add_argument("--bitscore_cutoff", type = int, help = "bitscore cutoff")
parser.add_argument("--protein_bucket", type = int, help = "size of bucket for protein grouping")



args = parser.parse_args()


cdd_process_domain_file = args.cdd_file
crispr_domain_file = args.crispr_file
index_file = args.index_file
pfam_file = args.pfam_tab_file
orf_fasta_file = args.fasta_file
summary_output_folder = args.output_folder
bit_score_cutoff = args.bitscore_cutoff
phylon = crispr_domain_file.split("/")[-5]
species = crispr_domain_file.split("/")[-4]
crispr_number = "CRISPR_" + crispr_domain_file.split("/")[-1].strip("_CRISPRdomain.asn")
assembly_id = crispr_domain_file.split("/")[-2]
protein_bucket = args.protein_bucket
is_crispr = False
crispr_type = ""
protein_info = {}
other_protein = {}
#Data is a dictionary with each ORF is an entry
#{ORF name:, CDD annotation 1: CRISPR annotation 1: CDD other annotation:, CRISPR_other_annotation:}
summary_dict = {}
system_summary = {}
system_summary[assembly_id] = {}
system_summary_out = ""


with open(cdd_process_domain_file, 'r+') as f:
    cdd_data = f.readlines()
with open(crispr_domain_file, 'r+') as f:
    crispr_data = f.readlines()
with open(index_file, 'r+') as f:
    file = f.readlines()
    index = {}
    for item in file:
        item = item.strip("\n").split(",")
        index[item[0]] = item[1]
with open(pfam_file, 'r+') as f:
    file = f.readlines()
    pfam_data = {}
    for item in file:
        item = item.strip("\n").split("\t")
        pfam_data[str(item[0])] = item[1] + '_' + item[2] + '_' + item[3]

orf_fasta_data = SeqIO.parse(orf_fasta_file, "fasta")
orf_record = []
for record in orf_fasta_data:
    orf_record.append(record.description)
orf_fasta_data = SeqIO.to_dict(SeqIO.parse(orf_fasta_file, "fasta"))
#print(orf_record)
#print(orf_fasta_data)

cdd_data = read_cdd_annotation(cdd_data, bit_score_cutoff)
crispr_data = read_crispr_annotation(crispr_data, index, pfam_data, bit_score_cutoff)
is_crispr, crispr_type, protein_info, other_protein = crispr_classification(cdd_data, crispr_data, orf_record, orf_fasta_data)
#print(is_crispr, crispr_type)
#print(cdd_data)
#print(crispr_data)


for orf in orf_record:
    summary_dict[orf] = {}
    summary_dict[orf]["ORF name"] = orf
    summary_dict[orf]["CDD annotation 1 and bit score"] = str(cdd_data[orf]["cdd1"] if orf in cdd_data else "")
    summary_dict[orf]["CRISPR annotation 1 and bit score"] = str(crispr_data[orf]["cdd1"] if orf in crispr_data else "")
    summary_dict[orf]["CDD other annotation"] = str(cdd_data[orf]["other_cdd"] if orf in cdd_data else "")
    summary_dict[orf]["CRISPR other annotation"] = str(crispr_data[orf]["other_cdd"] if orf in crispr_data else "")
    summary_dict[orf]["Protein seqeunce"] = str(orf_fasta_data[orf.split(" ")[0]].seq)

#print(summary_dict)
summary_data_frame = pd.DataFrame(summary_dict).transpose()
summary_data_frame["Assemble ID"] = assembly_id
summary_data_frame["CRISPR number"] = crispr_number
summary_data_frame["Species"] = species
summary_data_frame["Potential type"] = crispr_type
summary_data_frame["Is CRISPR"] = is_crispr


##Output the annotation information for each array
current_system_out = "/".join(crispr_domain_file.split("/")[0:-1]) + "/" + assembly_id + "_" + crispr_number +".csv"
outfile_exists = os.path.isfile(current_system_out)
with open(current_system_out, "a") as csvfile:
    header = ["", "CDD annotation 1 and bit score", "CDD other annotation", "CRISPR annotation 1 and bit score", "CRISPR other annotation", "ORF name", "Protein sequence", "Assemble ID", "CRISPR number", "Species", "Potential type", "Is CRISPR"]
    if not outfile_exists:
        headwriter = csv.DictWriter(csvfile, fieldnames=header)
        headwriter.writeheader()
    summary_data_frame.to_csv(csvfile, mode = "a", header = False)

#Summary of CRISPR array
if is_crispr:
    ##Append the Cas protein file both cvs sheet and fasta
    for protein in protein_info:
        if not os.path.exists(summary_output_folder + "crispr_array_crispr_protein_summary/"):
            os.makedirs(summary_output_folder + "crispr_array_crispr_protein_summary/")
        current_system_out = summary_output_folder + "crispr_array_crispr_protein_summary/" + protein + ".csv"
        outfile_exists = os.path.isfile(current_system_out)
        with open(current_system_out, "a") as csvfile:
            header = ["Species", "Assemble ID", "CRISPR number", "Sequence", "Annotation"]
            if not outfile_exists:
                headwriter = csv.DictWriter(csvfile, fieldnames=header)
                headwriter.writeheader()
            writer = csv.writer(csvfile)
            writer.writerow([species, assembly_id, crispr_number, protein_info[protein][0], protein_info[protein][1]])
        current_system_out = summary_output_folder + "crispr_array_crispr_protein_summary/" + protein + ".fasta"
        with open(current_system_out, "a") as fastafile:
            fastafile.write(">" + assembly_id + "_" + crispr_number + "_" + protein + "\n" + protein_info[protein][0] + "\n")

    ##Append other proteins to other protein file and fasta
    for protein in other_protein:
        size_range = len(other_protein[protein]) // protein_bucket * protein_bucket
        if not os.path.exists(summary_output_folder + "crispr_array_other_protein_summary/"):
            os.makedirs(summary_output_folder + "crispr_array_other_protein_summary/")
        current_system_out = summary_output_folder + "crispr_array_other_protein_summary/" + str(size_range) + "-" + str(size_range + protein_bucket) + ".csv"
        outfile_exists = os.path.isfile(current_system_out)
        with open(current_system_out, "a") as csvfile:
            header = ["Species", "Assemble ID", "CRISPR number", "Protein ID", "Sequence", "Annotation"]
            if not outfile_exists:
                headwriter = csv.DictWriter(csvfile, fieldnames=header)
                headwriter.writeheader()
            writer = csv.writer(csvfile)
            writer.writerow([species, assembly_id, crispr_number, protein, other_protein[protein][0], other_protein[protein][1]])
        current_system_out = summary_output_folder + "crispr_array_other_protein_summary/" + str(size_range) + "-" + str(size_range + protein_bucket) + ".fasta"
        with open(current_system_out, "a") as fastafile:
            fastafile.write(">" + assembly_id + "_" + crispr_number + "_" + protein + "\n" + other_protein[protein][0] + "\n")

    ##CRISPR array summary file, one file has all the information for every CRISPR array
    system_summary[assembly_id]["CRISPR_number"] = crispr_number
    system_summary[assembly_id]["Phylon"] = phylon
    system_summary[assembly_id]["Species"] = species
    system_summary[assembly_id]["Potential type"] = crispr_type
    system_summary[assembly_id]["Is CRISPR"] = is_crispr

    #print(system_summary)
    system_summary_df = pd.DataFrame(system_summary).transpose()
    #print(system_summary_df)
    system_summary_out = summary_output_folder + "CRISPR_array_summary.csv"
    outfile_exists = os.path.isfile(system_summary_out)
    with open(system_summary_out, "a") as csvfile:
        header = ["", "CRISPR_number", "Is CRISPR", "Phylon", "Potential type", "Species"]
        if not outfile_exists:
            headwriter = csv.DictWriter(csvfile, fieldnames=header)
            headwriter.writeheader()
        system_summary_df.to_csv(csvfile, mode = "a", header = False)

# Summary of none CRISPR array
if not is_crispr:
    ##Append the Cas protein file both cvs sheet and fasta
    for protein in protein_info:
        if not os.path.exists(summary_output_folder + "not_crispr_array_crispr_protein_summary/"):
            os.makedirs(summary_output_folder + "not_crispr_array_crispr_protein_summary/")
        current_system_out = summary_output_folder + "not_crispr_array_crispr_protein_summary/" + protein + ".csv"
        outfile_exists = os.path.isfile(current_system_out)
        with open(current_system_out, "a") as csvfile:
            header = ["Species", "Assemble ID", "CRISPR number", "Sequence", "Annotation"]
            if not outfile_exists:
                headwriter = csv.DictWriter(csvfile, fieldnames=header)
                headwriter.writeheader()
            writer = csv.writer(csvfile)
            writer.writerow([species, assembly_id, crispr_number, protein_info[protein][0], protein_info[protein][1]])
        current_system_out = summary_output_folder + "not_crispr_array_crispr_protein_summary/" + protein + ".fasta"
        with open(current_system_out, "a") as fastafile:
            fastafile.write(">" + assembly_id + "_" + crispr_number + "_" + protein + "\n" + protein_info[protein][0] + "\n")

    ##Append other proteins to other protein file and fasta
    for protein in other_protein:
        size_range = len(other_protein[protein]) // protein_bucket * protein_bucket
        if not os.path.exists(summary_output_folder + "not_crispr_array_other_protein_summary/"):
            os.makedirs(summary_output_folder + "not_crispr_array_other_protein_summary/")
        current_system_out = summary_output_folder + "not_crispr_array_other_protein_summary/" + str(
            size_range) + "-" + str(size_range + protein_bucket) + ".csv"
        outfile_exists = os.path.isfile(current_system_out)
        with open(current_system_out, "a") as csvfile:
            header = ["Species", "Assemble ID", "CRISPR number", "Protein ID", "Sequence", "Annotation"]
            if not outfile_exists:
                headwriter = csv.DictWriter(csvfile, fieldnames=header)
                headwriter.writeheader()
            writer = csv.writer(csvfile)
            writer.writerow([species, assembly_id, crispr_number, protein, other_protein[protein][0], other_protein[protein][1]])
        current_system_out = summary_output_folder + "not_crispr_array_other_protein_summary/" + str(
            size_range) + "-" + str(size_range + protein_bucket) + ".fasta"
        with open(current_system_out, "a") as fastafile:
            fastafile.write(">" + assembly_id + "_" + crispr_number + "_" + protein + "\n" + other_protein[protein][0] + "\n")

    ##CRISPR array summary file, one file has all the information for every CRISPR array
    system_summary[assembly_id]["CRISPR_number"] = crispr_number
    system_summary[assembly_id]["Phylon"] = phylon
    system_summary[assembly_id]["Species"] = species
    system_summary[assembly_id]["Potential type"] = crispr_type
    system_summary[assembly_id]["Is CRISPR"] = is_crispr

    # print(system_summary)
    system_summary_df = pd.DataFrame(system_summary).transpose()
    # print(system_summary_df)
    system_summary_out = summary_output_folder + "Not_CRISPR_array_summary.csv"
    outfile_exists = os.path.isfile(system_summary_out)
    with open(system_summary_out, "a") as csvfile:
        header = ["", "CRISPR_number", "Is CRISPR", "Phylon", "Potential type", "Species"]
        if not outfile_exists:
            headwriter = csv.DictWriter(csvfile, fieldnames=header)
            headwriter.writeheader()
        system_summary_df.to_csv(csvfile, mode="a", header=False)




