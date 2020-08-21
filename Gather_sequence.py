import os
import pandas as pd
import csv
import argparse

def get_minCED(data):
    out_data = []
    sequence_name = ""
    n = False
    location = []
    repeat_sequence = ""
    number_of_repeats = 0
    for line in data:
        if line.startswith("Sequence", 0, -1):
            sequence_name = line.split("'")[1]
        if line.startswith("CRISPR", 0, -1):
            location = line.rstrip("\n").split(" ")
        if n:
            if not line.startswith("Repeats", 0, -1):
                repeat_sequence = line.split("	")[2]
            else:
                number_of_repeats = line[(line.find("Repeats: ") + len("Repeats: ")):(line.find("	Average "))]
                out_data.append([sequence_name, location[-3], location[-1], "minCED", repeat_sequence, number_of_repeats])
                location = []
                repeat_sequence = ""
                number_of_repeats = 0
            n = False
        if line.startswith("--", 0, -1):
            n = True
    return out_data


def get_piler_cr(data):
    out_data = []
    sequence_name = ""
    n = False
    m = False
    for line in data:
        if line == "\n" and m:
            m = False
        if line.startswith("SUMMARY BY POSITION", 0, -1):
            n = True
        if n:
            if line.startswith(">", 0, -1):
                sequence_name = str(line.rstrip("\n").split(">")[-1])
        if n and m:
            sum_data = list(filter(None, line.rstrip("\n").split("  ")))
            out_data.append([sequence_name, sum_data[2], str(int(sum_data[2]) + int(sum_data[3])), "piler_cr", sum_data[-1], sum_data[4]])
        if line.startswith("===", 0, -1):
            m = True
    return out_data


def get_crispr_detect(data):
    out_data = []
    sequence_name = ""
    location = ""
    n = False
    m = False
    for line in data:
        if line.startswith("Array", 0, -1):
            location = line.split(" ")[2]
        if line.startswith(">", 0, -1):
            sequence_name = line.split("		")[0][1:]
        if n and m:
            sum_data = list(filter(None, line.split(" ")))
            out_data.append([sequence_name, location.split("-")[0], location.split("-")[1], "CRISPR_Detect",
                             sum_data[3].split("\t")[1], sum_data[0].rstrip("\t")])
            n = False
            m = False
        if n and line.startswith("=======", 0, -1):
            m = True
        if line.startswith("=======", 0, -1):
            n = True
    return out_data

def get_crisprfinder(data):
    out_data = []
    for line in data:
        if line != "\n" and not line.startswith("Strain", 0, -1):
            sum_data = line.split("\t")
            out_data.append([sum_data[1], sum_data[5], sum_data[6], "CRISPRfinder", sum_data[10], sum_data[14]])
    return out_data

#We will need input file directory, out put file directory, file type

parser = argparse.ArgumentParser(description = "Summarize different CRISPR array identification")
parser.add_argument("-i", "-input", type = str, help = "input file relative directory")
parser.add_argument("-o", "-output", type = str, help = "output file relative directory")
parser.add_argument("-t", "-type", type = str, help = "One of the following file type: minCED, piler_cr, CRISPR_Detect, CRISPRfinder")

args = parser.parse_args()
input_relative = args.i
output_relative = args.o
file_type = args.t

dirname = os.getcwd()

input_file = os.path.join(dirname, input_relative)
output_file = os.path.join(dirname, output_relative)
assert os.path.exists(input_file), "No file at " + str(input_file)
with open(input_file, 'r+') as f:
    data = f.readlines()

summary = []

if file_type == "minCED":
    summary = get_minCED(data)
elif file_type == "piler_cr":
    summary = get_piler_cr(data)
elif file_type == "CRISPR_Detect":
    summary = get_crispr_detect(data)
elif file_type == "CRISPRfinder":
    summary = get_crisprfinder(data)
else:
    print("Unexpected file type")



#The out put file should be a panda data frame
dataframe = pd.DataFrame(summary, columns = ["Seqence name", "CRISPR start", "CRISPR end", "Program", "Repeat sequence", "Number of repeats"])

outfile_exists = os.path.isfile(output_file)

with open(output_file, "a") as csvfile:
    header = ["Seqence name", "CRISPR start", "CRISPR end", "Program", "Repeat sequence", "Number of repeats"]
    if not outfile_exists:
        headwriter = csv.DictWriter(csvfile, fieldnames=header)
        headwriter.writeheader()
    dataframe.to_csv(csvfile, mode = "a", header = False, index = None)

#with open(output_file, "a") as csvfile:
#    csvfile.write("\n")







