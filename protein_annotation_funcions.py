import re

def read_cdd_annotation(cdd_data, bitscore_cutoff):
    #cdd output dictionary {ORF1:{cdd1:{name:bit-score}, other_cdd:{name2:bit-score2, name3:bit-score3}}
    in_region = False
    cdd_first = True
    output = {}
    name = ""
    for line in cdd_data:
        if line.startswith("QUERY", 0, -1):
            name = line.strip("\n").split("\t")[4]
        elif line.startswith("DOMAINS", 0, -1):
            in_region = True
        elif line.startswith("ENDDOMAINS", 0, -1):
            in_region = False
            cdd_first = True
        elif in_region == True and cdd_first == True:
            output[name] = {}
            line = line.strip("\n").split("\t")
            output[name]["cdd1"] = {}
            if float(line[7]) > bitscore_cutoff:
                output[name]["cdd1"][line[9]] = line[7]
                cdd_first = False
            output[name]["other_cdd"] = {}
        elif in_region == True and cdd_first == False:
            line = line.strip("\n").split("\t")
            if float(line[7]) > bitscore_cutoff:
                output[name]["other_cdd"][line[9]] = line[7]
    return output

def read_crispr_annotation(crispr_data, index, pfam_data, bitscore_cutoff):
    # index is dictionary with number corresponding to the protein type
    # crispr output dictionary {ORF1:{cdd1:{name:bit-score}, other_cdd:{name2:bit-score2, name3:bit-score3}}
    in_region = False
    cdd_first = True
    output = {}
    name = ""
    for line in crispr_data:
        if line.startswith("# Query: ", 0, -1):
            name = line.strip("\n").split(": ")[-1]
            continue
        if line.startswith("hits found", -11):
            count = int(line.split(" ")[1])
            if count > 0:
                in_region = True
                continue
            else:
                False
                continue
        if line.startswith("# RPSBLAST", 0, -1) or line.startswith("# BLAST processed", 0, -1):
            in_region = False
            cdd_first = True
        if in_region == True and cdd_first == True:
            output[name] = {}
            line = line.strip("\n").split("\t")
            output[name]["cdd1"] = {}
            #print(line)
            if float(line[-1]) > bitscore_cutoff:
                if index[line[1]] in pfam_data:
                    output[name]["cdd1"][pfam_data[index[line[1]]]] = line[-1]
                else:
                    output[name]["cdd1"][index[line[1]]] = line[-1]
                cdd_first = False
            output[name]["other_cdd"] = {}
            continue
        if in_region == True and cdd_first == False:
            line = line.strip("\n").split("\t")
            if float(line[-1]) > bitscore_cutoff:
                if index[line[1]] in pfam_data:
                    output[name]["other_cdd"][pfam_data[index[line[1]]]] = line[-1]
                else:
                    output[name]["other_cdd"][index[line[1]]] = line[-1]
                continue
    return output

def protein_to_crispr_type(protein_name, current_type):
    #print(protein_name)
    if re.search("Cas3", protein_name, re.IGNORECASE):
        if current_type == "" or current_type == "Type I":
            current_type = "Type I"
        elif current_type.find("Type I") == -1:
            current_type = current_type + "_" + "Type I"
    if re.search("Cas9", protein_name, re.IGNORECASE):
        if current_type == "" or current_type == "Type II":
            current_type = "Type II"
        elif current_type.find("Type II") == -1:
            current_type = current_type + "_" + "Type II"
    if re.search("Cas10", protein_name, re.IGNORECASE):
        if current_type == "" or current_type == "Type III":
            current_type = "Type III"
        elif current_type.find("Type III") == -1:
            current_type = current_type + "_" + "Type III"
    if re.search("Csf1", protein_name, re.IGNORECASE):
        if current_type == "" or current_type == "Type IV":
            current_type = "Type IV"
        elif current_type.find("Type IV") == -1:
            current_type = current_type + "_" + "Type IV"
    if re.search("Cpf1", protein_name, re.IGNORECASE) or re.search("Cas12", protein_name, re.IGNORECASE) or re.search("Cas14", protein_name, re.IGNORECASE):
        if current_type == "" or current_type == "Type V":
            current_type = "Type V"
        elif current_type.find("Type V") == -1:
            current_type = current_type + "_" + "Type V"
    if re.search("Cas13", protein_name, re.IGNORECASE):
        if current_type == "" or current_type == "Type VI":
            current_type = "Type VI"
        elif current_type.find("Type VI") == -1:
            current_type = current_type + "_" + "Type VI"
    return current_type


def crispr_classification(cdd_dict, crispr_dict, orf_record, orf_fasta):
    is_crispr = False
    crispr_type = ""
    protein_list = {}
    other_protein = {}
    record = orf_record[:]
    cas_protein_names = ["Cas", "Csm", "CRISPR", "Csn", "Cpf", "Csa", "Cse", "Csy", "Csm", "Csx", "Cmr", "Csx", "Csf"]
    if crispr_dict != {}:
        for annotation in crispr_dict:
            if crispr_dict[annotation]['cdd1'] != {}:
                is_crispr = True
    if is_crispr:
        for orf in crispr_dict:
            for annotation in crispr_dict[orf]['cdd1']:
                #print(annotation)
                crispr_type = protein_to_crispr_type(annotation, crispr_type)
                protein_list[annotation] = [str(orf_fasta[orf.split(" ")[0]].seq), crispr_dict[orf]['cdd1']]
                record.remove(orf)
    for orf in cdd_dict:
        for annotation in cdd_dict[orf]['cdd1']:
            #print(annotation)
            for name in cas_protein_names:
                if re.search(name, annotation, re.IGNORECASE):

                    is_crispr = True
                    crispr_type = protein_to_crispr_type(annotation, crispr_type)
                    if orf in record:
                        protein_list[annotation] = [str(orf_fasta[orf.split(" ")[0]].seq), cdd_dict[orf]['cdd1']]
                        record.remove(orf)
    if is_crispr:
        for name in record:
            other_protein[name] = [str(orf_fasta[name.split(" ")[0]].seq), [crispr_dict[name]['cdd1'] if name in crispr_dict else "", cdd_dict[name]['cdd1'] if name in cdd_dict else ""]]
    if not is_crispr:
        for name in record:
            other_protein[name] = [str(orf_fasta[name.split(" ")[0]].seq), [crispr_dict[name]['cdd1'] if name in crispr_dict else "", cdd_dict[name]['cdd1'] if name in cdd_dict else ""]]
    return is_crispr, crispr_type, protein_list, other_protein



