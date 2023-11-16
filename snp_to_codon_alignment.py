import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO
import re

def make_info_df(snp_into_xslx,sub_table):
    return pd.read_excel(snp_into_xslx,sub_table,header=1, index_col=0)
    

def headers_from_fasta(fasta_file):
    headers = []
    
    with open(fasta_file, 'r') as fasta_file:
        current_header = None
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                # If a line starts with '>', it's a header line
                current_header = line
                headers.append(current_header)
    
    return [line.replace('>', '') for line in headers]

def read_annotation_file(annotation_file):
    with open(annotation_file, 'r') as file:
        file_contents = file.readlines()
    return file_contents

def make_annotation_dict(annotation_file_list):
    Annotation_dict={}
    for i in annotation_file_list:
        match_geneid = re.search(r'gene_id "([A-Za-z0-9_]+)"',i)
        match_location = re.search(r'\tgene\t(\d+)\t(\d+)\t',i)
        match_direction = re.search(r'(?<=\t\.\t)([\-\+])(?=\t\.)',i)
        if match_geneid and match_location and match_direction:
            Annotation_dict[match_geneid.group(1)] = [int(match_location.group(1)),int(match_location.group(2)),match_direction.group(1),0,0]
    return Annotation_dict

def get_sequence_lengths(fasta_file,header):
    #print(header)
    #for record in SeqIO.parse(fasta_file, "fasta"):
    #    if record.id == header[0:len(record.id)]:
    #        seq_len = len(record.seq)
    return 4653728

def get_ORF(Annotation_dict,alignment_file,header):
    for key in Annotation_dict:
        if Annotation_dict[key][2] == "+":
            orf = Annotation_dict[key][0] % 3 + 1
            Annotation_dict[key][3]=orf
            
        elif Annotation_dict[key][2] == "-":
            orf = -1*(((get_sequence_lengths(alignment_file,header)-Annotation_dict[key][1]))%3+1)
            Annotation_dict[key][3]=orf
        #print(key, Annotation_dict[key])
    return Annotation_dict

def codon_finder( snp_df, Annotation_dict):
    for key, value in Annotation_dict.items():
        snp_indices = snp_df[snp_df.index.to_series().between(value[0], value[1])].index.tolist()
        value[4] = snp_indices
    return Annotation_dict

def find_bases(Annotation_dict):
    codons=[]
    for key,value in Annotation_dict.items():
        if len(value[-1])>0:
            # + reading frame
            if value[2] == "+":
                for i in value[-1]:
                    position_in_codon = (i-value[0])%3
                    if position_in_codon ==0:
                        codons.append([i,i+1,i+2])
                    if position_in_codon ==1:
                        codons.append([i-1,i,i+1])
                    if position_in_codon ==2:
                        codons.append([i-2,i-1,i])
            if value[2] == "-":
                for i in value[-1]:
                    position_in_codon = (i-value[1])%3
                    if position_in_codon ==0:
                        codons.append([i,i+1,i+2])
                    if position_in_codon ==1:
                        codons.append([i-1,i,i+1])
                    if position_in_codon ==2:
                        codons.append([i-2,i-1,i])
    return(codons)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find codons corresponding to snp in snp alignment file")
    parser.add_argument(
        "--snp_align_path",
        required=True,
        help="Path to input containing snp alignment .fasta file"
    )
    parser.add_argument(
        "--snp_info_path",
        required=True,
        help="Path to input containing snp info .xlsx file"
    )
    parser.add_argument(
        "--full_align_path",
        required=True,
        help="Path to full_alignment file"
    )

    parser.add_argument(
        "--annotation_path",
        required=True,
        help="Path to input directory containing annotation .gtf file"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Path to directory for output"
    )

    args = parser.parse_args()
    snp_align_path = args.snp_align_path
    full_align_path = args.full_align_path
    snp_info_path = args.snp_info_path
    annotation_path = args.annotation_path
    out_dir = args.outdir
    info_df = make_info_df(snp_info_path,"Supp. table 14" )
    headers = headers_from_fasta(snp_align_path)
    annotation_dict = make_annotation_dict(read_annotation_file(annotation_path))
    get_sequence_lengths(full_align_path,headers)
    annotation_dict=get_ORF(annotation_dict,full_align_path,headers)
    annotation_dict = codon_finder(info_df,annotation_dict)
    codon = find_bases(annotation_dict)
    print(codon)
    full_align_sequence = SeqIO.parse(full_align_path, 'fasta')
    # add refererence sequence so you dont have to iterate through all
    snp_align_sequence = SeqIO.parse(snp_align_path, 'fasta')
    extracted_sequences={}

    for seq_record in snp_align_sequence:# iterate throuh snp
        header = seq_record.id 
        print(header)
        if header in list(info_df.columns.to_list()):
            sequence = str(seq_record.seq)
            print(sequence)
        # Extract sequences based on the provided locations
            #print(codon)
            extracted_sequences[header] = ''.join(sequence[i - 1] for loc in codon for i in loc) #add base from ref to new codon file
    os.makedirs(args.outdir, exist_ok=True)
    try:
        file_name = "{}{}_codon.fasta".format(args.outdir,os.path.splitext(os.path.basename(annotation_path))[0])
        f_out = open(file_name, 'w')
    except IOError:
        print("Output file {} cannot be created".format(file_name))
        sys.exit(1)
    # Write the extracted sequences to a new FASTA file
    for header, sequence in extracted_sequences.items():
        print(type(header),type(sequence))
        f_out.write(">{header}\n{sequence}\n".format(header,sequence))
    f_out.close()

    '''
        for header in find_headers_in_fasta(data_dir):
            f_out.write('{}:{}\n'.format(header,''))
        f_out.close()
        '''