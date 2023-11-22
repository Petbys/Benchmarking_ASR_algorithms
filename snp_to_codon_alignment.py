import os
import sys
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
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

def get_sequence_lengths(sequence): # sequence parser
    #print(header)
    #for record in SeqIO.parse(fasta_file, "fasta"):
    #    if record.id == header[0:len(record.id)]:
    #        seq_len = len(record.seq)
    return len(record.seq)
#4653728
def get_ORF(Annotation_dict,alignment_file,header):
    for key in Annotation_dict:
        if Annotation_dict[key][2] == "+":
            orf = Annotation_dict[key][0] % 3 + 1
            Annotation_dict[key][3]=orf
            
        elif Annotation_dict[key][2] == "-":
            #orf = -1*(((get_sequence_lengths(alignment_file,header)-Annotation_dict[key][1]))%3+1)
            orf = -1*(((4653728-Annotation_dict[key][1]))%3+1) # find the complement 
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
            # forward reading frame
            if value[2] == "+":
                for i in value[-1]:
                    position_in_codon = (i-value[0])%3
                    if position_in_codon ==0:
                        codons.append([i,i+1,i+2,position_in_codon,value[2]])
                    if position_in_codon ==1:
                        codons.append([i-1,i,i+1,position_in_codon,value[2]])
                    if position_in_codon ==2:
                        codons.append([i-2,i-1,i,position_in_codon,value[2]])

            # reverse reading frame
            if value[2] == "-":
                for i in value[-1]:
                    position_in_codon = (value[1]-i)%3
                    if position_in_codon ==0:
                        codons.append([i-2,i-1,i,position_in_codon,value[2]])
                    if position_in_codon ==1:
                        codons.append([i-1,i,i+1,position_in_codon,value[2]])
                    if position_in_codon ==2:
                        codons.append([i,i+1,i+2,position_in_codon,value[2]])
    return codons

def complement_base(base):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement_dict.get(base, base)

def reverse_complement(sequence):
    complement_sequence = [complement_base(base) for base in sequence]
    return ''.join(complement_sequence)[::-1]  # Reversing the sequence

def translate_sequence(sequence):
    dna_seq = Seq(sequence)
    protein_seq = dna_seq.translate()
    return str(protein_seq)

def codon_to_aa(condon_seq):
    condon_seq= SeqIO.parse(condon_seq,"fasta")
    extracted_sequences={}
    for seq in condon_seq:
        header = seq.id 
        extracted_sequence=seq.seq
        print(len(extracted_sequence))
        extracted_sequences[header]=translate_sequence(extracted_sequence)
        print(len(extracted_sequences[header]))
    return extracted_sequences

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find codons corresponding to snp in snp alignment file")
    parser.add_argument(
        "--snp_align_path",
        required=True,
        help="Path to input containing snp alignment .fasta file"
    )
    parser.add_argument(
        "--info_df_path",
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
    info_df_path = args.info_df_path
    annotation_path = args.annotation_path
    out_dir = args.outdir
    info_df = make_info_df(info_df_path,"Supp. table 14" )
    info_df = info_df.apply(lambda col: np.where(col == '.', info_df['CO92 Ref.'], col))
    headers = headers_from_fasta(snp_align_path)
    annotation_dict = make_annotation_dict(read_annotation_file(annotation_path))
    #get_sequence_lengths(full_align_path,headers)
    annotation_dict=get_ORF(annotation_dict,full_align_path,headers)
    annotation_dict = codon_finder(info_df,annotation_dict)
    codon = find_bases(annotation_dict)

    

    #full_align_sequence = SeqIO.parse(full_align_path, 'fasta')
    #if re.search(pattern, input_string)
    # add refererence sequence so you dont have to iterate through all
    snp_align_sequence = SeqIO.parse(snp_align_path, 'fasta')
    extracted_sequences={}
    for record in SeqIO.parse(full_align_path, 'fasta'):
            header = record.id
            if bool(re.search(r'^Reference\w+', header)):
                sequence = str(record.seq)

    for seq_record in snp_align_sequence:# iterate throuh snp
        header = seq_record.id 
    #if header in list(info_df.columns.to_list()):
    #sequence = str(seq_record.seq)
        if header in list(info_df.columns.to_list()):
            extracted_sequence=""
            for i in codon:
                if i[4]=='+':
                    if i[3]==0:
                        extracted_sequence+=info_df.loc[i[0],header]+sequence[i[1]] +sequence[i[2]]
                    elif i[3]==1:
                        extracted_sequence+=sequence[i[0]] + info_df.loc[i[1],header]+sequence[i[2]]
                    elif i[3]==2:  
                        extracted_sequence+=sequence[i[0]] +sequence[i[1]]+ info_df.loc[i[2],header]
                elif i[4]=='-':
                    if i[3]==0:
                        extracted_sequence+=reverse_complement(sequence[i[0]] +sequence[i[1]]+ info_df.loc[i[2],header])
                    elif i[3]==1:
                        extracted_sequence+=reverse_complement(sequence[i[0]] + info_df.loc[i[1],header]+sequence[i[2]])
                    elif i[3]==2:
                        extracted_sequence+=reverse_complement(info_df.loc[i[0],header]+sequence[i[1]] +sequence[i[2]])
                extracted_sequences[header]= extracted_sequence
    os.makedirs(args.outdir, exist_ok=True)
    try:
        filename_codon = "{}{}_codon.fasta".format(args.outdir,os.path.splitext(os.path.basename(snp_align_path))[0])
        f_out = open(filename_codon, 'w')
    except IOError:
        print("Output file {} cannot be created".format(filename_codon))
        sys.exit(1)
    # Write the extracted sequences to a new FASTA file
    for header, sequence in extracted_sequences.items():
        #print(header,sequence)
        
        #f_out.write(">{header}\n{sequence}\n".format(header,sequence))
        f_out.write(f'>{header}\n{sequence}\n')
    f_out.close()

    try:
        filename_aa = "{}{}_aa.fasta".format(args.outdir,os.path.splitext(os.path.basename(snp_align_path))[0])
        f_out = open(filename_aa, 'w')
    except IOError:
        print("Output file {} cannot be created".format(filename_aa))
        sys.exit(1)
    # Write the extracted sequences to a new FASTA file
    extracted_aa_sequences = codon_to_aa(extracted_sequences)
    for header, sequence in extracted_aa_sequences.items():
        #print(header,sequence)
        #f_out.write(">{header}\n{sequence}\n".format(header,sequence))
        f_out.write(f'>{header}\n{sequence}\n')
    f_out.close()

    '''
        for header in find_headers_in_fasta(data_dir):
            f_out.write('{}:{}\n'.format(header,''))
        f_out.close()
        '''