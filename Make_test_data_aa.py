from Bio import SeqIO
import argparse
import os
import sys

def create_list_headers(headers):
    header_list = []
    with open(headers,"r") as headers:
        for header in headers:
            header_list.append(header.rstrip("\n"))
    return header_list
'''
def create_fasta_from_headers(header_list, input_fasta_file, output_fasta_file):
    sequences_to_write = []

    # Loop through the input FASTA file and extract sequences based on headers
    for record in SeqIO.parse(input_fasta_file, "fasta"):
        if record.id in header_list:
            sequences_to_write.append(record)

            # Stop if 10 sequences have been found
            if len(sequences_to_write) == 10:
                break
    print(header_list)
    # Write the sequences to a new FASTA file
    with open(output_fasta_file, "w") as output_handle:
        SeqIO.write(sequences_to_write, output_handle, "fasta")
'''
def extract_sequences(headers, fasta_file, output_file):
    extracted_sequences = {}
    current_header = None
    with open(fasta_file, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                current_header = line
                if current_header in headers:
                    extracted_sequences[current_header] = ''
            elif current_header in headers:
                extracted_sequences[current_header] += line
    print(extracted_sequences)
    with open(output_file, 'w') as f_out:
        for header, sequence in extracted_sequences.items():
            sequence_10_positions = sequence[:10]
            f_out.write(f">{header}\n{sequence_10_positions}\n")

if __name__ == '__main__':
        parser = argparse.ArgumentParser(
            description="testing pyasr")
        parser.add_argument(
            "--fasta",
            required=True,
        )
        parser.add_argument(
            "--headers",
            required=True,
        )
        parser.add_argument(
            "--output",
            required=True,
        )

        args = parser.parse_args()
        fasta_path= args.fasta
        headers_path = args.headers
        output_path = args.output

        header_list = create_list_headers(headers_path)
        print(header_list)
        extract_sequences(header_list, fasta_path, output_path)
