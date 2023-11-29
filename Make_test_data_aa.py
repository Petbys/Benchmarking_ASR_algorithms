from Bio import SeqIO
import argparse
def create_fasta_from_headers(header_list, input_fasta_file, output_fasta_file):
    sequences_to_write = []

    # Loop through the input FASTA file and extract sequences based on headers
    for record in SeqIO.parse(input_fasta_file, "fasta"):
        if record.id in header_list:
            sequences_to_write.append(record)

            # Stop if 10 sequences have been found
            if len(sequences_to_write) == 10:
                break

    # Write the sequences to a new FASTA file
    with open(output_fasta_file, "w") as output_handle:
        SeqIO.write(sequences_to_write, output_handle, "fasta")


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


    create_fasta_from_headers(headers_path, fasta_path, output_path)
