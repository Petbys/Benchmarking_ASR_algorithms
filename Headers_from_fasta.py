import os
import sys
import argparse

def find_headers_in_fasta(file_path):
    headers = []
    
    with open(file_path, 'r') as fasta_file:
        current_header = None
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                # If a line starts with '>', it's a header line
                current_header = line
                headers.append(current_header)
    
    return headers


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Get headers from fasta file")
    parser.add_argument(
        "--datadir",
        required=True,
        help="Path to input directory containing fasta file"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Path to directory where output"
    )

    args = parser.parse_args()
    data_dir = args.datadir
    out_dir = args.outdir

    os.makedirs(args.outdir, exist_ok=True)
    try:
        file_name = "{}/Headers.txt".format(args.outdir)
        f_out = open(file_name, 'w')
    except IOError:
        print("Output file {} cannot be created".format(file_name))
        sys.exit(1)
    for header in find_headers_in_fasta(data_dir):
        f_out.write('{}:{}\n'.format(header,''))
    f_out.close()