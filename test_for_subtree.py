def copy_fasta_with_10_positions(input_fasta, output_fasta):
    with open(input_fasta, 'r') as f_in, open(output_fasta, 'w') as f_out:
        header = None
        for line in f_in:
            if line.startswith('>'):
                header = line.rstrip()
                f_out.write(f"{header}\n")
            else:
                sequence_10_positions = line[:10]
                f_out.write(f"{sequence_10_positions}\n")

# Replace the file names accordingly
input_fasta_file = '10_first_test.fasta'
output_fasta_file = '10_first_subtree.fasta'

copy_fasta_with_10_positions(input_fasta_file, output_fasta_file)
