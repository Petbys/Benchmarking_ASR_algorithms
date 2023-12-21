# input file: alignments
from Bio import SeqIO


#for taxa in Alignment_file:
def get_metrics(Alignment_file):

    Sequences = [[i.seq,i.id] for i in Alignment_file]
    id = Sequences[0][1]
    alignment_length = len(Sequences[0][0])
    Ambiguous_x = 0
    mismatches = 0
    matches = 0
    gaps = 0
    for pos,base in enumerate(Sequences[0][0]):
        if base == 'X':
            Ambiguous_x +=1
        elif Sequences[1][0][pos] == '-':
            gaps +=1
        elif base != Sequences[1][0][pos]:
            mismatches +=1
        elif base == Sequences[1][0][pos]:
            matches +=1
    #print(alignment_length)
    #print(Ambiguous_x,mismatches,matches,gaps)
    #print(Ambiguous_x+mismatches+matches+gaps)

    Correctly_identified = matches/(alignment_length-Ambiguous_x)
    Ambiguous_sites = Ambiguous_x/alignment_length
    gaps = gaps / (alignment_length-Ambiguous_x)

    print('---------------------------------------------------------------')
    print(f'Metrics {id} ')
    print('---------------------------------------------------------------')
    print('\n')
    print(f' Amount of correctly identified characters: {round(Correctly_identified,3)*100}%')
    print(f' Amount of gaps: {round(gaps,4)*100}%')
    print(f' Amount of ambigous sites in true sequence: {round(Ambiguous_sites,3)*100}%')
    print('\n')
    return Correctly_identified,Ambiguous_sites,gaps


Alignment_file_paml_bolgar = SeqIO.parse('/Users/petterbystrom/Documents/Applied Bioinformatics/Aligned_bolgar.fasta', 'fasta')
Alignment_file_paml_bsk = SeqIO.parse('/Users/petterbystrom/Documents/Applied Bioinformatics/Aligned_bsk.fasta', 'fasta')
Alignment_file_arpip_bsk = SeqIO.parse('/Users/petterbystrom/Documents/Applied Bioinformatics/Aligned_bsk_arpip.fasta', 'fasta')
Alignment_file_arpip_bolgar = SeqIO.parse('/Users/petterbystrom/Documents/Applied Bioinformatics/Aligned_bolgar_arpip.fasta', 'fasta')
print('\n')
print('################## PAML ##################')
print('\n')
get_metrics(Alignment_file_paml_bolgar)
get_metrics(Alignment_file_paml_bsk)
print('\n')
print('################## ARPIP ##################')
print('\n')
get_metrics(Alignment_file_arpip_bolgar)
get_metrics(Alignment_file_arpip_bsk)