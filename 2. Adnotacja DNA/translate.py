import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

translation_dict = {
    "A": "GCU", "R": "CGU", "N": "AAU", "D": "GAU", "C": "UGU",
    "Q": "CAA", "E": "GAA", "G": "GGU", "H": "CAU", "I": "AUU",
    "L": "UUG", "K": "AAA", "M": "AUG", "F": "UUU", "P": "CCU",
    "S": "AGU", "T": "ACU", "W": "UGG", "Y": "UAU", "V": "GUU",
    "*": "UAA"
}

def retranslate_protein_to_mRNA(input_filename, output_filename):
    with open(input_filename) as input_handle:
        with open(output_filename, "w") as output_handle:
            for seq_record in SeqIO.parse(input_handle, "fasta"):
                translated_seq = ''.join(translation_dict.get(letter, 'NNN') for letter in seq_record.seq)
                mRNA_seq_record = SeqRecord(Seq(translated_seq), id=seq_record.id, description="Retranslated sequence")
                SeqIO.write(mRNA_seq_record, output_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='mbi2',
        description='Script to retranslate protein to mRNA. The program requires filenames (input and output) as arguments.'
    )
    parser.add_argument('input_filename', help='Input FASTA file containing protein sequences')
    parser.add_argument('output_filename', help='Output FASTA file to store retranslated mRNA sequences')
    args = parser.parse_args()
    
    retranslate_protein_to_mRNA(args.input_filename, args.output_filename)