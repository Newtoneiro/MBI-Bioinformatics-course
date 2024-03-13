import argparse
from Bio import SeqIO, SeqUtils


def gc_content(filename: str) -> float:
    """
    Calculate the GC content of a FASTA file.
    :param filename: str: Path to the FASTA file
    :return: float: GC content in percentage
    """
    with open(filename) as f:
        top, bottom = 0, 0
        for record in SeqIO.parse(f, "fasta"):
            top += SeqUtils.gc_fraction(record.seq) * len(record.seq)
            bottom += len(record.seq)

    return (top / bottom) * 100
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename', help='Input file in FASTA format', type=str, required=True)
    args = parser.parse_args()

    result = gc_content(args.filename)
    print(f"GC-Content: {result:0.2f} %")
