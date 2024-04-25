import pyranges as pr
import pandas as pd
import argparse
import os

DIRNAME = os.path.dirname(__file__)
REFFLAT_FILE = os.path.join(DIRNAME, "refFlat.txt")
VCF_FILE = os.path.join(DIRNAME, "coriell_chr1.vcf")


def get_range_refflat_and_vcf(refflat_file: str, vcf_file: str) -> tuple[pr.PyRanges, pr.PyRanges]:
    """
    Get PyRanges from refFlat.txt and vcf files
    :param refflat_file: path to the refFlat.txt file
    :param vcf_file: path to the vcf file
    :return: tuple of PyRanges from refFlat.txt and vcf files
    """
    range_refflat = pr.PyRanges(
        pd.read_csv(
            refflat_file,
            sep="\t",
            names=[
                "GenBrowser",
                "Name",
                "Chromosome",
                "Strand",
                "Start",
                "End",
                "cdsStart",
                "cdsEnd",
                "exonCount",
                "exonStarts",
                "exonEnds"
            ],
            header=None
        )
    )
    vcf = pd.read_csv(
        vcf_file,
        sep="\t",
        names=[
            "Chromosome",
            "Start",
            'ID',
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "_"
        ],
        header=None,
        comment="#"
    )
    vcf['End'] = vcf['Start']
    range_vcf = pr.PyRanges(vcf)

    return range_refflat, range_vcf


def variant_counter(refflat_file: str, chromosome: str, vcf_file: str, output_file: str, result_count: int):
    """
    Count the number of variants in each gene from the refFlat.txt file
    :param refflat_file: path to the refFlat.txt file
    :param vcf_file: path to the vcf file
    :param output_file: path to the output file
    """
    range_refflat, range_vcf = get_range_refflat_and_vcf(refflat_file, vcf_file)
    range_refflat = range_refflat[range_refflat.Chromosome == chromosome]
    unique_genes = range_refflat.GenBrowser.unique()
    
    genes_count = []
    for gene in unique_genes:
        rf_gene = range_refflat[range_refflat.GenBrowser == gene]
        intersects = range_vcf.intersect(rf_gene, how="containment")
        if intersects.empty:
            genes_count.append(0)
        else:
            genes_count.append(len(intersects.Start.unique()))

    result = pd.DataFrame(
        {
            "unique_genes": unique_genes,
            "genes_count": genes_count
        }
    )
    print(result.head(result_count))
    result.to_csv(output_file, sep=',', header=None, index=None)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--reffile', help='Path to the refFlat.txt file', type=str, default=REFFLAT_FILE)
    parser.add_argument('--chromosome', help='Chromosome number', type=str, default="chr1")
    parser.add_argument('--vcffile', help='Path to the vcf file', type=str, default=VCF_FILE)
    parser.add_argument('--output', help='Path to the output file', type=str, default="output.csv")
    parser.add_argument('--result_count', help='Number of results to display', type=int, default=5)
    args = parser.parse_args()

    variant_counter(args.reffile, args.chromosome, args.vcffile, args.output, args.result_count)