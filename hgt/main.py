import argparse
from .hgt import Hgt


def parse_args():
    parser = argparse.ArgumentParser(
        description='Detect HGTs in evolved bacterial samples. Tested with PacBio data.')
    parser.add_argument(
        'references', help='folder containg all reference genomes in fasta format.\
            Prefix of fasta file will be interpreted as strain name.'
    )
    parser.add_argument(
        'query_genome', help='query genome in genbank from which HGTs should be detected.')
    parser.add_argument(
        'strain', help='strain name for filtering output.')
    parser.add_argument('out_dir', help='output directory for storing files and plots')
    parser.add_argument('--plot', help='if this flag is added the alignment of every hgt and the annotation of the hgt is plotted. Its fast.', action='store_true')

    return parser.parse_args()

def main():
    args = parse_args()
    hgt = Hgt(args)
    hgt.chunk_assembly()
    hgt.mapper()
    hgt.get_mapping_stats()
    hgt.dump_origins()
    if args.plot:
        hgt.plot_hgt_annotations()
        hgt.plot_hgts()

main()