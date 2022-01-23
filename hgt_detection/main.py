import argparse
from .hgt import Hgt


def parse_args():
    parser = argparse.ArgumentParser(
        description='Detect HGTs in bacterial communities. Fast and simpple.')
    parser.add_argument(
        'references', help='folder containg all reference genomes in fasta format.\
            Prefix of fasta file will be interpreted as strain name.'
    )
    parser.add_argument(
        'query_genome', help='Query genome in genbank from which HGT should be detected.')
    parser.add_argument(
        'strain', help='Strain name for filtering input.')
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
        hgt.annotate_hgts()
        hgt.plot_hgts()
