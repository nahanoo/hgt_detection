import argparse
from hgt import Hgt
import yappi


def parse_args():
    parser = argparse.ArgumentParser(
        description='Detect HGTs in evolved bacterial samples. Tested with PacBio data.')
    parser.add_argument(
        'references', help='folder containg all reference genomes in fasta format.\
            Prefix of fasta file will be interpreted as strain name.'
    )
    parser.add_argument(
        'query_genome', help='query genome in genbank from which HGTs should be detected.')
    parser.add_argument('out_dir', help='output directory for storing files and plots')
    parser.add_argument('--plot', help='if this flag is added the alignment of every hgt and the annotation of the hgt is plotted. Its fast.', action='store_true')

    return parser.parse_args()

#@profile
def main():
    args = parse_args()
    hgt = Hgt(args)
    print('parsed')
    hgt.chunk_assembly()
    print('chunked')
    hgt.mapper()
    print('mapped')
    hgt.get_mapping_stats()
    print('analyzed')
    hgt.dump_origins()
    print('dumped')
    if args.plot:
        hgt.plot_hgt_annotations()
        hgt.plot_hgts()

main()