import argparse
from .hgt import Hgt


def parse_args():
    parser = argparse.ArgumentParser(
        description='Detect HGTs in bacteria evolved in co-cultures.')
    parser.add_argument(
        'mutant', help='genbank file of the mutant')
    parser.add_argument(
        'ancestor', help='Fasta file of the ancestor. Prefix of fasta file will be interpreted as strain name.')
    parser.add_argument(
        'references', help='Folder containing all genomes in fasta format of co-cultured bacteria.\
            Prefix of fasta file will be interpreted as strain name.'
    )
    parser.add_argument('out_dir', help='output directory')
    parser.add_argument(
        '--plot', help='plot alignments and annotations', action='store_true')

    return parser.parse_args()


def main():
    args = parse_args()
    hgt = Hgt(args)
    hgt.chunk_assembly()
    hgt.map_chunks()
    hgt.get_mapping_stats()
    hgt.dump_origins()
    if args.plot:
        hgt.plot_hgt_annotations()
        hgt.plot_hgts()
    hgt.clean()