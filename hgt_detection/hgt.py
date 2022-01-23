from os.path import join, split as path_split, exists
from os import listdir, mkdir
from pyexpat import features
from Bio import SeqIO
from subprocess import call, DEVNULL, STDOUT
import subprocess
import pysam
from Bio.SeqRecord import SeqRecord
import json
import argparse
import pandas as pd
from plotting import plot_alignment, plot_genbank


class Hgt:
    """This class chunks a genetic sequence (can be a read or
    an assembly) with a sliding window algorithm. Chunks are
    then mapped to all ancestreal genomes from an experiment
    and the contig name of the reference where the chunks align
    are returned."""

    def __init__(self, args):
        # Out direcotry to write files
        self.out_dir = args.out_dir
        # Dictionary with sample info from Sample class
        self.references = {fasta.split('.')[0]: join(
            args.references, fasta) for fasta in listdir(args.references)}
        # Fasta of sample assembly
        self.query = args.query_genome
        self.query_strain = path_split(self.query)[-1].split('.')[0]
        self.query_contigs = [contig for contig in SeqIO.parse(
            self.query, 'genbank')]
        with open(join(args.references, self.query_strain+'.fasta'),'w' ) as handle:
            SeqIO.write(self.query_contigs, handle, 'fasta')
        # Strain of query genome
        self.references[self.query_strain] = join(
            args.references, self.query_strain+'.fasta')
        # Contigs in list form
        self.query_contigs = [contig for contig in SeqIO.parse(
            self.query, 'genbank')]
        # Step size for sliding windows algorithm
        self.step = 10000
        self.window_size = 50000
        # Dictionary storing origins of sequence in assembly
        self.origins = {contig.id: dict() for contig in self.query_contigs}
        # Filtered origins
        self.filtered = {contig.id: dict() for contig in self.query_contigs}

        plots = join(self.out_dir, 'plots')
        if not exists(plots):
            mkdir(plots)

    def get_reference_names(self):
        """Returns dictionary. Reference_names with contig
        names as keys and strain as value."""
        reference_names = {}
        for strain, reference in self.references.items():
            for contig in SeqIO.parse(reference, 'fasta'):
                reference_names[contig.id] = strain
        return reference_names

    def chunker(self, seq, window_size, step):
        """Creates chunks of a sequence. window_size defines
        chunk size and step the amount of basepairs the window
        is moved forward."""
        # List which stores all chunks
        seqs = []
        seqlen = len(seq)
        for counter, i in enumerate(range(0, seqlen, step)):
            # Returns ether entire sequence or window depending on sequence length
            j = seqlen if i + window_size > seqlen else i + window_size
            chunk = seq[i:j]
            # Add chunk id to sequence id
            chunk.id = chunk.id + "." + str(counter)
            seqs.append(chunk)
            if j == seqlen:
                break
        return seqs

    def chunk_assembly(self):
        """Chunks an assembly of multiple contigs into different 
        chunks using a sliding window algorightm (see chunker function)."""
        assembly_chunks = []
        for contig in self.query_contigs:
            #record = SeqRecord(contig, id=name)
            # Creates chunks of every contig
            assembly_chunks += self.chunker(contig,
                                            self.window_size, self.step)
        target = join(self.out_dir,
                      "chunked_sequences.fasta")
        # Dumps chunks to fasta
        with open(target, "w") as handle:
            SeqIO.write(assembly_chunks, handle, "fasta")

    def mapper(self):
        """Maps chunked sequence to all ancesteral genomes
        experiment with minimap2.
        Minimap2 settings are set to accureate PacBio reads."""
        reads = join(self.out_dir, "chunked_sequences.fasta")
        for strain, reference in self.references.items():
            # Query sequences created with chunk_assembly()
            sam = join(self.out_dir, strain + ".sam")
            bam = join(self.out_dir, strain + ".sorted.bam")
            cmd = [
                "minimap2",
                "-ax",
                "asm5",
                reference,
                reads,
                ">",
                sam,
            ]
            # Calling minimap and surpressing stdout
            call(" ".join(cmd), shell=True, stdout=subprocess.DEVNULL,
                 stderr=STDOUT)
            cmd = ['samtools', 'sort', '-o', bam, sam]
            # Calling samtools and surpressing stdout
            call(" ".join(cmd), shell=True, stdout=DEVNULL,
                 stderr=STDOUT)
            cmd = ['samtools', 'index', bam]
            # Calling samtools and surpressing stdout
            call(" ".join(cmd), shell=True, stdout=DEVNULL,
                 stderr=STDOUT)

    def get_mapping_stats(self):
        """Checks all mapped sequences and returns contig name
        of reference."""
        reference_names = self.get_reference_names()
        for strain in self.references.keys():
            # Iterating over every sam file
            sam = join(
                self.out_dir, strain + ".sam"
            )
            a = pysam.AlignmentFile(sam, "rb")
            reads = []
            # Iterating over all reads
            # Read must be primary,mapped and have quality of 60
            for read in a:
                if (not read.is_unmapped) & (not read.is_secondary) & (read.mapq == 60):
                    reads.append(read)
            # Appends contig name of reference
            for read in reads:
                # Contig and position of query sequence form assembly are
                # stored in contig name
                # First part is contig name, second number is n step
                # step size is known and position therefore as well
                name = '.'.join(read.qname.split('.')[:-1])
                pos = int(read.qname.split('.')[-1])*self.step
                # Iteratign over aligned query sequence
                for j in range(pos+read.query_alignment_start, pos+read.query_alignment_end):
                    # Checking if position already in dictionary
                    if j not in self.origins[name].keys():
                        self.origins[name][j] = []
                        # Appending strain at given position
                        self.origins[name][j].append(
                            reference_names[read.reference_name])

                    else:
                        # Appending strain at given position
                        self.origins[name][j].append(
                            reference_names[read.reference_name])

    def concat_origins(self):
        """Concats origins. In input every position has strain info.
        This function concatenates positions following eachother and
        outputs the start position and the length of the foreing origin."""
        df = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'origins'])
        i = -1
        for name, contigs in self.filtered.items():
            # Previous pos
            prev = 0
            # Previous strain
            prev_strains = None
            for pos, strains in contigs.items():
                if (pos - 1 == prev) & (prev_strains == strains):
                    df.at[i, 'length'] += 1
                else:
                    i += 1
                    df.loc[i] = [name, pos, 1, strains]
                prev = pos
                prev_strains = strains
        return df

    def dump_origins(self):
        """Filters identified origins and dumps to json. Only positions
        with either two different origins or an origin which is not anceteral
        are of interest"""
        for name, contig in self.origins.items():
            for pos, strains in contig.items():
                # Applying filter
                if (len(set(strains)) > 1) or (list(set(strains))[0] != self.query_strain):
                    self.filtered[name][pos] = list(set(strains))
        df = self.concat_origins()
        df['origins'] = [', '.join(map(str, l)) for l in df['origins']]
        # Dumping to json
        j_f = json.dumps(self.filtered, indent=4)
        with open(join(self.out_dir, 'origins.json'), 'w') as handle:
            handle.write(j_f)
        df.to_csv(join(self.out_dir, 'origins.tsv'), sep='\t')

    def plot_hgts(self):
        out = join(self.out_dir, 'plots', 'hgt_alignments')
        if not exists(out):
            mkdir(out)
        df = pd.read_csv(join(self.out_dir, 'origins.tsv'), sep='\t')
        for i, row in df.iterrows():
            c = row['chromosome']
            p = row['position']
            origins = row['origins']
            origins = [element.lstrip() for element in origins.split(',')]
            origins = origins + [self.query_strain]
            for origin in set(origins):
                bam = join(self.out_dir, origin + '.sorted.bam')
                steps = int(p/self.step)
                read_names = ['.'.join([c, str(step)])
                              for step in range(steps-10, steps+10)]
                name = '.'.join([c, str(p), origin])
                plot_alignment(bam, read_names, name, out)

    def annotate_hgts(self):
        out = join(self.out_dir, 'plots', 'hgt_annotations')
        if not exists(out):
            mkdir(out)
        df = pd.read_csv(join(self.out_dir, 'origins.tsv'), sep='\t')
        for i, row in df.iterrows():
            c = row['chromosome']
            p = row['position']
            plot_genbank(self.query_contigs, c, p, out)


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
        'out_dir', help='output directory for storing files and plots')
    parser.add_argument(
        '--plot', help='if this flag is added the alignment of every hgt and the annotation of the hgt is plotted. Its fast.', action='store_true')

    return parser.parse_args()


args = parse_args()
hgt = Hgt(args)
