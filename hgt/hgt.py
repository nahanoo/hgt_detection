from os.path import join
from os import listdir
from Bio import SeqIO
from subprocess import call
import subprocess
import pysam
from Bio.SeqRecord import SeqRecord
import json
import argparse


class Hgt:
    """This class chunks a genetic sequence (can be a read or
    an assembly) with a sliding window algorithm. Chunks are
    then mapped to all ancestreal genomes from an experiment
    and the contig name of the reference where the chunks align
    are returned."""

    def __init__(self, args):
        # Dictionary with sample info from Sample class
        self.references = {fasta.split('.')[0]: join(
            args.references, fasta) for fasta in listdir(args.references)}
        # Fasta of sample assembly
        query = args.query_genome
        # Contigs in dict form
        self.contigs = self.get_contigs(query)
        # Out direcotry to write files
        self.out_dir = args.out_dir
        # Step size for sliding windows algorithm
        self.step = None
        # Dictionary storing origins of sequence in assembly
        self.origins = {key: dict() for key in self.contigs.keys()}
        # Filtered origins
        self.filtered = {key: dict() for key in self.contigs.keys()}

    def get_contigs(self, fasta):
        """Parses fastas and returns dictionary with contig name as
        key and sequence as value."""
        return {contig.id: contig.seq for contig in SeqIO.parse(fasta, "fasta")}


    def get_reference_names(self):
        """Returns dictionary. Reference_names with contig
        names as keys and strain as value."""
        reference_names = {}
        for strain, reference in self.references.items():
            for contig in SeqIO.parse(reference,'fasta'):
                reference_names[contig.id] = strain
        return reference_names

    def chunker(self, seq, window_size, step):
        """Creates chunks of a sequence. window_size defines
        chunk size and step the amount of basepairs the window
        is moved forward."""
        # List which stores all chunks
        seqs = []
        seqlen = len(seq)
        self.step = step
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
        window = 500
        self.step = 100
        for name, contig in self.contigs.items():
            record = SeqRecord(contig, id=name)
            # Creates chunks of every contig
            assembly_chunks += self.chunker(record, window, self.step)
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
        for strain,reference in self.references.items():
            # Query sequences created with chunk_assembly()
            sam = join(self.out_dir,strain + ".sam")
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
                 stderr=subprocess.STDOUT)

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
                print(name,self.step)
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

    def dump_origins(self):
        """Filters identified origins and dumps to json. Only positions
        with either two different origins or an origin which is not anceteral
        are of interest"""
        for name, contig in self.origins.items():
            for pos, strains in contig.items():
                # Applying filter
                if (len(set(strains)) > 1) or (strains[0] != self.sample['strain']):
                    self.filtered[name][pos] = list(set(strains))

        # Dumping to json
        j_f = json.dumps(self.filtered, indent=4)
        with open(join(self.sample['dir_name'], 'origins.json'), 'w') as handle:
            handle.write(j_f)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Detect HGT in bacterial communities. Fast and simpple.')
    parser.add_argument(
        'references', help='folder containg all reference genomes in fasta format.\
            Prefix of fasta file will be interpreted as strain name.'
    )
    parser.add_argument(
        'query_genome', help='query genome from which HGT should be detected.\
            Input fomrat needs to be fasta for now.')
    parser.add_argument('out_dir', help='output directory for dumping csv')

    return parser.parse_args()


args = parse_args()
hgt = Hgt(args)
