from os.path import join, split as path_split, exists
from os import listdir, mkdir
from Bio import SeqIO
from subprocess import call, DEVNULL, STDOUT
import subprocess
import pysam
import json
import pandas as pd
from .plotting import plot_alignment, plot_genbank


class Hgt:
    """This is a simple class for detecting HGTs in mutated
    bacterial strains. It chunks the query genome into smaller
    sequences using a sliding window algorithm. Those sequences
    are aligned to all parsed reference genomes.
    All positions with foreign sequences are outputted as a tsv.
    Additionally it annotates the foreign sequences and plots
    the alignemnt to all reference genomes."""

    def __init__(self, args):
        # Out direcotry to write files
        self.out_dir = args.out_dir
        # Strain name and path of references
        self.references = {fasta.split('.')[0]: join(
            args.references, fasta) for fasta in listdir(args.references)}
        # Genbank of sample assembly
        self.query = args.query_genome
        # Strain of query genome
        self.query_strain = args.strain
        # Contigs in list form
        self.query_contigs = [contig for contig in SeqIO.parse(
            self.query, 'genbank')]
        # Step size for sliding windows algorithm
        self.step = 100
        # Window size for sliding window algorithm
        self.window_size = 500
        # Dictionary storing origins of sequence in assembly
        self.origins = {contig.id: dict() for contig in self.query_contigs}
        # Dictionary for storing filtered origins
        self.filtered = {contig.id: dict() for contig in self.query_contigs}
        # Dataframe for anntoations
        self.annotated = pd.DataFrame(columns=['chromosome','position','length','product'])


        # If plots should be generated directory is created
        if args.plot:
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
            # Read must be primary,mapped
            for read in a:
                if (not read.is_unmapped):
                    reads.append(read)
            for read in reads:
                # Contig and position of query sequence form assembly are stored in contig name
                # First part is contig name, second number is n step
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
        """Concats origins. Outputs the start position and the length 
        of the foreing origin."""
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

    def annotate_hgts(self,df):
        genbank = {contig.id: {} for contig in self.query_contigs}
        for contig in self.query_contigs:
            for feature in contig.features:
                try:
                    start = feature.location.start
                    end = feature.location.end
                    product = feature.qualifiers['product']
                    genbank[contig.id][(start, end)] = product[0]
                except KeyError:
                    pass
        i = 0
        for counter,row in df.iterrows():
            c = row['chromosome']
            p = row['position']
            l = row['length']
            for j in range(p,p+l):
                for (start, end), product in genbank[c].items():
                    if (j >= start) & (j <= end):
                        self.annotated.loc[i] = [c, p, l, product]
                        i += 1
        self.annotated = self.annotated.drop_duplicates()



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
        # Creating string of origins stored as list
        df['origins'] = [', '.join(map(str, l)) for l in df['origins']]
        df.to_csv(join(self.out_dir, 'origins.tsv'), sep='\t', index=False)
        # Dumping to json and tsv
        self.annotate_hgts(df)
        j_f = json.dumps(self.filtered, indent=4)
        with open(join(self.out_dir, 'origins.json'), 'w') as handle:
            handle.write(j_f)
        self.annotated.to_csv(join(self.out_dir, 'origins.annotated.tsv'), sep='\t', index=False)

    def plot_hgts(self):
        # Create plot dir
        out = join(self.out_dir, 'plots', 'hgt_alignments')
        if not exists(out):
            mkdir(out)
        # Reads identified origins
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
                # Creating reads of interest
                read_names = ['.'.join([c, str(step)])
                              for step in range(steps-10, steps+10)]
                name = '.'.join([c, str(p), origin])
                plot_alignment(bam, read_names, name, out)

    def plot_hgt_annotations(self):
        # Create plot out dir
        out = join(self.out_dir, 'plots', 'hgt_annotations')
        if not exists(out):
            mkdir(out)
        df = pd.read_csv(join(self.out_dir, 'origins.tsv'), sep='\t')
        for i, row in df.iterrows():
            c = row['chromosome']
            p = row['position']
            plot_genbank(self.query_contigs, c, p, out)
