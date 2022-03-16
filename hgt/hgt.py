from os.path import join, split as path_split, exists
from os import listdir, mkdir, remove
from Bio import SeqIO
from subprocess import call, DEVNULL, STDOUT
import pysam
from Bio.SeqRecord import SeqRecord
import pandas as pd
from .plotting import plot_alignment, plot_genbank

class Hgt:
    """This is a simple class for detecting HGTs in mutated
    bacterial strains. It chunks the query genome into smaller
    sequences using a sliding window algorithm. Those sequences
    are aligned to all parsed reference genomes.
    All positions with foreign sequences are outputted as a tsv.
    Additionally it annotates the sequences and plots
    the alignments."""

    def __init__(self, args):
        # Out direcotry to write files
        self.out_dir = args.out_dir
        # Strain name and path of references
        self.references = {fasta.split('.')[0]: join(
            args.references, fasta) for fasta in listdir(args.references)}
        # Genbank of sample assembly
        self.query = args.mutant
        # Contigs in list form
        self.query_contigs = [contig for contig in SeqIO.parse(
            self.query, 'genbank')]
        # Features of mutant
        self.query_features = self.parse_genbank()
        # Ancestor fasta file
        self.ancestor = args.ancestor
        # Ancestor name
        self.ancestor_name = path_split(self.ancestor)[-1].split('.')[0]
        # Step size for sliding windows algorithm
        self.step = 100
        # Window size for sliding window algorithm
        self.window_size = 500
        # Dictionary storing origins of sequences in assembly
        self.origins = {contig.id: dict() for contig in self.query_contigs}
        # Dataframe for origins
        self.origins_df = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'origins'])
        # Dataframe for anntoations
        self.annotated = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'origins', 'product'])
        # List for items to deltet at the end
        self.trash = []

        # If plots should be generated directory is created
        if args.plot:
            plots = join(self.out_dir, 'plots')
            if not exists(plots):
                mkdir(plots)

    def parse_genbank(self):
        """Parses all features of a genbanka and stores locations
        and products in dictionary"""
        genbank = {contig.id: {} for contig in self.query_contigs}
        for contig in self.query_contigs:
            for feature in contig.features:
                # Some features don't have all desired keys
                try:
                    start = feature.location.start
                    end = feature.location.end
                    product = feature.qualifiers['product']
                    genbank[contig.id][(start, end)] = product[0]
                except KeyError:
                    pass
        return genbank

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
            chunk_id = seq.id
            chunk_seq = seq.seq[i:j]
            # Add chunk id to sequence id
            chunk = SeqRecord(seq=chunk_seq, id=chunk_id + "." + str(counter))
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
        self.chunks = join(self.out_dir,
                      "chunked_sequences.fasta")
        # Dumps chunks to fasta
        with open(self.chunks, "w") as handle:
            SeqIO.write(assembly_chunks, handle, "fasta")
        self.trash.append(self.chunks)

    def mapper(self, reference, reads, out):
        """Maps long accurate sequences to references with minimap2."""
        cmd = [
            "minimap2",
            "-ax",
            "asm5",
            reference,
            reads,
            ">",
            out,
        ]
        bam = out.replace('.sam', '.sorted.bam')
        # Calling minimap and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        cmd = ['samtools', 'sort', '-o', bam, out]
        # Calling samtools and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        cmd = ['samtools', 'index', bam]
        # Calling samtools and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        # Storing which files to delete at the end
        if out not in self.trash:
            self.trash.append(out)
            self.trash.append(bam)
            self.trash.append(bam+'.bai')

    def map_chunks(self):
        """Calls mapper to align mutant to references."""
        for strain, reference in self.references.items():
            out = join(self.out_dir, strain + ".sam")
            self.mapper(reference,self.chunks,out)

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
            # Read must be mapped
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

    def annotate_hgts(self):
        """Iterate over every detected hgt and check which product
        was inserted."""
        i = 0
        for counter, (c,p,l,o) in self.origins_df.iterrows():
            products = self.annotate_position(c, p, l)
            for product in products:
                self.annotated.loc[i] = [c, p, l, o, product]
                i += 1

    def annotate_position(self, c, p, l):
        """Returns products in a region in a genbank."""
        products = []
        for (start, end), product in self.query_features[c].items():
            if not set(range(start, end)).isdisjoint(range(p, p+l)):
                products.append(product)
        return products

    def concat_origins(self):
        """Concats origins. Outputs the start position and the length 
        of the foreing origin."""
        i = -1
        for name, contigs in self.origins.items():
            # Previous pos
            prev = 0
            # Previous strain
            prev_strains = None
            for pos, strains in contigs.items():
                strains = ', '.join(list(set(strains)))
                # Only concating when strains are identical
                if (pos - 1 == prev) & (prev_strains == strains):
                    self.origins_df.at[i, 'length'] += 1
                else:
                    i += 1
                    self.origins_df.loc[i] = [name, pos, 1, strains]
                prev = pos
                prev_strains = strains

    def check_sequence(self):
        """This checks if detected hgt was already present in ancestor."""
        for i, (chromosome, position, length, origins) in self.origins_df.iterrows():
            contigs = {contig.id: contig for contig in self.query_contigs}
            # Minimal sequence size which works well with minimap2 is around 500
            if length < 500:
                length = 500
            # Get correct padding
            if position + length > len(contigs[chromosome]):
                end = len(contigs[chromosome])
            else:
                end = length + position
            # Hgt sequence and sam
            query_seq = join(self.out_dir, 'query.fasta')
            query_out = join(self.out_dir, 'query.sam')
            # Dumping hgt sequence
            seq_id = contigs[chromosome].id
            seq = contigs[chromosome].seq[position:end]
            with open(query_seq, 'w') as handle:
                SeqIO.write(SeqRecord(seq=seq, id=seq_id), handle, 'fasta')
            # Mapping hgt sequence to ancestor
            self.mapper(self.ancestor,query_seq,query_out)
            a = pysam.AlignmentFile(query_out)
            unique = True
            # If sequence aligned to ancestor seq was already present
            for read in a:
                if (not read.is_unmapped):
                    unique = False
            # Adding ancestor to output if sequence was found
            if not unique:
                self.origins_df.at[i, 'origins'] += ', ' + \
                    self.ancestor_name
            # Deletes tmp seqs
            if query_seq not in self.trash:
                self.trash.append(query_seq)

    def dump_origins(self):
        """Filters identified origins and dumps to json. Only positions
        with either two different origins or an origin which is not anceteral
        are of interest"""
        self.concat_origins()
        self.check_sequence()
        self.origins_df.to_csv(
            join(self.out_dir, 'hgts.tsv'), sep='\t', index=False)
        self.annotate_hgts()
        self.annotated.to_csv(
            join(self.out_dir, 'hgts.annotated.tsv'), sep='\t', index=False)

    def plot_hgts(self):
        """Plots alignemtns to reference genomes and ancestor."""
        # Create plot dir
        out = join(self.out_dir, 'plots', 'hgt_alignments')
        if not exists(out):
            mkdir(out)
        # Reads identified origins
        sam = join(self.out_dir,self.ancestor_name+'.sam')
        self.mapper(self.ancestor,self.chunks,sam)
        for i, row in self.origins_df.iterrows():
            c = row['chromosome']
            p = row['position']
            origins = row['origins']
            origins = [element.lstrip() for element in origins.split(',')]
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
        for i, row in self.origins_df.iterrows():
            c = row['chromosome']
            p = row['position']
            # Plots features
            plot_genbank(self.query_contigs, c, p, out)

    def clean(self):
        """Cleans created files."""
        for item in self.trash:
            remove(item)

