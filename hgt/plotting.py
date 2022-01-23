from dna_features_viewer import GraphicFeature, GraphicRecord
import pysam
from os.path import join

def plot_alignment(bam, read_names, name, out):
        a = pysam.AlignmentFile(bam)
        contigs = dict()
        reads = []
        for read in a:
            if (not read.is_unmapped) & (not read.is_secondary) & (read.mapq == 60) & (read.qname in read_names):
                reads.append(read)
        for read in reads:
            if read.is_reverse:
                strand = -1
            else:
                strand = 1
            if read.reference_name in contigs.keys():
                contigs[read.reference_name].append((read.qname,
                                                     read.reference_start, read.reference_end, strand))
        
            else:
                contigs[read.reference_name] = []
                contigs[read.reference_name].append((read.qname,
                                                     read.reference_start, read.reference_end, strand))
        for contig, positions in contigs.items():
            starts = []
            ends = []
            features = []
            for pos in positions:
                starts.append(pos[1])
                ends.append(pos[2])
                gf = GraphicFeature(start=pos[1], end=pos[2], strand=pos[3],
                                    color="#ffd700", label=pos[0])
                features.append(gf)
            record = GraphicRecord(
                sequence_length=sorted(ends)[-1], features=features)
            record = record.crop((sorted(starts)[0], sorted(ends)[-1]))
            name = '.'.join(name.split('.') + [contig, '.pdf'])
            record.plot_on_multiple_pages(join(out, name),
                                          nucl_per_line=sorted(
                                              ends)[-1] - sorted(starts)[0],
                                          lines_per_page=10,
                                          plot_sequence=False
                                          )