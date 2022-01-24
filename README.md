# HGT

Detect HGTs in bacterial communities.

## Installation

```
git clone https://github.com/nahanoo/hgt.git
cd hgt
pip install .
```

This tool requires `samtools` and `minimap2` to be in your path.

## Usage

`detect_hgts --help`
```
usage: detect_hgts [-h] [--plot] references query_genome strain out_dir

Detect HGTs in bacterial communities. Fast and simpple.

positional arguments:
  references    folder containg all reference genomes in fasta format. Prefix
                of fasta file will be interpreted as strain name.
  query_genome  query genome in genbank from which HGTs should be detected.
  strain        strain name for filtering output.
  out_dir       output directory for storing files and plots

optional arguments:
  -h, --help    show this help message and exit
  --plot        if this flag is added the alignment of every hgt and the
                annotation of the hgt is plotted. Its fast.
```

The folder in references contains all the references in fasta format. Prefix will be interpeted as strain name.
Example: 
```
references/
    - agrobacterium_tumefaciens.fasta # strain name: agrobacterium_tumefaciens
    - comamonas_testosteroni.fasta # strain name: comamonas_testosteroni
    - microbacterium_saperdae.fasta # strain name: microbacterium_saperdae
    - ochrobactrum_anthropi.fasta # strain name: ochrobactrum_anthropi
```
The query_genome needs to be in genbank format. If your data is coming from an evolution experiment,
you likely have the fasta file of the ancetral strain in the references/ folder. By specifying the strain name used
for the ancetral genome in the references folder, output is filtered that only hits show up from the other
strains in references or if the sequence aligns to two differenct strains.

As an example we say that we want to analyze a strain which was evolved from the ancestral agrobacterium_tumefaciens.fasta sequence.
We can run `detect_hgts references/ agrobacterium_tumefaciens.evolved.genbank agrobacterium_tumefaciens ./`.
This will only ouput hits where sequences originate from other strain than `agrobacterium_tumefaciens.fasta` or where the sequence
aligns to `agrobacterium_tumefaciens.fasta` and at least one other strain.

## Output

Output is a simple tsv with chromosome, position, length and sequence origin.
Additionally the annotation of the hgt is outputted.
Furthermore every alignment of a hgt to all identified references is ouputted.

## Background

This packages chunks the query genome using a sliding window algorithm with a window size of 500 base-pairs and a step size of 100 base-pairs.
The position of the chunk is stored in the sequence name by concating the contig name and the step number. Therefore, we always know the original position
in the query.
Those sequences are then aligned to every reference genome. The aligned files are analyzed. Because the position in the query genome is stored in the sequence name we can trace back the position in the query and now know the sequence origin.