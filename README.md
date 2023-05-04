# Kable
A tool for sequence, annotation, and methylome comparison between multiple bacterial genomes

## Introduction
Various bioinformatics tools exist to annotate gene sequences in the bacterial genome and to find variants; however, linking variants to specific genes or identifying changes to gene boundaries require additional steps. Further, there are currently no tools to compare the location of base modifications between genomes. 

To address the shortcomings in current analysis tools, I present Kable. Kable utilizes a colored de Bruijn Graph to align multiple genomic sequences while simultaneously tracking associated feature information. Kable can handle gene annotations, base modification or any other data that can be described in a GFF or VCF format, including potential assembly errors identified by [al2var](https://github.com/jrhendrix/al2var). Kable provides three key features to the user. First, the user can search for a specific sequence or find a segment between two sequences and Kable will return any and all annotation data occurring within those sequences and display their exact locations in an alignment file. Second, Kable can identify variants between any number of genomes, outputting the locations of such variants to a VCF and displaying overlaying annotations in an alignment file. Finally, Kable offers the first method to directly compare the location of base modifications and record differences in the methylome to a specialized VCF file. 

## Getting Started

### Requirements
* Python v3+
* BioSeq
* Levenshtein python package


### Installation
Clone this repository and run kable.py

## Usage

### Build graph
As input, kable build requires a tab deliminated manifest file. Each line of the manifest file represents one sample and contains a list of files including one FASTA file and any number of feature files in GFF or TSV format that pertain to that sample. The line can also include an optional orientation indicator (‘-‘ or ‘+’). By default, Kable will input the FASTA file in the forward direcotry, but if the orientation indicator is equal to '-', then Kable will input the reverse complement of the sample's sequence. 

Kable build parses incoming genome sequences using a set k-mer length. By default, Kable uses a k size of 23, but the user can set a specific k-mer length using the ‘-k’ or ‘—kmer_size’ flag. Values as low as 3 are accepted, though the results may not be meaningful depending on the sequence length. If the user specifies a k-mer length below 3, Kable will automatically adjust the k-mer size to 3.

```
# Basic usage
python kable.py build -m manifest.tsv 

# Use a k-mer length of 33
python kable.py build -m manifest.tsv -k 33
```


### Add to graph
The add module requires a pre-existing kable-generated graph and feature file along with a new manifest file. The manifest should follow the same formatting pattern used by the build module. 

The add module will imoprt the graph, determine the k-mer size by the existing node length, then add any new samples included in the manifest file. New graph and feature files will be exported and labeled with 'add' in the savename.

```
# Basic usage
python kable.py add -g kable_output/kable_graph.gfa -f kable_output/kable_feats.tsv -m new_manifest.tsv
```

### Find sequence variants
The find_vars module identifies sequence variants between the genome sequences. This module exports a VCF file for each genome sequence. It also aligns the variable sequences between genomes and shows the location of gene annotations in relation to the variants. 


```
# Basic usage
python kable.py find_vars -g kable_output/kable_graph.gfa -f kable_output/kable_feats.tsv

# Search for reverse complement and set maximum search depth to 100
python kable.py find_vars -g kable_output/kable_graph.gfa -f kable_output/kable_feats.tsv -r -d 100

# Save sub graph and exclude annotaiton data
python kable.py find_vars -g kable_output/kable_graph.gfa -f kable_output/kable_feats.tsv -w -e
```


### Find base modifications

The kable mod_search module shows where base modifications occur in relation to gene annotations, summarizes the number and type of modifications in sequences, and counts the number and type of modifications in each genome. 

```
# Basic usage
python kable.py mod_search -g kable_output/kable_graph.gfa -f kable_output/kable_feats.tsv

# Search for reverse complement and set maximum search depth to 100
python kable.py mod_search -g kable_output/kable_graph.gfa -f kable_output/kable_feats.tsv -r -d 100

# Save sub graph and exclude annotaiton data
python kable.py mod_search -g kable_output/kable_graph.gfa -f kable_output/kable_feats.tsv -w -e
```

### Query for a specific sequence
The Kable query module searches for a provided sequence in the graph and returns a list of which genones and contigs contain the sequence. The module will also export a pseudo alignemnt file that shows the location of genome features such as gene annotations and base modifications in relation to the sequence.

```
# Basic usage
python kable.py query -g kable_output/kable_graph.gfa -f kable_output/kable_feats.tsv -q sequence

# Include base modifications in pseudo alignment file
python kable.py query -g kable_output/kable_graph.gfa -f kable_output/kable_feats.tsv -q sequence -m

# Search for the reverse complement sequence
python kable.py query -g kable_output/kable_graph.gfa -f kable_output/kable_feats.tsv -q sequence -r

```

### Find sequence between two points

### General usage
Regardless of the speicific module that is run, Kable will create an output direcotry within the current working directory called kable_output/. All log files and generated output will be placed within this directory and given the prefix kable along with a module dependent suffix. 

The location of this output directory can be changed for any module using the -p or --output-path flag as long as the path already exists. The name of the output directory can be changed for any module using the -o or --output_directory flag. Further, the prefix of the output files can be set using the -n or --savename flag. 





