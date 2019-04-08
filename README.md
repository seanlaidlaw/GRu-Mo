![GRu-Mo](img/logo.png)

# GRu-Mo: a school project pseudo-aligner based on sqlite

## Design
On launch, the given fasta file will be read line by line, being decomposed into multiple kmers and their respective hashes which are then used to populate an SQLite3 database, indexed on its hash column.

On the reading of a fastq file, the line is decomposed into every possible kmer for which the hash is calculated and looked up in the database. At this point the aligment score between the fasta and fastq is calculated based on penalties for mismatch, gap opening, and gap continuing, provided as arguments to the program. If the total penalty is lower than the cutoff specified as arugment then the quant column of the transcript corresponding to the the matched hash is incremented by one.

The idea of using sqlite was to be able to produce a working prototype in very few lines of code. Only having to write the functions for parsing and managing kmers means the entire project only took a couple of hundred lines of code. Additionally, storing only the relevant hashes as an indexed column in a SQLite3 means the memory usage usually associated with hash tables is not present.


## Benchmarking:

| Program  | Time for index and transcript quantification*     | Quantified Reads | Similarity Score\** |
|:--------:|:-------------------------------------------------:|:----------------:|--------:|
| Kallisto |                   20 417                          |     87 150       |  - |
| GRu-Mo   |                  215 418                          |     31 684       |  70.63%|

\* for an E. coli reads and genome

\** determined by number of bases that match precisely between the two strands

Compared to Kallisto, there is only a 10x time penalty.

## Usage
Minimal working example:

`./GRu_Mo --fasta example_data/e-coli-genome.fasta -fastq example_data/run.fastq`

However, to acheive better results, the penalty scores used to calculate the alignment, necessary to know what transcripts to quantify, can be explicitly specified:

`./GRu_Mo --kmer-size 14 --max-gap 3 --cutoff 20 --mismatch 3 --gap-pen 3 --gap-extend 1 --verbose 1 --fasta example_data/e-coli-genome.fasta -fastq example_data/run.fastq`
