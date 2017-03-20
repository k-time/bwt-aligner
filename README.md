# bwt-aligner
* Command line tool for error-tolerant short read mapping using the [Burrows-Wheeler transform (BWT)](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform), with FM-indexing, suffix tree search, and heuristics that prune branches at search time.
* See **[Bioinformatics Report.pdf](Bioinformatics%20Report.pdf)** and **[Presentation.pdf](Presentation.pdf)** for detailed project explanation.

Sequence alignment: https://en.wikibooks.org/wiki/Next_Generation_Sequencing_(NGS)/Alignment<br>
BWT for sequence alignment: http://bioinformatics.oxfordjournals.org/content/25/14/1754.long


## Files
We have included reference genome files and read files containing a random sample of 100 aligned reads for 3 different viruses. These genomes and reads are real data taken from the NCBI’s Sequence Read Archive.

The program files are described as follows:

* **bwt.py**: Implementation of BWT and exact search.
	
* **align_reads.py**: Takes a reference genome file and aligned reads file (and optional threshold level) as input. Uses our aligner to predict the positions of reads, and then compares them to the actual positions.  
	Usage: `python align_reads.py <genome file name> <read file name> [-t <threshold level>]`
* **search_bwt.py**: Implementation of inexact search algorithm. Can map a single read.  
	Usage: `python search_bwt.py [--no-indels] <reference_file> <read_file>`


## Recommended Usages
Here are our recommended usages to analyze our included data (with no -t flag, a default threshold of 3 is used):

* `python align_reads.py data/ebola_genome.fasta data/ebola_reads.fasta`
* `python align_reads.py data/coronavirus_genome.fasta data/coronavirus_reads.fasta`
* `python align_reads.py data/rsv_genome.fasta data/rsv_reads.fasta -t 7`
* `python search_bwt.py test` (test example of read mapping)


## Other Details
Entire read files tend to be rather large; they can be found by going to http://www.ncbi.nlm.nih.gov/sra, searching for an organism/virus, filtering by “DNA” and “aligned data”, and clicking on a result and downloading the aligned reads from the run. You can also can download the SRA Toolkit to download results through a command line interface.

Implemented using Python 2.7. Need to download library enum34 (run `pip install enum34`). We ran program on a MacBook; a more powerful machine will obviously speed up our results.
