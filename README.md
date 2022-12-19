# magicblast-meth (mbmeth.py)

Alignment of BS-Seq reads using Magic-BLAST. 

## Rationale 

There is NOT a tool availale for us to align BS-seq reads (DNA or RNA) to NCBI nt databse. This tool is developped to fill this gap. One application would be: if we want to check whether there is contamination or the contamination source of BS-seq reads. We can use magicblast-meth to align the reads to nt databse and get the accession number from the alignment results. Based on the accession ID, we can then extract the taxonomy information based on the ID to check whether there is any contamintation and/or what species are the contamintation come from. 

> Note: since Magic-BLAST support RNA-seq read alignment, magicblast-meth can be used to align `RNA Bisulfite-Seq` with the option `--spice T`. 

## Intro

This works for single-end reads and for **paired-end reads from the
directional protocol** (most common).

Uses the method employed by methylcoder and Bismark of *in silico*
conversion of all C's to T's in both reference and reads.

Recovers the original read (needed to tabulate methylation).

It allows us to align reads to nt database. 

## QuickStart

Without installation, you can use as `python ` with install, the
command is `bwameth.py`.

The commands:

```
mbmeth.py makeblastdb -r $REF
mbmeth.py magicblast -r $REF -r1 some_R1.fastq.gz -r2 some_R2.fastq.gz > some.output.sam
```
will create `some.output.bam` and `some.output.bam.bai`.
To align single end-reads, specify only 1 file: `-r1 some_read.fastq.gz`

## Installation

```
conda env create -n env4mbmeth --file environment.yaml python=3 
```

### Dependencies

```
Magic-BLAST
BLAST
samtools
```

## Acknowledgement

Special thanks to [bwa-meth](https://github.com/brentp/bwa-meth) because some of the codes were adapted from it. 

