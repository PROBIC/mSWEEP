# Example use case
This example is written for Themisto v3.0.0 or newer and mSWEEP v2.0.0 or newer.

## Downloading the data
Download a toy dataset from [Zenodo](https://zenodo.org/record/8113529).

## Indexing
Create a list containing paths to the input assemblies
```
ls -d $(pwd)"/assemblies/"*.fasta.gz > input_sequences.txt
```

Index the input data with [Themisto](https://github.com/algbio/themisto) (v3.0.0 or newer)
```
mkdir tmp
themisto build -k 31 -i input_sequences.txt -o themisto_index --temp-dir tmp -t 2 --mem-gigas 4
```
this will create the `themisto_index.tcolors` and `themisto_index.tdbg` files.

## (optional) Read correction
Correct errors in the reads using [fastp](https://github.com/opengene/fastp]
```
fastp --in1 215_1.fastq.gz --in2 215_2.fastq.gz --out1 corr_1.fastq.gz --out2 corr_2.fastq.gz -c --thread 2
```

## Pseudoalignment
Pseudoalign the reads with Themisto
```
themisto pseudoalign -q 215_1.fastq.gz -i themisto_index --temp-dir tmp -t 2 > 215_1.txt
themisto pseudoalign -q 215_2.fastq.gz -i themisto_index --temp-dir tmp -t 2 > 215_2.txt
```

### (optional) pseudoalign and compress alignment
Pseudoalign the reads as above and compress the alignment file with [alignment-writer](https://github.com/tmaklin/alignment-writer)
```
ntargets=$(wc -l clustering.txt | cut -f1 -d' ')
nreads=$((`gunzip -c 215_1.fastq.gz | wc -l` / 4 ))
themisto pseudoalign -q 215_1.fastq.gz -i themisto_index --temp-dir tmp -t 1 | alignment-writer -n $ntargets -r $nreads > 215_1.aln
themisto pseudoalign -q 215_1.fastq.gz -i themisto_index --temp-dir tmp -t 1 | alignment-writer -n $ntargets -r $nreads > 215_2.aln
```
this is particularly useful for very large alignments. mSWEEP v2.0.0 and newer can automatically detect the file format if the alignment files are compressed.

## Abundance estimation
Estimate abundances with
```
mSWEEP --themisto-1 215_1.txt --themisto-2 215_2.txt -i clustering.txt -t 2
```
this will print the abundances after the estimation is finished. To write the abundances to `215_abundances.txt`, run
```
mSWEEP --themisto-1 215_1.txt --themisto-2 215_2.txt -i clustering.txt -t 2 -o 215
```

## Binning and abundance estimation
Bin the reads by adding the `--bin-reads` toggle
```
mSWEEP --themisto-1 215_1.txt --themisto-2 215_2.txt -i clustering.txt -t 2 --bin-reads
```
which will create the `clust1.bin`, `clust2.bin`, `clust3.bin`, and `clust4.bin` files. These files can be used as input to [mGEMS](https://github.com/probic/mgems) to extract the reads from the input data to `themisto pseudoalign`.
