
### Experimental usage
#### Bootstrapping confidence intervals
mSWEEP can be used to produce confidence intervals for the abundance
estimates by bootstrapping the pseudoalignment counts and rerunning
the abundance estimation a number of times. This can be done
automatically by adding the '--iters' option to running mSWEEP:
```
> mSWEEP --themisto-1 215_1_alignment_1.txt --themisto-2 215_2_alignment_2.txt -i cluster_indicators.txt -t 2 -o 215 --iters 100
```

The bootstrapped abundance estimates will be appended to the output
file '215_abundances.txt' as new columns and can be used to calculate
confidence intervals for each of the abundance estimates.

Bootstrapping in mSWEEP by default resamples from all input pseudoalignments. If
you wish to resample fewer pseudoalignments, supply the
'--bootstrap-count <positive integer>' option with the number of resamples to perform
```
> mSWEEP --themisto-1 215_1_alignment_1.txt --themisto-2
215_2_alignment_2.txt -i cluster_indicators.txt -t 2 -o 215 --iters 100 --bootstrap-count 1000
```
This will result in wider confidence intervals because the bootstrap
iterations have less data available to them.

It is also possible to set the initial random seed with the '--seed
<positive integer>' option, which enables replicating the bootstrap
results across multiple runs.

#### Embedding colors in the Themisto index
Alternatively, it is possible to embed the grouping
information in Themisto's index, effectively treating any pseudoalignment in the
group as if the read aligned to all sequences in the group. This
approach will in most cases produce different results than the recommended one, but will
reduce the RAM, CPU, and disk space requirements for running Themisto
and mSWEEP.

__Themisto version v2.0.0 or newer__

```
## Build grouped Themisto index
mkdir themisto_grouped_index
mkdir tmp
themisto build -k 31 -i example.fasta -c clustering.txt -o themisto_grouped_index --temp-dir tmp

## Pseudoalign reads
themisto pseudoalign -query-file 215_1.fastq.gz -o 215_1_alignment.txt -i themisto_grouped_index --rc --temp-dir tmp --sort-output
themisto pseudoalign -query-file 215_2.fastq.gz -o 215_2_alignment.txt -i themisto_grouped_index --rc --temp-dir tmp --sort-output
```

__Themisto versions v0.1.1 to v1.2.0__

```
## Build grouped Themisto index
mkdir themisto_grouped_index
mkdir tmp
build_index --k 31 --input-file example.fasta --color-file clustering.txt --index-dir themisto_grouped_index --temp-dir tmp

## Pseudoalign reads
pseudoalign --query-file 215_1.fastq.gz --outfile 215_1_alignment.txt --rc --index-dir themisto_grouped_index --temp-dir tmp --sort-output
pseudoalign --query-file 215_2.fastq.gz --outfile 215_2_alignment.txt --rc --index-dir themisto_grouped_index --temp-dir tmp --sort-output
```

##### Using embedded colors with mSWEEP

```
## Extract unique cluster indicators from the clustering.txt file
awk '!seen[$0]++' clustering.txt > unique_clusters.txt

## Run mSWEEP with 2 threads
mSWEEP --themisto-1 215_1_alignment.txt --themisto-2 215_2_alignment.txt -i unique_clusters.txt -t 2
```
These example commands will construct the Themisto index with the
grouping indicators from 'clustering.txt' embedded in it. This means
that if you wish to change the grouping indicators, the index must be
reconstructed.
