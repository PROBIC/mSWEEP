## Toy data (Themisto)
Download the toy dataset supplied with mSWEEP v1.X.X from [Zenodo](https://zenodo.org/record/8101285).

Then, enter the toy data directory (mSWEEP-v1\_example/) and run the
build_index and pseudoalign commands from Themisto. The --mem-megas
option controls the maximum amount of RAM used in constructing the
index (if the limit is exceeded, Themisto will use disk storage), and
--n-threads the number of threads to use.

__Themisto version v2.0.0 or newer__

```
## Build the Themisto index without clustering
mkdir themisto_index
mkdir tmp

## Index using at most 2048 megabytes of memory and 2 threads.
themisto build -k 31 -i example.fasta -o themisto_index/index --temp-dir tmp --mem-megas 2048 --n-threads 2

## Pseudoalign reads using 2 threads
themisto pseudoalign -q 215_1.fastq.gz -o 215_1_alignment.txt -i themisto_index/index --temp-dir tmp --n-threads 2 --rc --sort-output
themisto pseudoalign -q 215_2.fastq.gz -o 215_2_alignment.txt -i themisto_index/index --temp-dir tmp --n-threads 2 --rc --sort-output
```

__Themisto versions v0.1.1 to v1.2.0__

```
## Build the Themisto index without clustering
mkdir themisto_index
mkdir tmp

## Index using at most 2048 megabytes of memory and 2 threads.
build_index --k 31 --input-file example.fasta --auto-colors --index-dir themisto_index --temp-dir tmp --mem-megas 2048 --n-threads 2

## Pseudoalign reads using 2 threads
pseudoalign --query-file 215_1.fastq.gz --outfile 215_1_alignment.txt --rc --index-dir themisto_index --temp-dir tmp --n-threads 2 --sort-output
pseudoalign --query-file 215_2.fastq.gz --outfile 215_2_alignment.txt --rc --index-dir themisto_index --temp-dir tmp --n-threads 2 --sort-output
```

Next, run mSWEEP on the alignment files (the -t option controls the
number of threads used in the estimation.)
```
## Run mSWEEP with 2 threads
mSWEEP --themisto-1 215_1_alignment.txt --themisto-2 215_2_alignment.txt -i clustering.txt -t 2 -o 215 --themisto-index themisto_index
```
This will write the relateve abundances to the "215_abundances.txt"
file in the folder mSWEEP was run in. If the '-o' option is not
specified, the abundances will print to cout.

Note that supplying the --themisto-index is optional but highly
recommended for Themisto versions v1.2.0 and older (running mSWEEP
without this option will *not* validate the input 'clustering.txt' and
may cause undefined behaviour).
