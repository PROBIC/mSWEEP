# mSWEEP
Fast and accurate bacterial community composition estimation on strain
level by using pseudoalignments and variational inference.

More about mSWEEP in the article [High-resolution sweep metagenomics
using fast probabilistic
inference](https://doi.org/10.12688/wellcomeopenres.15639.1) in
Wellcome Open Research (awaiting peer review)

If you use our method, please cite us as Mäklin T, Kallonen T, David S
 et al. High-resolution sweep metagenomics using fast probabilistic
 inference [version 1; peer review: awaiting peer review]. Wellcome
 Open Res 2020, 5:14 (https://doi.org/10.12688/wellcomeopenres.15639.1)

# Installation
mSWEEP can be obtained either in the form of a precompiled binary
* [Linux 64-bit binary](https://github.com/PROBIC/mSWEEP/releases/download/v1.3.2/mSWEEP_linux-v1.3.2.tar.gz)
* [macOS 64-bit binary](https://github.com/PROBIC/mSWEEP/releases/download/v1.3.2/mSWEEP_macOS-v1.3.2.tar.gz)
or by following the instructions below for compiling mSWEEP from source.

In addition to mSWEEP, you will need to install either [Themisto
(recommended)](https://github.com/jnalanko/Themisto) or
[kallisto](https://github.com/pachterlab/kallisto) for
pseudoalignment.

## Compiling from source
### Requirements
- C++11 compliant compiler.
- cmake

### Optional
- Compiler with OpenMP support.

If your compiler does not support OpenMP, only limited parts of mSWEEP
will benefit from parallellization. The prebuilt binaries are compiled
with OpenMP support and support parallellization.

### Compilation
Clone the mSWEEP repository (note the --recursive option in git clone!)
```
git clone --recursive https://github.com/PROBIC/mSWEEP.git
```
enter the directory and run
```
> mkdir build
> cd build
> cmake ..
> make
```
- This will compile the mSWEEP executable in build/bin/mSWEEP.

### Compilation tips for improving performance
1. If you intend to run mSWEEP on the machine used in compiling the
source code, you might want to add the '-march=native -mtune=native'
flags if compiling with GCC by running
```
> cmake -DCMAKE_CXX_FLAGS="-march=native -mtune=native" -DCMAKE_C_FLAGS="-march=native -mtune=native" ..
> make
```
Using these options significantly reduces the runtime of mSWEEP in
some environments (e.g. most HPC setups).

2. If the [Intel C++ compiler](https://software.intel.com/en-us/c-compilers) is
   available in your environment, you might want to use that to
   compile mSWEEP — especially if running on Intel hardware. The compiler can be specified by running
```
> cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc ..
> make
```

# Usage
## Reference data

A reference sequence collection and a grouping into clonal
complexes/sequence types is available in
[FigShare](https://figshare.com/articles/mSWEEP_reference_v1-0-0_tgz/8222636)
for the following species
- Campylobacter jejuni
- Escherichia coli
- Klebsiella pneumoniae
- Staphylococcus epidermidis

mSWEEP supports using a custom reference database. A typical workflow
for constructing the custom database might proceed as follows

1. Gather assembled sequences for the species of interest. Use
taxonomic profiling tools like [MetaFlow](https://doi.org/10.1007/978-3-319-31957-5_8) or
[MetaPhlAn2](https://doi.org/10.1038/nmeth.3589) to identify the
species in your sample if you are unsure
what sequences to include.

2. Provide a grouping for the assemblies (e.g. sequence types, clonal
   complexes, or the output of some clustering algorithm.)

3. If an assembly contains multiple contigs, merge them into a
   single contig. Do this for all assemblies.

4. Index the database with your pseudoalignment tool of choice and
   proceed with running the mSWEEP pipeline.

## Toy data (Themisto)
(Recommended) Enter the toy data directory (example/) and run the
build_index and pseudoalign commands from Themisto
```
## Build the Themisto index without clustering
mkdir themisto_index
mkdir tmp

build_index --k 31 --input-file example.fasta --auto-colors --index-dir themisto_index --temp-dir tmp

## Pseudoalign reads
pseudoalign --query-file 215_1.fastq.gz --outfile 215_1_alignment.txt --rc --index-dir themisto_index --temp-dir tmp
pseudoalign --query-file 215_2.fastq.gz --outfile 215_2_alignment.txt --rc --index-dir themisto_index --temp-dir tmp

## Run mSWEEP
mSWEEP --themisto-1 215_1_alignment.txt --themisto-2 215_2_alignment.txt -i clustering.txt
```

(Experimental) Alternatively, it is possible to embed the grouping
information in Themisto's index, effectively treating any pseudoalignment in the
group as if the read aligned to all sequences in the group. This
approach will in most cases produce different results than the recommended one, but will
reduce the RAM, CPU, and disk space requirements for running Themisto
and mSWEEP.
```
## Build grouped Themisto index
mkdir themisto_grouped_index
mkdir tmp
build_index --k 31 --input-file example.fasta --color-file clustering.txt --index-dir themisto_grouped_index --temp-dir tmp

## Pseudoalign reads
pseudoalign --query-file 215_1.fastq.gz --outfile 215_1_alignment.txt --rc --index-dir themisto_grouped_index --temp-dir tmp
pseudoalign --query-file 215_2.fastq.gz --outfile 215_2_alignment.txt --rc --index-dir themisto_grouped_index --temp-dir tmp

## Extract unique cluster indicators from the clustering.txt file
awk '!seen[$0]++' clustering.txt > unique_clusters.txt

## Run mSWEEP with 2 threads
mSWEEP --themisto-1 215_1_alignment.txt --themisto-2 215_2_alignment.txt -i unique_clusters.txt -t 2
```
These example commands will construct the Themisto index with the
grouping indicators from 'clustering.txt' embedded in it. This means
that if you wish to change the grouping indicators, the index must be
reconstructed.

## Toy data (kallisto)
```
kallisto index -i example_index example.fasta

kallisto pseudo -i example_index -o kallisto_out_folder 215_1.fastq.gz 215_2.fastq.gz
mSWEEP -f kallisto_out_folder -i clustering.txt
```
You should see that roughly 90% of the reads are assigned to group "clust2".

# General pipeline
## Preprocessing
- Obtain a set of reference sequences.
- Index the reference sequences for pseudoalignment with Themisto:
> build_index --k 31 --input-file reference_sequences.fasta --auto-colors --index-dir themisto_index --temp-dir tmp
- ... or with kallisto
> kallisto pseudo -i reference_kmi reference_sequences.fasta
- Define some grouping for the reference and save the grouping in a text file where each line contains the identifier of the grouping the corresponding reference sequence belongs to. For example with four sequences and two groups:
```
cluster1
cluster2
cluster2
cluster1
```
The grouping identifiers must be in the same order as their
corresponding sequences appear in the reference file. Alternatively,
you can use the 'matchfasta' utility to reorder the indicators.

### Reordering identifiers
If your grouping identifiers are not in the same order as in the fasta
file, you can use the 'matchfasta' utility tool supplied
with mSWEEP to reorder them based on the fasta file. matchfasta takes
as input a tab-separated table where the first column contains the name of the
sequence as it appears in the fasta file and the second column
contains the group:
```
seq_156	group2
seq_157	group1
seq_197	group4
seq_20	group3
seq_285	group4
seq_43	group2
seq_44	group1
seq_5	group3
```
and the fasta file. Running
> matchfasta --fasta sequences.fasta --groups groups_table.tsv > clustering_reordered.txt

will save the identifiers in the correct order to the
'clustering_reordered.txt' file. It is also possible to use a table
separated by a character different from tab, e. g. ',', by specifying the
separator with the '-d' argument:
> matchfasta --fasta sequences.fasta --groups groups_table.csv > clustering_reordered.txt -d ,

## Analysing reads (with Themisto)
- Pseudomap paired-end reads:
```
mkdir tmp
pseudoalign --index-dir themisto_index --query-file reads_1.fastq.gz --outfile reads_1_out.txt --temp-dir tmp --rc
pseudoalign --index-dir themisto_index --query-file reads_2.fastq.gz --outfile reads_2_out.txt --temp-dir tmp --rc
```

- Use mSWEEP to estimate cluster abundances from a single file with 2 threads:
> mSWEEP --themisto-1 215_1_alignment.txt --themisto-2 215_2_alignment.txt -i clustering.txt -t 2

Themisto can utilize multiple threads in the mapping phase. You can
run Themisto on multiple threads by specifying the number of threads
with e.g. the '--threads 8' flag.

## Analysing reads (with kallisto)
- Pseudomap paired-end reads:
> kallisto pseudo -i reference_kmi -o kallisto_out reads_1.fastq reads_2.fastq
- Or use [kallisto's batch mode](https://pachterlab.github.io/kallisto/manual) to analyze multiple files at once:
> kallisto pseudo -i reference_kmi -o kallisto_batch_out -b kallisto_batch.txt

- Use mSWEEP to estimate cluster abundances from a single file:
> mSWEEP -f kallisto_out_folder -i cluster_indicators.txt -o abundances.txt -t 2
- Or from a batch file and output to abundances/ folder:
- (**NOT RECOMMENDED** anymore as of 18 December 2019: kallisto's batch mode can be a bit unstable)
> mSWEEP -b kallisto_batch_out -i cluster_indicators.txt -o abundances

kallisto can utilize multiple threads in the mapping phase. You can run kallisto on multiple threads by specifying the number of threads with the '-t' flag.
When estimating samples submitted in a kallisto batch, mSWEEP can estimate multiple samples in parallel by specifying the number of threads with the '-t' flag.

## Bootstrapping confidence intervals (experimental)
mSWEEP can be used to produce confidence intervals for the abundance
estimates by bootstrapping the pseudoalignment counts and rerunning
the abundance estimation a number of times. This can be done
automatically by adding the '--iters' option to running mSWEEP:
```
> mSWEEP --themisto-1 alignment_1.txt --themisto-2 alignment_2.txt -i cluster_indicators.txt --iters 100 -o abundances.txt
```
or with kallisto
```
> mSWEEP -f kallisto_out_folder -i cluster_indicators.txt --iters 100 -o abundances.txt
```
The bootstrapped abundance estimates will be appended to the output
file as new columns and can be used to calculate confidence intervals
for each of the abundance estimates.

Bootstrapping can be performed on multiple threads by adding the '-t
4' option (NOTE: running the bootstrapping on multiple threads hasn't
been optimised and may use large amounts of memory).

# Running mSWEEP
mSWEEP accepts the following flags:

```
	--themisto-1 <themistoPseudoalignment1>
	Pseudoalignment results from Themisto for the 1st strand of paired-end reads.
	--themisto-2 <themistoPseudoalignment2>
	Pseudoalignment results from Themisto for the 2nd strand of paired-end reads.

	-f <pseudomappingFile>
	Pseudoalignment output file location from kallisto. Can't be used when -b is specified.
	-b <pseudomappingBatch>
	The kallisto batch matrix file location. Can't be used when -f is specified.

	-i <clusterIndicators>
	Group identifiers file. Must be supplied.
	-o <outputFile>
	Output file (folder when estimating from a batch) to write results in.
	-t <nrThreads>
	How many threads to use. (default: 1)

	--themisto-mode <PairedEndMergeMode>
	How to merge Themisto pseudoalignments for paired-end reads	(default: intersection).
	--iters <nrIterations>
	Number of times to rerun estimation with bootstrapped alignments (default: 1)

	--write-probs
	If specified, write the read equivalence class probabilities in a .csv matrix
	--print-probs
	Print the equivalence class probabilities rather than writing when using --write-probs
	--gzip-probs
	Gzip the .csv matrix output from --write-probs
	--help
	Print this message.

	ELBO optimization and modeling (these seldom need to be changed)
	--tol <tolerance>
	Optimization has converged when the bound changes less than the given tolerance.
	--max-iters
	Maximum number of iterations to run the gradient optimizer.
	-q <meanFraction>
	Fraction of the sequences in a group that the mean is set to. (default: 0.65)
	-e <dispersionTerm>
	Calibration term in the likelihood function. (default: 0.01)
```

# License
The source code from this project is subject to the terms of the MIT
license. A copy of the MIT license is supplied with the project, or
can be obtained at
[https://opensource.org/licenses/MIT](https://opensource.org/licenses/MIT).
