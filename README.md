# mSWEEP
Fast and accurate bacterial community composition estimation on strain
level by using pseudoalignments and variational inference.

More about mSWEEP in the article [High-resolution sweep metagenomics
using fast probabilistic inference](https://wellcomeopenresearch.org/articles/5-14)
in Wellcome Open Research.

If you use our method, please cite us as "Mäklin T, Kallonen T, David
S et al. High-resolution sweep metagenomics using fast probabilistic
inference [version 2; peer review: 2 approved]. Wellcome Open Res
2021, 5:14 (https://doi.org/10.12688/wellcomeopenres.15639.2)"

You can also watch a talk from [ISMB/ECCB
2019](https://www.iscb.org/ismbeccb2019) presenting the main results
from the mSWEEP paper [in YouTube](https://www.youtube.com/watch?v=VDfChoJwSKg).

# Installation
mSWEEP can be obtained either in the form of a precompiled binary
* [Linux 64-bit binary](https://github.com/PROBIC/mSWEEP/releases/download/v1.5.1/mSWEEP_linux-v1.5.1.tar.gz)
* [macOS 64-bit binary](https://github.com/PROBIC/mSWEEP/releases/download/v1.5.1/mSWEEP_macOS-v1.5.1.tar.gz)
or by following the instructions below for compiling mSWEEP from source.

In addition to mSWEEP, you will need to install either [Themisto
(recommended)](https://github.com/algbio/themisto) or
[kallisto](https://github.com/pachterlab/kallisto) for
pseudoalignment.

## Compiling from source
### Requirements
- C++11 compliant compiler.
- cmake (v3.0 or newer)

### Optional
- Compiler with OpenMP support.

If your compiler does not support OpenMP, mSWEEP can only be run in
single-threaded mode. The prebuilt binaries are compiled with OpenMP
support and allow parallellization.

### Compilation
Clone the mSWEEP repository
```
git clone https://github.com/PROBIC/mSWEEP.git
```
enter the directory and run
```
> mkdir build
> cd build
> cmake ..
> make
```
- This will compile the mSWEEP executable in build/bin/mSWEEP.

#### MPI support
mSWEEP can be compiled with MPI support, distributing the mixture
component estimation part of the program to several processes. To
compile with MPI support, set your environment appropriately and build
mSWEEP with the following commands (example case for Open MPI):
```
> mkdir build
> cd build
> module load mpi/openmpi
> cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
> make
```

The project should configure itself appropriately. To distribute the
computation to 4 processes after compiling, call mSWEEP with:
```
mpirun -np 4 mSWEEP --themisto-1 forward_aln.gz --themisto-2 reverse_aln.gz -i cluster_indicators.txt
```

In some cases it might be useful to use hybrid parallellization with
multiple threads per process. This can be accomplished through use of
the `-t` flag:
```
mpirun -np 2 mSWEEP --themisto-1 forward_aln.gz --themisto-2 reverse_aln.gz -i cluster_indicators.txt -t 2
```

which will distribute computation to two processes with two
threads. The optimal configuration will depend on the size and
structure of your data.

Note that when mSWEEP is called through MPI, the root process will
handle all read and write operations and only the estimation part is
distributed.

#### Compilation tips for improving performance
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
what sequences to include, or suspect that our samples contain
sequencing data from contaminating species.

2. a) Provide a grouping for the assemblies (e.g. sequence types,
   clonal complexes, or the output of some clustering algorithm.). If
   your species have an established multi-locus sequence typing
   scheme, the [mlst software](https://github.com/tseemann/mlst) is a
   good candidate for producing the grouping. Alternatively,
   [PopPUNK](https://github.com/johnlees/PopPUNK) can be used to
   cluster the sequences.

2. b) (Optional) Filter out reference sequences that cannot be reliably assigned
   to a group (eg. the sequence type cannot be determined) and perform
   other appropriate quality control measures.

3. If an assembly contains multiple contigs, merge them into a
   single contig. Do this for all assemblies.

4. Index the database with your pseudoalignment tool of choice and
   proceed with running the mSWEEP pipeline.

## Toy data (Themisto)
(Recommended) Enter the toy data directory (example/) and run the
build_index and pseudoalign commands from Themisto. The --mem-megas
option controls the maximum amount of RAM used in constructing the
index (if the limit is exceeded, Themisto will use disk storage), and
--n-threads the number of threads to use.
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
recommended (running mSWEEP without this option will *not* validate
the input 'clustering.txt' and may cause undefined behaviour).

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
```
## Build grouped Themisto index
mkdir themisto_grouped_index
mkdir tmp
build_index --k 31 --input-file example.fasta --color-file clustering.txt --index-dir themisto_grouped_index --temp-dir tmp

## Pseudoalign reads
pseudoalign --query-file 215_1.fastq.gz --outfile 215_1_alignment.txt --rc --index-dir themisto_grouped_index --temp-dir tmp --sort-output
pseudoalign --query-file 215_2.fastq.gz --outfile 215_2_alignment.txt --rc --index-dir themisto_grouped_index --temp-dir tmp --sort-output

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
- Obtain a set of reference sequences. (See the steps under Usage -> Reference data)
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
corresponding sequences appear in the reference file.

Multiple groupings can be supplied by adding them as columns to the
file containing the identifiers. When mSWEEP detects multiple
groupings, relative abundances will be estimated for all of them. When
outputting to file, the column index will be appended to the output
file name. For example with four sequences and two different groupings:
```
cluster1	clusterA
cluster2	clusterB
cluster2	clusterC
cluster1	clusterA
```
The default delimiter for the columns is a tab. The delimiter can be
changed to something else with the `--groups-delimiter` option.

Alternatively, you can use supply the reference sequences and a table
containing the groups to reorder the indicators.

### Reordering identifiers
If your grouping identifiers are not in the same order as in the fasta
file, you mSWEEP can reorder them based on the fasta file. In this
case, mSWEEP requires as input a tab-separated table where the first
column contains the name of the sequence as it appears in the fasta
file and the second column contains the group:
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
and the fasta file. If the above table is contained in the
'groups_list.tsv" file, running
```
mSWEEP --themisto-1 215_1_alignment.txt --themisto-2 215_2_alignment.txt -t 2 --fasta example.fasta --groups-list groups_list.tsv
```
will instruct mSWEEP to extract the identifiers from the table in the
correct order.

It is also possible to use a table
separated by a character different than tab, e. g. ',', by specifying the
separator with the '--groups-delimiter' argument:
```
mSWEEP --themisto-1 215_1_alignment.txt --themisto-2 215_2_alignment.txt -t 2 --fasta example.fasta --groups-list groups_list.tsv --groups-delimiter ,
```

Should you require the reordered identifiers written to a file
(e. g. for running mSWEEP on multiple input files), please refer to
the 'matchfasta' utility shipped alongside mSWEEP which performs the reordering.

Reordering the identifiers also supports providing several groupings at once:
```
seq_156	group2	groupA
seq_157	group1	groupB
seq_197	group4	groupC
seq_20	group3	groupD
seq_285	group4	groupE
seq_43	group2	groupA
seq_44	group1	groupB
seq_5	group3	groupF
```
which will instruct mSWEEP to estimate relative abundances for both groupings.

## Analysing reads (with Themisto)
- Pseudomap paired-end reads:
```
mkdir tmp
pseudoalign --index-dir themisto_index --query-file reads_1.fastq.gz --outfile reads_1_out.txt --temp-dir tmp --rc --sort-output
pseudoalign --index-dir themisto_index --query-file reads_2.fastq.gz --outfile reads_2_out.txt --temp-dir tmp --rc --sort-output
```

- Use mSWEEP to estimate cluster abundances from a single file with 2 threads:
```
mSWEEP --themisto-1 215_1_alignment.txt --themisto-2 215_2_alignment.txt -i clustering.txt -t 2
```

Themisto can also utilize multiple threads in the mapping phase. You can
run Themisto on multiple threads by specifying the number of threads
with the '--n-threads 8' flag.

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
	How to merge Themisto pseudoalignments for paired-end reads	(intersection or union, default: intersection).
	--themisto-index <ThemistoIndex>
	Path to the Themisto index the pseudoalignment was performed against (optional).

	--fasta <ReferenceSequences>
	Path to the reference sequences the pseudoalignment index was constructed from (optional)
	--groups-list <groupIndicatorsList>
	Table containing names of the reference sequences (1st column) and their group assignments (2nd column) (optional)
	--groups-delimiter <groupIndicatorsListDelimiter>
	Delimiter character for the --groups option (optional, default: tab)

	--iters <nrIterations>
	Number of times to rerun estimation with bootstrapped alignments (default: 1)
	--bootstrap-count <nrBootstrapCount>
	How many reads to resample when bootstrapping (integer, default: all)
	--seed <BootstrapSeed>
	Seed for the random generator used in bootstrapping (default: random)

	--write-probs
	If specified, write the read equivalence class probabilities in a .csv matrix
	--print-probs
	Print the read equivalence class probabilities to cout
	--write-likelihood
	Write the likelihood matrix to a file with "_likelihoods.txt" suffix if -o option is specified, print to cout if -o is not.
	--write-likelihood-bitseq
	Write the likelihoods in a format can be parsed by BitSeq's (https://github.com/bitseq/bitseq) functions.
	--gzip-probs
	Gzip the .csv matrix output from --write-probs and the likelihoods from --write-likelihood or --write-likelihood-bitseq.

	--read-likelihood
	Read in a likelihood matrix that has been written to file with the --write-likelihood toggle.

	--help
	Print this message.
	--version
	Print the version number.
	--cite
	Print citation information.

	ELBO optimization and modeling
	--no-fit-model
	Skip fitting the model entirely. Useful if only the likelihood matrix is required.
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
