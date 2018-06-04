# mSWEEP
Fast high-resolution sweep metagenomics using pseudoalignments and variational inference.

More about the method\: paper info

# Installation
## Requirements
- C++11 compliant compiler.
- cmake

## Installing the pipeline
- Install [kallisto](https://github.com/pachterlab/kallisto).
- NOTE: if you want to run kallisto in batch mode, install [the forked version](https://github.com/tmaklin/kallisto) which contains a fix to the pseudoalignment counts when running in batch mode.
- Clone the mSWEEP repository, enter the directory and run
```
> mkdir build
> cd build
> cmake ..
> make
```
- This will download and install mSWEEP.
- You can remove the build directory afterwards.

# Usage
# Toy data
There is a toy dataset included in the example/ folder. To run it, enter the directory and run the commands

```
kallisto index -i example_index example.fasta

kallisto pseudo -i example_index -o kallisto_out_folder 215_1.fastq.gz 215_2.fastq.gz
mSWEEP -f kallisto_out_folder -i clustering.txt

```
You should see that roughly 94.6% of the reads are assigned to group "clust2".

# General pipeline
## Preprocessing
- Obtain a set of reference sequences.
- Index the reference sequences for pseudoalignment with:
> kallisto pseudo -i reference_kmi reference_sequences.fasta
- Define some grouping for the reference and save the grouping in a text file where each line contains the identifier of the grouping the corresponding reference sequence belongs to. For example with four sequences and two groups:

```
cluster1
cluster2
cluster2
cluster1
```
The grouping identifiers must be in the same order as their corresponding sequences appear in the reference file.

## Analysing reads
- Pseudomap paired-end reads:
> kallisto pseudo -i reference_kmi -o kallisto_out reads_1.fastq reads_2.fastq
- Or use [kallisto's batch mode](https://pachterlab.github.io/kallisto/manual) to analyze multiple files at once:
> kallisto pseudo -i reference_kmi -o kallisto_batch_out -b kallisto_batch.txt

- Use mSWEEP to estimate cluster abundances from a single file:
> mSWEEP -f kallisto_out_folder -i cluster_indicators.txt -o abundances.txt
- Or from a batch file and output to abundances/ folder:
> mSWEEP -b kallisto_batch_out -i cluster_indicators.txt -o abundances

kallisto can utilize multiple threads in the mapping phase. You can run kallisto on multiple threads by specifying the number of threads with the '-t' flag.
When estimating samples submitted in a kallisto batch, mSWEEP can estimate multiple samples in parallel by specifying the number of threads with the '-t' flag.

mSWEEP accepts the following flags:

```
	-f <pseudomappingFile>
	The pseudoalignment output file location. Can't be used when -b is specified.
	-b <pseudomappingBatch>
	The kallisto batch matrix file location. Can't be used when -f is specified.
	-i <clusterIndicators>
	Group identifiers file. Must be supplied.
	-o <outputFile>
	Output file (folder when estimating from a batch) to write results in.
	-t <nrThreads>
	How many threads to use when processing a batch matrix (default: 1)

	--write-probs
	If specified, write the read equivalence class probabilities in a .csv matrix
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
