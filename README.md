# mSWEEP-GPU
Extension of mSWEEP with GPU acceleration for the abundance estimation algorithm. Adds two new options to the mSWEEP command line interface: `--algorithm` and `--emprecision`.

Algorithm descriptions:
- `rcggpu`: The original mSWEEP algorithm, implemented with LibTorch and can take advantage of a GPU. This algorithm is the fastest of the three, but uses ~25% more memory than the original algorithm. When using the GPU, memory usage is capped by the GPU memory (for example a Nvidia A100 may have 80GB of memory), so very large inputs may run into memory issues. With no GPU, this algorithm will run on the CPU, and is still marginally faster than the original algorithm.
- `emgpu`: A new algorithm that uses an expectation-maximization algorithm to estimate the relative abundances. This algorithm is slower than `rcggpu`, but uses ~40% less memory than the original algorithm. Inputs that run into memory issues on the GPU with `rcggpu` may still be able to take advantage of the GPU with this algorithm. Also, the precision of the algorithm can be changed with the `--emprecision` flag. The `float` precision is much faster and halves the memory usage, but loses some accuracy compared to the `double` precision and other algorithms. With no GPU, this algorithm will run on the CPU, and is faster than the original algorithm when using the `float` precision, but slower when using the `double` precision. **Cannot be used with `--bin-reads`.**
- `rcgcpu`: The original mSWEEP algorithm, running on the CPU. This algorithm is the slowest of the three, but can be improved through OpenMP support and supplying a large number of CPU's.

The default algorithm is set to `rcggpu`, which will use the GPU if available.

# mSWEEP
Fast and accurate bacterial community composition estimation on within-species
level by using pseudoalignments and variational inference.

More about mSWEEP in the article [High-resolution sweep metagenomics
using fast probabilistic inference](https://wellcomeopenresearch.org/articles/5-14)
in Wellcome Open Research.

# Installation
In addition to mSWEEP, you will need to install [Themisto](https://github.com/algbio/themisto) for pseudoalignment.

## Compiling from source
### Requirements
- C++17 compliant compiler.
- cmake (v3.0 or newer)
- [LibTorch](https://pytorch.org/get-started/locally/)

Notes on downloading LibTorch:
- If intended for use on the CSC or UH HPC clusters make sure to download the (Pre-cxx11 ABI) version for compatability with older glibc versions.
- If needing an older version of LibTorch due to not having the latest versions of ROCm (like on LUMI), these can be downloaded by changing the versions in the link.
  - For example: https://download.pytorch.org/libtorch/rocm6.0/libtorch-shared-with-deps-2.3.1%2Brocm6.0.zip to https://download.pytorch.org/libtorch/rocm5.2/libtorch-shared-with-deps-1.13.1%2Brocm5.2.zip

#### Optional
- CUDA Toolkit (if using LibTorch with CUDA support; version depending on downloaded LibTorch) or ROCm (if using LibTorch with ROCm support; version depending on downloaded LibTorch)
- Compiler with OpenMP support.

Without a CUDA or ROCm supported LibTorch and CUDA Toolkit or ROCm there will be no GPU acceleration.

If your compiler does not support OpenMP, the original mSWEEP can only be run in single-threaded mode.

### Compiling
Clone the mSWEEP repository
```
git clone https://github.com/Piketulus/mSWEEP-gpu.git
```
enter the directory and run
```
> mkdir build
> cd build
> cmake -DCMAKE_LIBTORCH_PATH=/absolute/path/to/libtorch ..
> cmake --build .
```
where `/absolute/path/to/libtorch` should be the absolute (!) path to the unzipped LibTorch distribution.

This will compile the mSWEEP executable in `build/bin/mSWEEP`.

For more info on compiling mSWEEP from source, please see the [documentation on compiling mSWEEP](/docs/compilation.md).

# Usage
More information about using mSWEEP is available in the [usage documentation](/docs/README.md).

## Abundance estimation
Estimate relative abundances of [Themisto](https://github.com/algbio/themisto) pseudoalignments `fwd.txt` and `rev.txt` using the lineages in `clustering.txt` with two threads by running 
```
mSWEEP --themisto-1 fwd.txt --themisto-2 rev.txt -i clustering.txt -t 2
```
or
```
mSWEEP --themisto fwd.txt,rev.txt -i clustering.txt -t 2
```

Both commands above will print the results. To write the result into a file called `result_abundances.txt` instead, run
```
mSWEEP --themisto-1 fwd.txt --themisto-2 rev.txt -i clustering.txt -t 2 -o result
```

## Binning reads
Estimate relative abundances and bin reads for all lineages with a relative abundance higher than 0.01 by running
```
mSWEEP --themisto-1 fwd.txt --themisto-2 rev.txt -i clustering.txt -t 2 --min-abundance 0.01
```

Extract reads only for reference lineages `lineage_1` and `lineage_2` by running
```
mSWEEP --themisto-1 fwd.txt --themisto-2 rev.txt -i clustering.txt -t 2 --target-groups lineage_1,lineage_2
```

Extract reads only for reference lineages `lineage_1` and `lineage_2` if their relative abundance is higher than 0.01 by running
```
mSWEEP --themisto-1 fwd.txt --themisto-2 rev.txt -i clustering.txt -t 2 --target-groups lineage_1,lineage_2 --min-abundance 0.01
```
i.e. the file format is automatically detected (alignment-writer v0.4.0 and newer).

## QC'ing binned reads
We recommend running [demix\_check](https://github.com/tmaklin/coreutils_demix_check) on the binned reads and/or [checkm](https://github.com/Ecogenomics/CheckM) on the bin-assembled genomes (BAGs) to evaluate the accuracy of the results.

## Working with large alignment files
### Compressing Themisto output files
For complex input data with many organisms, the pseudoalignment files from Themisto can get infeasibly large. In these cases, [alignment-writer](https://github.com/tmaklin/alignment-writer) can be used to compress the alignment files to <10% of the original size.

mSWEEP >=v2.0.0 can read the compressed alignments in directly by running
```
mSWEEP --themisto-1 fwd_compressed.aln --themisto-2 rev_compressed.aln -i clustering.txt -t 2

```

### Running estimation on large sparse alignments
If the target alignment is sparse, meaning that there are target groups which have few/no reads aligning against them in the whole sample, mSWEEP can be instructed to ignore these in the estimation by adding the `--min-hits 1` flag:
```
mSWEEP --themisto sparse.aln -i clustering.txt -t 2 --min-hits 1
```
This will reduce the runtime and memory use of the estimation proportional to how many target groups are removed. Using `--min-hits 1` does not affect the results beyond differences in computational accuracy.

The `--min-hits` flag also accepts values higher than 1 for pruning target groups with a small number of aligned reads. Using a value higher than 1 will change the resulting values.

## (experimental) Reliability of abundance estimates
Add the `--run-rate` flag to calculate a relative reliability value for each abundance estimate using a variation of the [RATE method](https://doi.org/10.1214/18-AOAS1222)
```
mSWEEP --themisto compressed.aln -i clustering.txt -t 2 --run-rate
```
This will append the RATE and KLD columns to the output. RATE values that exceed `1/(number of lineages in clustering.txt)` are considered reliably estimated.

If the reference contains many sequences that have zero or very few pseudoalignments, the denominator should be set to the number of lineages that have a nonzero abundance estimate instead of the total lineage count.

A reliably estimated value means that an abundance estimate from
mSWEEP has a large effect on how well the statistical model in mSWEEP
fits to the input alignment data. This translates to a high value in
the KLD column and the RATE columns, which is derived from the KLD
values by dividing each value by the sum of all KLDs.

__RATE as implemented in mSWEEP has not been tested thoroughly and is considered experimental.__ Consider using additional methods to verify the correctness of your results after filtering by RATE.

## More options
mSWEEP additionally accepts the following flags:

```
Usage: mSWEEP --themisto-1 <forwardPseudoalignments> --themisto-2 <reversePseudoalignments> -i <groupIndicatorsFile>
--verbose	Print status messages to cerr.
--version	Print mSWEEP version.
--cite	Print citation information.
--help	Print the help message.

Pseudoalignment files (required: -1 and -2, or only -x; will try to read from cin if none are given):
--themisto-1	Pseudoalignment results from Themisto for the 1st strand of paired-end reads.
--themisto-2	Pseudoalignment results from Themisto for the 2nd strand of paired-end reads.
--themisto	Single themisto alignment file or a comma separated list of several files.

Group indicators (required):
-i	Group indicators for the pseudoalignment reference.

Output prefix:
-o	Prefix for output files written from mSWEEP (default: print to cout).

Binning options:
--bin-reads	Run the mGEMS binning algorithm and write bins to the directory `-o` points to (default: false).
--target-groups	Only extract these groups, supply as comma separated list (default: extract all groups).
--min-abundance	Only extract groups that have a relative abundance higher than this value (default: 0).

Output options:
--write-probs	If specified, write the estimated read-to-group probabilities to a file with "_probs.tsv" suffix (default:false).
--print-probs	Print the read equivalence class probabilities to cout even if `-o` is given (default: false).
--write-likelihood	Write the internal likelihood matrix to a file with "_likelihoods.txt" suffix (default: false).
--write-likelihood-bitseq	Write the likelihoods in a format can be parsed by BitSeq's (https://github.com/bitseq/bitseq) functions (default: false).
--compress	Compress all output files using the given algorithm (one of z, bz2, lzma; default: don't compress).
--compression-level	Compression level (0-9; default: 6).

Input options:
--themisto-mode	How to merge pseudoalignments for paired-end reads (intersection, union, or unpaired; default: intersection).
--read-likelihood	Path to a precomputed likelihood file written with the --write-likelihood toggle. Can't be used with --bin-reads.

Estimation options:
-t	How many threads to use in abundance estimation (default: 1).
--no-fit-model	Do not estimate the abundances. Useful if only the likelihood matrix is required (default: false).
--max-iters	Maximum number of iterations to run the abundance estimation optimizer for (default: 5000).
--tol	Optimization terminates when the bound changes by less than the given tolerance (default: 0.000001).
--algorithm Which algorithm to use for abundance estimation (one of rcggpu, emgpu (can't be used with --bin-reads), rcgcpu (original mSWEEP); default: rcggpu).
--emprecision   Precision to use for the emgpu algorithm (one of float, double; default: double).

Bootstrapping options:
--iters	Number of times to rerun estimation with bootstrapped alignments (default: 0).
--seed	Seed for the random generator used in bootstrapping (default: random).
--bootstrap-count	How many pseudoalignments to resample when bootstrapping (default: number of reads).

Likelihood options:
-q	Mean for the beta-binomial component (default: 0.65).
-e	Dispersion term for the beta-binomial component (default: 0.01).
--alphas	Prior counts for the relative abundances, supply as comma-separated nonzero values (default: all 1.0).
--zero-inflation	Likelihood of an observation that contains 0 pseudoalignments against a reference group (default: 0.01).

Experimental options:
--run-rate	Calculate relative reliability for each abundance estimate using RATE (default: false).
--min-hits	Only consider target groups that have at least this many reads align to any sequence in them (default: 0).
```

# References
## Abundance estimation
If you use mSWEEP for abundance estimation, please cite us as
```
Mäklin T, Kallonen T, David S, Boinett CJ, Pascoe B, Méric G, Aanensen DM, Feil EJ, Baker S, Parkhill J, Sheppard SK, Corander J, and Honkela A
High-resolution sweep metagenomics using fast probabilistic inference [version 2; peer review: 2 approved]
Wellcome Open Resesearch 5:14 (2021)
https://doi.org/10.12688/wellcomeopenres.15639.2
```

## Binning
The binning algorithm enabled by the `--bin-reads` toggle is described in
```
Mäklin T, Kallonen T, Alanko J, Samuelsen Ø, Hegstad K, Mäkinen V, Corander J, Heinz E, and Honkela A
Bacterial genomic epidemiology with mixed samples
Microbial Genomics 7:11 (2021)
https://doi.org/10.1099/mgen.0.000691
```
Thee binning algorithm is also provided as the standalone software [mGEMS](https://github.com/probic/mgems).

## Specific versions
If you wish to cite a specific version of mSWEEP, visit the [releases
page](https://github.com/PROBIC/mSWEEP/releases) and find the doi for
the version of the program that you used. Then, cite the version as
```
Tommi Mäklin, and Antti Honkela. (2021).
PROBIC/mSWEEP: mSWEEP v1.5.0 (15 October 2021)
(v1.5.0). Zenodo. (https://doi.org/10.5281/zenodo.5571944)
```

# License
The source code from this project is subject to the terms of the MIT
license. A copy of the MIT license is supplied with the project, or
can be obtained at
[https://opensource.org/licenses/MIT](https://opensource.org/licenses/MIT).
