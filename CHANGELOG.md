# v2.0.0 (4 July 2023)
First major version increment of mSWEEP (breaks backwards-compatibility).

Output format changes:
- Add the total number of reads to the abundances file (resolves #21)
- Renamed `total_hits` to `num_aligned` in the abundances file (#21)

## New features
- Added an option to evaluate the [mGEMS binning](https://github.com/PROBIC/mGEMS) algorithm from the mSWEEP call with the `--bin-reads` toggle. (https://github.com/PROBIC/mSWEEP/commit/54004d8c2408764dbd970cd1a49c253ede92ce5d)
- Support reading alignments compressed with [alignment-writer](https://github.com/tmaklin/alignment-writer). (https://github.com/PROBIC/mSWEEP/commit/f169fccbf976a0994863e0ffca4a62718778594c)
- Read alignments from cin. (https://github.com/PROBIC/mSWEEP/commit/9878f8ef4d342b72877409001e153ffc0433318a)

## Removed features:
- Matching a fasta file to groups indicators is no longer supported (deprecated options `--fasta`, `--groups-list`, `--groups-delimiter`).

## Pseudoaligner support
- Removed kallisto support ([remove-kallisto-support](https://github.com/PROBIC/mSWEEP/tree/remove-kallisto-support))
- Removed support for Themisto v1.2.0 and older.  ([remove-kallisto-support](https://github.com/PROBIC/mSWEEP/tree/remove-kallisto-support))

## Installation
- Added a conda recipe and instructions on installing mSWEEP from bioconda. (#22)

## Build pipeline changes
- Require C++17 to build from source.
- Removed support for building zlib from source (https://github.com/PROBIC/mSWEEP/commit/5c94591b392932efd7b1d16cb6c3e3ddc688c8e1)
- Added the `CMAKE_BUILD_WITH_FLTO` flag for building with link-time optimization (https://github.com/PROBIC/mSWEEP/commit/ee1db015e9902eca05293a416576e89f1d54af7a)

## Internal changes
- Bump C++ standard to C++17.
- Rewrote most of the codebase.
- Fixed dependency versions to avoid conflicts.


# v1.6.3 (1 February 2023)
Fix build issues caused by an update in one of the dependencies.

# v1.6.2 (12 August 2022)
Updated dependencies and bug hunting.

## Bug fixes
- Fix interaction of --no-fit-model with WriteResults.
- Skip erroneously trying to write the probability matrix if --no-fit-model was toggled.

## Updated dependencies
### Use telescope v0.4.0
- About 10x speedup in reading pseudoalignments.
- May reduce the memory footprint on large input.

### Use rcgpar v1.0.2
- Fixes MPI estimation when the input data dimensions exceed the capacity of 32 bit signed integers.
- Enables compilation without MPI support even when MPI headers are present on the system.

## Internal changes
- Rename log.hpp -> msweep_log.hpp and correct the header guard to avoid some conflicts with dependencies.

## Build pipeline
- Use the CMAKE_ENABLE_MPI_SUPPORT flag to compile with or without support for MPI.


# v1.6.1 (5 May 2022)
## Changes:
- Updated dependency bxzstr to v1.1.0.
- Disabled zstd support from bxzstr by default.
- - Support can be enabled when compiling by changing -D ZSTD_FOUND=0 to -D ZSTD_FOUND=1 in config/CMakeLists-bxzstr.txt.in but requires also handling linking in the main CMakeLists.txt file.


# v1.5.2 (20 November 2021)
Added compatibility with the changes to Themisto's command line interface and new index file structure in Themisto v2.0.0.

## Changes
- Changed the --themisto-index argument so that the program will abort with an error telling the user to rerun mSWEEP without --themisto-index if Themisto v2.0.0 index format is detected.

## Documentation
- Updated documentation with usage instructions for both Themisto <=v1.2..0 and >=v2.0.0.


# v1.6.0 (15 November 2021)
November is sometimes in May edition.

## New features
- Added MPI support and instructions for using it.
- Added support for reading in likelihoods written with the --write-likelihoods toggle (resolves #12).

## Changes
- Many internal changes and code refactoring.
- Use a new implementation of the model fitting code from rcgpar, which contains tests, better multiprocessing support, and a distributable (MPI compatible) version of the model fitting code.

## Bugfixes
- Fixed --print-probs so that it always prints to cout like the documentation says.


# v1.5.1 (9 November 2021) 
Finally published edition.

New features
- New --version toggle prints the version of the program.
- New --cite toggle prints the citation information for the mSWEEP article in Wellcome Open Research.

Documentation
- Added info about the doi for specific versions of mSWEEP to the readme file.

Build pipeline changes
- Download dependencies that are used by mSWEEP and/or some other dependencies only once and reuse them.
- Download cxxio when building instead of shipping with mSWEEP.

Files restructuring
- Moved config files from the main folder into config/.

Code restructuring
- Renamed main.cpp to mSWEEP.cpp.
- Use functions from dependencies when available instead of copying them to the mSWEEP source code.


# v1.5.0 (15 October 2021) 
Fall foliage edition: code restructuring and new features.

## New Features
Options to extract the likelihood matrix that mSWEEP uses internally:
- --write-likelihood: output the likelihood matrix in tab separated matrix format. Will write to a file with the _likelihoods.txt suffix if -o is specified, otherwise the matrix will be emitted to cout.
- --write-likelihood-bitseq: same as above but the output will be in a format that is compatible with BitSeq's estimateExpression and estimateVBExpression programs. Files from this toggle will have the _bitseq_likelihoods.txt suffix.

Added --no-fit-model toggle to skip the relative abundance estimation part:
- --no-fit-model: skip estimating the relative abundances. Useful if only the likelihood matrix is needed.

Support supplying multiple groupings via the -i or --groups-list toggles:
- Several groupings can be supplied by appending them as columns to the argument given by either the -i or the --groups-list options.
- The column delimiter is defined by the --groups-delimiter argument (default: tab-separated.).
- If there are several groupings and output to file is requested, the output will be written to the file specified by the -o argument but with the column index appended. Otherwise the results from all runs will print to cout.

Bugfixes
- Removed the extra line at the end of output when running in bootstrap mode.

Internal changes
- Some code restructuring to make adding new features easier.
- Hopefully improved code readability and a bit of documentation.
- Renamed some variables and functions that used the old "bitfields" naming scheme.
- Resolved some compiler warnings that arose when compiling with -Wall -Wextra -Wpedantic.
- Made several integer types explicit with (u)int32_t style typing.
- The Grouping and Reference structs have been separated and made into proper classes.


# v1.4.0 (10 March 2020)
Beware the clichés of software naming edition.

## New features
- Support parallel processing through the '-t' flags with excellent scaling in larger problems.
- Add possibility to match the input grouping indicators to the fasta file through the '--fasta' and '--groups-list' options.
- Add the '--bootstrap-count' option which allows resampling fewer input alignments than the original sample contains.
- Add possibility to specify the initial random seed for bootstrapping through the '--seed' option.
- Support reading in files compressed with bz2 or lzma if compiled on a machine that supports them.

## Better error checking
- Validate that all input and output files exist and are accessible.
- Add possibility to validate the input grouping indicators when using Themisto pseudoalignments (resolves #4 ).
- Catch errors in several places that escaped in earlier versions.
- More informative error messages in the above-mentioned cases.

## More efficient resource usage
- Parallel proceessing in the RCG optimization using OpenMP.
- Memory usage reduced by ~40% and in large problems.
- Single core performance increased by ~10% in large problems.

## Better build pipeline
- Download dependencies when running cmake.
- Build without OpenMP if it is not supported.
- More aggressive compiler optimization flags.
- Support build and optimization with the Intel C compiler.

## Internal changes
- Improve code structure and legibility.
- Use an external library (telescope) to read in pseeudoalignments from both kallisto or Themisto.
- Better internal storage for the pseudoalignments.
- Change the (rareish) reset step in the RCG optimization to be computationally more expensive but consume significantly less memory.
- Separate bootstrap and regular sample processing classes.


# v1.3.2 (30 January 2020) 
Fix working with a grouped Themisto index.

- Add instructions how to use either a grouped or ungrouped index.
- mSWEEP will now not attempt to infer the grouping.
- Instead, everything should be handled by modfying the file supplied with -i.


# v1.3.1 (21 January 2020)
- Fix compilation issues on some systems.


# v1.2.2 (3 September 2019)
Quality-of-life improvements, including:

- Bootstrapping output format is now similar to estimation without.
- Add the number of bootstrap iterations to the output file.
- Print a status indicator when running bootstrapping.
- Internal changes to code structure.


# v1.1.0 (17 December 2018)
## Prepublication edition
This is the version that was used to run experiments in the [mSWEEP preprint (2019)](https://www.biorxiv.org/content/10.1101/332544v2), and the first release to print the version number when ran.
