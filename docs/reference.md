# Building a databases for mSWEEP
## Custom database
mSWEEP supports using a custom reference database. A typical workflow
for constructing the custom database might proceed as follows

1. Gather assembled sequences for the species of interest. If you are unsure of what species are present in your sample. use
taxonomic profiling tools like [MetaPhlAn2](https://doi.org/10.1038/nmeth.3589) to identify them.

2. a) Provide a grouping for the assemblies (e.g. sequence types,
   clonal complexes, or the output of some clustering
   algorithm.). [PopPUNK](https://github.com/johnlees/PopPUNK) is the
   recommended choice but if your species have an established multi-locus
   sequence typing scheme, that can also be used.

3. (Optional) Filter out reference sequences that cannot be
   reliably assigned to a group (eg. the sequence type cannot be
   determined) and perform other appropriate quality control
   measures. We recommend filtering assemblies using
   [checkm](https://github.com/Ecogenomics/CheckM) with a completeness
   threshold of >90% and contamination <5%.

4. (Optional) Set up the
   [demix\_check](https://github.com/tmaklin/coreutils_demix_check)
   index for additional QC checking of the bins produced by
   mSWEEP/mGEMS.

5. Index the database with Themisto.

## Prebuilt databases
For prebuilt reference databases, please refer to the supplementary datasets in the article [https://www.nature.com/articles/s41467-022-35178-5](https://www.nature.com/articles/s41467-022-35178-5).
