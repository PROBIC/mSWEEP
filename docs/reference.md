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
