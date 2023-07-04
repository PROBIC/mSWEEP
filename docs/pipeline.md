# General pipeline
This file describes how to run mSWEEP and mGEMS on a custom reference database following some best practices we have established. This assumes that you already have a set of reference assemblies to use. For more info on obtaining the reference assemblies, see steps under [Usage -> Reference data](/docs/reference.md)).

## Grouping reference sequences
Use [PopPUNK](https://github.com/bacpop/PopPUNK) to cluster the reference sequences
```
sed 's/^.*\///g' reference_assembly_paths.txt | sed 's/[.].*$//g' > reference_assembly_names.txt
paste reference_assembly_names.txt reference_assembly_paths.txt > poppunk_input_list.tsv

poppunk --create-db --r-files poppunk_input_list.tsv --output poppunk_fit

poppunk --fit-model bgmm --K 4 --ref-db poppunk_fit
```

Match the poppunk clusters with the assembly paths
```
echo -e "id\tcluster\tassembly" > ref_info.tsv
join -1 1 -2 1 <(sed '1d' poppunk_sketch/poppunk_sketch_clusters.csv | tr ',' '\t' | sort) <(sort poppunk_input_list.tsv) | tr ' ' '\t' >> ref_info.tsv
```

Check the [poppunk documentation](https://poppunk.readthedocs.io) for information on how to create a good clustering.

## Building a pseudoalignment index
Extract the reference assembly paths
```
cut -f3 ref_info.tsv | sed '1d' > ref_paths.txt
```

Index the reference sequences for pseudoalignment with [Themisto](https://github.com/algbio/themisto):
```
themisto build -k 31 -i ref_paths.txt -o themisto_index --temp-dir tmp
```

## Pseudoaligning reads
Error correct paired-end reads with [fastp](https://github.com/opengene/fastp)
```
fastp --in1 raw_fwd.fastq.gz --in2 raw_rev.fastq.gz --out1 corr_fwd.fastq.gz --out2 corr_rev.fastq.gz -c --thread 2
```

Pseudoalign the reads and compress the alignment file with [alignment-writer](https://github.com/tmaklin/alignment-writer)
```
ntargets=$(wc -l ref_paths.txt | cut -f1 -d' ')
nreads=$((`gunzip -c reads_fwd.fastq.gz | wc -l` / 4 ))
themisto pseudoalign -q corr_fwd.fastq.gz -i themisto_index --temp-dir tmp -t 1 | alignment-writer -n $ntargets -r $nreads > fwd.aln
themisto pseudoalign -q corr_rev.fastq.gz -i themisto_index --temp-dir tmp -t 1 | alignment-writer -n $ntargets -r $nreads > rev.aln
```

## Estimate abundances
Extract the cluster indicators
```
cut -f2 ref_info.tsv | sed '1d' > ref_clusters.txt
```

Estimate abundances with mSWEEP and bin reads for groups with abundance higher than 0.01 with
```
mSWEEP --themisto-1 fwd.aln --themisto-2 rev.aln -i ref_clusters.txt -t 2 --bin-reads --min-abundance 0.01 -o res
```

## Extract binned reads
Extract the reads with [mGEMS](https://github.com/probic/mgems)
```
bins=$(ls *.bin | tr '[[:space:]]' ',' | sed 's/,$//g')
mGEMS extract --bins $bins -r 215_1.fastq.gz,215_2.fastq.gz -o ./
```

## Setup demix\_check for quality checking
Setup the reference for [coreutils\_demix\_check](https://github.com/tmaklin/coreutils_demix_check) with
```
mkdir demix_ref
cd demix_ref
setup_reference.sh --ref_info ../ref_info.tsv --threads 2
cd ../
```

Check all reads that were binned with
```
for f in *.bin; do
	cluster=$(echo $f | sed 's/[.]bin$//g')
	fwd=$cluster"_1.fastq.gz"
	rev=$cluster"_2.fastq.gz"
	check_reads.sh --abundances res_abundances.txt --cluster $cluster --reference demix_ref --fwd $fwd --rev $rev --threads 2
	mv clu_score.tsv $cluster"_clu_score.tsv"
done
```
groups that receive a score of 1 or 2 are considered to pass the quality check (1 is better than 2) and can be used in assembly.
