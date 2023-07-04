# mSWEEP advanced usage
## Using a colored Themisto index
Instead of supplying the grouping to mSWEEP through the `-i` argument, it is possible to encode the grouping in the Themisto index. This is useful when performing hierarchical analyses with mSWEEP, where reads are first assigned to a species and then to some within-species groups.

Using a colored index results in a degenerated likelihood function in mSWEEP, which typically results in some false positives results. These are easily filtered out if a further within-species step is applied.

Embedding the grouping to the Themisto index results in effectively
treating any pseudoalignment in the group as if the read aligned to
all sequences in the group. This approach will produce different
results than the treating each sequence separately, but vastly reduces
the RAM, CPU, and disk space requirements for running Themisto and
mSWEEP.

The commands below assume you are using the example dataset from the [example use case](/docs/example.md)

### Embedding colors in the Themisto index
Create a list containing paths to the input assemblies
```
ls -d $(pwd)"/assemblies/"*.fasta.gz > input_sequences.txt
```

Concatenate the input assemblies according to their group assignments in the `clustering.txt` file
```
i=1
while read line; do
	group=$(sed -n "$i p" clustering.txt)
	gunzip -c $line >> $group"_concatenated.fasta"
	i=$((i + 1))
done < input_sequences.txt
```

Create a list of the concatenated files to use as input for `themisto build`
```
ls -d $(pwd)"/"*"_concatenated.fasta" > input_multifastas.txt
```

Build the colored index with Themisto
```
mkdir tmp
themisto build -k 31 -o colored_index --temp-dir tmp -t 2 --mem-gigas 4 -i input_multifastas.txt
```
this will create a colored index with 4 pseudoalignment targets.

### Pseudoalign reads
Pseudoalign reads against the colored index with
```
themisto pseudoalign -q 215_1.fastq.gz -i colored_index --temp-dir tmp -t 2 > 215_1.txt
```

### Using embedded colors with mSWEEP
Extract cluster indicators from the `input_multifastas.txt` file with awk
```
awk '!seen[$0]++' clustering.txt > unique_clusters.txt
sed 's/^.*\///g' input_multifastas.txt | sed 's/_concatenated[.]fasta$//g' > colored_clustering.txt
```

Run mSWEEP and bin reads for the subsequent within-species step
```
mSWEEP --themisto 215_1.txt -i colored_clustering.txt --bin-reads
```

## Bootstrapping confidence intervals
mSWEEP can be used to produce confidence intervals for the abundance
estimates by bootstrapping the pseudoalignment counts and rerunning
the abundance estimation a number of times. This can be done
automatically by adding the '--iters' option to running mSWEEP:
```
mSWEEP --themisto-1 215_1.txt --themisto-2 215_2.txt -i clustering.txt -t 2 -o 215 --iters 100
```

The bootstrapped abundance estimates will be appended to the output
file '215_abundances.txt' as new columns and can be used to calculate
confidence intervals for each of the abundance estimates.

Bootstrapping in mSWEEP by default resamples from all input pseudoalignments. If
you wish to resample fewer pseudoalignments, supply the
`--bootstrap-count <positive integer>` option with the number of resamples to perform
```
mSWEEP --themisto-1 215_1.txt --themisto-2 215_2.txt -i clustering.txt -t 2 -o 215 --iters 100 --bootstrap-count 1000
```
This will result in wider confidence intervals because the bootstrap
iterations have less data available to them.

It is also possible to set the initial random seed with the `--seed
<positive integer>` option, which enables replicating the bootstrap
results across multiple runs.
