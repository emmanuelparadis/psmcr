<h2>psmcr: Pairwise Sequentially Markovian Coalescent With R</h2>

psmcr is an R port of [`psmc`](https://github.com/lh3/psmc). Currently, `psmcr` includes:

- The function `psmc()` that runs the PSMC with (almost) all options as the original program. The bootstrap can be run from this function, with the possibility to run the replications in parallel.

- A `plot` function to display the result of `psmc` with various options.

- The function `vcf2DNAbin` that prepares a consensus genome from a VCF file with the possibility to choose the individual if the VCF file is from a population study.

- The function `seqBinning` to bin the sequences.

See below for some details.

<h3>Usage</h3>

There are several parameters to consider in order to obtain a sensible output from the PSMC: it is recommended to see the help page `?psmc` and the original documentation and publication accompanying [`psmc`](https://github.com/lh3/psmc).

The present R package runs in a reasonable time (i.e., a few minutes or less) with an input size of ~1 Mb and 20 iterations. An input size of 1 Mb is equivalent to a genome size of 100 Mb if the default binning size (100) is used. This makes possible to assess different parameter settings. Once these are found, it is possible to run the program as batch (i.e., non-interactively), for instance on a remote machine, with a small script like:

```r
x <- readRDS("data.rds")
o <- psmc(x, niters = 30, B = 100, .....)
saveRDS(o, "output_run_psmc.rds")
```

where the file data.rds contains the input (eventually binned) sequences saved with the R function `saveRDS`, `niters` is the number of iterations of the PSMC, `B` is the number of bootstrap replications, and `.....` are other parameter settings. Alternatively, if you prefer to save the input sequences as FASTA, the first line can be replaced with:

```r
library(ape)
x <- read.FASTA("data.fas")
```

The script can then be run in the usual way with (depending on the operating system):

```
R CMD BATCH script.R nohup &
```

Once the job is completed, the output can be analysed in R with:

```r
o <- readRDS("output_run_psmc.rds")
plot(o)
```

Of course, since the object `o` is a list, it can be analysed with standard R operations.

<h3>Getting Input From a VCF File</h3>

The input for `psmc` is the (eventually binned) consensus sequence of a diploid individual.  It is possible to build this sequence from:

- a reference genome in a FASTA file,
- a VCF file with genotypes from one or several individuals.

The function `vcf2DNAbin` takes as input the VCF file and the reference genome and outputs a consensus sequence. The option `individual` selects which individual to consider in the VCF file in case there are more than one. The reference genome can be specified in different ways:

- by default, `vcf2DNAbin` looks for a reference in the VCF file and opens it (after downloading if it is a remote file);
- a "DNAbin" object;
- a name of a FASTA file.

The function `seqBinning` can be used to bin the consensus sequences before calling `psmc`. The default bin size is 100, and the output is a set of sequences with K if there is at least one ambiguous base (representing a heterozygous site) within the bin, T otherwise.
