# Imputato

Imputato is an approach to phase and impute autopolyploid individuals. Almost incidentally, this also includes diploids.

The main idea is to decompose all samples into individual haploids. Haploid phasing and imputation is done repeatedly.
Within a sample, the *posterior* for the other haploids are then used to adjust the *priors* for those variants
where genotypes have been recorded. The basic idea is to adjust the priors in such a way that the eventual posteriors
are consistent with the recorded genotypes.

## Usage
The code itself is found in `imputato.cpp`. This can be compiled with a a modern C++ compiler. Flags such as
`-fopenmp -march=native` are highly recommended for parallelism and performance. Right now, the code is sort of hand tuned
for the number of threads, but it should be possible to set environment variables to override.

You also need to have Eigen installed (`libeigen3-dev` in Ubuntu). If it's a proper install on your setup, you might not need
to specify any further include directory.

At this point, the ploidy level (4 in this branch) and file names (`potato_chr1.map` and `potato_missing.gen`) are hard-coded
in the source.

## File formats
The map is supposed to simply be a one number per line. The first line is the number of markers, following lines are positions
in cM. A single Imputato run should just cover a single linkage group/chromosome.

The genotypes consists of a single line defining the number of individuals, and then one line per individual. Each such
line should contain the same number of columns as the number of markers in the map file. Values should be -1 (missing), or
one of the digits 0, 1, 2, 3, 4 for allele count, in the case of ploidy 4.

In the end, an output file `potato.out` is generated which is highly similar to the `potato_missing.gen` input, but with
imputed data in all positions. Even in cases of high uncertainty, a call will be made. There is also an auxiliary file
`potato.vcflike`, which is marker-centric and provides internal probabilities for the individuals, split in haploids.

## Citation
So far, Imputato has only been presented as a poster at the 7th International Conference of Quantiative Genetics in Vienna
July 22-26 2024. The code itself is a workbench for experimentation, rather than a final tool. However, it is covered by
an MIT sourc code license, so if it's useful in this form, please use it, but try to reference the eventual paper for giving
proper academic credit.

## Author
The author of Imputato is Carl Nettelblad, carl.nettelblad@it.uu.se.

