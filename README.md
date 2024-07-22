# Imputato

Imputato is an approach to phase and impute autopolyploid individuals. Almost incidentally, this also includes diploids.

The main idea is to decompose all samples into individual haploids. Haploid phasing and imputation is done repeatedly.
Within a sample, the *posterior* for the other haploids are then used to adjust the *priors* for those variants
where genotypes have been recorded. The basic idea is to adjust the priors in such a way that the eventual posteriors
are consistent with the recorded genotypes.

## Citation
So far, Imputato has only been presented as a poster at the 7th International Conference of Quantiative Genetics in Vienna
July 22-26 2024. The code itself is a workbench for experimentation, rather than a final tool. However, it is covered by
an MIT sourc code license, so if it's useful in this form, please use it, but try to reference the eventual paper for giving
proper academic credit.

## Author
The author of imputato is Carl Nettelblad, carl.nettelblad@it.uu.se.

