
<!-- README.md is generated from README.Rmd. Please edit that file -->

# superPop

<!-- badges: start -->

<!-- badges: end -->

When analysing large genetic datasets like UK Biobank, a common step in
quality control is to ensure that the data you’re analyzing come from an
ancestrally-homogenous group. Sophisticated software exists to estimate
global ancestry proportions on a genome-wide scale. For quick,
exploratory analyses, it is often useful to quickly identify individuals
from a particular ancestral super-population.

superPop helps you to quickly guess which ancestral superpopulation
someone is from. It’s designed for use with UK Biobank data.

It works by projecting individuals into principal component space using
a common set of SNPs between 1,000 genomes and UKB. It helps you scale
these PCs so they are comparable. It then applies a pre-trained random
forest classifier to group each individual into their likely parent
population.

This is meant for quick analysis, and clearly further results should be
inspected in more detail for population structure.

## Installation

You can install superPop from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("benjacobs123456/superPop")
```

## Example

Before starting, please perform some basic SNP quality control on your
data, and combine your SNP genotype data into a genome-wide PLINK1
binary fileset.

First, the function project\_pcs() projects UKB individuals into PC
space (a wrapper for plink2’s –score).

``` r
library(superPop)
project_pcs(bfile = "ukb_merged_chroms",path_to_plink = "plink.exe")
# This will produce an output file with the projected PC values called "projected_pcs.sscore"
```

Second, we scale the PCs so they are comparable to those calculated with
1kg.

``` r
new_pcs = transform_pcs("projected_pcs.sscore")
# This will produce an output file with the projected PC values called "projected_pcs.sscore"
```

Then we can predict which ancestral cluster each individual comes from
using:

``` r
predict_ancestry(new_pcs)
# which will produce an output file, ancestry_estimates.tsv, containing estimates for each individual. 
```

You can plot these projections, grouped by predicted cluster, with:

``` r
plot_predicted_ancestry(projected_pcs = new_pcs, projected_ancestry = "ancestry_estimates.tsv")
```
