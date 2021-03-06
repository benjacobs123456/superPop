% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{project_pcs}
\alias{project_pcs}
\title{Project UKB individual-level genotypes onto 1kg principal component space}
\usage{
project_pcs(
  bfile = "ukb_merged_chroms",
  out = "projected_pcs",
  n_threads = 1,
  path_to_plink_exe = "module unload plink; module load plink/2.0-20200328; plink2 "
)
}
\arguments{
\item{bfile}{a plink1 binary fileset containing genotypes for individuals whose ancestral superpopulation you want to estimate.}

\item{out}{a path to the output file}

\item{n_threads}{number of computing threads to use, defaults to 1}

\item{path_to_plink_exe}{the path to your plink2 executable}
}
\value{
a '.sscore' file with the projected PC values for individuals in your dataset in 1kg PC space.
}
\description{
Project UKB individual-level genotypes onto 1kg principal component space
}
\examples{
bfile = "ukb_merged_chroms"
path_to_plink = "plink.exe"
project_pcs()
}
