# PredLMM

## Data Description

The example data files provided are generated following the simumlation of section (3.2.1) of the main manuscript. In the folder named "Data", there are 

* a phenotype file: "example_pheno.csv"
* a covariate file: "example_covar.csv"
* PLINK Binary files: "example_geno"
* GCTA GRM files: "example_grm"

There are 5000 many individuals and 10,000 many SNPs. First two columns of the phenotype and covariate file have the family ID (FID) and individual ID (IID) of each individual. 
There are only one phenotype and one single covariate (intercept term). With the binary files, Genetic Relationship Matrix (GRM) files have been computed using GCTA. It is to be craefully noted that the order of the individuals in all the files (phenotype, covariate, GRM) have to be the same.


## Code Usage

In a general scenario, one should compute the GRM files first with the PLINK Binary files using GCTA software. Then, follow the jupyter notebook: PredLMM_notebook for estimating heribaility and variance of a phenotype adjusting the availble covariates. 



