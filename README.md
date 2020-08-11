# PredLMM

## Notebook Description

The jupyter notebook titled as PredLMM_notebook contains all the steps for implementing PredLMM on an example dataset provided in "Data" folder. The main module containing all the necessary python functions can be found inside the folder "PredLMM".

## Data Description

The example data files provided are generated following the simulation of section (3.2.1) of the main manuscript. In the folder named "Data", there are 

* a phenotype file: "example_pheno.csv"
* a covariate file: "example_covar.csv"
* PLINK Binary files: "example_geno"
* GCTA GRM files: "example_grm"

There are 5000 many individuals and 10,000 many SNPs. First two columns of the phenotype and covariate file have the family ID (FID) and individual ID (IID) of each individual. 
There are only one phenotype and one single covariate (intercept term). With the binary files, Genetic Relationship Matrix (GRM) files have been computed using GCTA. It is to be craefully noted that the order of the individuals in all the files (phenotype, covariate, GRM) have to be the same.


## Code Description

In a general scenario, one should compute the GRM files first with the PLINK Binary files using GCTA software. Then, follow the jupyter notebook: PredLMM_notebook for estimating heribaility and variance of a phenotype adjusting the availble covariates. 

In PredLMM_notebook, we estimate the heritability and variance of the phenotype twice:

* by fitting a LMM with a random subsample (of size subsample_size) and 
* by fitting the actual PredLMM algorithm treating the selected random subsample as the set of knots.

In our example, we considered the subsample_size (size of the set of knots) to be 500. One can change that depending upon the total population size. We have seen using 10% of the total population size as subsample_size works reliably. 

Along with the estimates, we also display the time taken by both the approaches.



