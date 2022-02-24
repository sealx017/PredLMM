# PredLMM

## Notebook Description

* The jupyter notebook titled as "PredLMM_notebook" contains all the steps for implementing PredLMM on an example dataset provided in "Data" folder. The notebook reads the phenotype, covariates, and the GRM files and estimates heritability based on those. 

* The PredLMM algorithm requires only a few particular blocks of the Genetic Relationship Matrix (GRM) and thus, computing the full GRM is not necessary. In the notebook named "PredLMM_notebook_with_GRM_computation_included", we explain how one can compute only the blocks of the GRM necessary for fitting the PredLMM algorithm and perform the subsequent analysis. We also allow incorporating user-specified SNP-weights e.g., LD-based weights computed using the software "LDAK5" into the estimation of the GRM-blocks.

* The main module containing all the necessary python functions can be found inside the folder named "PredLMM". 


## Data Description

The example data files provided are generated following the Section simulation study 1 of the main manuscript. In the folder named "Data", there are, 

* a phenotype file: "example_pheno.csv"
* a covariate file: "example_covar.csv"
* PLINK Binary files: "example_geno" (.bed, .bim, .fam)
* GCTA GRM files: "example_grm" (.bin)

There are 5000 individuals and 10,000 SNPs. First two columns of the phenotype and covariate file have the family ID (FID) and individual ID (IID) of each individual. There are only one phenotype and one single covariate (intercept term). With the binary files, Genetic Relationship Matrix (GRM) files have been computed using GCTA. It is to be craefully noted that the order of the individuals in all the files (phenotype, covariate, GRM) have to be the same.


## Code Description

In a general scenario, one should compute the GRM files first with the PLINK Binary files using GCTA software. Then, follow the jupyter notebook: PredLMM_notebook for estimating heribaility and variance of a phenotype adjusting the availble covariates. 

In PredLMM_notebook, we estimate the heritability and variance of the phenotype twice:

* by fitting a LMM with a random subsample (of size subsample_size) and 
* by fitting the actual PredLMM algorithm treating the selected random subsample as the set of knots.

In our example, we considered the subsample_size (size of the set of knots) to be 500. One can change that depending upon the total population size. We have seen using 10% of the total population size as subsample_size works reliably. 

Along with the estimates, we also display the time taken by both the approaches.

### References

1. Speed, D., Hemani, G., Johnson, M. R., & Balding, D. J. (2012). Improved heritability estimation from genome-wide SNPs. The American Journal of Human Genetics, 91(6), 1011-1021.

2. Yang, J., Lee, S. H., Goddard, M. E., & Visscher, P. M. (2011). GCTA: a tool for genome-wide complex trait analysis. The American Journal of Human Genetics, 88(1), 76-82.

