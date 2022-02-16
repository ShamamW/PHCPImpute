# PHCPImpute

An R package for imputation of pseudo-haplotypes data, based on the PHCP algorithm. \
The main function takes a pseudo-haplotype target and a set of reference genomes and returns the posterior probabilities for the possible genotypes in each SNP.
Full details described in ...

## Installation
The package was built under R version 4.1.0 \
To install the current packages, the devtools package should be installed. \
Additional required packages for usage: dplyr. 
```
# If devtools is not installed
install.packages("devtools") 
devtools::install_github("ShamamW/PHCPImpute")
library(PHCPImpute)
#if dplyr is not installed
install.packages("dplyr")
```


## Usage

```
phcp_impute(haps, ancient_tped, genetic_map, chr, choose_donors = F, variants)
```

## Input
### Required arguments
* haps - a data frame with phased reference genomes in .haps format, which is the output of SHAPEIT2. I.e., the first five columns are: \
&nbsp; 1. chromosome number without characters, e.g., 1 and not chr1 \
&nbsp; 2. variant ID \
&nbsp; 3. physical position in base-pairs \
&nbsp; 4. Allele 0 \
&nbsp; 5. Allele 1 \
The next columns are the phased haplotypes.

* ancient_tped - a data frame with the target pseudo-haplotype sequence in a Plink tped format. I.e, the first four columns are: \
&nbsp; 1. chromosome number without characters \
&nbsp; 2. variant ID \
&nbsp; 3. genetic distance \
&nbsp; 4. physical position in base-pairs \
The next two columns are two identical columns of the pseudo-haplotype sequence. 

* genetic_map: Genetic map for the reference genomes in Morgans. The file should have 3 columns: \
&nbsp; 1. chromosome name without characters\
&nbsp; 2. genetic distance (in Morgans) \
&nbsp; 3. position (in bp) 

* chr: integer. the number of chromosome to impute, without characters, e.g., 1 and not chr1

* genotyped: logical. Should the output include imputation results for SNPs that were genotyped in the pseudo-haplotype target? (defaults: FALSE).
* choose_donors: the number of donors' haplotypes to choose from the reference data. If set to FALSE, all donors will be used (default: FALSE).

### Optional arguments
* variants: a data frame with a list of variants from the reference panel to impute. \
The file should have the following columns with the following headers (additional columns will be ignored): \
&nbsp; 1. **chr** - chromosome number without characters \
&nbsp; 2. **bp** - physical position in base-pairs \
&nbsp; 3. **Ref** - reference allele \
&nbsp; 4. **Alt** - alternate allele 
* freq_donors: the maximal minor allele frequency in variants that are used for ranking and choosing the most informative donors. Variants with higher allele frequency will not be used for ranking the donors but will be used in the imputation (default: 0.5).

## Output
The output is a list of data frames with the posterior probabilities of the genotypes in the imputed SNPs. I.e., each SNP is assigned with the probabilities for Ref/Ref, Ref/Alt and Alt/Alt alleles. \
If the argument "genotyped" is set to TRUE, the list will contain a data frame called genotyped_variants with the results for the SNPs that were covered in the pseudo-haplotype target. \
If a data frame with variants to be imputed was supplied to the function ("variants" argument), the output list will include a data frame called chosen_variants with the posterior probabilities for these variants.


