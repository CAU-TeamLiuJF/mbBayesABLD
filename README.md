
=================

Version Updates
Version 1.03
Date: April 3, 2023

Description: This update enables the MB-BayesAS model to run based on the number of markers provided for each interval or based on a file that specifies the number of markers per interval.

1. Contact Information
liwn@cau.edu.cn

1. Installation
Download the source code and compile it using gcc, or directly run the pre-compiled version.

Dependencies: Requires the MKL library.

3. User Guide
(1) Parameter Explanation

```text
 Required parameters:
  -g, --raw=FILE             plink raw file, e.g. /path/plink.raw
  -l, --bfile=FILE           plink binary file prefix, e.g. /path/plink

 Optional parameters:
  -1, --VaSori=FILE          SNP (co)variance prior scale, phe/iden [phe]
  -2, --R0Sori=FILE          residual (co)variance prior scale, phe/iden [phe]
  -3, --VaSd=FILE            prior of (co)variance matrix of additive effect, wishart/invgamma [NULL]
  -4, --R0Sd=FILE            prior of (co)variance matrix of residual effect, wishart/invgamma [NULL]
  -5, --VaSSori=FILE         scale of SNP (co)variance prior scale, phe/iden [phe]
  -6, --R0SSori=FILE         scale of residual (co)variance prior scale, phe/iden [phe]
  -a, --id=INT               id position at phenotype file [1]
  -A, --varOut=FILE          posteriori mean of dispersion matrix [var.txt]
  -b, --binf=FILE            file contains number of SNPs per interval [NULL]
  -B, --burnin=INT           burnin steps [10]
  -c, --varC=INT             output first c predictions in MCMC chain [NULL]
  -d, --dfva=DOUBLE          degrees of freedom Va [npop + 1]
  -D, --dfve=DOUBLE          degrees of freedom R0 [npop + 1]
  -e, --region=CHR           Method for dividing snp into different groups, LD/fix [LD]
  -E, --effOut=FILE          posteriori mean of location effect [theta.txt]
  -f, --mapf=FILE            plink .map (.bim) use to define bins [NULL]
  -F, --fix=CHR              fixed effect position, e.g. '2 3' [NULL]
  -j, --bincol=INT           position of nSNPs in each bins at binf [1]
  -J, --maf=INT              SNPs with maf less than this value will be excluded [-0.01]
  -k, --varSeg=CHR           Method of allocating total genetic variance to bins, het/snp [het]
  -K, --geno=INT             Allow the missing proportion of the SNP [0.99]
  -L, --logf=FILE            redirect standard output to this file [NULL]
  -m, --miss=DOUBLE          missing value in phenotype [-99]
  -M, --mcmcOut=FILE         output predictions in MCMC chain [MCMC_var.txt]
  -N, --nsnp=FILE            Number of SNPs per interval [100]
  -o, --outdir=CHR           output directory [.]
  -p, --phef=FILE            phenotype file [NULL]
  -P, --mind=INT             Allow the missing genotype of the individual [0.99]
  -Q, --corf=FILE            correlation coefficient in a priori covariance matrix [NULL]
  -r, --report=INT           MCMC sampling report interval [100]
  -R, --ran=CHR              environmental random effect position [NULL]
  -s, --seed=LONG            random number seed [NULL]
  -S, --iter=INT             length of MCMC chain [100]
  -t, --thin=INT             thinning rate [10]
  -T, --thread=INT           number of threads for the OpenMP run-time library [1]
  -V, --gebvOut=FILE         breeding value output file [gebv.txt]
  -y, --phe=CHR              phenotypes position, e.g. '4 5' [last npop column(s)]

 Flag parameters:
  -C, --dontcenter           centering gene content matrix
  -H, --header               phenotype and bins files has header, [automatic detection]
  -n, --nocov                constrain the covariance between the residuals to zero
  -W, --varSamp              output sample of variances of marker effect

  -?, --help                 Give this help list
      --usage                Give a short usage message
      --version              Print program version
```

(2) Input Files
Genotype File: A PLINK raw format (*.raw) file. For details, refer to the PLINK file format description.

Phenotype File: A space-delimited plain text file containing individual IDs, fixed effects, and phenotypic values. Fixed effects can be represented as characters, with a maximum length of 20 characters. The order of individuals in this file should correspond to the first column of IDs in the raw file.

(3) Output Files
See parameter descriptions above for details.

(4) Example Code

```
BayesAS --bfile genotype --phef test.txt --nsnp 100
```
