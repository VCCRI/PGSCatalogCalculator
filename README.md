A tool that will grab PGS Catalog Data and allow for automatic scoring of samples and comparison against the EU Control Sample Set or Input Control Sample VCF

Required Software:
* R - 3.6.1
* Python - 2.7.16
* VT
* Bcftools
* PLINK 1.9

# Getting Started

Installation
```
python2 -m pip install -r https://raw.githubusercontent.com/VCCRI/PGSCatalogDownloader/master/requirements.txt
R
devtools::install_github("VCCRI/PGSCatalogDownloader")
```

## Required Files 

* Input VCF (Case - Required, Control - Optional, can provide your own)
* Reference Sequence Fasta File
* Filled in sample.yaml, file used to point package to relevant prerequisites

## Sample Run
```
require(PGSCatalogDownloader)
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
grabScoreId(inFile='sample.vcf.gz', inRef='/g/data/jb96/References_and_Databases/hs37d5.fa/hs37d5x.fa', inPGSID='PGS000073', inCL=cl)
# inFile = Input Case VCF
# inRef = Reference Sequence Fasta File
# inPGSID = PGS ID that you want to calculate the score for
# inPGSIDS = File that has newline seperated list of PGS IDs that you want to calculate the score for
# cl = Cluster that will be used to run the package
parallel::stopCluster(cl)
```

Please note that package looks for "sample.yaml" file in the current working directory to ensure that it references the correct packages and environment variables

Please find an example file in this repo: https://github.com/VCCRI/PGSCatalogDownloader/blob/master/sample.yaml

## Output

The tool will generate a boxplot, `quantile/boxplot.png`, and CSV that displays the relative risk of patient,`sample_out.csv`. 

These can be viewed concurrently by accessing `dashboard.Rmd`

These files demonstrate the stratified risk of samples against control samples and risk of each sample for the last condition respectively.
