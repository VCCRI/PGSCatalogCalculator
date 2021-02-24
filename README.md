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

* Input VCF (Case - Required)
* Reference Sequence Fasta File
* Filled in sample.yaml, file used to point package to relevant prerequisites
* Optional - Control VCF File

## Sample Run
```
require(PGSCatalogDownloader)
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
grabScoreId(inFile='sample.vcf.gz', inRef='/g/data/jb96/References_and_Databases/hs37d5.fa/hs37d5x.fa', inPGSID='PGS000073', inCL=cl, inControl='sample.vcf.gz')
parallel::stopCluster(cl)
```
## Input Parameters for grabScoreID

* inFile = Input Case VCF
* inRef = Reference Sequence FASTA File
* inPGSID = PGS ID that you want to calculate the score for
* inPGSIDS = File that has newline separated list of PGS IDs that you want to calculate the score for
* inCL = Cluster that will be used to run the package
* inControl = Control VCF (Optional)

Please note that package looks for "sample.yaml" file in the current working directory to ensure that it references the correct packages and environment variables

Please find an example file in this repo: https://github.com/VCCRI/PGSCatalogDownloader/blob/master/sample.yaml

Please define the output directory in the YAML file otherwise the tool will output files in the current R working directory

Please note that this tools creates intermediate files, error messages relating to these intermediate files should be ignored

## Output

The tool will generate a boxplot, `boxplot.png`, and CSV that displays the relative risk of patient and raw scores,`sample_out.csv`

These can be viewed concurrently by accessing `dashboard.Rmd`

These files demonstrate the stratified risk of samples against control samples and risk of each sample for the last condition respectively.
