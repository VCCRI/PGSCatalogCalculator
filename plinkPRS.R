plinkGRS <- function(inFile, inGRS){
  inPlink <- "/g/data/jb96/software/plink_1.9_linux_x86_64_20181202/plink"
  outFile <- gsub("_plink", "_score", inFile)
  ##TODO prettifry paste
  baseCommand <- paste(inPlink," --memory 2048 --bfile", inFile," --double-id  --allow-no-sex --score ", inGRS, "header sum --out", outFile)
  system(command=baseCommand)
  outFile <- paste0(outFile, ".profile")
  return(outFile)
}

prsiceGRS <- function(inFile, inGRS, outDir){
  # TODO Properly implement frequentist p value as it relies on GWAS or PRS Catalog having it
  baseDir <- "/g/data/jb96/sk3015/software/PRSice"
  baseCommand <- paste("Rscript", baseDir, "/PRSice.R --prsice", baseDir, "/PRSice_linux --target", inFile,
    "--base", inGRS, "--snp rsid --A1 effect_allele --A2 reference_allele --stat", effect_weight,  "--pvalue frequentist_add_pvalue ",
    ## TODO Change
    "--beta ",
    "--out", outDir, "PRS_Catalog",
    "--binary-target T")
  system(command=baseCommand)
  return(outFile)
}
