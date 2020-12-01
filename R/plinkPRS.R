plinkGRS <- function(inFile, inGRS, inYaml="sample.yaml"){
  #inPlink <- "/g/data/jb96/software/plink_1.9_linux_x86_64_20181202/plink"
  inPlink <- if(is.null(inYaml)) Sys.getenv("plink") else yaml::read_yaml(inYaml)$plink
  outFile <- gsub("_plink", "_score", inFile)
  ##TODO prettifry paste
  system2(command=inPlink, args=c("--memory","2048","--bfile", inFile,"--double-id", "--allow-no-sex", "--score", inGRS, "header", "sum", "--out", outFile),stdout=FALSE)
  outFile <- paste0(outFile, ".profile")
  return(outFile)
}

prsiceGRS <- function(inFile, inGRS, outDir){
  # TODO Properly implement frequentist p value as it relies on GWAS or PRS Catalog having it
  baseCommand <- paste("Rscript", baseDir, "/PRSice.R --prsice", baseDir, "/PRSice_linux --target", inFile,
    "--base", inGRS, "--snp rsid --A1 effect_allele --A2 reference_allele --stat", effect_weight,  "--pvalue frequentist_add_pvalue ",
    ## TODO Change
    "--beta ",
    "--out", outDir, "PRS_Catalog",
    "--binary-target T")
  system(command=baseCommand)
  return(outFile)
}