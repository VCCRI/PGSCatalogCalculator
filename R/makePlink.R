getMakePlink <- function(inVCF, inYaml="sample.yaml"){
  inOut <- setPlink(inVCF, inYaml="sample.yaml")
  if(is.null(inOut)){
    return(list(outFile=inOut, midSamples=NULL))
  } else if(grepl("\\.irem$", inOut)){
    midSamples <- getMindSamples(inOut)
    return(outFile=gsub("\\.irem$", "", inOut), midSamples=midSamples)
  } else {
    return(list(outFile=inOut, midSamples=NULL))
  }
}

getMindSamples <- function(inFile){
  inData <- data.table::fread(inFile, stringsAsFactors=FALSE, header=F)
  ## Get individual IID instead of family IDs
  return(inData$V2)
}

setPlink <- function(inVCF, inYaml="sample.yaml"){
  #inPlink <- "/g/data/jb96/software/plink_1.9_linux_x86_64_20181202/plink"
  inPlink <- if(is.null(inYaml)) Sys.getenv("plink") else yaml::read_yaml(inYaml)$plink
  ##ToDO prettifry paste
  if(is.null(inVCF)) return(NULL)
  if(grepl(".vcf.gz", inVCF)){
    inType <- paste("--vcf", inVCF, "--vcf-half-call h")
    outFile <- gsub("\\.vcf\\.gz", "_plink", inVCF)
  } else {
    inType <- paste("--bcf", inVCF)
    outFile <- gsub("\\.bcf", "_plink", inVCF)
 }
  #baseCommand <- paste(inPlink, inType, " --allow-extra-chr --maf 0.05 --mind 0.1 --geno 0.1 --hwe 1e-6 --vcf-filter --make-bed --chr 1-22 XY --memory 4096 --out", outFile)
  system2(command=inPlink, args=c(inType, "--allow-extra-chr", "--mind", "0.1", "--geno", "0.05",  "--vcf-filter", "--make-bed", "--chr", "1-22 XY", "--memory", "4096","--out", outFile), stdout=FALSE)
  if(!(file.exists(paste0(outFile, ".bim"))) & file.exists(paste0(paste0(outFile, ".irem")))){
    return(gsub("$", ".irem", outFile))
  } else if (!(file.exists(paste0(outFile, ".bim")))){
    return(NULL)
  } else {
    return(outFile)
  }
}

  

getMergePlink <- function(inControl, inDisease,inYaml="sample.yaml"){
  #inPlink <- "/g/data/jb96/software/plink_1.9_linux_x86_64_20181202/plink"
  inPlink <- if(is.null(inYaml)) Sys.getenv("plink") else yaml::read_yaml(inYaml)$plink
  outDir <- if(!(is.null(yaml::read_yaml(inYaml)$tempDir))){
    gsub("\\/$", "",yaml::read_yaml(inYaml)$tempDir)
  } else if(!(is.null(yaml::read_yaml(inYaml)$outputDir))){
    gsub("\\/$", "",yaml::read_yaml(inYaml)$outputDir)
  } else {
    dirname(inDisease)
  }
  outFile <- paste(outDir, gsub("_plink","_merge_plink", basename(inDisease)),sep="/")
  inDis <- paste0(inDisease, ".bed")
  inCont <- paste0(inControl, ".bed")
  if(!(file.exists(inDis))) inDisease <- gsub("_plink", "_plink-temporary", inDisease)
  if(!(file.exists(inCont))) inControl <- gsub("_plink", "_plink-temporary", inControl)
  if(!(file.exists(inDis)) & !(file.exists(inCont))){
    return(NULL)
  } else if(file.exists(inCont) & file.exists(inDis)){
    allCont <- unlist(lapply(c("bed", "bim", "fam"), function(x) paste(inControl, x, sep=".")))
    #baseCommand <- paste(inPlink, inType, " --allow-extra-chr --maf 0.05 --mind 0.1 --geno 0.1 --hwe 1e-6 --vcf-filter --make-bed --chr 1-22 XY --memory 4096 --out", outFile)
    system2(command=inPlink, args=c("--bfile", inDisease, "--bmerge", allCont, "--make-bed", "--allow-no-sex", "--memory", "4096","--out", outFile), stdout=FALSE)
    return(outFile)
  } else if(file.exists(inCont) | file.exists(inDis)){
    inOut <- c(file.exists(inCont), file.exists(inDis))
    return(c(inControl, inDisease)[which(inOut)])
     
  }

}


getRow <- function(inControl, inDisease){
  it <- itertools::ihasNext(product(inControl, inDisease))
  while (itertools::hasNext(it)) {
    x <- nextElem(it)
    makeFamFile(x[[1]]$inFile, x[[2]]$inFile)
    getMergePlink(x[[1]]$inFile, x[[2]]$inFile)
  }
}

checkPheno <- function(inFam){
  #Check if Pheno is Case Control
  return(all(c(1, 2) %in% as.numeric(inFam$V6)))
}

enforcePheno <- function(inFam){
  inFam$V3 <- 0
  inFam$V4 <- 0
  return(inFam)
}

writeFam <- function(inFam,inFamFile){
  fwrite(inFam, inFamFile, col.names =F)
}

#plinkFile <- getMakePlink("/g/data/jb96/sk3015/1000GDumpDecomp/europe_onechrom_DBN.vcf.gz")
#plinkFile <- getMakePlink("/g/data/ra5/sk3015/plinkMergeNoDup/plinkInput/2020-06-04/one_chrom_filt_plink.vcf.gz")
