mergeVCF <- function(inControl, inDisease){
  # TODO Write SNP Filter only
  outFile <- gsub("\\.vcf\\.gz", "_merged.vcf.gz", inDisease)
  inControlIndex <- file.exists(paste0(inControl, ".tbi"))
  inDiseaseIndex <- file.exists(paste0(inDisease, ".tbi"))
  if(!inDiseaseIndex) baseIndex(inDisease, inYaml)
  if(!inControlIndex) baseIndex(inControl, inYaml)
  baseCommand <- paste("bcftools merge --missing-to-ref -m none", inControl, inDisease, "-o", outFile, "-O z --threads 5")
  system(command=baseCommand)
  baseIndex(outFile)
  return(outFile)
}

makeFamFile <- function(inControl, inDisease, inOut){
  outFile <- gsub("(\\.vcf\\.gz)|\\.bcf", "_plink.fam", inOut)
  baseCommand <- system(paste("bcftools query -l",  inOut), intern=TRUE)
  # TODO Check if there is any overlap in control or diseased samples
  controlSamples <- system(paste("bcftools query -l", inControl), intern=TRUE)
  inDiseaseSamples <- system(paste("bcftools query -l", inDisease), intern=TRUE)
  inFam <- data.table::data.table(V1=baseCommand, V2=baseCommand, V3=0, V4=0, V5=0, V6=-9)
  inFam[V1 %in% controlSamples,"V6"] <- 1
  inFam[V1 %in% inDiseaseSamples,"V6"] <- 2
  return(list(inFam=inFam, outFile=outFile))
  #data.table::fwrite(inFam, outFile, col.names = F,sep=" ")
  #return(outFile)
}

writeFam <- function(inObj){
  if(is.null(inObjec)){
       return(FALSE)
  } else {
  data.table::fwrite(inObj$inFam, inObj$outFile, col.names = F,sep=" ")
  return(TRUE)
  }
}


filterMerged <- function(inFile, inName, inSNPs, inYaml="sample.yaml"){
  outFile <- gsub("\\.vcf\\.gz", paste0("_",inName, "_filt.vcf.gz"), inFile)
  bcftools <- if(is.null(inYaml)) Sys.getenv("bcftools") else yaml::read_yaml(inYaml)$bcftools
  if(!file.exists(outFile)){
    system2(command=bcftools, args=c("view", "-i", paste0("ID=@", inSNPs), inFile, "-o", outFile, "-O", "z","--threads", "1"), stdout=FALSE)
    print("Filtering BCF")
    #baseIndex(outFile)
    system2(command=bcftools, args=c("index", "-f", outFile),stdout=FALSE)
    print("Filtered Data")
    Sys.time()
    print(outFile)
  }
  return(outFile)
}


filterMergedReg <- function(inFile, inName, inRegion, inYaml="sample.yaml"){
  outFile <- gsub("\\.vcf\\.gz", paste0("_",inName, "_reg.vcf.gz"), inFile)
  bcftools <- if(is.null(inYaml)) Sys.getenv("bcftools") else yaml::read_yaml(inYaml)$bcftools
  system2(command=bcftools,args=c("view", "-R", inRegion, inFile, "-o", outFile, "-O", "z", "--threads", "1"))
  system2(command=bcftools, args=c("index", "-f", outFile),stdout=FALSE)
  return(outFile)
}
filterMergedPos <- function(inFile, inName, inSNPs, inYaml="sample.yaml"){
  ## TODO Implement using python script for the time being
  bcftools <- if(is.null(inYaml)) Sys.getenv("bcftools") else yaml::read_yaml(inYaml)$bcftools
  outFile <- gsub(paste0("(_", inName, "_filt", ")?\\.vcf\\.gz"), paste0("_",inName, ".bcf"), inFile)
  ##TODO Implement if file exists check
    print("PGS Catalog Start")
    #baseCommand <- paste0("python3 cyvf_pgs_catalog.py --score ", inSNPs, " --header --out ", outFile, " --vcf ", inFile)
    system2(command="python2", args=c(system.file("extdata", "vcf_pgs_catalog.py", package="PGSCatalogDownloader"),"--score", inSNPs, "--header", "--out", outFile, "--vcf", inFile), stdout=FALSE)
    print("PGS Catalog End")
    Sys.time()
    system2(command=bcftools, args=c("index", "-f", outFile),stdout=FALSE)
    print("Start Index")
    Sys.time()
  return(outFile)
}

## TODO Write QC Check for multivariate
#fil
#mergedFile <- mergeVCF(inControl="/g/data/jb96/sk3015/1000GDumpDecomp/europe_onechrom_DBN_snp.vcf.gz", inDisease="/g/data/ra5/sk3015/plinkMergeNoDup/plinkInput/2020-06-04/one_chrom_filt_snp.vcf.gz")

#filterMerged <- filterMerged(
#source('makePlink.R')
#plinkFile <- getMakePlink(mergedFile)

#makeFamFile(inControl="/g/data/jb96/sk3015/1000GDumpDecomp/europe_onechrom_DBN.vcf.gz", inDisease="/g/data/ra5/sk3015/plinkMergeNoDup/plinkInput/2020-06-04/one_chrom_filt.vcf.gz", inOut=mergedFile)
