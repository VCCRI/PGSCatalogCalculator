source('normScript.R')
mergeVCF <- function(inControl, inDisease){
  # TODO Write SNP Filter only
  outFile <- gsub("\\.vcf\\.gz", "_merged.vcf.gz", inDisease)
  inControlIndex <- file.exists(paste0(inControl, ".tbi"))
  inDiseaseIndex <- file.exists(paste0(inDisease, ".tbi"))
  if(!inDiseaseIndex) baseIndex(inDisease)
  if(!inControlIndex) baseIndex(inControl)
  baseCommand <- paste("bcftools merge --missing-to-ref -m none", inControl, inDisease, "-o", outFile, "-O z --threads 5")
  system(command=baseCommand)
  baseIndex(outFile)
  return(outFile)
}

makeFamFile <- function(inControl, inDisease, inOut){
  require(data.table)
  outFile <- gsub("(\\.vcf\\.gz)|\\.bcf", "_plink.fam", inOut)
  baseCommand <- system(paste("bcftools query -l",  inOut), intern=TRUE)
  # TODO Check if there is any overlap in control or diseased samples
  controlSamples <- system(paste("bcftools query -l", inControl), intern=TRUE)
  inDiseaseSamples <- system(paste("bcftools query -l", inDisease), intern=TRUE)
  inFam <- data.table(V1=baseCommand, V2=baseCommand, V3=0, V4=0, V5=0, V6=-9)
  inFam[V1 %in% controlSamples,"V6"] <- 1
  inFam[V1 %in% inDiseaseSamples,"V6"] <- 2
  fwrite(inFam, outFile, col.names = F,sep=" ")
  return(outFile)
}

filterMerged <- function(inFile, inName, inSNPs){
  outFile <- gsub("\\.vcf\\.gz", paste0("_",inName, "_filt.vcf.gz"), inFile)
  if(!file.exists(outFile)){
    baseCommand <- paste0("bcftools view -i ID=@", inSNPs, " ", inFile, " -o ", outFile, " -O z --threads 1")
    print("Filtering BCF")
    print(baseCommand)
    system(command=baseCommand)
    #baseIndex(outFile)
    indexCommand <- paste0("bcftools index -f ", outFile)
    print("Filtered Data")
    Sys.time()
    system(command=indexCommand)
    print(outFile)
  }
  return(outFile)
}


filterMergedReg <- function(inFile, inName, inRegion){
  outFile <- gsub("\\.vcf\\.gz", paste0("_",inName, "_reg.vcf.gz"), inFile)
  baseCommand <- paste("bcftools view -R", inRegion, inFile, "-o", outFile, "-O z --threads 1")
  system(command=baseCommand)
  indexCommand <- paste0("bcftools index -f ", outFile)
  system(command=indexCommand)
  return(outFile)
}
filterMergedPos <- function(inFile, inName, inSNPs){
  ## TODO Implement using python script for the time being
  outFile <- gsub(paste0("(_", inName, "_filt", ")?\\.vcf\\.gz"), paste0("_",inName, ".bcf"), inFile)
  ##TODO Implement if file exists check
    print("PGS Catalog Start")
    baseCommand <- paste0("python2 vcf_pgs_catalog.py --score ", inSNPs, " --header --out ", outFile, " --vcf ", inFile)
    print(baseCommand)
    #baseCommand <- paste0("python3 cyvf_pgs_catalog.py --score ", inSNPs, " --header --out ", outFile, " --vcf ", inFile)
    print(baseCommand)
    system(command=baseCommand)
    print("PGS Catalog End")
    Sys.time()
    indexCommand <- paste0("bcftools index -f ", outFile)
    system(command=indexCommand)
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
