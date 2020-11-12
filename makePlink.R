getMakePlink <- function(inVCF){
  inPlink <- "/g/data/jb96/software/plink_1.9_linux_x86_64_20181202/plink"
  ##ToDO prettifry paste
  if(grepl(".vcf.gz", inVCF)){
    inType <- paste("--vcf", inVCF, "--vcf-half-call h")
    outFile <- gsub("\\.vcf\\.gz", "_plink", inVCF)
  } else {
    inType <- paste("--bcf", inVCF)
    outFile <- gsub("\\.bcf", "_plink", inVCF)
 }
  #baseCommand <- paste(inPlink, inType, " --allow-extra-chr --maf 0.05 --mind 0.1 --geno 0.1 --hwe 1e-6 --vcf-filter --make-bed --chr 1-22 XY --memory 4096 --out", outFile)
  baseCommand <- paste(inPlink, inType, " --allow-extra-chr --mind 0.1 --geno 0.05 --hwe 1e-6 --vcf-filter --make-bed --chr 1-22 XY --memory 4096 --out", outFile)
  system(command=baseCommand)
  return(outFile)
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
