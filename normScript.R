

baseNorm <- function(inVCF, inRef, outFile){
  dir.create(outFile, showWarnings=FALSE)
  chrom <- c(as.character(1:22), "X", "Y")
  #chrom <- c(as.character(2:22), "X", "Y")
  #chrom <- c(as.character(1))
  require(doParallel)
  require(foreach)
  cl <- makeCluster(3)
  registerDoParallel(cl)
  inDBN <- foreach(inChrom=chrom, .export = c("baseIndex"), .packages="Rsamtools") %dopar% {
    filtSNP <- gsub("\\.[a-zA-Z']+(\\.gz)?$", paste0("_",inChrom, "_snp", ".vcf.gz"), inVCF)
    #filtSNP <- gsub("\\.bcf", paste0("_",inChrom, "_snp", ".vcf.gz"), inVCF)
    #filterSNP <- system(command=paste('bcftools view', inVCF, '-i \'TYPE="snp"\' -r', inChrom, '-O z -o', filtSNP, '--threads 2'))
    filterSNP <- system(command=paste('bcftools view', inVCF, '-r', inChrom, '-O z -o', filtSNP, '--threads 2'))
    #filterSNP <- system(command=paste('bcftools view', inVCF, '-r', inChrom, '-O z -o', filtSNP, '--threads 2'))
    baseCommand <- "/g/data3/jb96/software/vt/vt"
    for (x in c("Decom", "DecomBlock", "DBN")){
      assign(paste0("in", x), paste0(outFile, "_", inChrom, "_", x, ".vcf.gz")) 
    }
    system(command=paste(baseCommand, "decompose -o", inDecom, "-s", filtSNP), intern=TRUE)
    baseIndex(inDecom)
    system(command=paste(baseCommand, "decompose_blocksub -o", inDecomBlock, "-a", inDecom), intern=TRUE)
    baseIndex(inDecomBlock)
    print("Indexing Files 1")
    system(command=paste(baseCommand, "normalize -r", inRef, "-o", inDBN, inDecomBlock), intern=TRUE)
    print("Indexing Files 2")
    baseIndex(inDBN)
    
    return(inDBN)
  }
  stopCluster(cl)
  return(inDBN)
}


baseIndex <- function(inVCF){
  tryCatch({
    require(Rsamtools)
  }, warning = function(w){
    r = getOption("repos")
    r["CRAN"] = "https://cran.csiro.au/"
    options(repos = r)
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("Rsamtools")
    require(Rsamtools)
  }, error = function(e){
    r = getOption("repos")
    r["CRAN"] = "https://cran.csiro.au/"
    options(repos = r)
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("Rsamtools")
    require(Rsamtools)
  })
  print("Indexing File")
  system(command=paste('bcftools index', inVCF))
  #indexTabix(file=inVCF, format="vcf")
}
# Uncomment to add in main
#baseNorm(inVCF="/g/data/jb96/sk3015/1000GDumpTest/one_chrom.vcf.gz",
         #inRef="/g/data/jb96/References_and_Databases/hs37d5.fa/hs37d5x.fa",
         #outFile="/g/data/jb96/sk3015/1000GDumpDecomp/europe_onechrom")

# TODO Immplement interaction between norm script and pgs grabber

