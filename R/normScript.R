baseNorm <- function(inVCF, inRef, outFile, inYaml=NULL, inCL){
  dir.create(outFile, showWarnings=FALSE)
  chrom <- c(as.character(1:22), "X", "Y")
  #chrom <- c(as.character(2:22), "X", "Y")
  #chrom <- c(as.character(1))
  inDBN <- foreach::foreach(inChrom=chrom, .export = c("baseIndex")) %dopar% {
    filtSNP <- gsub("\\.[a-zA-Z']+(\\.gz)?$", paste0("_",inChrom, "_snp", ".vcf.gz"), inVCF)
    #filtSNP <- gsub("\\.bcf", paste0("_",inChrom, "_snp", ".vcf.gz"), inVCF)
    #filterSNP <- system(command=paste('bcftools view', inVCF, '-i \'TYPE="snp"\' -r', inChrom, '-O z -o', filtSNP, '--threads 2'))
    bcftools <- if(is.null(inYaml)) Sys.getenv("bcftools") else yaml::read_yaml(inYaml)$bcftools
    filterSNP <- system2(command=bcftools, args=c('view', inVCF, '-r', inChrom, '-O', 'z','-o', filtSNP, '--threads', '2'))
    #filterSNP <- system(command=paste('bcftools view', inVCF, '-r', inChrom, '-O z -o', filtSNP, '--threads 2'))
    #baseCommand <- "/g/data3/jb96/software/vt/vt"
    baseCommand <- if(is.null(inYaml)) Sys.getenv("vt") else yaml::read_yaml(inYaml)$vt 
    for (x in c("Decom", "DecomBlock", "DBN")){
      assign(paste0("in", x), paste0(outFile, "_", inChrom, "_", x, ".vcf.gz")) 
    }
    system2(command=baseCommand, args=c("decompose", "-o", inDecom, "-s", filtSNP), stdout=FALSE)
    baseIndex(inDecom, inYaml)
    system2(command=baseCommand, args=c("decompose_blocksub", "-o", inDecomBlock, "-a", inDecom), stdout=FALSE)
    baseIndex(inDecomBlock, inYaml)
    print("Indexing Files 1")
    system2(command=baseCommand, args=c("normalize", "-r", inRef, "-o", inDBN, inDecomBlock, "-q"), stdout=FALSE)
    print("Indexing Files 2")
    baseIndex(inDBN, inYaml)
    
    return(inDBN)
  }
  return(inDBN)
}


baseIndex <- function(inVCF, inYaml=NULL){
  print("Indexing File")

  bcftools <- if(is.null(inYaml)) Sys.getenv("bcftools") else yaml::read_yaml(inYaml)$bcftools
  system2(command=bcftools, args=c('index', "-f", inVCF), stdout=FALSE)
  #indexTabix(file=inVCF, format="vcf")
}
# Uncomment to add in main
#baseNorm(inVCF="/g/data/jb96/sk3015/1000GDumpTest/one_chrom.vcf.gz",
         #inRef="/g/data/jb96/References_and_Databases/hs37d5.fa/hs37d5x.fa",
         #outFile="/g/data/jb96/sk3015/1000GDumpDecomp/europe_onechrom")

# TODO Immplement interaction between norm script and pgs grabber

