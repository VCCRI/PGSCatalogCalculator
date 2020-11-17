require(data.table)

getPGSEFO <- function(inFile){
  tryCatch({
     require(readxl)
  }, warning = function(w){
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.csiro.au/"
    install.packages("readxl")
    require(readxl)
  }, error = function(e){
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.csiro.au/"
    install.packages("readxl")
    require(readxl)
  })
  inData <- read_xlsx(inFile, sheet=7)
  inData <- setDT(inData)
  require(assertthat)
  needCols <- c("Ontology Trait ID", "Ontology Trait Label",	"Ontology Trait Description", "Ontology URL")
  assert_that(ncol(inData) == 4)
  assert_that(all.equal(colnames(inData), needCols))
  return(inData)
}

getEUBuild <- function(inFile){
  tryCatch({
     require(readxl)
  }, warning = function(w){
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.csiro.au/"
    install.packages("readxl")
    require(readxl)
  }, error = function(e){
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.csiro.au/"
    install.packages("readxl")
    require(readxl)
  })
  inData <- read_xlsx(inFile, sheet=5)
  inData <- setDT(inData)
  require(assertthat)
  #needCols <- c("Ontology Trait ID", "Ontology Trait Label",	"Ontology Trait Description", "Ontology URL")
  assert_that(ncol(inData) == 17)
  inIds <- inData[grepl("European", `Broad Ancestry Category`),`Polygenic Score (PGS) ID`]
  return(inIds)
}

getPGSMeta <- function(inFile){
  tryCatch({
     require(readxl)
  }, warning = function(w){
    r = getOption("repos")
    r["CRAN"] = "https://cran.csiro.au/"
    install.packages("readxl")
    require(readxl)
  }, error = function(e){
    r = getOption("repos")
    r["CRAN"] = "https://cran.csiro.au/"
    install.packages("readxl")
    require(readxl)
  })
  inData <- read_xlsx(inFile, sheet=2)
  inData <- setDT(inData)
  require(assertthat)
  print("Test")
  #Mapped Traits is a comma seperated list
  needCols <- c(
  "Polygenic Score (PGS) ID"
  ,"PGS Name"
  ,"Reported Trait"
  ,"Mapped Trait(s) (EFO label)"
  ,"Mapped Trait(s) (EFO ID)"
  ,"PGS Development Method"
  ,"PGS Development Details/Relevant Parameters"
  ,"Original Genome Build"
  ,"Number of Variants"
  ,"Number of Interaction Terms"
  ,"PGS Publication (PGP) ID"
  ,"Publication (PMID)"
  ,"Publication (doi)"
  ,"Score and results match the original publication"
  )
  assert_that(ncol(inData) == 14)
  assert_that(all.equal(colnames(inData), needCols))
  ##Generate FTP URL
  inData[,ftp_link := paste0("http://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/", `Polygenic Score (PGS) ID`, "/ScoringFiles/", `Polygenic Score (PGS) ID`, ".txt.gz")]
  return(inData)
}

getRTmpFiles <- function(){
  URL <- c("ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata.xlsx")
  #x <- getURL(URL)
  #return(rawConnection(x))
  tempFile <- tempfile()
  inTemp <- download.file(URL, destfile=tempFile, mode='wb', method="libcurl")
  print("Downloaded File")
  return(tempFile)
}

getFTPFiles <- function(inTerm, inData){
  needCols <- c(
  "Polygenic Score (PGS) ID"
  ,"ftp_link"
  ,"Number of Variants"
  )
  return(inData[grepl(tolower(inTerm), tolower(`Mapped Trait(s) (EFO label)`)), ..needCols])
}

getPGId <- function(inData){
  require(assertthat)
  assert_that(all.equal(unique(inData$`Polygenic Score (PGS) ID`), inData$`Polygenic Score (PGS) ID`))
  return(unique(inData$`Polygenic Score (PGS) ID`))
}

getFileSize <- function(x){
  require(RCurl)
  res <- url.exists(x, .header=TRUE)
  fileSize <- as.numeric(res['Content-Length'])
  require(assertthat)
  assert_that(!is.na(as.numeric(fileSize)))
  return(fileSize)
}


isSmall <- function(inData){
  require(assertthat)
  assert_that(length(inData) == 3)
  needCols <- c(
  "Polygenic Score (PGS) ID"
  ,"ftp_link"
  ,"Number of Variants"
  )
  assert_that(all.equal(names(inData), needCols))
  inData[,size_file := lapply(ftp_link, getFileSize)]
  inData[,is_small := `Number of Variants` <= 100 & as.numeric(size_file) <= 1024*5,]
  if(nrow(inData[is_small == TRUE,]) == 0) print("No Data Matches")
  return(inData)
}
getPGSFiles <- function(inData){
  require(assertthat)
  require(httr)
  print(inData)
  print(class(inData))
  #assert_that("data.table" %in% class(inData))
  #assert_that(is.recursive(inData) == FALSE)
  inFileURL <- getElement(inData, "ftp_link")
  pgsId <- getElement(inData, "Polygenic Score (PGS) ID")
  response <- httr::HEAD(inFileURL)
  assert_that(length(inFileURL) == 1)
  tempFile <- tempfile()
  print("Start File")
  inTemp <- download.file(inFileURL, destfile=tempFile, mode='wb', quiet=T)
  print("Downloading File")
  return(list(file=tempFile, pgsId=pgsId))
}

setClassPGSGRS <- function(){
  setClassUnion("nullOrDatatable", c("NULL", "data.table"))
  setClassUnion("nullOrCharacter", c("NULL", "character"))
  setClass("pgsInput", slots=(list(pgsInput="nullOrDatatable", pgsType="nullOrCharacter", pgsId="character")))
}
pgsGRS <- function(inFile){
  require(data.table)
  require(assertthat)
  ## TODO Recognise when there multiple columns with the same column name
  ## Check if variant is A,C,T,G handle edge case where variant is 
  ## Check DR3/DR4-DQ8 PGS000021
  inputFile <- getElement(inFile, "file")
  pgsId <- getElement(inFile, "pgsId")
  inData <- fread(cmd=paste("zcat -f", inputFile) , stringsAsFactors=F)
  needCols <- c("effect_allele","effect_weight")
  if(any(!(needCols %in% colnames(inData)))){
    return(new("pgsInput", pgsInput=NULL, pgsType=NULL, pgsId = pgsId))
  }
  assert_that(all(needCols %in% colnames(inData)))
  if(c("rsID") %in% colnames(inData)){
    pgsType <- "rsID"
  } else {
    pgsType <- "chr_pos"
    inData$id <- paste0(inData$chr, inData$chr_pos, inData$effect_allele)
  }
  return(new("pgsInput", pgsInput=inData, pgsType=pgsType, pgsId = pgsId))
}

newPGSGRS <- function(inRow){
  inFile <- getPGSFiles(inRow)
  inClass <- pgsGRS(inFile)
}

getFiles <- function(){
  inTemp <- getRTmpFiles()
  metaFile <- getPGSMeta(inTemp)
  #ftpLinks <- getFTPFiles(args[3], metaFile)
  ftpLinks <- getFTPFiles('*', metaFile)
  ## TODO Faster method to check smaller file
  #checkOne <- ftpLinks[2,]
  #checkOne <- head(isSmall(ftpLinks)[is_small == FALSE,],2:2)
  ftpLinks <- ftpLinks[grepl("PGS000006",`Polygenic Score (PGS) ID`),]
  #ftpLinks <- ftpLinks[grepl("PGS000039",`Polygenic Score (PGS) ID`),]
  #ftpLinks <- ftpLinks[,]
  require(doParallel)
  require(foreach)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  ### TODO Remove apply to increase performance
  #inFile <- getPGSFiles(ftpLinks[1,])
  #setClassPGSGRS()
  #inClass <- pgsGRS(inFile)
  #retFiles <- inClass
  retFiles <- foreach(d=iter(ftpLinks, by='row'), .export = c("getPGSFiles", "pgsGRS", "setClassPGSGRS"), .packages=c("data.table", "httr", "assertthat")) %dopar% {
    inFile <- getPGSFiles(d)
    setClassPGSGRS()
    inClass <- pgsGRS(inFile)
  }
  ftpFiles <-  do.call("rbind", apply(ftpLinks, 1, getPGSFiles))
  stopCluster(cl)
  print('Finished Download')
  return(retFiles)
}

getEUFiles <- function(){
  inTemp <- getRTmpFiles()
  metaFile <- getPGSMeta(inTemp)
  euBuilds <- getEUBuild(inTemp)
  euMeta <- metaFile[(!grepl("PGS000006",`Polygenic Score (PGS) ID`)) & (`Polygenic Score (PGS) ID` %in% euBuilds) & (`Original Genome Build` %in% c("GRCh37", "hg19")),]
  ftpLinks <- getFTPFiles('*', euMeta)
  require(doParallel)
  require(foreach)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  retFiles <- foreach(d=iter(ftpLinks, by='row'), .export = c("getPGSFiles", "pgsGRS", "setClassPGSGRS"), .packages=c("data.table", "httr", "assertthat")) %dopar% {
    inFile <- getPGSFiles(d)
    setClassPGSGRS()
    inClass <- pgsGRS(inFile)
  }
  ftpFiles <-  do.call("rbind", apply(ftpLinks, 1, getPGSFiles))
  stopCluster(cl)
  print('Finished Download')
  return(retFiles)
}

checkFiles <- function(){
  setClassPGSGRS()
  pgsObj <- pgsGRS(list(file='/home/114/sk3015/Analysis/PGS_PRSice/mockCheckSN.tsv', pgsId='testing'))
}


getRSIds <- function(inPgsInput){
  require(assertthat)
  assert_that(inPgsInput@pgsType == "rsID")
  inGRSInput <- inPgsInput@pgsInput
  tempFile <- tempfile()
  needCols <- c("rsID")
  fwrite(inGRSInput[,..needCols], tempFile,row.names=F, col.names=F, sep="\t")
  #return(unique((inPgsInput@pgsInput)$rsID))
  return(tempFile)
}

helperDwQ <- function(inData){
  download.file(inData$links, destfile=inData$fileName, method="libcurl", quiet=T)
  if(grepl("(pgen\\.zst)|(pvar\\.zst)", inData$fileName, perl=T)){
    decomCom <- paste(plink2, "--zst-decompress", inData$fileName, ">", gsub("\\.zst$", "", inData$fileName))
    system(command=decomCom)
  }
}

get1000Genomes <- function(inFile, plink2="/g/data/jb96/sk3015/software/plink2"){
  ## Cannot use PLINK 1.9 because corresponding VCF GT are incorrect
  ## Done by manually checking GT against base 1000G ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  ## TODO Multithread downloads
  require(doParallel)
  require(data.table)
  inZst <- gsub("\\.vcf\\.gz", ".pgen.zst", inFile)
  inPvar <- gsub("\\.vcf\\.gz", ".pvar.zst", inFile)
  inPSam <- gsub("\\.vcf\\.gz", ".psam", inFile)
  links <- c("https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1",
             "https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1",
             "https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1")
  allFiles <- data.table(fileName=c(inZst, inPvar, inPSam), links=links)
  inDir <- dirname(inFile)
  if(!(dir.exists(inDir))){
    dir.create(inDir)
  }
  cl <- makeCluster(3)
  registerDoParallel(cl)
  totalDowCodes <- foreach(b=iter(allFiles, by='row'), .export = "helperDwQ") %dopar% helperDwQ(b)
  stopCluster(cl)
  baseFile <- gsub("\\.pgen", "", gsub("\\.zst", "", inZst))
  system(command=paste(plink2, "--export vcf --keep-if SuperPop==EUR --out /g/data/jb96/sk3015/1000GTest/sampletesting --pfile", baseFile))
  return(paste0(baseFile,".vcf.gz"))
}

  #inPgen <- gsub("\\.vcf\\.gz", ".pgen.zst", inFile)
  #dwPgen <- download.file("https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1", destfile=inPgen, method="wget")
  #iPVar <- gsub("\\.vcf\\.gz", ".pvar.zst", inFile)
  #dwPVar<- download.file("https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1", destfile=iPVar, method="wget")
  


getChrPos <- function(inPgsInput){
  require(assertthat)
  assert_that(inPgsInput@pgsType == "chr_pos")
  inGRSInput <- inPgsInput@pgsInput
  inGRSInput$chr_position_end <- inGRSInput$chr_position
  inGRSInput$chr_position <- as.numeric(inGRSInput$chr_position-1)
  tempFile <- tempfile()
  tempFile <- gsub("$", ".bed", tempFile)
  needCols <- c("chr_name", "chr_position", "chr_position_end")
  fwrite(inGRSInput[,..needCols], tempFile,row.names=F, col.names=F, sep="\t")
  #return(unique((inPgsInput@pgsInput)$rsID))
  return(tempFile)
}

getPGSId <- function(inPgsInput){
  #inGRSInput <- inPgsInput@pgsInput
  return(inPgsInput@pgsId)
}

getPGSType <- function(inPgsInput){
  #inGRSInput <- inPgsInput@pgsInput
  return(inPgsInput@pgsType)
}

getGRSInputChrPos <- function(inPgsInput){
  require(assertthat)
  assert_that(inPgsInput@pgsType == "chr_pos")
  inGRSInput <- inPgsInput@pgsInput
  tempFile <- tempfile()
  inGRSInput$id <- paste(inGRSInput$chr_name, inGRSInput$chr_position, inGRSInput$effect_allele, sep=':')
  needCols <- c("id", "effect_allele", "effect_weight", "chr_name", "chr_position")
  fwrite(inGRSInput[,..needCols], tempFile,row.names=F, sep=" ")
  return(tempFile)
}

getGRSInputRsID <- function(inPgsInput){
  require(assertthat)
  assert_that(inPgsInput@pgsType == "rsID")
  inGRSInput <- inPgsInput@pgsInput
  inGRSInput$rsID <- paste(inGRSInput$rsID, inGRSInput$effect_allele, sep=":")
  tempFile <- tempfile()
  needCols <- c("rsID", "effect_allele", "effect_weight")
  browser()
  fwrite(inGRSInput[,..needCols], tempFile,row.names=F, sep=" ")
  return(tempFile)
}

getNormControl <- function(){
  normControl <- baseNorm(inVCF=inControl,
       inRef="/g/data/jb96/References_and_Databases/hs37d5.fa/hs37d5x.fa",
       outFile=gsub("\\.vcf\\.gz", paste0("_norm", ".vcf.gz"), inContol))
  return(normControl)
}


 
runGRSCalc <- function(inObjec, inFile){
  grsFile <- getGRSInputRsID(inObjec)
  rsIDFile <- getRSIds(inObjec)
  pgsID <- getPGSId(inObjec)
  print(inFile)
  filterRsid <- filterMerged(inFile=inFile, inName=pgsID, inSNPs=rsIDFile)
  filterMerged <- filterMergedPos(inFile=filterRsid, inName=pgsID, inSNPs=grsFile)
  source('makePlink.R')
  plinkFile <- getMakePlink(filterMerged)
  source('plinkPRS.R')
  outScore <- plinkGRS(inFile= plinkFile, inGRS=grsFile)
  system(command=paste("rm", filterMerged))
  system(command=paste("rm", filterRsid))
  return(outScore)
}

runGRSCalcChrPos <- function(inObjec, inFile){
  grsFile <- getGRSInputChrPos(inObjec)
  regionId <- getChrPos(inObjec)
  pgsID <- getPGSId(inObjec)
  filterMerged <- filterMergedPos(inFile=inFile, inName=pgsID, inSNPs=grsFile)
  source('makePlink.R')
  plinkFile <- getMakePlink(filterMerged)
  source('plinkPRS.R')
  outScore <- plinkGRS(inFile= plinkFile, inGRS=grsFile)
  system(command=paste("rm", filterMerged))
  return(outScore)
}

getScore <- function(inFile){
  filtFile <- gsub(paste0(inDir, "/"), "", inFile)
  filtFile <- gsub("_score\\.profile", "", filtFile)
  filtSample <- gsub(".+_DBN_", "", filtFile)
  inData <- fread(inFile, stringsAsFactors=FALSE)
  inData[,PGS_RECORD_ID := filtSample]
  return(inData)
}


getAllScores <- function(inDir){
  inFiles <- list.files(path=inDir, ".*_DBN_PGS\\d+_score.profile", full.names=TRUE)
  totalSamples <- rbindlist(lapply(inFiles, getFile))

  aggSamples <- totalSamples[,sum(as.numeric(SCORESUM)), by=.(IID, PGS_RECORD_ID)]

  setnames(aggSamples, "V1", "PRS")
  return(aggSamples)
}

packages <- function(inPackage){
  tryCatch({
     require(inPackage, character.only=T)
  }, warning = function(w){
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.csiro.au/"
    options(repos = r)
    install.packages(inPackage)
    require(inPackage, character.only=T)
  }, error = function(e){
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.csiro.au/"
    options(repos = r)
    install.packages(inPackage)
    require(inPackage, character.only=T)
  })
}
main <- function(){
  lapply(list("assertthat","cowplot","data.table","doParallel","foreach","ggplot2","glm2","httr","optparse","promises","questionr","RCurl", "readxl","Rsamtools", "future"), packages)
  source('normScript.R')
  source('mergeVCF.R')
  require(optparse)
  library("optparse")

  option_list <- list(
    make_option(c("-f", "--file"), type="character", default=NULL,
                help="Input VCF filename", metavar="character"),
      make_option(c("-o", "--out"), type="character", default="sample_norm/",
                help="output directory", metavar="character"),
      make_option(c("-r", "--ref"), type="character", default="index.fa",
                help="reference sequence file", metavar="character"),
      make_option(c("-p", "--pgs-id"), type="character", default="PGS000073",
                help="Single PGS ID to score file", metavar="character"),
      make_option(c("-P", "--pgs-id-file"), type="character", default="sample.csv",
                help="PGS ID that has list of Files", metavar="character")
  )

  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
  if (is.null(opt$file)){
    print_help(opt_parser)
    stop("At least one argument must be supplied", call.=FALSE)
  }
  require(future)
  require(promises)
  plan("multiprocess", workers=6)
  inspecFile <- readRDS(url("https://pgscatalogscraper.s3-us-west-2.amazonaws.com/eu_metadata.RDS", "rb"))

  normFile <- future({
    source('normScript.R')
    baseNorm(inVCF=opt$file,
           inRef=opt$ref,
           outFile=gsub("\\.[a-zA-Z']+(\\.gz)?$","_norm", opt$file))
  }, packages=c("Rsamtools"))

  #while(!(resolved(inspecFile) & resolved(normFile))) {}
  while(!(resolved(normFile))) {}

  if(resolved(normFile)){
     ## Need to fix ID
     normFile <- value(normFile)
     inputFile <- inspecFile
     scoreFiles <- lapply(inputFile, function(inObjec){
      scoreFiles <- data.table() 
      if(inObjec@`pgsId` == opt$`pgs-id` | is.na(opt$`pgs-id`)){
        require(doParallel)
        require(foreach)
        cl <- makeCluster(10)
        registerDoParallel(cl)
        source('makePlink.R')
        source('plinkPRS.R')
        scoreFiles <- foreach(d=iter(unlist(normFile)), .export = c("getPGSType","runGRSCalc", "runGRSCalcChrPos", "getChrPos", "getGRSInputChrPos", "getGRSInputRsID", "getGRSInputRsID", "getRSIds", "getPGSId", "filterMerged", "filterMergedPos", "filterMergedReg", "getMakePlink", "plinkGRS", "baseIndex"), .packages=c("data.table"), .combine=rbind) %dopar% {
          if(getPGSType(inObjec) == "rsID"){
           scoreFile <- runGRSCalc(inObjec, d)
           return(data.table(inFile=scoreFile,inRecord=getPGSId(inObjec)))
          } else {
           scoreFile <- runGRSCalcChrPos(inObjec, d)
           return(data.table(inFile=scoreFile,inRecord=getPGSId(inObjec)))
          }
      }
        stopCluster(cl)
      }
      return(scoreFiles)
    })
 }
 source("calcScores.R")
 source("getCharts.R")
 setPlots(rbindlist(scoreFiles))
 
}
