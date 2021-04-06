getPGSEFO <- function(inFile){
  inData <- readxl::read_xlsx(inFile, sheet=7)
  inData <- data.table::setDT(inData)
  needCols <- c("Ontology Trait ID", "Ontology Trait Label",	"Ontology Trait Description", "Ontology URL")
  assertthat::assert_that(ncol(inData) == 4)
  assertthat:assert_that(all.equal(colnames(inData), needCols))
  return(inData)
}

getEUBuild <- function(inFile){
  inData <- readxl::read_xlsx(inFile, sheet=5)
  inData <- data.table::setDT(inData)
  #needCols <- c("Ontology Trait ID", "Ontology Trait Label",	"Ontology Trait Description", "Ontology URL")
  assertthat::assert_that(ncol(inData) == 17)
  inIds <- inData[grepl("European", `Broad Ancestry Category`),`Polygenic Score (PGS) ID`]
  return(inIds)
}

getPGSMeta <- function(inFile){
  inData <- readxl::read_xlsx(inFile, sheet=2)
  inData <- data.table::setDT(inData)
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
  ,"FTP link"
  ,"License/Terms of Use"
  )
  assertthat::assert_that(ncol(inData) == 16)
  assertthat::assert_that(all.equal(colnames(inData), needCols))
  ##Generate FTP URL
  inData[,ftp_link := `FTP link`]
  return(inData)
}

getRTmpFiles <- function(){
  URL <- c("ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata.xlsx")
  #x <- getURL(URL)
  #return(rawConnection(x))
  tempFile <- tempfile()
  inTemp <- download.file(URL, destfile=tempFile, mode='wb', method="wget")
  print("Downloaded File")
  return(tempFile)
}

getFTPFiles <- function(inTerm, inData){
  needCols <- c(
  "Polygenic Score (PGS) ID"
  ,"ftp_link"
  ,"Number of Variants"
  )
  return(inData[grepl(tolower(inTerm), tolower(`Mapped Trait(s) (EFO label)`)), needCols, with=FALSE])
}

getPGId <- function(inData){
  assertthat::assert_that(all.equal(unique(inData$`Polygenic Score (PGS) ID`), inData$`Polygenic Score (PGS) ID`))
  return(unique(inData$`Polygenic Score (PGS) ID`))
}

getFileSize <- function(x){
  res <- RCurl::url.exists(x, .header=TRUE)
  fileSize <- as.numeric(res['Content-Length'])
  assertthat::assert_that(!is.na(as.numeric(fileSize)))
  return(fileSize)
}


isSmall <- function(inData){
  assertthat::assert_that(length(inData) == 3)
  needCols <- c(
  "Polygenic Score (PGS) ID"
  ,"ftp_link"
  ,"Number of Variants"
  )
  assertthat::assert_that(all.equal(names(inData), needCols))
  inData[,size_file := lapply(ftp_link, getFileSize)]
  inData[,is_small := `Number of Variants` <= 100 & as.numeric(size_file) <= 1024*5,]
  if(nrow(inData[is_small == TRUE,]) == 0) print("No Data Matches")
  return(inData)
}
getPGSFiles <- function(inData){
  #assertthat::assert_that("data.table" %in% class(inData))
  #assertthat::assert_that(is.recursive(inData) == FALSE)
  inFileURL <- getElement(inData, "ftp_link")
  pgsId <- NULL
  pgsId <- getElement(inData, "Polygenic Score (PGS) ID")
  response <- httr::HEAD(inFileURL)
  assertthat::assert_that(length(inFileURL) == 1)
  tempFile <- NULL
  tempFile <- tempfile()
  print("Start File")
  inTemp <- download.file(inFileURL, destfile=tempFile, mode='wb', quiet=T,method="wget")
  print("Downloading File")
  return(list(file=tempFile, pgsId=pgsId))
}

setClassPGSGRS <- function(){
  setClassUnion("nullOrDatatable", c("NULL", "data.table"))
  setClassUnion("nullOrCharacter", c("NULL", "character"))
  setClass("pgsInput", slots=(list(pgsInput="nullOrDatatable", pgsType="nullOrCharacter", pgsId="character")))
}
pgsGRS <- function(inFile){
  ## TODO Recognise when there multiple columns with the same column name
  ## Check if variant is A,C,T,G handle edge case where variant is 
  ## Check DR3/DR4-DQ8 PGS000021
  inputFile <- getElement(inFile, "file")
  pgsId <- getElement(inFile, "pgsId")
  inData <- data.table::fread(cmd=paste("zcat -f", inputFile) , stringsAsFactors=F)
  needCols <- c("effect_allele","effect_weight")
  if(any(!(needCols %in% colnames(inData)))){
    return(new("pgsInput", pgsInput=NULL, pgsType=NULL, pgsId = pgsId))
  }
  assertthat::assert_that(all(needCols %in% colnames(inData)))
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
  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl)
  ### TODO Remove apply to increase performance
  #inFile <- getPGSFiles(ftpLinks[1,])
  #setClassPGSGRS()
  #inClass <- pgsGRS(inFile)
  #retFiles <- inClass
  retFiles <- foreeach::foreach(d=iterators::iter(ftpLinks, by='row'), .export = c("getPGSFiles", "pgsGRS", "setClassPGSGRS"), .packages=c("data.table", "httr", "assertthat")) %dopar% {
    inFile <- getPGSFiles(d)
    setClassPGSGRS()
    inClass <- pgsGRS(inFile)
  }
  ftpFiles <-  do.call("rbind", apply(ftpLinks, 1, getPGSFiles))
  parallel::stopCluster(cl)
  print('Finished Download')
  return(retFiles)
}

getEUFiles <- function(){
  inTemp <- getRTmpFiles()
  metaFile <- getPGSMeta(inTemp)
  euBuilds <- getEUBuild(inTemp)
  euMeta <- metaFile[(!grepl("PGS000006",`Polygenic Score (PGS) ID`)) & (`Polygenic Score (PGS) ID` %in% euBuilds) & (`Original Genome Build` %in% c("GRCh37", "hg19")),]
  ftpLinks <- getFTPFiles('*', euMeta)
  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl)
  retFiles <- foreach::foreach(d=iterators::iter(ftpLinks, by='row'), .export = c("getPGSFiles", "pgsGRS", "setClassPGSGRS"), .packages=c("data.table", "httr", "assertthat")) %dopar% {
  #retFiles <- apply(ftpLinks, 1, function(d){
    inFile <- getPGSFiles(d)
    setClassPGSGRS()
    inClass <- pgsGRS(inFile)
  #})
  }
  ftpFiles <-  do.call("rbind", apply(ftpLinks, 1, getPGSFiles))
  parallel::stopCluster(cl)
  print('Finished Download')
  return(retFiles)
}

checkFiles <- function(){
  setClassPGSGRS()
  pgsObj <- pgsGRS(list(file='/home/114/sk3015/Analysis/PGS_PRSice/mockCheckSN.tsv', pgsId='testing'))
}


getRSIds <- function(inPgsInput){
  assertthat::assert_that(inPgsInput@pgsType == "rsID")
  inGRSInput <- inPgsInput@pgsInput
  tempFile <- tempfile()
  needCols <- c("rsID")
  data.table::fwrite(inGRSInput[,needCols, with=FALSE], tempFile,row.names=F, col.names=F, sep="\t")
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
  cl <- parallel::makeCluster(3)
  doParallel::registerDoParallel(cl)
  totalDowCodes <- foreach::foreach(b=iterators::iter(allFiles, by='row'), .export = "helperDwQ") %dopar% helperDwQ(b)
  parallel::stopCluster(cl)
  baseFile <- gsub("\\.pgen", "", gsub("\\.zst", "", inZst))
  system(command=paste(plink2, "--export vcf --keep-if SuperPop==EUR --out /g/data/jb96/sk3015/1000GTest/sampletesting --pfile", baseFile))
  return(paste0(baseFile,".vcf.gz"))
}

  #inPgen <- gsub("\\.vcf\\.gz", ".pgen.zst", inFile)
  #dwPgen <- download.file("https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1", destfile=inPgen, method="wget")
  #iPVar <- gsub("\\.vcf\\.gz", ".pvar.zst", inFile)
  #dwPVar<- download.file("https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1", destfile=iPVar, method="wget")
  


getChrPos <- function(inPgsInput){
  assertthat::assert_that(inPgsInput@pgsType == "chr_pos")
  inGRSInput <- inPgsInput@pgsInput
  inGRSInput$chr_position_end <- inGRSInput$chr_position
  inGRSInput$chr_position <- as.numeric(inGRSInput$chr_position-1)
  tempFile <- tempfile()
  tempFile <- gsub("$", ".bed", tempFile)
  needCols <- c("chr_name", "chr_position", "chr_position_end")
  data.table::fwrite(inGRSInput[,needCols, with=FALSE], tempFile,row.names=F, col.names=F, sep="\t")
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
  assertthat::assert_that(inPgsInput@pgsType == "chr_pos")
  inGRSInput <- inPgsInput@pgsInput
  tempFile <- tempfile()
  inGRSInput$id <- paste(inGRSInput$chr_name, inGRSInput$chr_position, inGRSInput$effect_allele, sep=':')
  needCols <- c("id", "effect_allele", "effect_weight", "chr_name", "chr_position")
  data.table::fwrite(inGRSInput[,needCols,with=FALSE], tempFile,row.names=F, sep=" ")
  return(tempFile)
}

getGRSInputRsID <- function(inPgsInput){
  assertthat::assert_that(inPgsInput@pgsType == "rsID")
  inGRSInput <- inPgsInput@pgsInput
  inGRSInput$rsID <- paste(inGRSInput$rsID, inGRSInput$effect_allele, sep=":")
  tempFile <- tempfile()
  needCols <- c("rsID", "effect_allele", "effect_weight")
  data.table::fwrite(inGRSInput[,needCols, with=FALSE], tempFile,row.names=F, sep=" ")
  return(tempFile)
}

getNormControl <- function(){
  normControl <- baseNorm(inVCF=inControl,
       inRef="/g/data/jb96/References_and_Databases/hs37d5.fa/hs37d5x.fa",
       outFile=gsub("\\.vcf\\.gz", paste0("_norm", ".vcf.gz"), inContol))
  return(normControl)
}

runGRSCalc <- function(inObjec, inDis,inFam,inCont=NULL){
  grsFile <- getGRSInputRsID(inObjec)
  rsIDFile <- getRSIds(inObjec)
  pgsID <- getPGSId(inObjec)
  filterRsidDis <- filterMerged(inFile=inDis, inName=pgsID, inSNPs=rsIDFile)
  filterMergedDis <- filterMergedPos(inFile=filterRsidDis, inName=pgsID, inSNPs=grsFile)
  plinkFileDis <- getMakePlink(inVCF=filterMergedDis)
  if(!(is.null(inCont))) {
    filterRsidCont <- filterMerged(inFile=inCont, inName=pgsID, inSNPs=rsIDFile)
    filterMergedCont <- filterMergedPos(inFile=filterRsidCont, inName=pgsID, inSNPs=grsFile)
    plinkFileCont <- getMakePlink(inVCF=filterMergedCont)
    if(file.exists(paste0(plinkFileCont, ".fam")) &  file.exists(paste0(plinkFileDis, ".fam"))) makeFamFile(inControl=paste0(plinkFileCont, ".fam"), inDisease=paste0(plinkFileDis, ".fam"), inOut)
    plinkFile <- getMergePlink(inControl=plinkFileCont, inDisease=plinkFileDis)
  } else {
    plinkFile <- plinkFileDis
  }

  isControl <- any(inFam)
  #inCont <- writeFam(inFam, paste0(plinkFile, ".fam"))
  outScore <- plinkGRS(inFile= plinkFile, inGRS=grsFile, inControl=isControl)
  ##system(command=paste0("rm ", filterMerged, "*"))
  #system(command=paste0("rm ", inFile, "*"))
  ##system(command=paste0("rm ", filterRsid,"*"))
  #system(command=paste0("rm ", plinkFile,"*"))
  return(outScore)
}

runGRSCalcChrPos <- function(inObjec, inDis, inFam,inCont=NULL){
  grsFile <- getGRSInputChrPos(inObjec)
  regionId <- getChrPos(inObjec)
  pgsID <- getPGSId(inObjec)
  filterMergedDis <- filterMergedPos(inFile=inDis, inName=pgsID, inSNPs=grsFile)
  plinkFileDis <- getMakePlink(inVCF=filterMergedDis)
  if(!(is.null(inCont))) {
    filterMergedCont <- filterMergedPos(inFile=inCont, inName=pgsID, inSNPs=grsFile)
    plinkFileCont <- getMakePlink(inVCF=filterMergedCont)
    if(file.exists(paste0(plinkFileCont, ".fam")) &  file.exists(paste0(plinkFileDis, ".fam"))) makeFamFile(inControl=paste0(plinkFileCont, ".fam"), inDisease=paste0(plinkFileDis, ".fam"), inOut)
    plinkFile <- getMergePlink(inControl=plinkFileCont, inDisease=plinkFileDis)
  } else {
     plinkFile <- plinkFileDis
  }
  isControl <- if(is.null(inCont)) FALSE else TRUE
  #inCont <- writeFam(inFam, paste0(plinkFile, ".fam"))
  outScore <- plinkGRS(inFile= plinkFile, inGRS=grsFile, inControl=isControl)
  #inCont <- writeFam(inFam, paste0(plinkFile, ".fam"))
  #outScore <- plinkGRS(inFile= plinkFile, inGRS=grsFile, inControl=inCont)
  #system(command=paste0("rm ", filterMerged, "*"))
  #system(command=paste0("rm ", inFile, "*"))
  #system(command=paste0("rm ", plinkFile,"*"))
  return(outScore)
}

getScore <- function(inFile){
  filtFile <- gsub(paste0(inDir, "/"), "", inFile)
  filtFile <- gsub("_score\\.profile", "", filtFile)
  filtSample <- gsub(".+_DBN_", "", filtFile)
  inData <- data.table::fread(inFile, stringsAsFactors=FALSE)
  inData[,PGS_RECORD_ID := filtSample]
  return(inData)
}


getAllScores <- function(inDir){
  inFiles <- list.files(path=inDir, ".*_DBN_PGS\\d+_score.profile", full.names=TRUE)
  totalSamples <- data.table::rbindlist(lapply(inFiles, getFile))
  aggSamples <- totalSamples[,sum(as.numeric(SCORESUM)), by=.(IID, PGS_RECORD_ID)]
  data.table::setnames(aggSamples, "V1", "PRS")
  return(aggSamples)
}


#' @export
getLatestMeta <- function(inFile=tempfile()){
  if(is.null(inFile))  inFile <- tempfile()
  if(file.exists(inFile)) return(inFile)
  inspecD <- download.file("https://pgscatalogscraper.s3-us-west-2.amazonaws.com/eu_metadata.RDS",destfile=inFile, method="wget", mode="rb")
  return(inFile)
}
  


checkFilesExist <- function(inputFile=NULL, inputRef=NULL, yamlFile = "sample.yaml", controlVCF=NULL){
  inCheck <- list(inputFile, inputRef, yamlFile)
  inFiles <- unlist(lapply(inCheck, function(x){
           return(!file.exists(x))
  }))
  if(!(is.null(controlVCF))){
       if(!(file.exists(controlVCF))) stop("Control VCF not found")
  }
  if(length(inCheck[inFiles]) >= 1){
    stop(paste(inCheck[inFiles], "does not exist"))
  }
  if(Sys.which(yaml::read_yaml(yamlFile)$bcftools) == ""){
    stop("Please specify bcftools location in sample.yaml file in current directory")
  }

  if(Sys.which(yaml::read_yaml(yamlFile)$vt) == ""){
    stop("Please specify vt location in sample.yaml file in current directory")
  }

  if(Sys.which(yaml::read_yaml(yamlFile)$plink) == ""){
    stop("Please specify plink location in sample.yaml file in current directory")
  }
  if(is.null(yaml::read_yaml(yamlFile)$outputDir)){
    print("No Output Dir Specified")
  }
  if(is.null(yaml::read_yaml(yamlFile)$tempDir)){
    print("No Temp Dir Specified")
  }
}
  
#' @export
grabScoreId <- function(inFile=NULL, inPGSID=NULL, inPGSIDS=NULL, inRef=NULL, inYamlFile="sample.yaml", inCL=NULL, inControl=NULL, inMeta=NULL){
  checkFilesExist(inputFile=inFile, inputRef=inRef, yamlFile=inYamlFile, controlVCF=inControl)
  if (!(is.null(inPGSIDS))){
    # This is first because of weird bug where pgs-id value replicates pgs-id-file even though it is null
    file <- data.table::fread(inPGSIDS, stringsAsFactors=F, header=F)
    file <- file$V1
  } else if(!(is.null(inPGSID))) {
    file <- inPGSID
  } else {
    print("Calculating all Scores")
    file <- NA
  }
  inspecTemp <- getLatestMeta(inMeta)
  inspecFile <- readRDS(inspecTemp)
  outDir <- if(!(is.null(yaml::read_yaml(inYamlFile)$tempDir))){
    gsub("\\/$", "",yaml::read_yaml(inYamlFile)$tempDir)
  } else if(!(is.null(yaml::read_yaml(inYamlFile)$outputDir))){
    gsub("\\/$", "",yaml::read_yaml(inYamlFile)$outputDir)
  } else {
    dirname(inFile)
  }
  famFile <- FALSE
  disSample <- NULL
  normFile <- 
    baseNorm(inVCF=inFile,
           inRef=inRef,
           outFile=paste(outDir, gsub("\\.[a-zA-Z']+(\\.gz)?$","_norm", basename(inFile)),sep="/"),
           inYaml=inYamlFile)
  if(!(is.null(inControl))) {
    controlFile <- 
      baseNorm(inVCF=inControl,
             inRef=inRef,
             outFile=paste(outDir, gsub("\\.[a-zA-Z']+(\\.gz)?$","_norm", basename(inControl)),sep="/"),
             inYaml=inYamlFile)
    disSample <- getDisSample(inControl=controlFile[[1]], inDisease=normFile[[1]])
    famFile <- TRUE
  }
  ##while(!(resolved(inspecFile) & resolved(normFile))) {}
     ## Need to fix ID
     #normFile <- future::value(normFile)
     inputFile <- inspecFile
     scoreFiles <- lapply(inputFile, function(inObjec){
      scoreFiles <- data.table::data.table() 
      if(inObjec@`pgsId` %in% file | is.na(file)){
        if(!(is.null(inControl))) {
          combFrame <- data.table::setDT(list(normFile=unlist(normFile), controlFile=unlist(controlFile), inFamFile=rep(famFile, times=length(unlist(normFile)))))
        } else {
          combFrame <- data.table::setDT(list(normFile=unlist(normFile), inFamFile=rep(famFile, times=length(unlist(normFile)))))
        }
        scoreFiles <- foreach::foreach(d=iterators::iter(unique(combFrame),by='row'), .export = c("getPGSType","runGRSCalc", "runGRSCalcChrPos", "getChrPos", "getGRSInputChrPos", "getGRSInputRsID", "getGRSInputRsID", "getRSIds", "getPGSId", "filterMerged", "filterMergedPos", "filterMergedReg", "getMakePlink", "plinkGRS", "baseIndex", "makeFamFile", "getMergePlink", "inControl"), .packages=c("data.table", "yaml", "assertthat"), .combine=rbind) %dopar% {
          if(getPGSType(inObjec) == "rsID"){
            if(any(("controlFile" %in% names(d)))){
             scoreFile <- runGRSCalc(inObjec=inObjec, inDis=d$normFile, inCont=d$controlFile,inFam=d$inFamFile)
           } else {
             scoreFile <- runGRSCalc(inObjec=inObjec, inDis=d$normFile, inFam=d$inFamFile)
           }
           if(is.null(scoreFile)) return(NULL)
           return(data.table::data.table(inFile=scoreFile,inRecord=getPGSId(inObjec)))
          } else {
            if(any(("controlFile" %in% names(d)))){
             scoreFile <- runGRSCalcChrPos(inObjec=inObjec, inDis=d$normFile, inCont=d$controlFile,inFam=d$inFamFile)
            } else {
             scoreFile <- runGRSCalcChrPos(inObjec=inObjec, inDis=d$normFile,inFam=d$inFamFile)
            }
           if(is.null(scoreFile)) return(NULL)
           return(data.table::data.table(inFile=scoreFile,inRecord=getPGSId(inObjec)))
          }
      }
      }
      return(scoreFiles)
    })

  scoreFiles <- unique(rbindlist(scoreFiles))
  scoreDat <- getAggDf(scoreFiles$inFile, scoreFiles$inRecord)
  #if(!(is.null(inControl))){
    #inControl <- grabScoreControl(inPGSID=inPGSID, inPGSIDS=inPGSIDS, inRef=inRef, inYamlFile=inYamlFile, inCL=inCL, inControl=inControl,inRDS=inspecFile)
  #}
  #mergeData <- mapply(getMergePlink, , SIMPLIFY=FALSE)
  #controlDIter <- iterators::iter(inControl, by='row')
  #scoreDIter <- iterators::iter(scoreFiles, by='row')
  #mergeData <- getRow(controlDIter, scoreDIter)
  #mergeData <- lapply(mergeData, function(x) getMergePlink(inMerge))

  #getMergePlink(scoreFiles, inControl)
  if(!(is.null(disSample))){
    inControl <- scoreDat[!(IID %in% as.character(disSample)),]
    scoreDat <- scoreDat[IID %in% as.character(disSample),]
  }
  setPlots(scoreDat, inControl, inOutDir=outDir)
 
}
