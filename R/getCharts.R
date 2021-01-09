readScore <- function(inData, inControl, inNumBins=4){
	## Function that reads in data and fam file and will split
	## data into 4 equal parts unless otherwise specified
	#inData <- fread(inFile, stringsAsFactors=F)
  inData <- inData
  inControl <- inControl
  #inFam <- fread(inFamFile, stringsAsFactors=F, header=F)
	## Assert that data is what we expect
	# Num of cols are 4
	# No duplicated samples and no na PRS scores
	assertthat::assert_that(ncol(inData) == 3)
	assertthat::assert_that(any(duplicated(inData$IID)) == FALSE)
	assertthat::assert_that(any(is.na(as.numeric(inData$PRS))) == FALSE)
  # If there is overlap between the base file and target file
  # use controls to create deciles else use target set
  setkey(inData, IID)
  setkey(inControl, IID)
  needCols <- c("V2", "V6")
	assertthat::assert_that(ncol(inControl) == 3)
	assertthat::assert_that(any(duplicated(inControl$IID)) == FALSE)
	assertthat::assert_that(any(is.na(as.numeric(inControl$PRS))) == FALSE)
  #inFam <- inFam[,needCols,with=FALSE]
  #inFam$V6 <- as.numeric(inFam$V6)
  #assertthat::assert_that(any(is.na(as.numeric(inFam$V6))) == FALSE)
  nrowInit <- nrow(inData)
  #inData <- inData[inFam, nomatch=0]
  inControl[,subject_type := "control"]
  inData[,subject_type := "case"]
  inTotal <- rbind(inData, inControl)
  #assert_that(nrow(inData) == nrowInit)
  ## Need to be careful with ifelse
  inTotal$PRS <- as.numeric(inTotal$PRS)
  assertthat::assert_that(any(is.na(as.numeric(inTotal$PRS))) == FALSE)
  assertthat::assert_that(length(unique(inTotal$subject_type)) == 2)
  assertthat::assert_that(all.equal(sort(unique(inTotal$subject_type)), c("case", "control")))
  inTotal <- tryCatch({
    inTotal$bin_PRS <- -1
    breaksD <- as.numeric(as.character(quantile(as.numeric(inTotal[subject_type == "control",PRS]), probs = seq(0.25, 1, by = 1/inNumBins))))
    seqC <- as.numeric(seq(0.25,1-1/inNumBins, by=1/inNumBins))
    inTotal[,bin_PRS := as.numeric(as.character(cut(as.numeric(PRS), breaks=breaksD, labels=seqC, include.lowest = TRUE, right=TRUE)))]
    inTotal[,bin_PRS := ifelse(is.na(bin_PRS), 1, bin_PRS)]
    assertthat::assert_that(nrow(inTotal[bin_PRS == -1,]) == 0)
    inTotal
  }, error = function(e){
    ##Default behaviour if inNumBins is not defined
    inTotal[,bin_PRS := as.numeric(as.character(cut(inTotal$PRS, breaks=c(quantile(inTotal$PRS, probs = seq(0, 1, by = 0.1))), labels=seq(0,0.9, by=0.1), include.lowest = FALSE, right=TRUE)))]
  })
  minBin <- min(as.numeric(inTotal$bin_PRS))
  inTotal[,is_bott := ifelse(bin_PRS == minBin, TRUE,FALSE)]
  assertthat::assert_that(length(unique(inTotal[is_bott == TRUE, bin_PRS])) == 1)
  assertthat::assert_that(length(unique(inTotal[is_bott == FALSE, bin_PRS])) == (inNumBins-1))
  assertthat::assert_that(!(c("0.25") %in% as.character(unique(inTotal[is_bott == FALSE, bin_PRS]))))
  assertthat::assert_that(c("0.25") %in% as.character(unique(inTotal[is_bott == TRUE, bin_PRS])))
  #assert_that(nrow(unique(inTotal)) == nrow(inTotal))
  inTotal <- unique(inTotal)
  #fwrite(inTotal[bin_PRS == 0.75& subject_type=="case",],paste0("~/quantilePlot/",unique(inTotal$PGS_RECORD_ID), ".csv"))
	return(inTotal)
}

getORGLM <- function(inGLMMod){
  inGLM <- inGLMMod$logitModel
  inQuant <- inGLMMod$quant
  if(inQuant != 1){
    retFrame <- tryCatch({
      oddsR <- questionr::odds.ratio(inGLM, level=0.95)
      assertthat::assert_that(length(oddsR$OR) == 2)
      assertthat::assert_that(is.numeric(inQuant))
      retFrame <- cbind(splitNo=inQuant, lowerBound=oddsR$`2.5 %`[2], oddsRatio=oddsR$OR[2], upperBound=oddsR$`97.5 %`[2])
    }, error = function(e){
      inORs <- exp(summary(inGLM)$coefficients["PRS",1] +
       qnorm(c(0.025,0.5,0.975)) * summary(inGLM)$coefficients["PRS",2])
      retFrame <- cbind(splitNo=inQuant, lowerBound=inORs[1], oddsRatio=inORs[2], upperBound=inORs[3])
    })
  } else {
    retFrame <- cbind(splitNo=inQuant, lowerBound=1, oddsRatio=1, upperBound=1)
  }
  # No DT because of doubles
  return(retFrame)
}

createTable <- function(inQuant, inData, inNumBins=4){
  inData$bin_data <- inData$bin_PRS*inNumBins
  botData <- inData[bin_data == 1 ,uniqueN(IID), by=.(subject_type, bin_data)]
  botData$category <- "z_bottom"
  quantData <- inData[bin_data == inQuant ,uniqueN(IID), by=.(subject_type, bin_data)]
  quantData$category <- "quantile"
  contTable <- rbind(botData, quantData)
  if(nrow(contTable[subject_type == "case",])  == 0){
    nullRow <- data.table::data.table(subject_type = c("case", "case"), bin_data
= rep(unique(contTable$bin_data),2), V1=c(0, 0),  category=c("z_bottom", "quantile"))
    contTable <- rbind(contTable, nullRow)
    contTable$V1 <- as.numeric(contTable$V1)
  }
  return(list(contTable=data.table::dcast(contTable, subject_type ~ category, value.var="V1", fun.aggregate = sum), quant=inQuant))
}

getORContT <- function(inData){
  contTable <- inData$contTable
  contTable$subject_type <- NULL
  contTable <- as.matrix(contTable)
  inQuant <- inData$quant
  if(inQuant != 1){
    oddsR <- questionr::odds.ratio(contTable, level=0.95)
    assertthat::assert_that(length(oddsR$OR) == 1)
    assertthat::assert_that(is.numeric(inQuant))
    retFrame <- cbind(splitNo=inQuant, lowerBound=oddsR$`2.5 %`[1], oddsRatio=oddsR$OR[1], upperBound=oddsR$`97.5 %`[1])
  } else {
    retFrame <- cbind(splitNo=inQuant, lowerBound=1, oddsRatio=1, upperBound=1)
  }
  return(retFrame)
}

getQuantilePlot <- function(inData){
  inDate <- "2020-02-17"
  inDir <- paste0("/home/114/sk3015/Analysis/makeGRS/quantilePlot/",inDate)
  #dir.create(inDir, showWarnings=FALSE)
  inFile <- paste0(inDir,"/boxplot.png")
  ggplot2::theme_set(cowplot::theme_cowplot())
  png(filename=inFile, width=1920, height=1080)
  ggplot2::ggplot(data.frame(inData, stringsAsFactors=F), ggplot2::aes(x=splitNo, y=oddsRatio)) +
  ggplot2::geom_line()+
  ggplot2::geom_pointrange(ggplot2::aes(ymin=lowerBound, ymax=upperBound))
  dev.off()
}

getPlots <- function(inFiles, inRecord, inControl){
  disData <- getAggDf(inFiles, inRecord)
  if(is.null(inControl)){
    controlData <- data.table::fread(system.file("extdata", "control_samples.csv", package="PGSCatalogDownloader"), stringsAsFactors=F)
  } else {
    controlData <- inControl
  }
  lPLots <- lapply(list(unique(disData$PGS_RECORD_ID)), function(x){
    scoreData <- readScore(disData[PGS_RECORD_ID == x,], controlData[PGS_RECORD_ID == x,])
    scoreData[,risk := ifelse(as.numeric(bin_PRS) <= 0.25, 'Low Risk', ifelse(as.numeric(bin_PRS) > 0.75, "High Risk", "Medium Risk"))]
    aggScore <- scoreData 
    aggScore <- scoreData[subject_type == "case",]
    needCols <- c("IID", "PRS","risk")
    if(is.null(inControl)){
      needScore <- aggScore[,needCols, with=FALSE]
    } else {
      needCols <- c("IID", "subject_type","PRS","risk")
      needScore <- scoreData[,needCols, with=FALSE]
      colnames(needScore) <- c("Sample", "Subject Type", "PRS","Risk")
    }
    fwrite(needScore, "sample_out.csv")
    inContTable <- lapply(1:4, function(x)createTable(inQuant=x, inData=scoreData))
    orDF <- do.call("rbind", lapply(inContTable, getORContT))
    orDF <- setDT(data.frame(orDF, stringsAsFactors=F))
    ##ifelse not used in case of factor
    orDF[,upperBound := ifelse(is.infinite(upperBound), as.numeric(lowerBound), as.numeric(upperBound))]
    ggplot2::theme_set(cowplot::theme_cowplot())
    return(ggplot2::ggplot(data.frame(orDF, stringsAsFactors=F), ggplot2::aes(x=splitNo, y=oddsRatio)) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin=lowerBound, ymax=upperBound)) + ggplot2::ggtitle(x))
  })
  return(lPLots)
}

setPlots <- function(inFiles, inControl=NULL){
  allPlots <- getPlots(inFiles$inFile, inFiles$inRecord, inControl)
  inTime <- Sys.time()
  inDir <- paste0("quantilePlot/")
  dir.create(inDir, showWarnings=FALSE)
  inFile <- paste0(inDir, "boxplot.png")
  print(inFile)
  #png(filename=inFile, width=1920, height=1080)
  ggplot2::ggsave(filename=inFile, plot=cowplot::plot_grid(plotlist=allPlots))
  #dev.off()
}
