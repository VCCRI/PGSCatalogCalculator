require(data.table)
require(assertthat)


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
  #inFam <- inFam[,..needCols]
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
  print(unique(inTotal$PGS_RECORD_ID))
  assert_that(length(unique(inTotal$subject_type)) == 2)
  assert_that(all.equal(sort(unique(inTotal$subject_type)), c("case", "control")))
  inTotal <- tryCatch({
    inTotal$bin_PRS <- -1
    breaksD <- as.numeric(as.character(quantile(as.numeric(inTotal[subject_type == "control",PRS]), probs = seq(0.25, 1, by = 1/inNumBins))))
    seqC <- as.numeric(seq(0.25,1-1/inNumBins, by=1/inNumBins))
    inTotal[,bin_PRS := as.numeric(as.character(cut(as.numeric(PRS), breaks=breaksD, labels=seqC, include.lowest = TRUE, right=TRUE)))]
    inTotal[,bin_PRS := ifelse(is.na(bin_PRS), 1, bin_PRS)]
    assert_that(nrow(inTotal[bin_PRS == -1,]) == 0)
    inTotal
  }, error = function(e){
    ##Default behaviour if inNumBins is not defined
    inTotal[,bin_PRS := as.numeric(as.character(cut(inTotal$PRS, breaks=c(quantile(inTotal$PRS, probs = seq(0, 1, by = 0.1))), labels=seq(0,0.9, by=0.1), include.lowest = FALSE, right=TRUE)))]
  })
  minBin <- min(as.numeric(inTotal$bin_PRS))
  inTotal[,is_bott := ifelse(bin_PRS == minBin, TRUE,FALSE)]
  assert_that(length(unique(inTotal[is_bott == TRUE, bin_PRS])) == 1)
  assert_that(length(unique(inTotal[is_bott == FALSE, bin_PRS])) == (inNumBins-1))
  assert_that(!(c("0.25") %in% as.character(unique(inTotal[is_bott == FALSE, bin_PRS]))))
  assert_that(c("0.25") %in% as.character(unique(inTotal[is_bott == TRUE, bin_PRS])))
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

createGLM <- function(inQuant, inData, inNumBins=4){
  ##Useless fking function
  ##Retarded, what was I thinking
  require(glm2)
  inData$bin_data <- inData$bin_PRS*inNumBins
  assertthat::assert_that(unique(c("is_bott", "PRS") %in% colnames(inData)))
  assertthat::assert_that(length(unique(inData$is_bott)) == 2, msg=paste("Unexpected Values for bottom quantiles with Quantile", inQuant))
  if(inQuant == 1){
    assertthat::assert_that(length(unique(inData[bin_data == inQuant | bin_data == 1,bin_data])) == 1)
  } else {
    assertthat::assert_that(length(unique(inData[bin_data == inQuant | bin_data == 1,bin_data])) == 2)
  }
  print(summary(glm(factor(is_bott) ~ PRS, data=inData[bin_data == inQuant | bin_data == 1,], family="binomial")))
  print(summary(glm(factor(is_bott) ~ 0+bin_PRS, data=inData, family="binomial")))
  #browser()
  #return(list(logitModel = glm(factor(is_bott) ~ PRS, data=inData[bin_data == inQuant | bin_data == 1,], family="binomial"), quant=inQuant))
}

createTable <- function(inQuant, inData, inNumBins=4){
  inData$bin_data <- inData$bin_PRS*inNumBins
  botData <- inData[bin_data == 1 ,uniqueN(IID), by=.(subject_type, bin_data)]
  botData$category <- "z_bottom"
  quantData <- inData[bin_data == inQuant ,uniqueN(IID), by=.(subject_type, bin_data)]
  quantData$category <- "quantile"
  contTable <- rbind(botData, quantData)
  if(nrow(contTable[subject_type == "case",])  == 0){
    nullRow <- data.table(subject_type = c("case", "case"), bin_data
= rep(unique(contTable$bin_data),2), V1=c(0, 0),  category=c("z_bottom", "quantile"))
    contTable <- rbind(contTable, nullRow)
    contTable$V1 <- as.numeric(contTable$V1)
  }
  return(list(contTable=dcast(contTable, subject_type ~ category, value.var="V1", fun.aggregate = sum), quant=inQuant))
}

getORContT <- function(inData){
  require(questionr)
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
  require(ggplot2)
  require(cowplot)
  inDate <- "2020-02-17"
  inDir <- paste0("/home/114/sk3015/Analysis/makeGRS/quantilePlot/",inDate)
  #dir.create(inDir, showWarnings=FALSE)
  inFile <- paste0(inDir,"/boxplot.png")
  theme_set(theme_cowplot())
  png(filename=inFile, width=1920, height=1080)
  ggplot(data.frame(inData, stringsAsFactors=F), aes(x=splitNo, y=oddsRatio)) +
  geom_line()+
  geom_pointrange(aes(ymin=lowerBound, ymax=upperBound))
  dev.off()
}

totalData <- fread("diseased.csv", stringsAsFactors=F)
controlData <- fread("control_samples.csv", stringsAsFactors=F)
lapply(unique(totalData$PGS_RECORD_ID), function(x){
  #testData <- readScore("~/total_1.csv", "~/plink_missvar.fam")
  #if(x %in% c("PGS000005", "PGS000007", "PGS000008", "PGS000009", "PGS000016", "PGS000035")){
  if(x %in% c("PGS000009")){
  testData <- readScore(totalData[PGS_RECORD_ID == x,], controlData[PGS_RECORD_ID == x,])
  #print(testData)
  #inGLMModel <- lapply(1:4, function(x)createGLM(inQuant=x, inData=testData))
  inContTable <- lapply(1:4, function(x)createTable(inQuant=x, inData=testData))
  browser()
  orDF <- do.call("rbind", lapply(inContTable, getORContT))
  print(orDF)
  #getQuantilePlot(orDF)
  #orDF <- do.call("rbind", lapply(inGLMModel, getORGLM))
  ## Putting outside because R
  #print(head(orDF))
  require(ggplot2)
  require(cowplot)
  theme_set(theme_cowplot())
  inDate <- "2020-08-10"
  inDir <- paste0("~/quantilePlot/",inDate)
  #dir.create(inDir, showWarnings=FALSE)
  inFile <- paste0(inDir,"/", x, "-boxplot.png")
  print(inFile)
  browser()
  #theme_set(theme_cowplot())
  #png(filename=inFile, width=1920, height=1080)
 #orDF[1,4] <- 1
 #orDF[2,4] <- 1
 #orDF[3,4] <- 1
 #orDF[4,4] <- 31.33891
 #orDF[4,3] <- 31.33891
  ggplot(data.frame(orDF, stringsAsFactors=F), aes(x=splitNo, y=oddsRatio)) +
  geom_pointrange(aes(ymin=lowerBound, ymax=upperBound))
  #dev.off()
  }
})
