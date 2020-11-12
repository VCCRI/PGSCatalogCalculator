getFile <- function(inFile, inRecord){
  require(data.table)
  if(length(inFile) != 0 & file.exists(inFile)){
    inData <- fread(inFile, stringsAsFactors=FALSE)
    inData[,PGS_RECORD_ID := inRecord]
    #inData[,subject_type := "diseased"]
    return(inData)
  } else {
    return(NULL)
  }
}

getAggDf <- function(inFiles, inRecords){
  require(data.table)

  totalSamples <- rbindlist(mapply(getFile, inFiles, inRecords))

  aggSamples <- totalSamples[,sum(as.numeric(SCORESUM)), by=.(IID, PGS_RECORD_ID)]

  setnames(aggSamples, "V1", "PRS")
  return(aggSamples)
}
