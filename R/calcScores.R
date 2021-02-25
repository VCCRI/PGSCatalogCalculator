getFile <- function(inFile, inRecord){
  if(length(inFile) != 0 & file.exists(inFile)){
    inData <- data.table::fread(inFile, stringsAsFactors=FALSE)
    inData[,PGS_RECORD_ID := inRecord]
    #inData[,subject_type := "diseased"]
    return(inData)
  } else {
    return(NULL)
  }
}

getAggDf <- function(inFiles, inRecords){
  totalSamples <- data.table::rbindlist(mapply(getFile, inFiles, inRecords))
  aggSamples <- totalSamples[,sum(as.numeric(SCORESUM)), by=.(IID, PGS_RECORD_ID)]

  setnames(aggSamples, "V1", "PRS")
  return(aggSamples)
}
