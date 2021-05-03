getFile <- function(inFile, inRecord, midSamples){
  if(length(inFile) != 0 & file.exists(inFile)){
    inData <- data.table::fread(inFile, stringsAsFactors=FALSE)
    inData[,PGS_RECORD_ID := inRecord]
    inData <- inData[!(IID %in% midSamples),]
    #inData[,subject_type := "diseased"]
    return(data.table(inData))
  } else {
    return(NULL)
  }
}

getAggDf <- function(inFiles, inRecords,midSamples){
  totalSamples <- data.table::rbindlist(mapply(getFile, inFiles, inRecords, midSamples,SIMPLIFY=FALSE))
  aggSamples <- totalSamples[,sum(as.numeric(SCORESUM)), by=.(IID, PGS_RECORD_ID)]

  setnames(aggSamples, "V1", "PRS")
  return(aggSamples)
}
