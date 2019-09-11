#' GapFilter
#'
#' A series of filters which removes a large subset of peaks from non-targeted analysis based on user input
#' Specifically focused on removing exposure compounds and their metabolites from toxicological experiments
#'
#' @param batchCorr If the user is planning on using batchCorr after GapFilter they can setup the document with injections included already and turn 'batchCorr' to TRUE
#' @param samplesToGrab Identifier in sample name of the files (e.g. "QC") which GapFilter will use to determine whether a feature should be kept in or not
#' @param ratio Percentage of the sample type used for which need to have gap status "No gap" for a feature to be kept
#' @param
#'
#' @return Outputs a result file with all the peaks and samples split into two worksheets depending on whether they were gap-filtered or not
#'
#' @import openxlsx
#' @importFrom grDevices graphics utils
#'
#' @export
#'
GapFilter<-function(batchCorr=FALSE, samplesToGrab="", ratio=100){

  library(openxlsx)
  ###Selecting the input-file
  cat("Inpupt the raw data file structured according to the example file connected to the package:\n\n")
  filename<-file.choose()
  rawData=read.xlsx(xlsxFile = filename, sheet=1)

  ###Structuring samples
  samps<-rawData[-1,1]

  if(batchCorr==TRUE){
    inj<-as.numeric(rawData[-1,2])
    peakAreas<-rawData[-1,-c(1:2)]
  } else {
    peakAreas<-rawData[-1,-1]
  }

  ###Incorporating RTs in to feature names in a way recognizible easily by batchCorr
  for(i in 1:length(peakAreas[1,])){
    colnames(peakAreas)[i]<-paste(colnames(peakAreas)[i],"_",rawData[1,(i+2)],"-",i)
  }

  ####Gap status filter prior to batchCorr####
  gapFilt<-read.xlsx(xlsxFile=filename, sheet=2, rowNames = TRUE)

  gapFilt<-gapFilt[grepl(samplesToGrab, rownames(gapFilt)),]
  peaksToRemove<-vector()

  ###IF interested in having different settings for how to filter, need to be introduced here
  ###Removes a feature which has more fully missing features than the ratio set by the user
  for(i in 1:length(gapFilt[1,])){
    if((sum(grepl("Full gap", gapFilt[,i]))/length(gapFilt[,i])*100) >= ratio){
      peaksToRemove<-c(peaksToRemove, i)
    }
  }

  #Output showing how many features were removed. Also, create .xlsx-document with
  peaksRemoved<-peakAreas[,peaksToRemove]
  cat((length(peaksRemoved[1,])/length(peakAreas[1,])*100)," % of features removed due to not being found in ", ratio,"% of user-specified samples.")
  peakAreas<-peakAreas[,-peaksToRemove]

  #Formatting data for print for features to keep
  if(batchCorr==TRUE){
  notFiltered<-cbind(samps, inj,peakAreas)
  } else {
  notFiltered<-cbind(samps,peakAreas)
  }

  #Formatting data for print for features that have been filtered
  if(batchCorr==TRUE){
    filtered<-cbind(samps, inj,peaksRemoved)
  } else {
    filtered<-cbind(samps, peaksRemoved)
  }


  #Creating workbook and saving to drive
  finalDocument <- createWorkbook()
  addWorksheet(finalDocument, sheetName="FeaturesKept")
  addWorksheet(finalDocument, sheetName="GapFilteredFeatures")

  writeData(finalDocument, sheet=1, x=notFiltered)
  writeData(finalDocument, sheet=2, x=filtered)

  saveWorkbook(finalDocument,"GapFilterResults.xlsx", overwrite = TRUE)
}
