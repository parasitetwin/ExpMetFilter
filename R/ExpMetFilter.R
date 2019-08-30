#' ExpMetFilter
#'
#' A series of filters which removes a large subset of peaks from non-targeted analysis based on user input
#' Specifically focused on removing exposure compounds and their metabolites from toxicological experiments
#'
#' @param File The path to, including the name and file-ending of, where the result file should be created
#' @param report When set to 'True' will output a PDF file with how many features were removed due to each filtering process
#'
#' @return Outputs a result file with all the peaks and sample split into whether they were filtered out or not
#'
#' @import openxlsx
#' @import RAMClustR
#' @importFrom grDevices graphics utils
#' @export
#'
ExpMetFilter<-function(File=NULL, report=FALSE){
  
  #File=NULL
  # KeepNegVal=FALSE
  # RemoveMets=TRUE
  # MergeAllSets=TRUE
  # report=FALSE
  # PseudoMS2filter=TRUE
  
  #Metacode
  #Input document should have one outlining document with the number of filters, what they are, the internal standards etc.
  
  #Loop for the number of worksheets present in document
  #Applying thershold filter to each dataset according to outlining document        (Filters filled in for all worksheets           #Collected information: Removed peaks + threshSummary
  #Removing all features which contain any negative value                           (Switch Yes/No)                                 #Collected information: Removed peaks + negSummary
  #Splitting each dataset with the IS designated to each dataset acc. out. doc.     (IS filled in for all worksheets)
  
  #Removing probable metabolites of exposure compound from each dataset             (List of probable exposure compound-metabolites for each dataset)
  #Correlating exposure compound at certain RT with plausible in-source fragments   (List of compounds that have gone through batchCorr with RTs & exact masses)
  
  #Output in format ready for MUVR + printing a report for all the stuff performed above (Number of features removed in each step etc.)
  
  #requireNamespace()
  library(openxlsx)
  library(RAMClustR)
  
  cat("Choose the unified ExpMISFFilter-file containing the datasets you want to clean:\n\n")
  
  ####Choosing file and establishing number of datasets
  if(is.null(File)==TRUE){
    fileName<-file.choose()
  } else if (is.character(File)==TRUE){
    fileName<-File
  } else if (is.list(File)==TRUE){
    
  }
  
  options(warn=-1)
  
  #Calculating number of datasets
  nDataset<-(length(getSheetNames(fileName))-1)/2
  
  #Reading and storing metainfo for later purposes
  metaInfo<-read.xlsx(fileName, sheet=1)
  datasetNames<-as.character(colnames(metaInfo[-c(1,4)]))
  completeDataset<-list()
  allRemovedFeatures<-list()
  
  ISmass<-metaInfo[1,-c(1,4)]
  ThresholdFilter<-metaInfo[2,-c(1,4)]
  removeMets<-metaInfo[3,-c(1,4)]
  InSourceFragRemove<-metaInfo[4,-c(1,4)]
  sampleFilter<-metaInfo[5,-c(1,4)]
  dPPM <- metaInfo[6,-c(1,4)]
  expMass<- metaInfo[7,-c(1,4)]
  MergeAllSets <- metaInfo[8,-c(1,4)]
  KeepNegVal <- metaInfo[9,-c(1,4)]
  
  #Checking if number of samples are the same, if merge datasets has been indicated by the user
  # if(is.na(MergeAllSets[1]) == FALSE){
  #   lengths<-NA
  #   for(n in 1:nDataset){
  #     tempSampleNames<-read.xlsx(fileName, sheet=(n*2))[,1]
  #     lengths[n]<-length(tempSampleNames[-which(grepl(sampleFilter[1,n], tempSampleNames))])
  #   }
  #   if((nDataset-sum(duplicated(lengths))-1)==0){
  #     cat("The number of samples left in each dataset after applying respective filters are equal in numbers.\n\n")
  #   } else {
  #     stop("The number of samples in datasets, post-filtering, are not equal. Aborting ExpMISFFilter.\n\n")
  #   }
  # }
  
  ####Beginning of datafiltering loop####
  for(n in 1:nDataset){
    #Reading the datasets from the xlsx-file
    #Also collecting relevant threshold filters and IS
    cat("###Processing the", colnames(metaInfo)[n+1],"dataset###\n\n")
    
    otu<-read.xlsx(fileName, sheet=(n*2))
    samps<-otu[,1]
    if(any(is.na(samps))==TRUE){
      samps[is.na(samps)]<-"RT"
    }
    Xotu<-otu[,-1]
    
    #Turning lists of characters into data.frames with numeric
    if(typeof(Xotu)=='list'){
      Xotu<-lapply(Xotu, function(x) as.numeric(as.character(x)))
      Xotu<-as.data.frame(Xotu)
    }
    
    ISmass<-as.character(round(as.numeric(metaInfo[1,(n+1)]),5))
    ThresholdFilter<-as.numeric(metaInfo[2,(n+1)])
    removeMets<-metaInfo[3,(n+1)]
    InSourceFragRemove<-metaInfo[4,(n+1)]
    sampleFilter<-metaInfo[5,(n+1)]
    dPPM<-as.double(metaInfo[6,(n+1)])
    expMass<-metaInfo[7,(n+1)]
    MergeAllSets<-metaInfo[8,(n+1)]
    KeepNegVal<-metaInfo[9,(n+1)]
    blankFilter<-metaInfo[10,(n+1)]
    
    #Checking for several sample types to be removed
    if(grepl(",",sampleFilter)){
      sampleFilter<-strsplit(sampleFilter,",")
      temp<-NA
      for(i in 1:length(sampleFilter[[1]])){
        temp<-c(temp,sampleFilter[[1]][i])
      }
      sampleFilter<-temp[-1]
    }
    
    #Building meta-information document for correlation to exposure compound removal
    if(is.na(InSourceFragRemove)==FALSE){
      experiment<-list()
      experiment$`design`<-as.data.frame(metaInfo[c(13:17),(n+1)])
      rownames(experiment$`design`)<-metaInfo[c(13:17),1]
      if(metaInfo[17,2]=="GC-MS"){
        experiment$instrument<-as.data.frame(metaInfo[c(19:30),(n+1)])
        rownames(experiment$instrument)<-metaInfo[c(19:30),1]
      } else if (metaInfo[17,2]=="LC-MS") {
        experiment$instrument<-as.data.frame(metaInfo[c(32:45),(n+1)])
        rownames(experiment$instrument)<-metaInfo[c(32:45),1]
      } else {
        PseudoMS2filter = FALSE
      }
    }
    
    #####Filter based on feature intensity####
    #If the maximum intensity in any sample is bigger than the filter then check if it is the first compound to fulfill that requirement
    #If so:  Create first Xotu.filt-input
    #If not: Cbind at the end of the data-frame
    #Also collecting all the discarded peaks to be printed in report
    
    if(is.na(ThresholdFilter)==TRUE){
      cat("User is not applying any feature intensity filter.\n\n")
      strA<-paste0("User is not applying any feature intensity filter.\n\n")
    }else{
      Xotu.1<-NULL
      thresholdFeatures<-NA
      
      for(l in 1:length(Xotu[1,])){
        if(max(Xotu[,l])>ThresholdFilter){
          
          if(is.null(Xotu.1)==TRUE){
            Xotu.1<-as.data.frame(Xotu[,l])
            colnames(Xotu.1)[1]<-colnames(Xotu)[l]
          } else {
            Xotu.1<-cbind(Xotu.1,Xotu[,l])
            colnames(Xotu.1)[length(Xotu.1[1,])]<-colnames(Xotu)[l]
          }
          
        } else {
          thresholdFeatures<-as.data.frame(cbind(thresholdFeatures, Xotu[,l]))
          colnames(thresholdFeatures)[length(thresholdFeatures[1,])]<-colnames(Xotu)[l]
        }
      }
      cat("Summary of threshold-filter:\n")
      cat((100-(length(Xotu.1[1,])/length(Xotu[1,])*100)),"% of the features have been filtered.\n\n")
      strA<-paste0("Summary of threshold-filter:\n",(100-(length(Xotu.1[1,])/length(Xotu[1,])*100)),"% of the features have been filtered.\n\n")
      threshSummary<-(100-(length(Xotu.1[1,])/length(Xotu[1,])*100))
      Xotu<-Xotu.1
    }
    
    #Report-pieces gathered: Removed features + threshSummary
    #Xotu = Outgoing
    
    ####Splitting all masses with internal standard intensities####
    #If user has set 'ISCorr' to TRUE, then all intensities split with IS
    #IS-mass determined in metaInfo and collected for each dataset from there
    #If user has IS for one dataset but not another, the one without should be signified with a '0' in the document
    
    if(is.na(ISmass)==TRUE){
      cat("User is not employing any IS to correct for evaporation or uneven number of injections. No correction applied.\n\n")
      strB<-paste0("User is not employing any IS to correct for evaporation or uneven number of injections. No correction applied.\n\n")
    }else{
      RTs<-Xotu[1,]
      Xotu<-Xotu[-1,]
      Xotu.2<-as.data.frame(Xotu[,1])
      
      #Avoiding colnames turning into X-versions of themselves
      if(is.integer(which(grepl(ISmass, colnames(Xotu))==TRUE)) && length(which(grepl(ISmass, colnames(Xotu))==TRUE)) == 0L){
        oldNames<-colnames(Xotu)
        colnames(Xotu)<-substring(colnames(Xotu), first=2)
        colnames(Xotu)<-round(as.double(colnames(Xotu)),5)
        
        ISint<-Xotu[,which(grepl(ISmass, colnames(Xotu))==TRUE)]
      } else {
        ISint<-Xotu[,which(grepl(ISmass, colnames(Xotu))==TRUE)]
      }
      
      
      for(l in 1:length(Xotu[1,])){
        Xotu.2<-cbind(Xotu.2,(Xotu[,l]/ISint))
        colnames(Xotu.2)[length(Xotu.2[1,])]<-colnames(Xotu)[l]
      }
      
      #Started with a bogus-column, so has to be removed
      Xotu.2<-as.data.frame(Xotu.2[,-1])
      colnames(Xotu.2)<-colnames(RTs)
      Xotu.2<-rbind(RTs,Xotu.2)
      Xotu<-Xotu.2
      cat("Summary of IS-correction:\nAll feature intensities split through the intensities of mass:", ISmass,"\n\n")
      strB<-paste0("Summary of IS-correction:\nAll feature intensities split through the intensities of mass:", ISmass,"\n\n")
    }
    #Report pieces gathered features: None
    #Xotu = Outgoing data frame
    
    #####PseudoMS2 filter#####
    if(is.na(InSourceFragRemove)==TRUE){
      cat("User is not interested in using the PseudoMS2 filter\n\n")
      removedCorPeaks<-NA
      strB<-paste0("User is not interested in using the PseudoMS2 filter\n\n")
    } else {
      
      #colnames(Xotu)<-gsub(".1$", "", colnames(Xotu))
      #colnames(Xotu)<-gsub(".2$", "", colnames(Xotu))
      colnames(Xotu)<-paste(colnames(Xotu),"_",Xotu[1,], sep="")
      rownames(Xotu)<-samps
      Xotu<-Xotu[-1,]
      
      value<-integer(0)
      
      #experiment2 <- defineExperiment(csv = TRUE) #Check out the format of this and save useful information in my own initializing document. Can avoid this step that way
      write.csv(Xotu, "ramclustR.csv", row.names = TRUE)
      RC <- ramclustR(ms="ramclustR.csv", ExpDes=experiment, st=0.25, sr=0.5, mspout = FALSE)
      
      #Finding exposure compound, figuring out which cluster it's in and removing everything that clusters with it
      dPPMHigh<-as.double(expMass)+(as.double(expMass)*(as.double(dPPM)/1000000))
      dPPMLow<-as.double(expMass)-(as.double(expMass)*(as.double(dPPM)/1000000))
      
      matches<-which(RC$fmz>dPPMLow & RC$fmz<dPPMHigh)
      clustMatches<-RC$featclus[matches]
      
      ####If list of metabolites is present, also remove suspected in-source fragments of metabolites####
      if(is.na(removeMets)==FALSE){
        expCompMet<-read.xlsx(fileName, sheet=((n*2)+1))
        expCompMet<-unique(expCompMet)
        
        for(k in 1:length(expCompMet[,1])){
          #Finding exposure compound, figuring out which cluster it's in and removing everything that clusters with it
          dPPMHigh<-as.double(expCompMet[k,1])+(as.double(expCompMet[k,1])*(as.double(dPPM)/1000000))
          dPPMLow<-as.double(expCompMet[k,1])-(as.double(expCompMet[k,1])*(as.double(dPPM)/1000000))
          
          matches<-c(matches,which(RC$fmz>dPPMLow & RC$fmz<dPPMHigh))
          clustMatches<-c(clustMatches,RC$featclus[matches])
          clustMatches<-unique(clustMatches)
        }
        
      }
      
      #Putting RTs back as a row instead of part as column-names
      RTs<-strsplit(colnames(Xotu),"_")[[1]][2]
      colnames(Xotu)[1]<-strsplit(colnames(Xotu),"_")[[1]][1]
      
      for(i in 2:length(Xotu[1,])){
        RTs<-c(RTs,strsplit(colnames(Xotu),"_")[[i]][2])
        colnames(Xotu)[i]<-strsplit(colnames(Xotu),"_")[[i]][1]
      }
      
      Xotu<-rbind(RTs,Xotu)
      #Xotu.4<-NA
      
      #Binding removedCorPeaks together and removing peaks that corelated to exposure compound or its' metabolites
      removedCorPeaks<-NA
      
      if(length(clustMatches)>0){
        for(i in 1:length(clustMatches)){
          clustMembers<-RC$fmz[which(RC$featclus==clustMatches[i])]
          if(identical(which(colnames(Xotu)%in%clustMembers),value)==TRUE){
            next()
          } else {
            removedCorPeaks<-cbind(removedCorPeaks, Xotu[,which(colnames(Xotu)%in%clustMembers)])
            Xotu<-Xotu[,-which(colnames(Xotu)%in%clustMembers)]
          }
        }
        
        removedCorPeaks<-as.data.frame(removedCorPeaks[,-1])
        cat("\n\nSummary of in-source fragment correlation-filter:\n",length(removedCorPeaks[1,]),"features were removed due to clustering with the exposure compound.\n\n")
        strC<-paste0("\n\nSummary of in-source fragment correlation-filter:\n",length(removedCorPeaks[1,]),"features were removed due to clustering with the exposure compound.\n\n")
      } else if(length(clustMatches)==0){
        Xotu<-Xotu
        cat("\n\nNo features were removed due to clustering of in-source fragments\n\n")
        strC<-paste0("\n\nNo features were removed due to clustering of in-source fragments\n\n")
      }
    }
    #Report pieces gathered: peaks clustered within RT and with high correlation
    #Xotu = Outgoing data frame
    
    #####Removing theoretical metabolites and the parental mass from featureset####
    if(is.na(removeMets)==FALSE){
      expCompMet<-read.xlsx(fileName, sheet=((n*2)+1))
      
      #if(is.na(removedCorPeaks)==FALSE){
      #  expCompMet[c(rep("CorPeak",length(colnames(removedCorPeaks)))),]<-c(as.double(colnames(removedCorPeaks)))
      #}
      
      
      dPPMactual<-vector()
      metPeaksRemoved<-data.frame()
      
      #preparing masses for comparisons
      
      whichToCorrect<-which(is.na(as.double(colnames(Xotu))))
      #colnames(Xotu)[whichToCorrect]<-gsub(".1$", "", colnames(Xotu)[whichToCorrect])
      #colnames(Xotu)[whichToCorrect]<-gsub(".2$", "", colnames(Xotu)[whichToCorrect])
      colnames(Xotu)<-substring(colnames(Xotu)[whichToCorrect],2)
      orbiPeaksMZ<-as.double(colnames(Xotu))
      
      #Checking every peak gathered through orbitrap-analysis against the list of potential metabolites of the exposure compound
      #Making a list over the masses, the dPPM between the orbiPeak and the expCompMet and noting the number of the orbiPeak to be removed
      for(k in 1:length(orbiPeaksMZ)){
        dPPMHigh<-orbiPeaksMZ[k]+(orbiPeaksMZ[k]*(as.double(dPPM)/1000000))
        dPPMLow<-orbiPeaksMZ[k]-(orbiPeaksMZ[k]*(as.double(dPPM)/1000000))
        
        matches<-subset(expCompMet, expCompMet[,1]>dPPMLow & expCompMet[,1]<dPPMHigh)
        
        if(length(matches[,1]>0)){
          for(i in 1:length(matches[,1])){
            #Calculating actual dPPM for all matches within dPPM range
            dPPMactual[i]<-((abs(orbiPeaksMZ[k]-matches[i,1])/orbiPeaksMZ[k])*1000000)
          }
          
          metPeaksRemoved<-as.data.frame(rbind(metPeaksRemoved, cbind(orbiPeaksMZ[k],min(dPPMactual), k)))
        }
      }
      
      #Removing the m/z's which have been matched to potential expCompMets and saving xlsx without them. 3rd column = position in Xotu.1.done
      if(dim(metPeaksRemoved)!=c(0,0)){
        removedPeaks<-Xotu[,(metPeaksRemoved[,3])]
        Xotu<-Xotu[,-(metPeaksRemoved[,3])]
        
        for(k in 1:length(metPeaksRemoved[,1])){
          colnames(removedPeaks)[k]<-paste(colnames(removedPeaks)[k],"_dPPM:", round(metPeaksRemoved[k,2],5))
        }
        
        cat("Summary of theoretical exposure-compound metabolite removal:\n", length(metPeaksRemoved[,1]), "features were removed due to being within ", as.double(dPPM), "ppm m/z of theoretical compound metabolites.\n\n")
        strD<-paste0("Summary of theoretical exposure-compound metabolite removal:\n", length(metPeaksRemoved[,1]), "features were removed due to being within", as.double(dPPM), "ppm m/z of theoretical compound metabolites.\n\n")
        
      } else {
        removedPeaks <- NA
        
        cat("No theoretical exposure compound metabolites were found or removed.\n\n")
        strD<-paste0("No theoretical exposure compound metabolites were found or removed.\n\n")
      }
      
      
      
      
      ####Removing within specified PPM and RT of the pseudo fragments
    } else {
      cat("User is not interested in using a metabolite-filter\n\n")
      strD<-paste0("User is not interested in using a metabolite-filter\n\n")
    }
    #Report pieces gathered features: peaks removed and their dPPMs
    #Xotu = Outgoing data frame
    
    #####Filter out features with negative values####
    #If use has set 'KeepNegVal' to TRUE, then nothing is done
    #Else every feature which conctains negative value will be removed
    #Could probably implement imputation-methods here
    
    if(is.na(KeepNegVal)==TRUE){
      cat("User has decided to keep features with negative values. No features filtered.\n\n")
      strE<-paste0("User has decided to keep features with negative values. No features filtered.\n\n")
    } else {
      Xotu[Xotu<0]<-NA
      Xotu.5<-NA
      negFeatures<-NA
      
      for(l in 1:length(Xotu[1,])){
        if(sum(is.na(Xotu[,l]))==0){
          Xotu.5<-cbind(Xotu.5,Xotu[,l])
          colnames(Xotu.5)[length(Xotu.5[1,])]<-colnames(Xotu)[l]
        } else {
          negFeatures<-as.data.frame(cbind(negFeatures, as.double(Xotu[,l])))
          colnames(negFeatures)[length(negFeatures[1,])]<-colnames(Xotu[l])
        }
      }
      
      if(is.null(dim(Xotu.5))==FALSE){
        Xotu.5<-Xotu.5[,-1]
        cat("Summary of negative value filter:\n")
        cat((100-(length(Xotu.5[1,])/length(Xotu[1,])*100)),"% of the features have been filtered due to containing negative values.\n\n")
        strE<-paste0((100-(length(Xotu.5[1,])/length(Xotu[1,])*100)),"% of the features have been filtered due to containing negative values.\n\n")
        negSummary<-(100-(length(Xotu.5[1,])/length(Xotu[1,])*100))
        Xotu<-Xotu.5
      } else {
        cat("None of the features have been filtered due to containing negative values.\n\n")
        strE<-paste0("None of the features have been filtered due to containing negative values.\n\n")
      }
      
      #negFeatures[is.na(negFeatures)]<-0
      #negFeatures<-as.double(negFeatures)
    }
    #Report pieces gathered: Removed features + negSummary
    #Outgoing: Xotu
    #####Blank filter#####
    if(is.na(blankFilter)==FALSE){
      blanks<-Xotu[which(grepl("Blank",samps)),]
      tempSamps<-samps[-which(grepl("Blank", samps))]
      blankSamps<-samps[which(grepl("Blank",samps))]
      rownames(blanks)<-blankSamps
      QCblanks<-blanks[which(grepl("QC",blankSamps)),]
      blanks<-blanks[-which(grepl("QC",blankSamps)),]
      Xotu.6<-Xotu[-which(grepl("Blank",samps)),]
      
      QCrat<-NA
      SampRat<-NA
      blankPeaks<-NA
      Xotu.keep<-NA
      
      for(bl in 1:length(Xotu.6[1,])){
        QCRat<-mean(as.double(QCblanks[,bl]))/max(as.double(Xotu.6[which(grepl("QC",tempSamps)),bl]))
        SampRat<-mean(as.double(blanks[,bl]))/max(as.double(Xotu.6[-c(which(grepl("QC",tempSamps)),1),bl]))
        if(QCRat>0.4 || SampRat>0.4){
          blankPeaks<-cbind(blankPeaks,Xotu.6[,bl])
          colnames(blankPeaks)[length(blankPeaks[1,])]<-colnames(Xotu.6)[bl]
        } else {
          Xotu.keep<-cbind(Xotu.keep,Xotu.6[,bl])
          colnames(Xotu.keep)[length(Xotu.keep[1,])]<-colnames(Xotu.6)[bl]
        }
      }
      
      
      
      if(is.null(dim(Xotu.keep))==FALSE){
        cat("Blank filter summary:\nA total of ", (length(blankPeaks[1,])/length(Xotu[1,])*100), "% of peaks removed due to presence in blanks.")
        strF<-paste0("Blank filter summary:\nA total of ", (length(blankPeaks[1,])/length(Xotu[1,])*100), "% of peaks removed due to presence in blanks.")
        Xotu<-Xotu.keep[,-1]
        blankPeaks<-blankPeaks[,-1]
        blankPeaks<-rbind(blankPeaks,blanks[,which(colnames(blankPeaks)%in%colnames(blanks))],QCblanks[,which(colnames(blankPeaks)%in%colnames(blanks))])
      } else {
        cat("Blank filter summary:\nNo peaks removed due to presence in blanks.")
        strF<-paste0("Blank filter summary:\nNo peaks removed due to presence in blanks.")
      }
    
    } else {
      cat("User not interested in applying blank-filter to data")
      strF<-paste0("User not interested in applying blank-filter to data")
    }
    
    
    #####Removing unwanted samples from the dataset#####
    
    samps.preFilter<-samps
    if(is.na(sampleFilter)==FALSE){
      for(i in 1:length(sampleFilter)){
        Xotu<-Xotu[-which(grepl(sampleFilter[i], samps)),]
        samps<-samps[-which(grepl(sampleFilter[i], samps))]
      }
    }
    
    #Removing blank samples from sample name list
    if(is.na(blankFilter)==FALSE){
      samps<-samps[-which(grepl("Blank",samps))]
    }
    
    #####Storing the current dataset in the list completeDataset####
    
    completeDataset[[n]]<-Xotu
    
    #Compiling and storing removed features
    
    if(is.na(ThresholdFilter)==FALSE){
      thresholdFeatures<-rbind(rep("threshold", times=length(thresholdFeatures[1,])),thresholdFeatures)
      thresholdFeatures<-thresholdFeatures[,-1]
    } else {
      thresholdFeatures<-NA
    }
    
    if(is.na(KeepNegVal)==FALSE && negSummary!=0){
      negFeatures<-rbind(rep(as.character("neg_value"), times=length(negFeatures[1,])), negFeatures)
      negFeatures<-negFeatures[,-1]
    } else {
      negFeatures <- NA
    }
    
    if(is.na(removeMets)==FALSE && is.null(removedPeaks)==FALSE){
      removedPeaks<-rbind(rep("metPeaks", times=length(removedPeaks[1,])), removedPeaks)
      removedPeaks<-removedPeaks[,-1]
    } else {
      removedPeaks <- NA
    }
    
    if(is.na(blankFilter)==FALSE && is.null(blankPeaks)==FALSE){
      blankPeaks<-rbind(rep("blankPeaks", times=length(blankPeaks[1,])), blankPeaks)
      blankPeaks<-blankPeaks[,-1]
    } else {
      blankPeaks <- NA
    }
    
    #Testing if corpeaks is just NA or not
    if(is.na(InSourceFragRemove)==FALSE){
      if(sum(is.na(removedCorPeaks))==length(removedCorPeaks)){
        removedCorPeaks<-rbind(rep("corPeaks", times=1), removedCorPeaks)
        for(i in 2:length(otu[,1])){
          removedCorPeaks<-rbind(removedCorPeaks,NA)
        }
        removedCorPeaks<-as.character(removedCorPeaks)
      } else {
        removedCorPeaks<-rbind(rep("corPeaks", times=length(removedCorPeaks[1,])), removedCorPeaks)
        removedCorPeaks<-removedCorPeaks[,-1]
      }
    } else {
      removedCorPeaks<-NA
    }
    
    
    allRemovedFeatures[[n]]<-cbind(thresholdFeatures, negFeatures, removedPeaks, removedCorPeaks, blankPeaks)
    rownames(allRemovedFeatures[[n]])[2:length(allRemovedFeatures[[n]][,1])]<-samps.preFilter
    
    if(is.na(sampleFilter)==FALSE){
      for(i in 1:length(sampleFilter)){
        allRemovedFeatures[[n]]<-allRemovedFeatures[[n]][-which(grepl(sampleFilter[i], rownames(allRemovedFeatures[[n]]))),]
      }
      #Removing blank samples
      if(is.na(blankFilter)==FALSE){
        allRemovedFeatures[[n]]<-allRemovedFeatures[[n]][-which(grepl("Blank",rownames(allRemovedFeatures[[n]]))),]
      }
    }
    
    
    ###Printing report for the dataset
    if(report==TRUE){
      pdf(paste0(colnames(metaInfo)[n+1],' Report.pdf'), width=12, height=14)
      plot(NA, xlim=c(0,8), ylim=c(0,8), bty='n',
           xaxt='n', yaxt='n', xlab='', ylab='')
      text(1,7,paste0("Report for ",colnames(metaInfo)[n+1]," dataset"))
      text(1,6,strA, pos=4)
      text(1,5,strB, pos=4)
      text(1,4,strC, pos=4)
      text(1,3,strD, pos=4)
      text(1,2,strE, pos=4)
      text(1,1,strF, pos=4)
      dev.off()
      
      #addWorksheet(finalDocument, sheetName="Report") #Consider making this into a pdf for professionalism and win
      #writeData(finalDocument, sheet=(length(finalDocument$sheet_names)), reportData)
    }
    
    cat("\n\n\n\n")
  } #Dataset filter loop ending asdfsadfsadfsaddfafasdfasdfasdfasdfasdfasdfasdfsadfasdfsadfsadfasdfasdfasdfdsafasdfas
  
  colnames(Xotu)
  
  #####Creating and saving report + output from filtering and corrections
  finalDocument <- createWorkbook()
  
  #Avoiding colnames turning into X-versions of themselves
    oldNames<-colnames(Xotu)
    colnames(Xotu)<-substring(colnames(Xotu), first=2)
  
  ####Printing to document####
  if(is.na(MergeAllSets)==FALSE){
    
    ##Merging the final datasets
    mergedDatasets<-NA
    addWorksheet(finalDocument, sheetName="CorrectedFeatures")
    
    for(i in 1:nDataset){
      completeDataset[[i]]<-as.data.frame(rbind(rep(as.character(datasetNames[i]),length(completeDataset[[i]][1,])),completeDataset[[i]]), stringsAsFactors=FALSE)
      mergedDatasets<-as.data.frame(cbind(mergedDatasets,completeDataset[[i]]))
    }
    
    mergedDatasets<-cbind(rbind(NA,as.data.frame(as.character(samps))),mergedDatasets)
    mergedDatasets<-mergedDatasets[,-2]
    
    writeData(finalDocument, sheet=1, x=mergedDatasets)
    
    ##Mergin the removed features
    mergedRemovedFeatures<-NA
    addWorksheet(finalDocument, sheetName="RemovedFeatures")
    
    for(i in 1:nDataset){
      allRemovedFeatures[[i]]<-as.data.frame(rbind(rep(datasetNames[i],length(allRemovedFeatures[[i]][1,])),allRemovedFeatures[[i]]), stringsAsFactors=FALSE)
      mergedRemovedFeatures<-as.data.frame(cbind(mergedRemovedFeatures,allRemovedFeatures[[i]]), stringsAsFactors=FALSE)
    }
    
    mergedRemovedFeatures<-cbind(rbind(NA, NA,as.data.frame(as.character(samps))),mergedRemovedFeatures)
    mergedRemovedFeatures<-mergedRemovedFeatures[,-2]
    
    writeData(finalDocument, sheet=2, x=mergedRemovedFeatures)
    
  } else {
    for(i in 1:nDataset){
      
      ##Adding one worksheet / dataset
      addWorksheet(finalDocument, sheetName=paste("CorrectedFeatures_", i))
      completeDataset[[i]]<-cbind(rbind(as.data.frame(as.character(samps))),completeDataset[[i]])
      writeData(finalDocument, sheet=i, completeDataset[[i]])
      
      ##Adding one worksheet / removed featureset
      addWorksheet(finalDocument, sheetName=paste("RemovedFeatures_", i))
      allRemovedFeatures[[i]]<-cbind(rbind(NA,as.data.frame(as.character(samps))),allRemovedFeatures[[i]])
      writeData(finalDocument, sheet=(i+1), allRemovedFeatures[[i]])
    }
  }
  
  saveWorkbook(finalDocument,"ExpMISFFilter_Results.xlsx", overwrite = TRUE)
  #threshSummary, negSummary. Report-data to be printed
}#Function ending
