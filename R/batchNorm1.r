setwd("C:/R/QCData")                    # Specify work directory
rm(list=ls())

library(xcms)

peakInfo=function(PT) {
  peakInfo=matrix(unlist(strsplit(colnames(PT),'@')),ncol=2,byrow=TRUE)
  peakInfo[,1]=substr(peakInfo[,1],3,max(nchar(peakInfo[,1])))
  peakInfo=matrix(as.numeric(peakInfo),ncol=2)
  colnames(peakInfo)=c('mz','rt')
  rownames(peakInfo)=paste('feature',1:nrow(peakInfo),sep='_')
  return(peakInfo)
}

CV=function(mat) {
  mean=apply(mat,2,function(x) mean(x,na.rm=TRUE))
  sd=apply(mat,2,function(x) sd(x,na.rm=TRUE))
  cv=sd/mean
}

refOut=function(PT,meta,grpType='R',CVlimit=0.3) {
  batch=meta[,1]
  grp=meta[,2]
  uniqBatch=unique(batch)
  CVMat=matrix(nrow=length(uniqBatch),ncol=ncol(PT))
  rownames(CVMat)=uniqBatch
  aveIntMat=matrix(nrow=length(uniqBatch),ncol=ncol(PT))
  rownames(aveIntMat)=uniqBatch
  for (b in 1:length(uniqBatch)) {
    bat=uniqBatch[b]
    PTbatch=PT[batch==bat & grp==grpType,]
    CVMat[b,]=ifelse(CV(PTbatch)<=CVlimit,TRUE,FALSE)
    aveIntMat[b,]=apply(PTbatch,2,mean)
  }
  return(list(CV=CVMat,aveInt=aveIntMat))
}

refCorr=function(PT,meta,refs,FCLimit=5) {
  batch=meta[,1]
  PTcorr=PT
  cvm=refs$CV  # Take out CV flags
  nBatch=nrow(cvm)
  uniqBatch=rownames(cvm)
  nFeat=ncol(cvm)
  refCorrMat=matrix(FALSE,nrow=nBatch,ncol=nFeat)  # Allocate matrix for "features normalised by ref samples"
  cvFlags=apply(cvm,2,function(x) sum(x,na.rm=TRUE))
  aveInt=refs$aveInt
  meanIntRat=matrix(1,nrow=nBatch,ncol=nBatch) # Calcluate average feature intensity ratios across batches
  for (b in 1:(nBatch-1)) {
    for (bb in (b+1):nBatch) {
      meanIntRat[bb,b]=mean(aveInt[bb,])/mean(aveInt[b,])
      meanIntRat[b,bb]=1/meanIntRat[bb,b]
    }
  }
  # Find features candidates where Ref normalisation may be possible (where CV<CVLimit)
  whichFeatsCV=which(cvFlags>1)
  lenCV=length(whichFeatsCV)
  # Step through candidates
  for (lc in 1:lenCV) {
    feat=whichFeatsCV[lc]
    featIntRat=matrix(1,nrow=nBatch,ncol=nBatch) # Calcluate specific feature intensity ratios across batches
    for (b in 1:(nBatch-1)) {
      for (bb in (b+1):nBatch) {
        featIntRat[bb,b]=aveInt[bb,feat]/aveInt[b,feat]
        featIntRat[b,bb]=1/featIntRat[bb,b]
      }
    }
    featFlags=abs(log(featIntRat/meanIntRat))<=log(FCLimit) # Find candidates for ref norm acc to FC<FCLimit
    for (b in 1:nBatch) {
      if (cvm[b,feat]==FALSE | is.na(cvm[b,feat])) { # Combine CV and FC criteria
        featFlags[,b]=featFlags[b,]=FALSE
      }
    }
    if (any(colSums(featFlags)>1)) {
      ## Find reference batch
      refBatch=min(which(colSums(featFlags)==max(colSums(featFlags))))
      refCorr=featFlags[,refBatch]
      WhichRefCorr=which(refCorr==TRUE)
      refInt=aveInt[refBatch,feat]
      for (b in which(refCorr==TRUE)[-1]) {
        corrFact=refInt/aveInt[b,feat]
        PTcorr[batch==uniqBatch[b],feat]=PT[batch==uniqBatch[b],feat]*corrFact
      }
      refCorrMat[,feat]=refCorr 
    }
  }
  return(list(PTRef=PTcorr,refCorr=refCorrMat,PTOrg=PT))
}

featIntRat=aveInt[,c]/aveInt[1,c]


PT=bAQR$PTalign
peakInfo=peakInfo(PT)
refs=refOut(PT,meta)
RC=refCorr(PT,meta,refs,FCLimit=5)
