setwd("C:/R/QCData")  # Specify work directory
rm(list=ls())         # Clear memory

library(xcms)

## Read samples
QC1=xcmsSet('OrgData',method='centWave',prefilter=c(3,440),peakwidth=c(5,76),snthresh=6,mzdiff=0.0045,ppm=15)
QC2=group(QC1,bw=10,minfrac=0.75,minsamp=1,mzwid=0.015,sleep=0)
QC3=retcor(QC2, family="s", span=0.2)
QC_nofill=group(QC3,bw=1.5,mzwid=0.015,minfrac=.75)
QC_fill=fillPeaks(QC_nofill,method='chrom')

QCB=grabAlign(QC_nofill,batch='Batch_B',grp='QB')
RefB=grabAlign(QC_nofill,batch='Batch_B',grp='Ref')
QCF=grabAlign(QC_nofill,batch='Batch_F',grp='QF')
RefF=grabAlign(QC_nofill,batch='Batch_F',grp='Ref')
QCH=grabAlign(QC_nofill,batch='Batch_H',grp='QH')
RefH=grabAlign(QC_nofill,batch='Batch_H',grp='Ref')
PTnofill=rbind(QCB,RefB,QCF,RefF,QCH,RefH)
# dim(PTnofill) # 90 x 11667

QCB=grabAlign(QC_fill,batch='Batch_B',grp='QB')
RefB=grabAlign(QC_fill,batch='Batch_B',grp='Ref')
QCF=grabAlign(QC_fill,batch='Batch_F',grp='QF')
RefF=grabAlign(QC_fill,batch='Batch_F',grp='Ref')
QCH=grabAlign(QC_fill,batch='Batch_H',grp='QH')
RefH=grabAlign(QC_fill,batch='Batch_H',grp='Ref')
PTfill=rbind(QCB,RefB,QCF,RefF,QCH,RefH)
# dim(PTfill) # 90 x 11667

batch=c(rep('B',nrow(QCB)+nrow(RefB)),rep('F',nrow(QCF)+nrow(RefF)),rep('H',nrow(QCH)+nrow(RefH)))
grp=c(rep('Q',nrow(QCB)),rep('R',nrow(RefB)),rep('Q',nrow(QCF)),rep('R',nrow(RefF)),rep('Q',nrow(QCH)),rep('R',nrow(RefH)))
meta=cbind(batch,grp)
# dim(meta) # 90 x 2
# save(PTfill,file='PT_fill.RData')
# save(PTnofill,file='PT_nofill.RData')
# save(meta,file='meta.RData')
