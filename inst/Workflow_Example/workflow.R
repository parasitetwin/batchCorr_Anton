#######################
## Read in raw data

library(xcms)

filepath="C:/R/QCData/OrgData"  # Specify directory of centroid files

## Read samples
QC1=xcmsSet(filepath,method='centWave',prefilter=c(3,440),peakwidth=c(5,76),snthresh=6,mzdiff=0.0045,ppm=15)
QC2=group(QC1,bw=10,minfrac=0.75,minsamp=1,mzwid=0.015,sleep=0)
QC3=retcor(QC2, family="s", span=0.2)
QC_nofill=group(QC3,bw=1,mzwid=0.015,minfrac=.75)
QC_fill=fillPeaks(QC_nofill,method='chrom')

## Organise into peaktable with missing data
QCB=grabAlign(QC_nofill,batch='Batch_B',grp='QB')
RefB=grabAlign(QC_nofill,batch='Batch_B',grp='Ref')
QCF=grabAlign(QC_nofill,batch='Batch_F',grp='QF')
RefF=grabAlign(QC_nofill,batch='Batch_F',grp='Ref')
QCH=grabAlign(QC_nofill,batch='Batch_H',grp='QH')
RefH=grabAlign(QC_nofill,batch='Batch_H',grp='Ref')
PTnofill=rbind(QCB,RefB,QCF,RefF,QCH,RefH)
# dim(PTnofill) # 90 x 11667

## Organise into peaktable without missing data
QCB=grabAlign(QC_fill,batch='Batch_B',grp='QB')
RefB=grabAlign(QC_fill,batch='Batch_B',grp='Ref')
QCF=grabAlign(QC_fill,batch='Batch_F',grp='QF')
RefF=grabAlign(QC_fill,batch='Batch_F',grp='Ref')
QCH=grabAlign(QC_fill,batch='Batch_H',grp='QH')
RefH=grabAlign(QC_fill,batch='Batch_H',grp='Ref')
PTfill=rbind(QCB,RefB,QCF,RefF,QCH,RefH)
# dim(PTfill) # 90 x 11667

## Set up metadata (Quick'n'Dirty approach)
batch=c(rep('B',nrow(QCB)+nrow(RefB)),rep('F',nrow(QCF)+nrow(RefF)),rep('H',nrow(QCH)+nrow(RefH)))
grp=c(rep('Q',nrow(QCB)),rep('R',nrow(RefB)),rep('Q',nrow(QCF)),rep('R',nrow(RefF)),rep('Q',nrow(QCH)),rep('R',nrow(RefH)))
meta=cbind(batch,grp)
# dim(meta) # 90 x 2

## Save relevant data files
# save(PTfill,file='PT_fill.RData')
# save(PTnofill,file='PT_nofill.RData')
# save(meta,file='meta.RData')


##########################
## Perform batch alignment

# Extract peakinfo (i.e. m/z and rt of features)
peakIn=peakInfo(PTnofill)
# Flag presence/missingness on batch level
bF=batchFlag(PTnofill,meta,peakIn)
# Find possible alignment candidates per sample type
aIQ=alignIndex(bF,grpType='Q',mzdiff=0.002,rtdiff=15,report=T,reportName='../splits_aIQ')
aIR=alignIndex(bF,grpType='R',mzdiff=0.002,rtdiff=15,report=T,reportName='../splits_aIR')
# Plot achieved alignments
plotAlign(bF,aIQ,plotType='pdf',reportName='../clustPlots_aIQ')
plotAlign(bF,aIR,plotType='pdf',reportName='../clustPlots_aIR')
# Aggregate alignments across sample types
aQR=aggregateIndex(aIQ,aIR)
# Perform alignment -> Peaktable
bA=batchAlign(bF,aQR,PTfill,meta)
# Extract new peak table
PT=bA$PTalign


#########################################
## Perform within batch drift corrections


################################################
## Perform between batch intensity normalisation

# Extract aligned peak table
PT=bAQR$PTalign
peakInfo=peakInfo(PT)
# Aggregate ref sample info on batch level
refs=refOut(PT,meta)
# Total number of batch features
allBFeats=prod(dim(RC$refCorr))
# Check different fold change limit-settings
fcs=seq(1.2,10,by=0.2)
propCorr=numeric(length(fcs))
c=0
for (fc in fcs) {
  c=c+1
  cat(c)
  RC=refCorr(PT,meta,refs,FCLimit=fc)
  # Find proportion of batch features normalized by reference samples under FClimit constraint.
  propCorr[c]=sum(RC$refCorr)/allBFeats
}
# Plot supplementary figure for FC limit determination
plot(fcs,propCorr,xlab='Fold change limit', ylab='Proportion normalized by reference samples')
lines(fcs,propCorr)
# Perform ref normalization using FClimit=5
RC=refCorr(PT,meta,refs,FCLimit=fc)
# Extract reference normalised peak table
PT_RefNorm=RC$PTRef

# NB: Final normalization is not complete!!!
# Batches where Ref heuristic criterion is not passed still need to be population normalized...
# Info on which batches need to be population normalized is found in RC$refCorr
