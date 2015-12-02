## Batch align

# Extract peakinfo (i.e. m/z and rt of features)
peakIn=peakInfo(PTnofill)
#
bF=batchFlag(PTnofill,meta,peakIn)
aIQ=alignIndex(bF,grpType='Q',mzdiff=0.0025,rtdiff = 15,report=T,reportName='../splits_aIQ')
aIR=alignIndex(bF,grpType='R',mzdiff=0.0025,rtdiff = 15,report=T,reportName='../splits_aIR')
plotAlign(bF,aIQ,plotType='pdf',reportName='../clustPlots_aIQ')
plotAlign(bF,aIR,plotType='pdf',reportName='../clustPlots_aIR')
aQR=aggregateIndex(aIQ,aIR)
bA=batchAlign(bF,aQR,PTfill,meta)

PT=bA$PTalign
