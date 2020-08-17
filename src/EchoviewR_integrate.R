

# Choose your ship: 'MS', 'KPH', 'KJH', or 'CDH'
vessel<-'CDH'

###Path to rawdata
if (vessel == 'CDH') {
	rawdata<-'D:\\Working\\KRILL\\KRILL2019\\cruise_folders\\CDH2019828\\ACOUSTIC\\EK80\\EK80_RAW_REDUCED'
} else if (vessel == 'KPH') {
    rawdata<-'D:\\Working\\KRILL\\KRILL2019\\cruise_folders\\S2019701_PKRONPRINSHAAKON_9566\\ACOUSTIC\\EK80\\EK80_ORIGINALRAWDATA\\EK80DK'
} else if (vessel == 'MS') {
    rawdata<-'D:\\Working\\KRILL\\KRILL2019\\cruise_folders\\MS2019999\\ACOUSTIC\\ES80\\ES80_RAWDATA_REDUCED'
} else if (vessel == 'KJH') {
    rawdata<-'D:\\Working\\KRILL\\KRILL2019\\cruise_folders\\KJH2019999\\ACOUSTIC\\EK60\\EK60_ORIGINALRAWDATA'
}

wd<-paste0('D:\\Working\\KRILL\\KRILL2019\\data\\echo-integration-2020\\', vessel, '\\')

###Read reference file (Start and stop of transects)
refFile<-read.csv(paste0(wd,'transects.csv'),sep=',', header=TRUE)

###List .raw files
raw<-list.files(rawdata,pattern='.raw$',full.names=TRUE)

###Select the .raw-files associated with a transect
START<-grep('On',as.character(refFile$Event))
STOP<-grep('Off',as.character(refFile$Event))

for (t in 1:length(START)) {
  # Turn time into something suitable for a filename
  time = sub("\\..*", "", gsub(":", "", refFile$Time[START[t]]))
  # Make up the date/time in a form suitable for a filename
  timestamp<-paste0(gsub("-", "", refFile$Date[START[t]]), 'T', time)
  
  print(paste0('Processing event: ', refFile$Event[START[t]], ' (', timestamp, ')'))
  library(stringr)
  transectId<-word(refFile$Event[START[t]], -1)
  
  STARTFILE<-grep(refFile$Associated_rawfile[START[t]], raw)
  STOPFILE<-grep(refFile$Associated_rawfile[STOP[t]], raw)

  rawFiles<-raw[STARTFILE:STOPFILE]

  ###Read calibration file
  calFile<- paste0(wd,vessel,'_calibration.ecs')

  ###Read EV-template specific for Cabo
  EVtemplate<-paste0(wd,'CCAMLR_SWARM_120kHz_only_', vessel, '.EV')
  
  ###Read region definitions if present (for us, these come from LSSS).
  evrFile<-paste0(wd, vessel, '.evr')
  
  library(EchoviewR,quietly=TRUE)
  
  EVAppObj=COMCreate('EchoviewCom.EvApplication') 
  
  outputEVFile<-paste0(wd, timestamp, '_', vessel, '.ev')
  
  EVFile=EVCreateNew(EVAppObj=EVAppObj,
                     templateFn=EVtemplate,
                     EVFileName=outputEVFile,
                     filesetName="fisheries",
                     dataFiles=rawFiles,
                     CloseOnSave = FALSE)$EVFile
  
  EVAddCalibrationFile(EVFile=EVFile, filesetName='fisheries', calibrationFile=calFile)
  
  if (file.exists(evrFile))
    EVImportRegionDef(EVFile=EVFile, evrFile=evrFile)
  
  EVSaveFile(EVFile=EVFile)
  
  swarmDetResults=EVSchoolsDetect(EVFile = EVFile,
                                  acoVarName='120 Dilation filter 3x3',
                                  outputRegionClassName = 'aggregation',
                                  deleteExistingRegions = TRUE,
                                  distanceMode = "GPS distance",
                                  maximumHorizontalLink = 15,
                                  maximumVerticalLink = 5,
                                  minimumCandidateHeight = 3, 
                                  minimumCandidateLength = 15,
                                  minimumSchoolHeight = 3,
                                  minimumSchoolLength = 15,
                                  dataThreshold = -70)
                                  
  exportFileName=paste(wd,'krillNASC_', transectId, '_', timestamp, '_', vessel,'.csv',sep='') 
  
  EVExportUnderlying(EVFile=EVFile, 
                     variableName="Krill NASC from mean Sv (export here for NASC values)", 
                     pingRange = c(-1, -1), 
                     filePath=exportFileName)
  
  EVAppObj$Quit()
}

