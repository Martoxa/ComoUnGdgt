source("Data_manipulation.R")
source("Control_checks.R")
source("Indices.R")
source("Calibrations.R")
source("bayRmbt_model.R")
source("bayRmbt_forward.R")
source("bayRmbt_predict.R")
source("bayspR_tex.R")
source("bayspR_tex_analog.R")
source("RthChordDistances.R")

## Data Manipulation

dataArea<-read.csv(file("test_ds_missing.csv"),head=TRUE,sep=',') #Load a dataset with peak areas and missing IIIc and IIIc'

BrFractionalAbundances<-FracA(dataArea,out = "br",coerce = TRUE) #Calculate the fractional abundances of only brGDGTs. coerce=TRUE will add missing compounds as 0s, this is FALSE by default
IsoFractionalAbundances<-FracA(dataArea,out = "iso") #Calculate the FA of only isoGDGTs. Other options include "all", which is default, which outputs the FA of br and iso, calculated separately, and "both" which calculates the FA based on the sum of ALL GDGTs.

## Indices

BrMBT5<-MBT5(dataArea) #Can calculate MBT'5Me directly from the peak area data. It does however calculate FA beforehand.
IsoTEX<-TEX86(IsoFractionalAbundances,data = "fa") #Or the indices can be calculated from FA data.

MBT5(IsoFractionalAbundances,data = "fa") #Function cannot run by default if at least one of the GDGTs from the original equation is missing
CBTp(IsoFractionalAbundances,data = "fa",complete = FALSE) # complete=FALSE will allow the calculation to run with missing compounds, and will show which ones were excluded

BrFAMissing<-BrFractionalAbundances[,1:14] #Remove Ic from the dataset
CBTp(BrFAMissing,data = "fa") #Again, CBT' won't be calculated unless all compounds are present by default
CBTp(BrFAMissing,data = "fa",complete=FALSE) #complete=FALSE will calculate it with missing data

## Calibrations brGDGTs

#Linear calibrations

RussellMAAT<-linearCalib(dataArea,env="Temperature",calibration = "Russell") #Environmental calibrations can be generated directly from peak areas. Published calibrations ara available and can be access specifying the environmental parameter (T, pH), and the calibration.
DeJongepH<-linearCalib(BrFAMissing,env="pH",calibration = "DeJonge",data = "fa",complete = FALSE) #Calibrations can be generated from FA too. And just like before, missing data can be circumvented with complete=FALSE

#Other calibrations

MartinezMAF<-bayrmbt_predict(dataArea,10,10,Tmodel = "T0",Type="lake") #Other types of calibrations can be adapted too. BayMBT now runs directly from peak areas.

## Calibrations isoGDGTs

#BAYSPAR
Marine<-dataArea[dataArea$Type=="M",]
Marine<-Marine[5:9,] #Just a handful samples from a similar area
TierneyBayspar<-bayspR_tex(Marine,mean(Marine$Longitude),mean(Marine$Latitude),mean(Marine$MAAT),"SST")