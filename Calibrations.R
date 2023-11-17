## Calibrations

#Russell et al., 2018
AfricaRussellMAAT<-alist(MBT5=,m=32.42,b=-1.21,(MBT5*m)+b)
AfricaRussellpH<-alist(CBTp=,m=2.65,b=8.95,(CBTp*m)+b)

#Martinez-Sosa et al., 2022
GlobalMartinezMAF<-alist(MBT5=,m=1/0.03,b=0.075/0.03,(MBT5*m)+b)

#DeJonge et al., 2014
GlobalDeJongeMAAT<-alist(MBT5=,m=31.45,b=-8.57,(MBT5*m)+b)
GlobalDeJongepH<-alist(CBTp=,m=1.59,b=7.15,(CBTp*m)+b)

Temperature<-list(Russell=AfricaRussellMAAT,Martinez=GlobalMartinezMAF,DeJonge=GlobalDeJongeMAAT)

pH<-list(Russell=AfricaRussellpH,DeJonge=GlobalDeJongepH)

Calibrations<-list(Temperature=Temperature,pH=pH)

## Aplication

linearCalib<-function(raw,env,calibration,data="area",out="all",na.ignore=FALSE,coerce=FALSE,complete=TRUE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore=na.ignore,coerce=coerce)
  } else if (data == "fa"){fa<-raw}
  if(env=="Temperature"){prediction<-as.function(Calibrations[[env]][[calibration]])(MBT5(fa,data="fa",complete = complete,na.ignore=na.ignore))
  } else if(env=="pH"){prediction<-as.function(Calibrations[[env]][[calibration]])(CBTp(fa,data="fa",complete = complete,na.ignore=na.ignore))
  } else stop("Select one of the available environmental parameters: Temperature or pH")
  prediction
}