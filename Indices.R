MBT5<-function(raw){
  (raw$Ia+raw$Ib+raw$Ic)/(raw$Ia+raw$Ib+raw$Ic+raw$IIa5+raw$IIb5+raw$IIc5+raw$IIIa5)
} #checked Baxter et al., 2023

MBT6<-function(raw){
  (raw$Ia+raw$Ib+raw$Ic)/(raw$Ia+raw$Ib+raw$Ic+raw$IIa6+raw$IIb6+raw$IIIa6)
}

fC<-function(raw){
  ((raw$IIb5+raw$Ib+raw$IIb6)+2*(raw$Ic))/((raw$IIIa5+raw$IIIa6+raw$IIa5+raw$IIa6+raw$Ia)+(raw$IIb5+raw$Ib+raw$IIb6)+(raw$Ic))*0.5
}
CBTp<-function(raw){
  log10((raw$Ic+raw$IIa6+raw$IIb6+raw$IIc6+raw$IIIa6+raw$IIIb6+raw$IIIc6)/(raw$Ia+raw$IIa5+raw$IIIa5))
}

IR<-function(raw){
  (raw$IIa6+raw$IIb6+raw$IIIa6)/(raw$IIa5+raw$IIb5+raw$IIIa5+raw$IIa6+raw$IIb6+raw$IIIa6)
}
IBT<-function(raw){
  -log10((raw$IIa6+raw$IIIa6)/(raw$IIa5+raw$IIIa5))
}
CI<-function(raw){
  (raw$Ia+raw$Ib)/(raw$Ia+raw$IIa5+raw$IIIa5+raw$IIIa6)
}
TEX86<-function(raw){
  (raw$GDGT2+raw$GDGT3+raw$Crenp)/(raw$GDGT1+raw$GDGT2+raw$GDGT3+raw$Crenp)
}
BIT<-function(raw,na.ignore=FALSE){
  if(na.ignore==TRUE){raw[is.na(raw)]<-0}
  (raw$Ia+raw$IIa5+raw$IIa6+raw$IIIa5+raw$IIIa6)/(raw$Ia+raw$IIa5+raw$IIa6+raw$IIIa5+raw$IIIa6+raw$Crena)
}  #checked Baxter et al., 2023

MI<-function(raw){
  (raw$GDGT1 + raw$GDGT2 + raw$GDGT3)/(raw$GDGT1 + raw$GDGT2 + raw$GDGT3 + raw$Crena + raw$Crenp)
}
RI<-function(raw){
  (raw$Ib+2*(raw$Ic))/(raw$Ia+raw$Ib+raw$Ic)
}

pGDGT0<-function(raw){
  raw$GDGT0/(raw$GDGT0+raw$Crena)
} #Checked Zander, et al. in prep