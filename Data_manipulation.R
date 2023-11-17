FracA<-function(raw,out="all",coerce=FALSE,na.ignore=FALSE){
  GDGTn<-c("GDGT0","GDGT1","GDGT2","GDGT3","Crena","Crenp","IIIa5","IIIa6","IIIb5","IIIb6","IIIc5","IIIc6","IIa5","IIa6","IIb5","IIb6","IIc5","IIc6","Ia","Ib","Ic")
  if(na.ignore==TRUE){raw[is.na(raw)]<-0}
  if(coerce==TRUE){
    noGDGT<-GDGTn[is.na(match(GDGTn,colnames(raw)))]
    nullGDGT<-data.frame(matrix(0,ncol = length(noGDGT),nrow=dim(raw)[1]))
    colnames(nullGDGT)<-noGDGT
    raw<-cbind(raw,nullGDGT)
    message(paste("Coerced ",paste(noGDGT,collapse = " ")))
  }
  if(out=="both"){
    out<-raw[,na.omit(match(GDGTn,colnames(raw)))]/apply(raw[,na.omit(match(GDGTn,colnames(raw)))],1,sum)
    out[is.na(out)]<-0
    out
  } else if (out == "iso"){
    out<-raw[,na.omit(match(GDGTn[1:6],colnames(raw)))]/apply(raw[,na.omit(match(GDGTn[1:6],colnames(raw)))],1,sum)
    out[is.na(out)]<-0
    out
  } else if (out == "br"){
    out<-raw[,na.omit(match(GDGTn[7:21],colnames(raw)))]/apply(raw[,na.omit(match(GDGTn[7:21],colnames(raw)))],1,sum)
    out[is.na(out)]<-0
    out
  } else if (out == "all") {
    out<-cbind(raw[,na.omit(match(GDGTn[1:6],colnames(raw)))]/apply(raw[,na.omit(match(GDGTn[1:6],colnames(raw)))],1,sum),raw[,na.omit(match(GDGTn[7:21],colnames(raw)))]/apply(raw[,na.omit(match(GDGTn[7:21],colnames(raw)))],1,sum))
    out[is.na(out)]<-0
    out
  }

}