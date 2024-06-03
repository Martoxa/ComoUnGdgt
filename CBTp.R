CBTp<-function(raw,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){raw[is.na(raw)]<-0}
  formula<-"log10(( raw$Ic + raw$IIa6 + raw$IIb6 + raw$IIc6 + raw$IIIa6 + raw$IIIb6 + raw$IIIc6 )/( raw$Ia + raw$IIa5 + raw$IIIa5 ))"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(raw))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(raw,formula)}
  out<-eval(parse(text = formula))
  out[-10000 > out]<-NA
  out[out>10000]<-NA
  out
}