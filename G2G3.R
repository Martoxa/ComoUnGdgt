G2G3<-function(raw,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){raw[is.na(raw)]<-0}
  formula<-"raw$GDGT2 / raw$GDGT3"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(raw))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(raw,formula)}
  out<-eval(parse(text = formula))
  out[out>10000]<-NA
  out[is.nan(out)]<-NA
  out
}