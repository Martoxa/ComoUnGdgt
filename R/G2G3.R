G2G3<-function(gdgt,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"gdgt$GDGT2 / gdgt$GDGT3"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  out<-eval(parse(text = formula))
  out[out>10000]<-NA
  out[is.nan(out)]<-NA
  out
}
