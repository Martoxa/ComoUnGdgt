BITx<-function(raw,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){raw[is.na(raw)]<-0}
  formula<-"( raw$Ia + raw$IIa5 + raw$IIa6 + raw$IIIa5 + raw$IIIa6 )/( raw$Ia + raw$IIa5 + raw$IIa6 + raw$IIIa5 + raw$IIIa6 + raw$Cren )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(raw))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(raw,formula)}
  eval(parse(text = formula))
}  #checked Baxter et al., 2023