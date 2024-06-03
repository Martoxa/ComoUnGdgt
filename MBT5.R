MBT5<-function(raw,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){raw[is.na(raw)]<-0}
  formula<-"( raw$Ia + raw$Ib + raw$Ic )/( raw$Ia + raw$Ib + raw$Ic + raw$IIa5 + raw$IIb5 + raw$IIc5 + raw$IIIa5 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(raw))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(raw,formula)}
  eval(parse(text = formula))
} #checked Baxter et al., 2023