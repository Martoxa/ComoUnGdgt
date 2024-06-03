IsoRat<-function(raw,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){raw[is.na(raw)]<-0}
  formula<-"( raw$IIa6 + raw$IIb6 + raw$IIIa6 )/( raw$IIa5 + raw$IIb5 + raw$IIIa5 + raw$IIa6 + raw$IIb6 + raw$IIIa6 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(raw))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(raw,formula)}
  eval(parse(text = formula))
}