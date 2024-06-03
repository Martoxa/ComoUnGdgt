RinIdx<-function(raw,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){raw[is.na(raw)]<-0}
  formula<-"( raw$Ib + 2*( raw$Ic ))/( raw$Ia + raw$Ib + raw$Ic )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(raw))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(raw,formula)}
  eval(parse(text = formula))
}