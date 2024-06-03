MetIdx<-function(raw,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){raw[is.na(raw)]<-0}
  formula<-"( raw$GDGT1 + raw$GDGT2 + raw$GDGT3 )/( raw$GDGT1 + raw$GDGT2 + raw$GDGT3 + raw$Cren + raw$Creni )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(raw))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(raw,formula)}
  eval(parse(text = formula))
}