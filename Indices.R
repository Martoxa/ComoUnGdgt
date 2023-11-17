##brGDGTs only

MBT5<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"( fa$Ia + fa$Ib + fa$Ic )/( fa$Ia + fa$Ib + fa$Ic + fa$IIa5 + fa$IIb5 + fa$IIc5 + fa$IIIa5 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  eval(parse(text = formula))
} #checked Baxter et al., 2023

MBT6<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"( fa$Ia + fa$Ib + fa$Ic )/( fa$Ia + fa$Ib + fa$Ic + fa$IIa6 + fa$IIb6 + fa$IIIa6 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  eval(parse(text = formula))
}

fC<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"(( fa$IIb5 + fa$Ib + fa$IIb6 )+2*( fa$Ic ))/(( fa$IIIa5 + fa$IIIa6 + fa$IIa5 + fa$IIa6 + fa$Ia )+( fa$IIb5 + fa$Ib + fa$IIb6 )+( fa$Ic ))*0.5"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  eval(parse(text = formula))
}
CBTp<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"log10(( fa$Ic + fa$IIa6 + fa$IIb6 + fa$IIc6 + fa$IIIa6 + fa$IIIb6 + fa$IIIc6 )/( fa$Ia + fa$IIa5 + fa$IIIa5 ))"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  out<-eval(parse(text = formula))
  out[-10000 > out]<-NA
  out[out>10000]<-NA
  out
}
IR<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"( fa$IIa6 + fa$IIb6 + fa$IIIa6 )/( fa$IIa5 + fa$IIb5 + fa$IIIa5 + fa$IIa6 + fa$IIb6 + fa$IIIa6 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  eval(parse(text = formula))
}
IBT<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"-log10(( fa$IIa6 + fa$IIIa6 )/( fa$IIa5 + fa$IIIa5 ))"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  out<-eval(parse(text = formula))
  out[-10000 > out]<-NA
  out[out>10000]<-NA
  out
}
CI<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"( fa$Ia + fa$Ib )/( fa$Ia + fa$IIa5 + fa$IIIa5 + fa$IIIa6 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  eval(parse(text = formula))
}
RI<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"( fa$Ib + 2*( fa$Ic ))/( fa$Ia + fa$Ib + fa$Ic )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  eval(parse(text = formula))
}

## isoGDGTs only
TEX86<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"( fa$GDGT2 + fa$GDGT3 + fa$Crenp )/( fa$GDGT1 + fa$GDGT2 + fa$GDGT3 + fa$Crenp )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  eval(parse(text = formula))
}
MI<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"( fa$GDGT1 + fa$GDGT2 + fa$GDGT3 )/( fa$GDGT1 + fa$GDGT2 + fa$GDGT3 + fa$Crena + fa$Crenp )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  eval(parse(text = formula))
}
pGDGT0<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"fa$GDGT0 /( fa$GDGT0 + fa$Crena )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  eval(parse(text = formula))
} #Checked Zander, et al. in prep
G2G3<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"fa$GDGT2 / fa$GDGT3"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  out<-eval(parse(text = formula))
  out[out>10000]<-NA
  out[is.nan(out)]<-NA
  out
}

## br and isoGDGTs
BIT<-function(raw,complete=TRUE,data="area",out="all",na.ignore=FALSE){
  if(data=="area"){fa<-FracA(raw,out=out,na.ignore = na.ignore)
  } else if (data == "fa"){fa<-raw}
  if(na.ignore==TRUE){fa[is.na(fa)]<-0}
  formula<-"( fa$Ia + fa$IIa5 + fa$IIa6 + fa$IIIa5 + fa$IIIa6 )/( fa$Ia + fa$IIa5 + fa$IIa6 + fa$IIIa5 + fa$IIIa6 + fa$Crena )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(fa))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(fa,formula)}
  eval(parse(text = formula))
}  #checked Baxter et al., 2023
