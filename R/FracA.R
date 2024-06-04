#' Calculate Fractional Abundance
#'
#' Calculates the fractional abundance of the selected GDGTs presented as peak areas.
#' @param gdgt Data frame with the peak areas of GDGTs. Each row should correspond to a sample and each column to each GDGT with the appropriated name format.
#' @param group GDGTs to be transformed into fractional abundances. Options are "iso" (isoGDGTs), "br" (brGDGTs), "br_extended" (5/6 isomer and 7-methyl), "isoGMGTs", or "OHGDGTs".
#' @param coerce Logic variable. If FALSE (default) the function will use the provided data. If TRUE, the function will incorporate missing GDGTs from each group as 0s for the calculations.
#' @param how Select how the fractional abundance will be calculated. Options are "each" (default), which calculates the fractional abundance of each group independently; and "all", which calculates the fractional abundance of each GDGT over the total sum of selected GDGTs.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the fractional abundances of all GDGTs selected for all samples provided.

FracA<-function(gdgt,group,coerce=FALSE,how="each",na.ignore=FALSE){
  GDGTs<-list("br"=c("IIIa5","IIIa6","IIIb5","IIIb6","IIIc5","IIIc6",
                     "IIa5","IIa6","IIb5","IIb6","IIc5","IIc6",
                     "Ia","Ib","Ic"),
              "iso"=c("GDGT0","GDGT1","GDGT2","GDGT3","Cren","Creni"),
              "br_extended"=c("IIIa56","IIIa7","IIIb7","IIa7"),
              "isoGMGTs"=c("GMGT0","GMGT1","GMGT2","GMGT3","GMGT4"),
              "OHGDGTs"=c("GDGTX0","GDGTX1","GDGTX2","GDGTX3")
  )
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  if(coerce==TRUE){
    noGDGT<-unlist(GDGTs[match(group,names(GDGTs))])[is.na(match(unlist(GDGTs[match(group,names(GDGTs))]),colnames(gdgt)))]
    nullGDGT<-data.frame(matrix(0,ncol = length(noGDGT),nrow=dim(gdgt)[1]))
    colnames(nullGDGT)<-noGDGT
    gdgt<-cbind(gdgt,nullGDGT)
    message(paste("Coerced ",paste(noGDGT,collapse = " ")))
  }
  if(how=="each"){
    out<-as.data.frame(matrix(NA,ncol = length(unlist(GDGTs[match(group,names(GDGTs))])),nrow = dim(gdgt)[1]))
    colnames(out)<-unlist(GDGTs[match(group,names(GDGTs))])
    for (i in 1:length(group)) {
    out[,match(GDGTs[[group[i]]],colnames(out))]<-gdgt[,na.omit(match(GDGTs[[group[i]]],colnames(gdgt)))]/apply(gdgt[,na.omit(match(GDGTs[[group[i]]],colnames(gdgt)))],1,sum)
    }
    out[is.na(out)]<-0
  } else if(how=="all"){
    out<-as.data.frame(matrix(NA,ncol = length(unlist(GDGTs[match(group,names(GDGTs))])),nrow = dim(gdgt)[1]))
    colnames(out)<-unlist(GDGTs[match(group,names(GDGTs))])
    for (i in 1:length(group)) {
      out[,match(GDGTs[[group[i]]],colnames(out))]<-gdgt[,na.omit(match(GDGTs[[group[i]]],colnames(gdgt)))]/apply(gdgt[,na.omit(match(unlist(GDGTs[match(group,names(GDGTs))]),colnames(gdgt)))],1,sum)
    }
    out[is.na(out)]<-0
  }
  out
}









