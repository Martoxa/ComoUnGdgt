#' Ratio of GDGT2 to GDGT3
#'
#' Calculates the ratio of GDGT2 to GDGT3 fractional abundances of isoGDGTs. Each row should correspond to a sample and each column to each isoGDGT with the appropriated name format.
#' @param gdgt Data frame with the peak areas or fractional abundances of isoGDGTs. Each row should correspond to a sample and each column to each brGDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all isoGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the GDGT2/GDGT3 value for each sample provided.
#' @export

G2G3<-function(gdgt,complete=TRUE,na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"gdgt$GDGT2 / gdgt$GDGT3"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  out<-eval(parse(text = formula))
  out[out>10000]<-NA
  out[is.nan(out)]<-NA
  out
}
