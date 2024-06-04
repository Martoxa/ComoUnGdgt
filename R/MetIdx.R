#' Methane Index
#'
#' Calculates the Methane index (Zhang, 2011) for isoGDGTs presented as either peak areas or fractional abundances.
#' @param gdgt Data frame with the peak areas or fractional abundances of isoGDGTs. Each row should correspond to a sample and each column to each isoGDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all isoGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the MI value for each sample provided.
#' @references Zhang, Y. G., Zhang, C. L., Liu, X. L., Li, L., Hinrichs, K. U., & Noakes, J. E. (2011). Methane Index: A tetraether archaeal lipid biomarker indicator for detecting the instability of marine gas hydrates. Earth and Planetary Science Letters, 307(3-4), 525-534.


MetIdx<-function(gdgt,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"( gdgt$GDGT1 + gdgt$GDGT2 + gdgt$GDGT3 )/( gdgt$GDGT1 + gdgt$GDGT2 + gdgt$GDGT3 + gdgt$Cren + gdgt$Creni )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  eval(parse(text = formula))
}
