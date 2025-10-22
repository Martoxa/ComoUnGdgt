#' Ri/b index
#'
#' Calculates Ri/b (Xie et al., 2012) for GDGTs presented as either peak areas or fractional abundances.
#' @param gdgt Data frame with the peak areas or fractional abundances of GDGTs. Each row should correspond to a sample and each column to each GDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all GDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the Ri/b value for each sample provided.
#' @references Xie, S., Pancost, R. D., Chen, L., Evershed, R. P., Yang, H., Zhang, K., ... & Xu, Y. (2012). Microbial lipid records of highly alkaline deposits and enhanced aridity associated with significant uplift of the Tibetan Plateau in the Late Miocene. Geology, 40(4), 291-294.
#' @export


Rib<-function(gdgt,complete=TRUE,na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"( gdgt$GDGT0 + gdgt$GDGT1 + gdgt$GDGT2 + gdgt$GDGT3 + gdgt$Cren + gdgt$Creni )/( gdgt$Ia + gdgt$Ib + gdgt$Ic + gdgt$IIa5 + gdgt$IIa6 + gdgt$IIb5 + gdgt$IIb6 + gdgt$IIc5 + gdgt$IIc6 + gdgt$IIIa5 + gdgt$IIIa6 + gdgt$IIIb5 + gdgt$IIIb6 + gdgt$IIIc5 + gdgt$IIIc6 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  eval(parse(text = formula))
}
