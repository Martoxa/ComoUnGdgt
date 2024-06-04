#' TEX86 index
#'
#' Calculates the TEX86 index (Schouten, 2002) for isoGDGTs presented as either peak areas or fractional abundances.
#' @param gdgt Data frame with the peak areas or fractional abundances of isoGDGTs. Each row should correspond to a sample and each column to each brGDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all 6 isoGDGTs are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the TEX86 value for each sample provided.
#' @references Schouten, S., Hopmans, E. C., Schefuß, E., & Damste, J. S. S. (2002). Distributional variations in marine crenarchaeotal membrane lipids: a new tool for reconstructing ancient sea water temperatures?. Earth and Planetary Science Letters, 204(1-2), 265-274.

TEX86<-function(gdgt,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"( gdgt$GDGT2 + gdgt$GDGT3 + gdgt$Creni )/( gdgt$GDGT1 + gdgt$GDGT2 + gdgt$GDGT3 + gdgt$Creni )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  eval(parse(text = formula))
}
