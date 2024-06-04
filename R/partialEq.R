#' Generates partial equations for GDGT indices
#'
#' Compares the provided list of GDGTs to the required list of GDGTs for the index that one wants to calculate. Any GDGTs not provided in the dataset are removed from the equation and a warning is generated.
#' @param gdgt Data frame of peak areas or fractional abundances of GDGTs where each row is a sample and each column a GDGT with the appropriate name format.
#' @param eq Equation for the index that wants to be calculated. This is provided by the index function as a text string.
#' @return The function returns the modified equation for the index.

partialEq<-function(gdgt,eq){
  out<-unlist(strsplit(eq," "))
  if(TRUE %in% is.na(match(out[grep("gdgt",out)],paste0("gdgt$",colnames(gdgt))))){
    warning(
      paste("Removed variables:",paste(
        sub("gdgt.","",out[grep("gdgt",out)[is.na(
          match(out[grep("gdgt",out)],paste0("gdgt$",colnames(gdgt)))
        )]]),collapse=","))
    )
  }
  out[grep("gdgt",out)[is.na(match(out[grep("gdgt",out)],paste0("gdgt$",colnames(gdgt))))]]<-0
  out<-paste(out,collapse=" ")
  out
}
