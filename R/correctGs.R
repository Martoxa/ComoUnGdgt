#' Generates equations for GDGT indices
#'
#' Compares the provided list of GDGTs to the required list of GDGTs for the index that one wants to calculate. Any GDGTs not provided in the dataset are removed from the equation and a warning is generated.
#' @param eq Equation for the index that wants to be calculated. This is provided by the index function as a text string.
#' @return The function returns the equation for the index.

correctGs<-function(eq){
  sub("gdgt.","",unlist(strsplit(eq," "))[grep("gdgt",unlist(strsplit(eq," ")))])
}

