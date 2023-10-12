correctGs<-function(eq){
  sub("fa.","",unlist(strsplit(eq," "))[grep("fa",unlist(strsplit(eq," ")))])
}

partialEq<-function(fa,eq){
  out<-unlist(strsplit(eq," "))
  if(TRUE %in% is.na(match(out[grep("fa",out)],paste0("fa$",colnames(fa))))){
    warning(
      paste("Removed variables:",paste(
          sub("fa.","",out[grep("fa",out)[is.na(
                match(out[grep("fa",out)],paste0("fa$",colnames(fa)))
                )]]),collapse=","))
      )
  }
  out[grep("fa",out)[is.na(match(out[grep("fa",out)],paste0("fa$",colnames(fa))))]]<-0
  out<-paste(out,collapse=" ")
  out
}