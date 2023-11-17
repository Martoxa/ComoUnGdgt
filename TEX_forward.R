TEX_forward<-function(lat,lon,t,type="SST",mode="standard",tolerance=NA){
  if(mode=="analog"){
    if(is.na(tolerance)){stop("To use analog mode, enter a search tolerance in TEX units")}
  }
  
  if(type=="SST"){
    if(mode=="standard"){Params<-readRDS(file="SST_param_std.rds")
    } else if (mode=="analog"){
      Params<-readRDS(file="SST_param_ana.rds")
      DataIn<-readRDS(file = "SST_DataIn.rds")
    } else stop("Calibration mode must be either 'standard' or 'analog'")
  } else if(type=="subT"){
    if(mode=="standard"){
      Params<-readRDS(file="subT_param_std.rds")
      Data_In<-readRDS(file = "subT_DataIn.rds")
    } else if (mode=="analog"){Params<-readRDS(file="subT_param_ana.rds")
    } else stop("Calibration mode must be either 'standard' or 'analog'")
  } else stop("Calibration type must be either 'SST' or 'subT'")
  
  grid_half_space<-10
  
  if(mode=="standard"){
    Nloc<-length(lon)
    inder_g<-matrix(unlist(lapply(1:Nloc,function(a){which(abs(Params$Locs.Comp[,1]-lon[a])<=grid_half_space & abs(Params$Locs.Comp[,2]-lat[a])<= grid_half_space)})),nrow = Nloc,1)
    alpha_samples<-Params$alpha.samples.comp[inder_g,]
    beta_samples<-Params$beta.samples.comp[inder_g,]
  } else if(mode=="analog"){ 
    N_bg<-dim(matrix(unlist(DataIn[[1]][2]),ncol=2))[1]
    spatialMean<-matrix(unlist(lapply(1:N_bg,function(a){mean(unlist(DataIn[[1]][8])[unlist(DataIn[[1]][5])==a])})),nrow=N_bg,ncol=1)
    inder_g<- spatialMean >= mean(t)-tolerance & spatialMean <= mean(t)+tolerance
    if(sum(inder_g)==0){stop("Your search tolerance is too narrow")}
    alpha_samples<-as.data.frame(alpha_samples)[inder_g,]
    beta_samples<-as.data.frame(beta_samples)[inder_g,]
    tau2_samples<-kronecker(matrix(1,dim(alpha_samples)[1],1),Params$tau2.samples)
    alpha_samples<-matrix(as.matrix(alpha_samples),nrow=1,ncol=dim(alpha_samples)[1]*dim(alpha_samples)[2])
    beta_samples<-matrix(as.matrix(beta_samples),nrow=1,ncol=dim(beta_samples)[1]*dim(beta_samples)[2])
  }  
    tex86<-matrix(rnorm(n=dim(alpha_samples)[1]*dim(alpha_samples)[2],mean=c(t*beta_samples+alpha_samples),sd=c(kronecker(matrix(1,length(t),1),sqrt(Params$tau2.samples)))),nrow=dim(alpha_samples)[1],ncol = dim(alpha_samples)[2])
    
}