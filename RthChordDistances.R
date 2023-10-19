RthChordDistances<-function(llPoints1,llPoints2){
  RR<-6378.137 #Radius of the Earth in km
  N<-dim(llPoints1)[1]
  M<-dim(llPoints2)[1]
  
  Pts_paired_Vec<-cbind(kronecker(llPoints1,matrix(1,nrow = M,ncol=1)),kronecker(matrix(1,nrow = N,ncol=1),llPoints2))
  Half_Angles_AsVec<-asin(
    sqrt(
      sin((Pts_paired_Vec[,2]-Pts_paired_Vec[,4])*pi/180/2)^2 +
      cos(Pts_paired_Vec[,2]*pi/180)*cos(Pts_paired_Vec[,4]*pi/180)*sin(abs(Pts_paired_Vec[,1]-Pts_paired_Vec[,3])*pi/180/2)^2
      )
  )
  Chords_as_vec<-2*RR*sin(Half_Angles_AsVec)
  Chord_Dists_Mat<-t(matrix(Chords_as_vec,nrow = M,ncol = N))
  Chord_Dists_Mat
}