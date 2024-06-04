#' Earth Chord Distances
#'
#' Calculates the chordal distances (in km) between points 1 and 2.
#' @param llPoints1 Nx2 matrix where each row is a longitude-latitude pair (-180 to 180, and -90 to 90, respectively).
#' @param llPoints2 Mx2 matrix where each row is a longitude-latitude pair (-180 to 180, and -90 to 90, respectively).
#' @return A NxM matrix of chordal distance between the ith row of llPoint1 and the jth row of llPoint2. Radius of the earth is taken as 6378.137
#' @author Jessica Tierney https://github.com/jesstierney/BAYSPAR

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
