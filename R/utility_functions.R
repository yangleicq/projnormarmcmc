#'
#'
#'Takes a vector x and an integer n and return n by x matrix
#'@param x a vector
#'@param n an integer
#'@return n by x matrix
#'@examples
#' rep_row(c(1,2,3), 3)
#'@export


rep_row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#'
#'
#'Takes a vector x and an integer n and return x by n matrix
#'@param x a vector
#'@param n an integer
#'@return x by n matrix
#'@examples
#' rep_col(c(1,2,3), 3)
#'@export


rep_col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


#'
#'
#'Calculate the cross product of two 3-dimensional vectors v1, v2
#'@param v1 3-dimensional vector
#'@param v2 3-dimensional vector
#'@return cross product of v1,v2, a 3-dimensional vector
#'@examples
#' cross_product(c(1,2,3),c(4,5,6))
#'@export


cross_product <- function(v1,v2){
  v1_cross_v2 <- rbind(v1,v2)
  u <- c()
  u[1] <- det(v1_cross_v2[,2:3])
  u[2] <- -det(v1_cross_v2[,c(1,3)])
  u[3] <- det(v1_cross_v2[,1:2])
  return(u)
}

#'
#'
#'Calculate the dihedral angel of 4 consecutive points in three dimensional space. The formula is
#'from wikipedia https://en.wikipedia.org/wiki/Dihedral_angle
#'@param p1 point 1
#'@param p2 point 2
#'@param p3 point 3
#'@param p4 point 4
#'@return dihedral angel of 4 consecutive points
#'@examples
#' dihedral(c(1,2,3),c(4,5,6),c(7,8,9),c(10,11,12))
#'@export

dihedral <- function(p1,p2,p3,p4){
  b1 <- p2-p1
  b2 <- p3-p2
  b3 <- p4-p3
  n1 <- cross_product(b1,b2)
  n1 <- n1/sqrt(sum(n1^2))
  n2 <- cross_product(b2,b3)
  n2 <- n2/sqrt(sum(n2^2))
  b2.unit <- b2/sqrt(sum(b2^2))
  x <- sum(n1*n2)
  y <- sum(cross_product(n1,n2)*b2.unit)
  angle <- atan2(y,x)
  return(angle)
}

#'
#'
#'Calculate the bond angel of 3 consecutive points in three dimensional space.
#'@param p1 point 1
#'@param p2 point 2
#'@param p3 point 3
#'@return bond angel of 3 consecutive points
#'@examples
#' bond(c(1,2,3),c(4,5,6),c(7,8,9))
#'@export

bond <- function(p1,p2,p3){
  b1 <- p1-p2
  b2 <- p3-p2
  angle <- acos(sum(b1*b2)/sqrt(sum(b1^2)*sum(b2^2)))
  return(angle)
}

