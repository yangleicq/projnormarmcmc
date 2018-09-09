#' generate Covariance matrix from conditional covariance matrix and \eqn{\beta}
#'
#' @param sigmamat_conditional (k-1) by (k-1) conditional convariance matrix
#' @param gamma vector of length (k-1)

#' @return Covariance matrix.
#' An k by k non-negative definite matrix
#' @examples
#' sigmagen(diag(1,1),c(-0.5,0.5))
#' @export


sigmagen<-function(sigmamat_conditional, gamma){
  return(rbind( cbind(sigmamat_conditional + gamma%*%t(gamma),gamma) , cbind(t(gamma),1) ))
}
