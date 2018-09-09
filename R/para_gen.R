#' Generate a set of true parameters
#'
#' @return An object of "para_true".
#' An object of class "para_true" is a list containing at least the following components:
#' @examples
#' para_gen()
#' @export

para_gen <- function(){
  k=3
  H=3
  rotation <- t(matrix(c(0, 0.8, 0.6, 0.8, -0.36, 0.48, 0.6, 0.48, -0.64),3,3))
  phi_true <- array(0,c(k,k,H))

  phi_true[,,1] <- rotation%*%diag(c(0.1,0.2,0.3))%*%t(rotation)
  phi_true[,,2] <- rotation%*%diag(c(0.2,-0.3,-0.2))%*%t(rotation)
  phi_true[,,3] <- rotation%*%diag(c(0.4,0.6,0.4))%*%t(rotation)

  mu_true <- matrix(0,k,H)
  mu_true[,1] <- c(1.26,-0.62,0.23)
  mu_true[,2] <- c(2,1,3)
  mu_true[,3] <- c(0.1,0.2,0.3)
  sigmamat_conditional_true=array(0,c((k-1),(k-1),H))
  sigmamat_conditional_true[,,1]=diag(2)
  sigmamat_conditional_true[,,2]=diag(2)
  sigmamat_conditional_true[,,3]=diag(2)
  sigmamat_true <- array(0,c(k,k,H))
  gamma_true <- matrix(c(0.3,-0.2,-0.5,-0.3,0.4,0.5),(k-1),H)


  sigmamat_true[,,1] <- sigmagen(sigmamat_conditional_true[,,1],gamma_true[,1])
  sigmamat_true[,,2] <- sigmagen(sigmamat_conditional_true[,,2],gamma_true[,2])
  sigmamat_true[,,3] <- sigmagen(sigmamat_conditional_true[,,3],gamma_true[,3])

  para_true <- list(mu_true=mu_true,
                    gamma_true = gamma_true,
                    sigmamat_conditional_true = sigmamat_conditional_true,
                    phi_true=phi_true)
  return(para_true)
}
