#' Generate mcmc chain based on the data
#'
#' @param sphericaldata An n by 3 matrix
#' @param b_mat An n by regime number matrix
#' @param iterations iteration number
#' @param seed needed for random number generation
#' @param para_true true parameters
#' @param r_true the true lengths, unknown in real spherical data case
#' @return An object of "mcmc_chain".
#' An object of class "mcmc_chain" is a list containing at least the following components:
#' @importFrom  tensorA to.tensor mul.tensor
#' @importFrom MASS mvrnorm
#' @examples
#' mcmc_chain(sphericaldata, b_mat, iterations=5000, seed=1234, para_true=NULL,r_true=NULL)
#' @export


mcmc_chain <- function(sphericaldata, b_mat, iterations=5000, seed=1234, para_true=NULL,r_true=NULL){
  ##dimension
  k <- ncol(sphericaldata)
  ##regime number
  H <- ncol(b_mat)
  ##sample size
  n <- nrow(sphericaldata)
  ##right now only consider AR(1)
  p <- 1

  ####values used for simulation
  n_vec <- colSums(b_mat)
  ind <- rowSums(b_mat*rep_row(1:H,n))
  ind1 <- ind
  ind1[1] <- 999
  b_pre <- array(0,c(k,k*H,n))
  for(nn in 1:n){
    b_pre[1:k,((ind[nn]-1)*k+(1:k)),nn] <- diag(k)
  }

  b_pre <- to.tensor(b_pre)
  b_op <- to.tensor(array(0,c(k,k*H,n)))


  #### initial values
  ### mu is mean vector of length (k*H)
  mu_present <- numeric(k*H)
  mu_present_mat <- matrix(0,k,H)
  for(h in 1:H){
    mu_present_mat[,h] <- mu_present[1:k+(h-1)*k]
  }
  ### gamma_mat is a (k-1) by H matrix with [,h] as gamma_present for regime h
  gamma_mat <- matrix(0, (k-1), H)
  ### sigma_conditional_arr is (k-1) by (k-1) by H array with [,,h]
  ### as sigma_conditional_present for regime h
  sigmamat_conditional_arr <- array(0, c((k-1),(k-1),H))
  for(h in 1:H)
    sigmamat_conditional_arr[,,h] <- diag(k-1)
  ### phi_arr is k by k by H array with [,,h]
  ### as phi_present for regime h
  phi_arr <- array(0,c(k,k,H))


  ### vector lengths
  r_present  <-  rep(1,n)
  ###compute the precision_full_arr from sigmamat_conditional_arr and gamma_mat
  precision_full_arr <- array(0,c(k,k,H))
  for(h in 1:H){
    sigmamat_full <- rbind( cbind(sigmamat_conditional_arr[,,h] + gamma_mat[,h]%*%t(gamma_mat[,h]),gamma_mat[,h]) , cbind(t(gamma_mat[,h]),1) )
    precision_full_arr[,,h] <- solve(sigmamat_full)
  }

  phi_time <- phi_arr[, , ind]
  augmenteddata_raw <- sphericaldata*rep_col(r_present,k)



  ####containers for sampling variables
  mu_all <- matrix(0, iterations, k*H)
  gamma_mat_all <- array(0,c((k-1), H, iterations))
  sigmamat_conditional_arr_all <- array(0, c((k-1),(k-1),H,iterations))
  phi_arr_all <- array(0,c(k,k,H,iterations))


  ####priors
  ###we use the same prior for all H regimes for simplicity
  ##mu~Normal
  mu_prior_mu <- numeric(k*H)
  mu_prior_sigmamat <- diag(k*H)*1e5
  mu_prior_precision <- solve(mu_prior_sigmamat)
  ##gamma~Normal
  gamma_prior_mu <- numeric(k-1)
  gamma_prior_sigmamat <- diag(k-1)*1e5
  gamma_prior_precision <- solve(gamma_prior_sigmamat)
  ##sigma~Wishart
  sigmamat_conditional_prior_df <- rep(4,H)
  sigmamat_conditional_prior_mat <- diag(k-1)
  ##phi~Normal
  phi_prior_mu <- numeric(k)
  phi_prior_sigmamat <- diag(k)*1e5
  phi_prior_precision <- solve(phi_prior_sigmamat)

  set.seed(seed)

  ####
  for(iter in 1:iterations){

    augmenteddata_raw_ten <- to.tensor(augmenteddata_raw)
    names(augmenteddata_raw_ten) <- c("n","k")

    phi_time <- to.tensor(phi_arr[, , ind])
    names(phi_time) <- c("k0","k","n")


    ###### getting mu
    precision_full_time <- to.tensor(precision_full_arr[, , ind])
    names(precision_full_time) <- c("k","k1","n")
    names(b_pre) <- c('k','kH','n')
    b_op <- to.tensor(array(0,c(k,k*H,n)))
    names(b_op) <- c('k','kH','n')
    b_op[,,1] <- b_pre[,,1]
    ##phi_jB_{j-1}
    lag1b_op <- mul.tensor(phi_time[,,2:n], i=c("k"), b_pre[,,1:(n-1)], by=c('n'))
    names(lag1b_op) <- c('k','kH','n')
    ##b_op <- B_{j}-phi_jB_{j-1}
    b_op[,,2:n] <- b_pre[,,2:n] - lag1b_op
    b_op2 <- b_op
    names(b_op2) <- c('k1','kH1','n')

    left <- mul.tensor(b_op,i=c("k"), precision_full_time, by = c('n'))

    mu_posterior_precision <- apply( mul.tensor( left,i=c("k1"), b_op2, by = c('n')), c(1,2), sum) +
      mu_prior_precision
    mu_posterior_sigmamat <- solve (mu_posterior_precision)

    augmenteddata <- augmenteddata_raw
    augmenteddata[2:n,] <- augmenteddata_raw[2:n,]-
      t(unclass(mul.tensor(phi_time[,,2:n], i=c("k"), augmenteddata_raw_ten[1:(n-1),], by=c('n'))))

    augmenteddata_ten <- to.tensor(augmenteddata)
    names(augmenteddata_ten) <- c("n","k1")

    mu_posterior_mu <- mu_posterior_sigmamat%*%(rowSums(mul.tensor(left, i = c("k1"), augmenteddata_ten, by = c("n")))+
                                                  mu_prior_precision%*%mu_prior_mu)


    mu_present <- mvrnorm(1,mu=mu_posterior_mu,Sigma=mu_posterior_sigmamat)
    mu_all[iter,] <- mu_present


    ###getting gamma_h
    mu_present_ten <- to.tensor(mu_present)
    names(mu_present_ten) <- "kH"

    augmenteddatacentered <- augmenteddata-
      t(unclass(mul.tensor(b_op, i=c("kH"), mu_present_ten)))

    for(h in 1:H){
      sigmamat_conditional_present <- sigmamat_conditional_arr[,,h]
      precision_conditional_present <- solve(sigmamat_conditional_present)
      gamma_posterior_precision <- sum(augmenteddatacentered[ind==h,k]^2)*precision_conditional_present + gamma_prior_precision
      gamma_posterior_sigmamat <- solve(gamma_posterior_precision)
      gamma_posterior_mu <- drop(gamma_posterior_sigmamat%*%(precision_conditional_present%*%(t(augmenteddatacentered[ind==h,1:(k-1)])%*%
                                                                                                augmenteddatacentered[ind==h,k])+
                                                               gamma_prior_precision%*%gamma_prior_mu))
      gamma_present <- mvrnorm(1,mu=gamma_posterior_mu,Sigma=gamma_posterior_sigmamat)

      gamma_mat[,h]<- gamma_present

    }

    gamma_mat_all[,,iter] <- gamma_mat


    ####getting sigmamat_conditional
    gamma_time <- t(gamma_mat[,ind])
    difference <- augmenteddatacentered[,1:(k-1)]-rep_col(augmenteddatacentered[,k],2)*gamma_time

    for(h in 1:H){
      sigmamat_contitional_posterior_mat <- sigmamat_conditional_prior_mat+t(difference[ind==h,])%*%difference[ind==h,]
      sigmamat_conditional_inverse <- rWishart(1,(n_vec[h]+sigmamat_conditional_prior_df[h]), Sigma=solve(sigmamat_contitional_posterior_mat))[,,1]
      sigmamat_conditional_arr[,,h] <- solve(sigmamat_conditional_inverse)
    }
    sigmamat_conditional_arr_all[,,,iter] <- sigmamat_conditional_arr



    for(h in 1:H){
      sigmamat_full <- rbind( cbind(sigmamat_conditional_arr[,,h] + gamma_mat[,h]%*%t(gamma_mat[,h]),gamma_mat[,h]) , cbind(t(gamma_mat[,h]),1) )
      precision_full_arr[,,h] <- solve(sigmamat_full)
    }
    ####getting phi

    for(h in 1:H){
      mu_present_mat[,h] <- mu_present[1:k+(h-1)*k]
    }
    mu_present_time <- mu_present_mat[,ind]
    w <- augmenteddata_raw - t(mu_present_time)

    for(kk in 1:k)
      for(h in 1:H)
      {
        phi_posterior_sigmamat <- solve(sum(w[(which(ind1==h)-1),kk]^2)*precision_full_arr[,,h] + phi_prior_precision)
        phi_posterior_mu <- phi_posterior_sigmamat %*% (rowSums(rep_row(w[(which(ind1==h)-1),kk],3) * (precision_full_arr[,,h]%*%(
          t(w[ind1==h,]-w[(which(ind1==h)-1),-kk]%*%t(phi_arr[,-kk,h])))))+
            drop(phi_prior_precision%*%phi_prior_mu))

        phi_present <- mvrnorm(1,mu = phi_posterior_mu,Sigma = phi_posterior_sigmamat)
        phi_arr[,kk,h] <- phi_present
      }
    phi_arr_all[,,,iter] <- phi_arr

    if(is.null(r_true)){####getting r_n
      nn=1
      mu_con <- mu_present_mat[,ind[nn]]
      a <- drop(sphericaldata[nn,]%*%precision_full_arr[,,ind[nn]]%*%sphericaldata[nn,])+
        drop(sphericaldata[nn,]%*%t(phi_arr[,,ind[nn+1]])%*%precision_full_arr[,,ind[nn+1]]%*%phi_arr[,,ind[nn+1]]%*%sphericaldata[nn,])
      b <- drop(sphericaldata[nn,]%*%precision_full_arr[,,ind[nn]]%*%mu_con)+
        drop(sphericaldata[nn,]%*%t(phi_arr[,,ind[nn+1]])%*%precision_full_arr[,,ind[nn+1]]%*%(
          phi_arr[,,ind[nn+1]]%*%mu_present_mat[,ind[nn]]-mu_present_mat[,ind[nn+1]]+augmenteddata_raw[(nn+1),]))
      bovera <- b/a
      uniform1 <- runif(1)
      uniform2 <- runif(1)
      v <- uniform1*exp(-0.5*a*(r_present[nn]-bovera)^2)
      root <- sqrt(-2*log(v)/a)
      upper <- bovera + root
      lower <- bovera + max(-bovera,-root)
      r_present[nn] <- ((upper^k-lower^k)*uniform2+lower^k)^(1/k)
      augmenteddata_raw[nn,] <- r_present[nn]*sphericaldata[nn,]

      for(nn in 2:(n-1)){
          mu_con <- mu_present_mat[,ind[nn]] +  phi_arr[,,ind[nn]]%*%(augmenteddata_raw[(nn-1),]-mu_present_mat[,ind[nn-1]])
          a <- drop(sphericaldata[nn,]%*%precision_full_arr[,,ind[nn]]%*%sphericaldata[nn,])+
            drop(sphericaldata[nn,]%*%t(phi_arr[,,ind[nn+1]])%*%precision_full_arr[,,ind[nn+1]]%*%phi_arr[,,ind[nn+1]]%*%sphericaldata[nn,])
          b <- drop(sphericaldata[nn,]%*%precision_full_arr[,,ind[nn]]%*%mu_con)+
            drop(sphericaldata[nn,]%*%t(phi_arr[,,ind[nn+1]])%*%precision_full_arr[,,ind[nn+1]]%*%(
              phi_arr[,,ind[nn+1]]%*%mu_present_mat[,ind[nn]]-mu_present_mat[,ind[nn+1]]+augmenteddata_raw[(nn+1),]))
          bovera <- b/a
          uniform1 <- runif(1)
          uniform2 <- runif(1)
          v <- uniform1*exp(-0.5*a*(r_present[nn]-bovera)^2)
          root <- sqrt(-2*log(v)/a)
          upper <- bovera + root
          lower <- bovera + max(-bovera,-root)
          r_present[nn] <- ((upper^k-lower^k)*uniform2+lower^k)^(1/k)
          augmenteddata_raw[nn,] <- r_present[nn]*sphericaldata[nn,]}
      nn=n
      mu_con <- mu_present_mat[,ind[nn]] +  phi_arr[,,ind[nn]]%*%(augmenteddata_raw[(nn-1),]-mu_present_mat[,ind[nn-1]])
      a <- drop(sphericaldata[nn,]%*%precision_full_arr[,,ind[nn]]%*%sphericaldata[nn,])
      b <- drop(sphericaldata[nn,]%*%precision_full_arr[,,ind[nn]]%*%mu_con)
      bovera <- b/a
      uniform1 <- runif(1)
      uniform2 <- runif(1)
      v <- uniform1*exp(-0.5*a*(r_present[nn]-bovera)^2)
      root <- sqrt(-2*log(v)/a)
      upper <- bovera + root
      lower <- bovera + max(-bovera,-root)
      r_present[nn] <- ((upper^k-lower^k)*uniform2+lower^k)^(1/k)
      augmenteddata_raw[nn,] <- r_present[nn]*sphericaldata[nn,]

    }else {augmenteddata_raw <- rep_col(r_true,3)*sphericaldata}

    print(iter)

  }


  result_mcmc <- list("para_true"=para_true,
                      "b_mat"=b_mat,
                      "mu_all"=mu_all,
                      "gamma_mat_all"=gamma_mat_all,
                      "sigmamat_conditional_arr_all"=sigmamat_conditional_arr_all,
                      "phi_arr_all"=phi_arr_all)

  return(result_mcmc)
}
