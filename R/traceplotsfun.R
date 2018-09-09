#' Generate traceplot for mu
#'
#' @param trun number where the truncation is performed
#' @param iterations iteration number
#' @param reind regime index
#' @param result_mcmc results from mcmc_chain or mcmc_chain_vectorr, An object of "mcmc_chain".
#' @param mu_true the true parameter mu, NULL in real case
#' @param pdf_plot if true, save plot to pdf file
#' @param ... argument to be passed to pdf_plot
#' @importFrom  coda HPDinterval mcmc
#' @examples
#' plotmu(9001,10000,1,result_mcmc,mu_true,pdf_plot = T, width= 7.5,height= 7.5)
#' @export



plotmu <- function(trun = 9000,
                   iterations = 10000,
                   reind = 1,
                   result_mcmc,
                   mu_true=NULL,plot0=TRUE, pdf_plot=FALSE,...){
  if(pdf_plot) pdf(paste0("mu_",reind,".pdf"),...)
  mu_estimate <- matrix(0,3,3)
  colnames(mu_estimate)<-c("mean","lower","upper")
  rownames(mu_estimate)<-c("mu_1","mu_2","mu_3")
  if(plot0) par(mfrow=c(3,1))
  for(i in 1:3){
    dat <- result_mcmc$mu_all[trun:iterations,i+(reind-1)*3]
    mu_estimate[i,1]<-mean(dat)
    ci <- HPDinterval(mcmc(dat))
    mu_estimate[i,2]<-ci[1]
    mu_estimate[i,3]<-ci[2]
    if(plot0){
      plot(dat,type="l",ylab = substitute(mu[paste(i,'  ',reind)],list(i=i,reind=reind)))
      if(!is.null(mu_true))
        abline(h=mu_true[i,reind],col="red")
      abline(h=mu_estimate[i,1],col="green")
      abline(h=mu_estimate[i,2],col="green")
      abline(h=mu_estimate[i,3],col="green")
    }
  }
  if(pdf_plot) dev.off()
  return(mu_estimate)
}


#' Generate traceplot for gamma
#'
#' @param trun number where the truncation is performed
#' @param iterations iteration number
#' @param reind regime index
#' @param result_mcmc results from mcmc_chain or mcmc_chain_vectorr, An object of "mcmc_chain".
#' @param gamma_true the true parameter gamma, NULL in real case
#' @param pdf_plot if true, save plot to pdf file
#' @param ... argument to be passed to pdf_plot
#' @importFrom  coda HPDinterval mcmc
#' @examples
#' plotgamma(9001,10000,1,result_mcmc,gamma_true,pdf_plot = T, width= 7.5,height= 6)
#' @export


plotgamma <- function(trun = 9000,
                      iterations = 10000,
                      reind = 1,
                      result_mcmc,
                      gamma_true=NULL,plot0=TRUE,pdf_plot=FALSE,...){
  if(pdf_plot) pdf(paste0("gamma_",reind,".pdf"),...)
  gamma_estimate <- matrix(0,2,3)
  colnames(gamma_estimate)<-c("mean","lower","upper")
  rownames(gamma_estimate)<-c("gamma_1","gamma_2")
  if(plot0) par(mfrow=c(2,1))
  for(i in 1:2){
    dat <- result_mcmc$gamma_mat_all[i,reind,trun:iterations]
    gamma_estimate[i,1]<-mean(dat)
    ci <- HPDinterval(mcmc(dat))
    gamma_estimate[i,2]<-ci[1]
    gamma_estimate[i,3]<-ci[2]
    if(plot0){
      plot(dat,type="l",ylab = substitute(gamma[paste(i,'  ',reind)],list(i=i,reind=reind)))
      if(!is.null(gamma_true))
        abline(h=gamma_true[i,reind],col="red")
      abline(h=gamma_estimate[i,1],col="green")
      abline(h=gamma_estimate[i,2],col="green")
      abline(h=gamma_estimate[i,3],col="green")
    }
  }

  if(pdf_plot) dev.off()
  return(gamma_estimate)
}


#' Generate traceplot for sigmamat_conditional
#'
#' @param trun number where the truncation is performed
#' @param iterations iteration number
#' @param reind regime index
#' @param result_mcmc results from mcmc_chain or mcmc_chain_vectorr, An object of "mcmc_chain".
#' @param sigmamat_conditional_true the true parameter sigmamat_conditional_true, NULL in real case
#' @param pdf_plot if true, save plot to pdf file
#' @param ... argument to be passed to pdf_plot
#' @importFrom  coda HPDinterval mcmc
#' @examples
#' plotsigmamat_conditional(9001,10000,1,result_mcmc,sigmamat_conditional_true,pdf_plot = T, width= 9,height= 6)
#' @export


plotsigmamat_conditional <- function(trun = 9000,
                                     iterations = 10000,
                                     reind = 1,
                                     result_mcmc,
                                     sigmamat_conditional_true=NULL,plot0=TRUE, pdf_plot=FALSE,...){
  if(pdf_plot) pdf(paste0("sigmamat_conditional_",reind,".pdf"),...)
  sigcon_estimate <- matrix(0,4,3)
  colnames(sigcon_estimate)<-c("mean","lower","upper")
  rownames(sigcon_estimate)<-c("sigcon_1_1","sigcon_1_2","sigcon_2_1","sigcon_2_2")
  if(plot0) par(mfrow=c(2,2))
  for(i in 1:2)
    for(j in 1:2){
      dat <- result_mcmc$sigmamat_conditional_arr_all[i,j,reind,trun:iterations]
      sigcon_estimate[(2*i+j-2),1]<-mean(dat)
      ci <- HPDinterval(mcmc(dat))
      sigcon_estimate[(2*i+j-2),2]<-ci[1]
      sigcon_estimate[(2*i+j-2),3]<-ci[2]
      if(plot0){
        plot(dat,type="l",ylab = substitute(sigma[paste(i,j,'  ',reind)],list(i=i,j=j,reind=reind)))
        if(!is.null(sigmamat_conditional_true))
          abline(h=sigmamat_conditional_true[i,j,reind],col="red")

        abline(h=sigcon_estimate[(2*i+j-2),1],col="green")
        abline(h=sigcon_estimate[(2*i+j-2),2],col="green")
        abline(h=sigcon_estimate[(2*i+j-2),3],col="green")
      }
    }

  if(pdf_plot) dev.off()
  return(sigcon_estimate)
}




#' Generate traceplot for phi
#'
#' @param trun number where the truncation is performed
#' @param iterations iteration number
#' @param reind regime index
#' @param result_mcmc results from mcmc_chain or mcmc_chain_vectorr, An object of "mcmc_chain".
#' @param phi_true the true parameter phi_true, NULL in real case
#' @param pdf_plot if true, save plot to pdf file
#' @param ... argument to be passed to pdf_plot
#' @importFrom  coda HPDinterval mcmc
#' @examples
#' plotphi(9001,10000,1,result_mcmc,phi_true,pdf_plot = T, width= 9,height= 9)
#' @export


plotphi <- function(trun = 9000,
                    iterations = 10000,
                    reind = 1,
                    result_mcmc,
                    phi_true=NULL,plot0=TRUE, pdf_plot=FALSE,...){
  if(pdf_plot) pdf(paste0("phi_",reind,".pdf"),...)
  phi_estimate <- matrix(0,9,3)
  colnames(phi_estimate)<-c("mean","lower","upper")
  rownames(phi_estimate)<-c("phi_1_1","phi_1_2","phi_1_3",
                            "phi_2_1","phi_2_2","phi_2_3",
                            "phi_3_1","phi_3_2","phi_3_3")
  if(plot0) par(mfrow=c(3,3))
  for(i in 1:3)
    for(j in 1:3){
      dat <- result_mcmc$phi_arr_all[i,j,reind,trun:iterations]
      phi_estimate[(3*i+j-3),1]<-mean(dat)
      ci <- HPDinterval(mcmc(dat))
      phi_estimate[(3*i+j-3),2]<-ci[1]
      phi_estimate[(3*i+j-3),3]<-ci[2]
      if(plot0){
        plot(dat,type="l",ylab = substitute(phi[paste(i,j,'  ',reind)],list(i=i,j=j,reind=reind)))
        if(!is.null(phi_true))
          abline(h=phi_true[i,j,reind],col="red")

        abline(h=phi_estimate[(3*i+j-3),1],col="green")
        abline(h=phi_estimate[(3*i+j-3),2],col="green")
        abline(h=phi_estimate[(3*i+j-3),3],col="green")
      }

    }

  if(pdf_plot) dev.off()
  return(phi_estimate)
}

#' Generate traceplots for all parameters&regimes, save them to pdfs
#'
#' @param trun number where the truncation is performed
#' @param iterations iteration number
#' @param reall vector of all regime numbers
#' @param result_mcmc results from mcmc_chain or mcmc_chain_vectorr, An object of "mcmc_chain".
#' @param para_true list of true parameters, NULL in the real case
#' @importFrom  coda HPDinterval mcmc
#' @examples
#' allplot(9001,10000,c(1,2),result_mcmc, para_true)
#' @export


allplot <- function(trun, iterations, reall,result_mcmc,plot0=TRUE,pdf_plot=TRUE ,para_true=NULL){
  para_estimate <- list()
  for(reind in reall){
    para_estimate[[paste0("mu",reind)]]<-
      plotmu(trun,iterations,reind,result_mcmc,para_true$mu_true,plot0=plot0,pdf_plot = pdf_plot, width= 7.5,height= 7.5)
    para_estimate[[paste0("gamma",reind)]]<-
      plotgamma(trun,iterations,reind,result_mcmc,para_true$gamma_true,plot0=plot0,pdf_plot = pdf_plot, width= 7.5,height= 6)
    para_estimate[[paste0("sigcon",reind)]]<-
      plotsigmamat_conditional(trun,iterations,reind,result_mcmc,para_true$sigmamat_conditional_true,plot0=plot0,pdf_plot = pdf_plot, width= 9,height= 6)
    para_estimate[[paste0("phi",reind)]]<-
      plotphi(trun,iterations,reind,result_mcmc,para_true$phi_true,plot0=plot0,pdf_plot = pdf_plot, width= 9,height= 9)

  }
  return(para_estimate)
}
