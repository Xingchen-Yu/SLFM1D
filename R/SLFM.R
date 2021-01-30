#'@title SLFM1D
#'
#'@description Fits a one-dimensional Spherical Latent Factor Model (Circular Factor Model) on a dataset consists of binary responses
#'
#'@param ymat A binary matrix containing 0 and 1 or missing values(NA)
#'@param n_pos Number of posterior samplers to collect
#'@param burnin Burnin/warmup period for the MCMC
#'@param thin Thining of the posterior chain, the default is 1
#'@param hyperparams Hyperparameters
#'@param initial_values Initial values for input
#'@param core Number of cpu cores to use for parallel computing
#'@param cluster_seed Seed number for reproducibility when 2 or more cores are used.
#'@return A list containing the posterior samples for each model parameter
#'@examples
#'data(ymat)
#'set.seed(2021)
#'posterior_chain = SLFM(ymat = ymat,n_pos=200,burnin=50,core=1)
#'@importFrom stats "rgamma"
#'@importFrom stats "runif"
#'@import "snow"
#'@import "snowfall"
#'@import "Rcpp"
#'@import "RcppArmadillo"
#'@import "rlecuyer"
#'@useDynLib SLFM1D
#' @importFrom Rcpp "sourceCpp"
#'@export
#'
SLFM = function(ymat,n_pos=1000,burnin=500,thin = 1, hyperparams=list(a = 1, b = 1/10, ccc_a = 1, ccc_b=25,kappa_a = 1,omega_sd=0.1,kappa_sd=0.5,
                                                                      j_epi = 0.02, i_epi = 0.02, j_leap = 5, i_leap = 5,skip = 50,jitter = T),
                initial_values=NULL,core=2,cluster_seed=1234){
  ### checking and installing required packages###
  # required_package = c('Rcpp','snowfall','RcppArmadillo')
  # check_package = sum(unlist(lapply(required_package, require, character.only = TRUE)))==3
  # if(check_package ==F){
  #   install.packages(required_package,repos = "http://cran.us.r-project.org")
  #   lapply(required_package, require, character.only = TRUE)
  # }
  #######################################################
  # print('Compiling Source Code')
  # sourceCpp(code=RcppCode)
  # print('Compiling done')
  sfInit(parallel=TRUE,cpus=core)
  sfClusterSetupRNG( type="RNGstream",seed=cluster_seed)
  sfLibrary("Rcpp", character.only=TRUE)
  sfLibrary("RcppArmadillo", character.only=TRUE)
  sfLibrary('SLFM1D', character.only=TRUE)
  # sfExport(list=c('RcppCode'),namespace = 'SLFM1D')
  #######################################################

  iter = n_pos * thin + burnin
  jitter = hyperparams[['jitter']]

  mu = 0
  a = hyperparams[['a']]
  b = hyperparams[['b']]
  ccc_a = hyperparams[['ccc_a']]
  ccc_b = hyperparams[['ccc_b']]
  kappa_a = hyperparams[['kappa_a']]
  omega_sd = hyperparams[['omega_sd']]

  nr = nrow(ymat)
  nc = ncol(ymat)

  t_sig = rep(hyperparams[['kappa_sd']],nc)
  omega = a/b
  ccc = ccc_a/ccc_b
  skip_temp = hyperparams[['skip']]
  skip = ifelse(skip_temp>iter,iter,skip_temp)

  if(length(initial_values)>0){

    print('loading initial values from user inputs')
    kappa = initial_values[['kappa']]
    beta = initial_values[['beta']]
    tau_no = initial_values[['tau_no']]
    tau_yes = initial_values[['tau_yes']]

  }else{

    kappa = rep(0.1,nc)
    beta = runif(nr,-pi/2,pi/2)
    tau_no = runif(nc,-pi,pi)
    tau_yes = runif(nc,-pi,pi)

  }

  if(jitter == T){

    b_range = c(hyperparams[['i_epi']]/2,hyperparams[['i_epi']]*2)
    yn_range = c(hyperparams[['j_epi']]/2,hyperparams[['j_epi']]*2)

    l_range_tau = c(1,2 * hyperparams[['j_leap']])
    l_range = c(1,2*hyperparams[['i_leap']])

  }else{

    delta = rep(hyperparams[['i_epi']],nr)
    delta_yes = delta_no = rep(hyperparams[['j_epi']],nc)

    delta2 = delta/2
    delta2_yes = delta_yes/2
    delta2_no = delta_no/2
  }

  lol = c(sin(mu),cos(mu))
  cbeta_prior = lol * omega
  c_alpha = kappa_a*nc+ccc_a

  ######parallel miscellaneous############################
  nr_par = round(seq(0,nr,length.out = core+1))
  nc_par = round(seq(0,nc,length.out = core+1))
  sfExport(list=ls())
  na = which(is.na(ymat==T)) - 1
  len_na = length(na)
  na.position = which(is.na(ymat)==T, arr.ind = T)
  impute = ifelse(nrow(na.position)==0,F,T)
  i_index = as.numeric(na.position[,1]) - 1
  j_index = as.numeric(na.position[,2]) - 1

  # pos_pred = pos_pred2 = pos_pred3 = matrix(0,nr,nc)
  no_na = which(is.na(ymat)==F)
  yes_master = no_master = kappa_master = matrix(0,nc,n_pos)
  beta_master = matrix(0,nr,n_pos)
  omega_master = ccc_master = rep(0,n_pos)
  likeli_chain = rep(0,iter)

  omega_ratio = 0
  beta_accept_rs_all = rep(0,nr)
  kappa_accept_rs = rep(0,nc)
  kappa_accept_rs_all = no_accept_rs_all = yes_accept_rs_all = rep(0,nc)
  #######################################
  # if(core>1){
  #   print('loading Rcpp code to all cores')
  # }
  # sfClusterEval(sourceCpp(code = RcppCode))
  # sfClusterEval(sourceCpp(file="./src/functions.cpp"))
  # sfClusterEval(sourceCpp(file="./src/RcppExports.cpp"))
  # sfClusterEval(sourceCpp(file="./src/functions_header.h"))
  core_1 = core - 1
  node = 0:core_1
  j = 1

  print(paste0('Running ',iter," iterations"))
  for(i in 1:iter){
    cat("\rProgress: ",i,"/",iter)
    if(i %in% seq(1,iter,skip)){
      #### step size and leap steps jittering
      if(jitter == T){

        leap = sample(l_range[1]:l_range[2],nr,replace=T)
        leap_tau = sample(l_range_tau[1]:l_range_tau[2],nc,replace=T)

        delta = runif(nr,b_range[1],b_range[2])
        delta_yes = delta_no = runif(nc,yn_range[1],yn_range[2])

        delta2 = delta/2
        delta2_yes = delta_yes/2
        delta2_no = delta_no/2
      }

      #######################################
      ### tuning the proposal variance of the scale parameter \kappa during the initial 2000 iterations
      if(i<burnin){
        ks = kappa_accept_rs/skip
        kappa_skip = min(ks)
        ### targeting acceptance ratio between 0.3 and 0.6#####
        out = update_tsig(0.6,0.3,t_sig,ks,nc)
        t_sig = out[[1]]
        kappa_mod = out[[2]]

        kappa_accept_rs = rep(0,nc)
        # print(paste0('percent kappa changed ', kappa_mod))
        # print(paste0('min kappa acceptance ', kappa_skip))
      }
      sfExport("delta",'delta2','delta_yes','delta2_yes','delta_no','delta2_no','t_sig','leap','leap_tau')

      if(impute==T){
        ### impute the missing values
        ymat = impute_NA(na, i_index, j_index, ymat, tau_yes, tau_no, beta, kappa, len_na)
        sfExport('ymat')
      }
    }
    ################################################################################
    ### Update scale parameter \kappa
    out = sfLapply(node,wrapper_kappa,nc_par,nr,beta,tau_yes,tau_no,kappa,ymat,kappa_a,ccc,t_sig)
    kappa = unlist(lapply(out,"[[",1))
    sfExport("kappa")
    ###update hyperprior parameter for \kappa
    ccc = rgamma(1,c_alpha,ccc_b+sum(kappa))
    sfExport("ccc")

    kappa_ratio = sum(sapply(out,"[[",2))/nc
    haha = unlist(lapply(out,"[[",3))
    kappa_accept_rs = kappa_accept_rs  + haha
    kappa_accept_rs_all = kappa_accept_rs_all  + haha
    ################################################################################
    ### update ideal points \beta_i's
    out = sfLapply(node,wrapper_beta,nr_par,delta,delta2,leap,nc,omega,cbeta_prior,beta,tau_yes,tau_no,kappa,ymat)
    beta = unlist(lapply(out,"[[",1))
    beta_ratio = sum(sapply(out,"[[",2))/nr

    haha = unlist(lapply(out,"[[",3))
    beta_accept_rs_all = beta_accept_rs_all  + haha
    sfExport("beta")
    ################################################################################
    ### update \psi_j's
    out = sfLapply(node,wrapper_yes,nc_par,delta_yes,delta2_yes,leap_tau,nr,beta,tau_yes,tau_no,kappa,ymat)
    tau_yes = unlist(lapply(out,"[[",1))
    yes_ratio = sum(sapply(out,"[[",2))/nc

    haha = unlist(lapply(out,"[[",3))
    yes_accept_rs_all = yes_accept_rs_all  + haha
    sfExport("tau_yes")
    ################################################################################
    ### update \zeta_j's
    out = sfLapply(node,wrapper_no,nc_par,delta_no,delta2_no,leap_tau,nr,beta,tau_yes,tau_no,kappa,ymat)
    tau_no = unlist(lapply(out,"[[",1))
    no_ratio = sum(sapply(out,"[[",2))/nc

    haha = unlist(lapply(out,"[[",3))
    no_accept_rs_all = no_accept_rs_all  + haha
    sfExport("tau_no")
    ################################################################################
    ### update hyperparameter \omega for \beta_i's
    out = update_omega(omega,beta,nr,a,b,omega_sd)
    omega = out[[1]]
    omega_ratio = omega_ratio + out[[2]]
    cbeta_prior = lol * omega
    sfExport('omega',"cbeta_prior")

    ### compute joint loglikelihood
    waic_out = sfLapply(node,wrapper_waic,nc_par,nr,beta,tau_yes,tau_no,kappa,ymat)
    temp = do.call("cbind",lapply(waic_out,"[[",1))
    sum_temp = sum(temp)
    likeli_chain[i] = sum_temp

    ##maybe adding some warning features?
    ## output average acceptance ratio and joint loglikelihood every 100 iterations
    # if(i %in% seq(1,iter,100)){
    #   cat("\rProgress: ",i,"/",iter)
    #   print(paste('beta acceptance ratio is',round(beta_ratio,digits = 2)))
    #   print(paste('yes acceptance ratio is',round(yes_ratio,digits = 2)))
    #   print(paste('no acceptance ratio is',round(no_ratio,digits = 2)))
    #   print(paste('kappa acceptance ratio is',round(kappa_ratio,digits = 2)))
    #   print(paste('omega acceptance ratio is',round(omega_ratio/i,digits = 2)))
    #   # print(paste('loglikelihood is ',round(sum_temp,digits = 0)))
    # }
    ### record paratmer after burnin and compute waic using running sums
    if(i>burnin & (i %in% seq(burnin+1,iter,thin))){
      beta_master[,j] = beta
      yes_master[,j] = tau_yes
      no_master[,j] = tau_no
      omega_master[j] = omega
      kappa_master[,j] = kappa
      # pos_pred = pos_pred + do.call("cbind",lapply(waic_out,"[[",2))
      # pos_pred2 = pos_pred2 + temp
      # pos_pred3 = pos_pred3 + do.call("cbind",lapply(waic_out,"[[",3))
      j = j + 1
    }
  }
  sfStop()
  return(list(beta=beta_master,psi=yes_master,zeta=no_master,omega=omega_master,kappa=kappa_master,likeli=likeli_chain))
}


