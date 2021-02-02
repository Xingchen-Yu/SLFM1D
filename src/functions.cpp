#include "functions_header.h"
// using namespace "functions_header";

//#'@export
// [[Rcpp::export]]
List update_beta(int t, arma::vec nr_par, arma::vec epsilon_vec, arma::vec epsilon2_vec, arma::vec leap_vec,int nc, double omega,
                 arma::vec cbeta_prior, arma::vec beta,arma::vec tau_yes,arma::vec tau_no,arma::vec kappa,arma::mat ymat){
  arma::vec nu = arma::zeros(p);
  arma::vec x_temp = arma::zeros(p);
  arma::vec x = arma::zeros(p);
  double count = 0;
  double alpha, ae, cosat, sinat, h, h_new, accept,beta_new,beta_prev;
  int start = nr_par(t);
  int end = nr_par(t+1);
  int e_s = end - start;
  arma::vec beta_out = arma::zeros(e_s);
  arma::vec accept_chain = arma::ones(e_s);

  for(int i = start; i< end; i++,count++ ){
    double epsilon = epsilon_vec(i);
    double epsilon2 = epsilon2_vec(i);
    int leap = leap_vec(i);
    arma::rowvec ymat_row = ymat.row(i);
    beta_prev = beta_new = beta(i);
    x = s_to_e(beta_prev,x);
    nu = arma::randn(p);
    nu = (dp - x * x.t()) * nu;

    h = likeli_beta(omega,beta_prev,tau_yes,tau_no,kappa,nc,ymat_row) - as_scalar(0.5 * nu.t() * nu);
    for(int j = 0; j < leap; j++) {
      nu = nu + epsilon2 * gradient_beta(cbeta_prior, beta_new, tau_yes, tau_no,kappa, nc, ymat_row);
      nu = (dp - x * x.t()) * nu;

      alpha = norm(nu);
      ae = alpha * epsilon;
      cosat = cos(ae);
      sinat = sin(ae);
      x_temp = x;
      //geodesic flow//
      x = x * cosat + nu/alpha * sinat;
      nu = nu * cosat - alpha * x_temp * sinat;
      //
      beta_new = atan2(x(0),x(1));
      //
      nu = nu + epsilon2 * gradient_beta(cbeta_prior, beta_new, tau_yes, tau_no,kappa ,nc,ymat_row);
      nu = nu = (dp - x * x.t()) * nu;

    }
    h_new = likeli_beta(omega,beta_new,tau_yes,tau_no,kappa,nc,ymat_row) - as_scalar(0.5 * nu.t() * nu);
    accept = exp(h_new - h);
    if(accept < as_scalar(arma::randu(1))){
      beta_new = beta_prev;
      e_s = e_s - 1;
      accept_chain(count) = 0;
    }
    beta_out(count) = beta_new;
  }
  return List::create(beta_out,e_s,accept_chain);
}
//#'@export
// [[Rcpp::export]]
List update_tau_yes(int t, arma::vec nc_par,arma::vec epsilon_vec, arma::vec epsilon2_vec, arma::vec leap_vec, int nr,
                    arma::vec beta,arma::vec tau_yes,arma::vec tau_no,arma::vec kappa,arma::mat ymat){
  int start = nc_par(t);
  int end = nc_par(t+1);
  int e_s = end - start;
  arma::vec yes_out = arma::zeros(e_s);
  arma::vec accept_chain = arma::ones(e_s);

  arma::vec nu = arma::zeros(p);
  arma::vec x_temp = arma::zeros(p);
  arma::vec x = arma::zeros(p);
  double alpha, ae, cosat, sinat, h, h_new, accept,tau_new,tau_prev;
  double count = 0;

  for(int j = start; j< end; j++,count++ ){
    double kappa_j = kappa(j);
    double epsilon = epsilon_vec(j);
    double epsilon2 = epsilon2_vec(j);
    int leap = leap_vec(j);
    //reuse
    arma::vec temp_no = pow(acos(cos(tau_no(j)-beta)),2);
    //
    arma::vec ymat_col = ymat.col(j);
    tau_prev = tau_new = tau_yes(j);
    x = s_to_e(tau_prev,x);
    nu = arma::randn(p);
    nu = (dp - x * x.t()) * nu;

    h = likeli_tau_yes(tau_prev,beta,kappa_j,nr,ymat_col,temp_no) - as_scalar(0.5 * nu.t() * nu);
    for(int i = 0; i < leap; i++) {
      nu = nu + epsilon2 * gradient_tau_yes(tau_new, beta,kappa_j,nr,ymat_col,temp_no);
      nu = (dp - x * x.t()) * nu;

      alpha = norm(nu);
      ae = alpha * epsilon;
      cosat = cos(ae);
      sinat = sin(ae);
      x_temp = x;
      //geodesic flow//
      x = x * cosat + nu/alpha * sinat;
      nu = nu * cosat - alpha * x_temp * sinat;
      //
      tau_new = atan2(x(0),x(1));
      //
      nu = nu + epsilon2 * gradient_tau_yes(tau_new,beta,kappa_j,nr,ymat_col,temp_no);
      nu = nu = (dp - x * x.t()) * nu;

    }
    h_new = likeli_tau_yes(tau_new,beta,kappa_j,nr,ymat_col,temp_no) - as_scalar(0.5 * nu.t() * nu);
    accept = exp(h_new - h);
    if(accept < as_scalar(arma::randu(1))){
      e_s -= 1;
      tau_new = tau_prev;
      accept_chain(count) = 0;
    }
    yes_out(count) = tau_new;
  }
  return List::create(yes_out,e_s,accept_chain);
}
//#'@export
// [[Rcpp::export]]
List update_tau_no(int t,arma::vec nc_par, arma::vec epsilon_vec, arma::vec epsilon2_vec, arma::vec leap_vec, int nr,
                   arma::vec beta,arma::vec yes_out,arma::vec tau_no,arma::vec kappa,arma::mat ymat){
  int start = nc_par(t);
  int end = nc_par(t+1);
  int e_s = end - start;
  arma::vec no_out = arma::zeros(e_s);
  arma::vec accept_chain = arma::ones(e_s);

  arma::vec nu = arma::zeros(p);
  arma::vec x_temp = arma::zeros(p);
  arma::vec x = arma::zeros(p);
  double alpha, ae, cosat, sinat, h, h_new, accept,tau_new,tau_prev;
  double count = 0;

  for(int j = start; j< end; j++,count++ ){
    double kappa_j = kappa(j);
    double epsilon = epsilon_vec(j);
    double epsilon2 = epsilon2_vec(j);
    int leap = leap_vec(j);
    //reuse
    arma::vec temp_yes = pow(acos(cos(yes_out(j)-beta)),2);
    //
    arma::vec ymat_col = ymat.col(j);
    tau_prev = tau_new = tau_no(j);
    x = s_to_e(tau_prev,x);
    nu = arma::randn(p);
    nu = (dp - x * x.t()) * nu;

    h = likeli_tau_no(tau_prev,beta,kappa_j,nr,ymat_col,temp_yes) - as_scalar(0.5 * nu.t() * nu);
    for(int i = 0; i < leap; i++) {
      nu = nu + epsilon2 * gradient_tau_no(tau_new, beta,kappa_j,nr,ymat_col,temp_yes);
      nu = (dp - x * x.t()) * nu;

      alpha = norm(nu);
      ae = alpha * epsilon;
      cosat = cos(ae);
      sinat = sin(ae);
      x_temp = x;
      //geodesic flow//
      x = x * cosat + nu/alpha * sinat;
      nu = nu * cosat - alpha * x_temp * sinat;
      //
      tau_new = atan2(x(0),x(1));
      //
      nu = nu + epsilon2 * gradient_tau_no(tau_new,beta,kappa_j,nr,ymat_col,temp_yes);
      nu = nu = (dp - x * x.t()) * nu;

    }
    h_new = likeli_tau_no(tau_new,beta,kappa_j,nr,ymat_col,temp_yes) - as_scalar(0.5 * nu.t() * nu);
    accept = exp(h_new - h);
    if(accept < as_scalar(arma::randu(1))){
      tau_new = tau_prev;
      e_s -= 1;
      accept_chain(count) = 0;
    }
    no_out(count) = tau_new;
  }
  return List::create(no_out,e_s,accept_chain);
}
//#'@export
// [[Rcpp::export]]
List update_omega(double omega, arma::vec beta, int nr, double a, double b,double omega_sd){
  double eta = log(omega);
  double eta_new = omega_sd * as_scalar(arma::randn(1)) + eta;
  double omega_new = exp(eta_new);
  double accept;
  int accept_out = 1;

  accept = exp(likeli_omega(omega_new,beta,nr,a,b) + eta_new - likeli_omega(omega,beta,nr,a,b) - eta);
  if(accept < as_scalar(arma::randu(1))){
    omega_new = omega;
    accept_out = 0;
  }
  return List::create(omega_new,accept_out);
}
//#'@export
// [[Rcpp::export]]
List update_kappa(int t,arma::vec nc_par, int nr,arma::vec beta, arma::vec yes_out, arma::vec no_out,
                  arma::vec kappa,arma::mat ymat, double kappa_a, double ccc, arma::vec t_sig_vec){

  int start = nc_par(t);
  int end = nc_par(t+1);
  int e_s = end - start;
  arma::vec kappa_out = arma::zeros(e_s);
  arma::vec accept_chain = arma::ones(e_s);

  double count = 0;

  for(int j = start; j< end; j++,count++ ){
    arma::vec asset = sb_c_vec(pow(acos(cos(no_out(j)-beta)),2) - pow(acos(cos(yes_out(j)-beta)),2));
    double kappa_last = kappa(j);
    double t_sig = t_sig_vec(j);
    double tt = log(kappa_last);
    double tt_new = t_sig * as_scalar(arma::randn(1)) + tt;
    double kappa_new = exp(tt_new);
    double accept;
    arma::vec ymat_col = ymat.col(j);

    accept = exp(likeli_kappa(nr,asset,kappa_new,kappa_a,ccc,ymat_col) + tt_new -
      likeli_kappa(nr,asset,kappa_last,kappa_a,ccc,ymat_col) - tt);

    if(accept < as_scalar(arma::randu(1))){
      kappa_new = kappa_last;
      e_s -= 1;
      accept_chain(count) = 0;
    }
    kappa_out(count) = kappa_new;
  }
  return List::create(kappa_out,e_s,accept_chain);
}
//#'@export
// [[Rcpp::export]]
List waic_cpp(int t, arma::vec nc_par,int nr,arma::vec beta, arma::vec tau_yes, arma::vec tau_no, arma::vec kappa,arma::mat ymat){
  int start = nc_par(t);
  int end = nc_par(t+1);
  arma::mat amat = arma::zeros(nr,end-start);
  arma::mat amat_exp = amat;
  arma::mat amat_2 = amat;
  double count = 0;
  double temp = 0;
  for(int j = start; j< end; j++,count++){
    double kappa_j = kappa(j);
    arma::vec ymat_col = ymat.col(j);
    arma::vec asset = sb_c_vec(square(acos(cos(tau_no(j)-beta))) - square(acos(cos(tau_yes(j)-beta))));

    for(int i = 0 ;i < nr; i++){
      if(ymat_col(i) == 1){
        amat(i,count) = temp = betaInc_log(kappa_j,kappa_j,asset(i));
      }else{
        amat(i,count) = temp = betaInc_log_lower(kappa_j,kappa_j,asset(i));
      }
      amat_exp(i,count) = exp(temp);
      amat_2(i,count) = pow(temp,2);
    }
  }
  return List::create(amat,amat_exp,amat_2);
}
//#'@export
// [[Rcpp::export]]
arma::mat predict(int t, arma::vec nc_par,int nr,arma::vec beta, arma::vec tau_yes, arma::vec tau_no, arma::vec kappa,arma::mat ymat){
  int start = nc_par(t);
  int end = nc_par(t+1);
  arma::mat amat = arma::zeros(nr,end-start);
  arma::mat amat_exp = amat;
  arma::mat amat_2 = amat;
  double count = 0;
  double temp = 0;
  for(int j = start; j< end; j++,count++){
    double kappa_j = kappa(j);
    arma::vec ymat_col = ymat.col(j);
    arma::vec asset = sb_c_vec(square(acos(cos(tau_no(j)-beta))) - square(acos(cos(tau_yes(j)-beta))));

    for(int i = 0 ;i < nr; i++){
      temp = betaInc(kappa_j,kappa_j,asset(i));
      if(temp >= 0.5){
        amat(i,count) = 1;
      }
    }
  }
  return amat;
}


