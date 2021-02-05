#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]

//[[Rcpp::depends(RcppArmadillo)]]
const double pi2 = pow(M_PI,2);
const double pi22 = 2 * pow(M_PI,2);
const double pi22_i = 1/pi22;
const double log_pi = log(M_PI);
const int MAXIT = 1000;
const double EPS = 3e-12;
const double FPMIN = 1e-30;
const double p = 2;
const arma::mat dp = arma::eye(p,p);
const double mu = 0;
const double log2_p_pi = log(2) + log(M_PI);

arma::vec s_to_e(double angle, arma::vec x){
  x(0) = sin(angle);
  x(1) = cos(angle);
  return x;
}

double betaContFrac(double a, double b, double x) {

  double qab = a + b;
  double qap = a + 1;
  double qam = a - 1;
  double c = 1;
  double d = 1 - qab * x / qap;
  if (fabs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  double h = d;
  int m;
  for (m = 1; m <= MAXIT; m++) {
    int m2 = 2 * m;
    double aa = m * (b-m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    h *= (d * c);
    aa = -(a+m) * (qab+m) * x / ((a+m2) * (qap+m2));
    d = 1 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    double del = d*c;
    h *= del;
    if (fabs(del - 1) < EPS) break;
  }
  return h;
}

double betaInc(double a, double b, double x) {
  if (x == 0)
    return 0;
  else if (x == 1)
    return 1;
  else {
    double logBeta = lgamma(a+b) - lgamma(a) - lgamma(b)+ a * log(x) + b * log(1-x);
    if (x < (a+1) / (a+b+2))
      return exp(logBeta) * betaContFrac(a, b, x) / a;
    else
      return 1 - exp(logBeta) * betaContFrac(b, a, 1-x) / b;
  }
}

double betaInc_log(double a, double b, double x) {
  if (x == 1)
    return 0;
  else {
    double logBeta = lgamma(a+b) - lgamma(a) - lgamma(b)+ a * log(x) + b * log(1-x);
    if (x < (a+1) / (a+b+2))
      return logBeta + log(betaContFrac(a,b,x))-log(a);
    else
      return log(1 - exp(logBeta) * betaContFrac(b, a, 1-x) / b);
  }
}

double betaInc_log_lower(double a, double b, double x) {
  if (x == 0)
    return 0;
  else {
    double logBeta = lgamma(a+b) - lgamma(a) - lgamma(b)+ a * log(x) + b * log(1-x);
    if (1-x < (b+1) / (a+b+2))
      return logBeta + log(betaContFrac(b,a,1-x))-log(b);
    else
      return log(1 - exp(logBeta) * betaContFrac(a, b, x) / a);
  }
}

double bessi0_exp(double x)
{
  double ax,ans;
  double y;
  if ((ax=fabs(x)) < 3.75) {
    y=x/3.75;
    y*=y;
    ans=(1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
                                            +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))))*exp(-ax);
  } else {
    y=3.75/ax;
    ans=(1/sqrt(ax))*(0.39894228+y*(0.1328592e-1
                                      +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
                                                                           +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
                                                                           +y*0.392377e-2))))))));
  }
  return ans;
}

double dgamma_log (double x, double a, double b){
  return a * log(b) - lgamma(a) + (a - 1) * log(x) - b * x;
}

double likeli_omega(double omega, arma::vec beta,int nr,double a, double b){
  return  - nr * (log2_p_pi + log(bessi0_exp(omega))) +
    omega * sum(cos(beta)-1) + dgamma_log (omega, a, b);
}

double sb_c(double x){
  return (x + pi2)/pi22;
}
arma::vec sb_c_vec(arma::vec x){
  return (x + pi2)/pi22;
}


arma::vec Rf_rbinom_vec(arma::vec x,arma::vec shape1,const int n){
  arma::vec kobe = arma::zeros(n);
  double temp,shape;
  for(int i=0 ;i < n; i++){
    shape = shape1(i);
    temp = betaInc(shape,shape,x(i));
    //kobe(i) = Rf_rbinom(1,temp);
    if(as_scalar(arma::randu(1))<= temp){
      kobe(i) = 1;
    }else{
      kobe(i) = 0;
    }
  }
  return kobe;
}

double likeli_beta (double omega, double beta, arma::vec tau_yes,
                    arma::vec tau_no,arma::vec kappa ,int nc,arma::rowvec ymat_row){

  const arma::vec asset = sb_c_vec(square(acos(cos(tau_no-beta))) - square(acos(cos(tau_yes-beta))));
  arma::vec avec = arma::zeros(nc);
  double shape;

  for(int i = 0 ;i < nc; i++){
    shape = kappa(i);
    if(ymat_row(i) == 1){
      avec(i) = betaInc_log(shape,shape,asset(i));
    }else{
      avec(i) = betaInc_log_lower(shape,shape,asset(i));
    }
  }
  return sum(avec) + omega * cos(beta-mu);
}


arma::vec gradient_beta(arma::vec cbeta_prior, double beta, arma::vec tau_yes, arma::vec tau_no,arma::vec kappa,
                        const int nc, arma::rowvec ymat_row){

  arma::vec temp_x = arma::zeros(p);
  arma::vec please = arma::zeros(p);

  double temp,shape,shape2,logbeta,check,beta_tau_yes,beta_tau_no,avec,buffer_yes,buffer_no,yes_temp,no_temp;
  for(int i = 0 ;i < nc; i++){
    shape = kappa(i);
    yes_temp = tau_yes(i);
    no_temp = tau_no(i);
    beta_tau_yes = yes_temp - beta;
    beta_tau_no = no_temp - beta;
    buffer_yes = acos(cos(beta_tau_yes));
    buffer_no = acos(cos(beta_tau_no));
    temp = sb_c(pow(buffer_no,2)-pow(buffer_yes,2));


    shape2 = 2 * shape;
    check = (shape + 1) / (shape2 +2);
    if(ymat_row(i) == 1){
      if (temp < check){
        avec = shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec = pi22_i/( temp * (1-temp) *(exp(logbeta)- betaContFrac(shape,shape,1-temp)/shape));
      }
    }else{
      if (1 - temp < check){
        avec = -shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,1-temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec =  pi22_i/( temp * (1-temp) *(betaContFrac(shape,shape,temp)/shape - exp(logbeta)));
      }
    }
    please += 2 * avec * (buffer_yes/fabs(sin(beta_tau_yes)) * s_to_e(yes_temp,temp_x) -
      buffer_no/fabs(sin(beta_tau_no)) * s_to_e(no_temp,temp_x));
  }

  return please + cbeta_prior;
}


double likeli_tau_yes(double tau_yes,arma::vec beta,double shape,const int nr,
                      arma::vec ymat_col,arma::vec temp_no){

  const arma::vec asset = sb_c_vec(temp_no - square(acos(cos(tau_yes-beta))));
  arma::vec avec = arma::zeros(nr);

  for(int i = 0 ;i < nr; i++){
    if(ymat_col(i) == 1){
      avec(i) = betaInc_log(shape,shape,asset(i));
    }else{
      avec(i) = betaInc_log_lower(shape,shape,asset(i));
    }
  }
  return sum(avec);
}


arma::vec gradient_tau_yes(double tau_yes,arma::vec beta,double shape,const int nr,
                           arma::vec ymat_col,arma::vec temp_no){

  arma::vec temp_x = arma::zeros(p);
  arma::vec please = arma::zeros(p);

  double temp,shape2,logbeta,check,beta_tau_yes,avec,buffer_yes,beta_temp;
  for(int i = 0 ;i < nr; i++){
    beta_temp = beta(i);
    beta_tau_yes = tau_yes - beta_temp;
    buffer_yes = acos(cos(beta_tau_yes));
    temp = sb_c(temp_no(i) - pow(buffer_yes,2));

    shape2 = 2 * shape;
    check = (shape + 1) / (shape2 +2);
    if(ymat_col(i) == 1){
      if (temp < check){
        avec = shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec = pi22_i/( temp * (1-temp) *(exp(logbeta)- betaContFrac(shape,shape,1-temp)/shape));
      }
    }else{
      if (1 - temp < check){
        avec = -shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,1-temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec =  pi22_i/( temp * (1-temp) *(betaContFrac(shape,shape,temp)/shape - exp(logbeta)));
      }
    }
    please += 2 * avec * buffer_yes/fabs(sin(beta_tau_yes)) * s_to_e(beta_temp,temp_x);
  }

  return please;
}

double likeli_kappa(int nr, arma::vec asset,double shape,double kappa_a, double ccc,arma::vec ymat_col){
  arma::vec avec = arma::zeros(nr);

  for(int i = 0 ;i < nr; i++){
    if(ymat_col(i) == 1){
      avec(i) = betaInc_log(shape,shape,asset(i));
    }else{
      avec(i) = betaInc_log_lower(shape,shape,asset(i));
    }
  }
  return sum(avec) + dgamma_log(shape,kappa_a,ccc);
}

double likeli_tau_no(double tau_no,arma::vec beta,double shape,const int nr,
                     arma::vec ymat_col,arma::vec temp_yes){

  const arma::vec asset = sb_c_vec(square(acos(cos(tau_no-beta))) - temp_yes);
  arma::vec avec = arma::zeros(nr);

  for(int i = 0 ;i < nr; i++){
    if(ymat_col(i) == 1){
      avec(i) = betaInc_log(shape,shape,asset(i));
    }else{
      avec(i) = betaInc_log_lower(shape,shape,asset(i));
    }
  }
  return sum(avec);
}


arma::vec gradient_tau_no(double tau_no,arma::vec beta,double shape,const int nr,
                          arma::vec ymat_col,arma::vec temp_yes){

  arma::vec temp_x = arma::zeros(p);
  arma::vec please = arma::zeros(p);

  double temp,shape2,logbeta,check,beta_tau_no,avec,buffer_no,beta_temp;
  for(int i = 0 ;i < nr; i++){
    beta_temp = beta(i);
    beta_tau_no = tau_no - beta_temp;
    buffer_no = acos(cos(beta_tau_no));
    temp = sb_c(pow(buffer_no,2) - temp_yes(i));

    shape2 = 2 * shape;
    check = (shape + 1) / (shape2 +2);
    if(ymat_col(i) == 1){
      if (temp < check){
        avec = shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec = pi22_i/( temp * (1-temp) *(exp(logbeta)- betaContFrac(shape,shape,1-temp)/shape));
      }
    }else{
      if (1 - temp < check){
        avec = -shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,1-temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec =  pi22_i/( temp * (1-temp) *(betaContFrac(shape,shape,temp)/shape - exp(logbeta)));
      }
    }
    please += - 2 * avec * buffer_no/fabs(sin(beta_tau_no)) * s_to_e(beta_temp,temp_x);
  }

  return please;
}
