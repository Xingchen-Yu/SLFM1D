#
# wrapper_beta = function(t){
#   update_beta(t,nr_par,delta,delta2,leap,nc,omega,cbeta_prior,beta,tau_yes,tau_no,kappa,ymat)
# }
#
# wrapper_yes = function(t){
#   update_tau_yes(t,nc_par,delta_yes,delta2_yes,leap_tau,nr,beta,tau_yes,tau_no,kappa,ymat)
# }
#
# wrapper_no = function(t){
#   update_tau_no(t,nc_par,delta_no,delta2_no,leap_tau,nr,beta,tau_yes,tau_no,kappa,ymat)
# }
#
# wrapper_kappa = function(t){
#   update_kappa(t,nc_par,nr,beta,tau_yes,tau_no,kappa,ymat,kappa_a,ccc,t_sig)
# }
#
# wrapper_waic = function(t){
#   waic_cpp(t,nc_par,nr,beta,tau_yes,tau_no,kappa,ymat)
# }
#
# wrapper_predict = function(t){
#   predict(t,nc_par,nr,beta,tau_yes,tau_no,kappa,ymat)
# }

update_tsig = function(upper,lower,t_sig,check,nnn){
  l_s = which(check<=lower)
  u_s = which(check>=upper)

  t_sig[l_s] = t_sig[l_s] * 0.95
  t_sig[u_s] = t_sig[u_s] * 1.05
  return(list(t_sig,length(c(l_s,u_s))/nnn))
}
