

update_tsig = function(upper,lower,t_sig,check,nnn){
  l_s = which(check<=lower)
  u_s = which(check>=upper)

  t_sig[l_s] = t_sig[l_s] * 0.95
  t_sig[u_s] = t_sig[u_s] * 1.05
  return(list(t_sig,length(c(l_s,u_s))/nnn))
}
