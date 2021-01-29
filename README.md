  This is an R package which implements the the circular factor model proposed in "Spatial voting models in circular spaces: A case study of the U.S. House of Representatives". 
  # Install Package
  devtools::install_github('Xingchen-Yu/SLFM1D')\
  library(SLFM1D)\
  data(ymat)\
  set.seed(2021)\
  posterior_chain = SLFM(ymat = ymat,n_pos=200,burnin=50,core=1)
