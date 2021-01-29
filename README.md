  This is an R package which implements the the circular factor model proposed in "Spatial voting models in circular spaces: A case study of the U.S. House of Representatives". 
  The roll call voting data can be found at voteview.com. The input data is an IxJ matrix where I represents the number of legislators (subjects) and J stands for the number of bills (items).
  # Install Package
  devtools::install_github('Xingchen-Yu/SLFM1D')\
  library(SLFM1D)\
  data(ymat)\
  set.seed(2021)\
  posterior_chain = SLFM(ymat = ymat,n_pos=200,burnin=50,core=1)
