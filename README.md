  # Overview
  This is an R package that implements the the circular factor model proposed in "Spatial voting models in circular spaces: A case study of the U.S. House of Representatives". 
  The roll call voting data can be found at [here](http://voteview.com). The input data is an IxJ matrix where I represents the number of legislators (subjects) and J stands for the number of bills (items). For roll call voting data, it is recommended to initialize democrates between (-$\pi$/2,0) and republicans between (0,$\pi/2$) for faster convergence and mixing.
  # Install Package
  devtools::install_github('Xingchen-Yu/SLFM1D')\
  library(SLFM1D)\
  data(ymat)\
  set.seed(2021)\
  posterior_chain = SLFM(ymat = ymat,n_pos=200,burnin=50,core=1)
