# Hext : Hidden Markov model with Gaussian Mixture capability

 Gaussian mixture hidden Markov model for multiple sequence data with mixed-type extensions. Gaussian mixture emission is added. Forward-backward algorithm implemented in log-domain. This is a fork of the mhsmm package.

#Install
To install directly from github, open a terminal, type R, then

    devtools::install_github('htso/Hext')

#Dependencies
You need these packages on CRAN,    

    install.packages("mvtnorm", "MCMCpack", "MASS", "corpcor", "Rsolnp")

#Datasets
I include a couple of generic datasets, which are different flavor of gaussian mixture from low to high dimension. To load a dataset, just type

    data(simdat2single)

#Run
I provide a demo in the /demo subfolder. To run it, 

    demo("hext-demo", package="Hext")

#WARNING
The demo script may take a long time to finish.

#Platforms
Tested it on Linux (ubuntu 14.04).





