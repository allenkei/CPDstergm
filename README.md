# install the package
library(devtools)\
install_github("allenkei/CPDstergm")\
library(CPDstergm)

# replicate the experiment
set.seed(1)\
SBM_list <- sim_SBM_batch(num_nets=10, n=50, rho=0)\
network_stats=c("edges", "mutual")\
sim1_result1 <- CPD_STERGM_batch(SBM_list, network_stats)\
colMeans(sim1_result1)
