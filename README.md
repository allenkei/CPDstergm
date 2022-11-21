# install the package
library(devtools)\
install_github("allenkei/CPDstergm")\
library(CPDstergm)

# reproduce the experiment
num_node <- c(50,100,500)\
network_stats=c("edges", "mutual")\
result <- matrix(NA, nrow=3, ncol=4)

for(i in 1:3){\
&nbsp;&nbsp;&nbsp; set.seed(1)\
&nbsp;&nbsp;&nbsp; SBM_list <- sim_SBM_batch(num_nets=10, n=num_node[i], rho=0.9)\
&nbsp;&nbsp;&nbsp; sim_result <- CPD_STERGM_batch(SBM_list, network_stats)\
&nbsp;&nbsp;&nbsp; result[i,] <- colMeans(sim_result)\
}



