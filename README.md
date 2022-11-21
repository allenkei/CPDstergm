# install the package
library(devtools)\
install_github("allenkei/CPDstergm")\
library(CPDstergm)

# reproduce the experiment
## SBM

num_node <- c(50,100,500)\
network_stats=c("edges", "mutual")\
result <- matrix(NA, nrow=3, ncol=4)

for(i in 1:3){\
&nbsp;&nbsp;&nbsp; set.seed(1)\
&nbsp;&nbsp;&nbsp; SBM_list <- sim_SBM_batch(num_nets=10, n=num_node[i], rho=0.0)\
&nbsp;&nbsp;&nbsp; sim_result <- CPD_STERGM_batch(SBM_list, network_stats)\
&nbsp;&nbsp;&nbsp; result[i,] <- colMeans(sim_result)\
}

## STERGM

num_node <- c(50,100,500)\
y1_target <- c(250,500,2500)\
network_stats=c("edges", "mutual")\
coefs_pos <- matrix(c(-1, -1, -1, -1, -2,  1, -2,  1), nrow=2, ncol=4, byrow = T)\
coefs_neg <- matrix(c( -1,  -1, -1,  -1, -2,  -1, -2,  -1), nrow=2, ncol=4, byrow = T)\
result <- matrix(NA, nrow=3, ncol=4)

for(i in 1:3){\
&nbsp;&nbsp;&nbsp; set.seed(1)\
&nbsp;&nbsp;&nbsp; STERGM_list <- sim_STERGM_batch(num_nets=10, n=num_node[i], network_stats, coefs_pos, coefs_neg, y1_stats=y1_target[i], node_attr=NA)\
&nbsp;&nbsp;&nbsp; sim_result <- CPD_STERGM_batch(STERGM_list, network_stats)\
&nbsp;&nbsp;&nbsp; esult[i,] <- colMeans(sim_result)\
}






