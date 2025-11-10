# R package: CPDstergm
This package detects multiple change points in time series of graphs, using Separable Temporal Exponential-family Random Graph Model (STERGM). The optimization problem with Group Fused Lasso regularization on the model parameters is solved by Alternating Direction Method of Multipliers (ADMM).

# The paper
Change Point Detection on A Separable Model for Dynamic Networks\
Yik Lun Kei\*, Hangjian Li\*, Yanzhen Chen, Oscar Hernan Madrid Padilla\
Transactions on Machine Learning Research (TMLR) 2025\
[PDF](https://arxiv.org/pdf/2303.17642.pdf)\
Funded by NSF DMS-2015489\
\* denotes equal contribution

# Install the package in R
```r
library(devtools)
install_github("allenkei/CPDstergm")
library(CPDstergm)
```

# Demonstration
```r
# using the DJIA data from ecp package
library(ecp)
data(DJIA)
market <- DJIA$market
date_vec <- DJIA$dates[1:1138]
rownames(market) <- date_vec

# construct the network data
start <- which(date_vec == '2007-01-01')
end <- which(date_vec == '2010-01-04')
date_range <- start:end
mydata <- array(NA,c(29, 29, length(date_range)))
df <- list()

for(i in 1:length(date_range)){
  temp <- market[date_range[i]:(date_range[i]+3),]
  temp <- ifelse(cor(temp)< 0, 1, 0)
  diag(temp) <- 0
  df[[i]] <- temp
}

# detect change points
result <- CPD_STERGM(df, directed=FALSE, network_stats=c("edges", "triangles"), list_of_lambda=10^c(2:5))

seq_date <- rownames(market[start:end,])
seq_date[result[[1]]] # detected change points (dates)
```


# Results in paper
The code to reproduce the results in the paper can be found at:
https://github.com/allenkei/CPDstergm_demo
