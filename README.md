# R package: CPDstergm
This package detects multiple change points in time series of graphs, using Separable Temporal Exponential-family Random Graph Model (STERGM). The optimization problem with Group Fused Lasso regularization on the model parameters is solved by Alternating Direction Method of Multipliers (ADMM).

# The paper
Change Point Detection on A Separable Model for Dynamic Networks\
Yik Lun Kei\*, Hangjian Li\*, Yanzhen Chen, Oscar Hernan Madrid Padilla\
Transactions on Machine Learning Research (TMLR) 2025\
[PDF](https://arxiv.org/pdf/2303.17642.pdf)\
\* denotes equal contribution

# Install the package
library(devtools)\
install_github("allenkei/CPDstergm")\
library(CPDstergm)

# Demonstration
The code to reproduce the results in the paper can be found at:
https://github.com/allenkei/CPDstergm_demo
