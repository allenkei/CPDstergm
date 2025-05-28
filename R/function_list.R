## usethis namespace: start
#' @useDynLib CPDstergm, .registration = TRUE
## usethis namespace: end
NULL


#' @import network ndtv dplyr tergm network networkDynamic ecp
#' @importFrom stats as.formula median qnorm rbinom sd simulate complete.cases
#' @importFrom ergm ergmMPLE san
#' @importFrom Rcpp sourceCpp
NULL











#' Simulation of dynamic networks from the Stochastic Block Model
#' @description This function simulates multiple sequences of dynamic networks, from the Stochastic Block Model (SBM).
#' @param num_seq The number of sequences of dynamic networks.
#' @param n The number of node.
#' @param rho The correlation coefficient between consecutive time points for a dyad.
#'
#' @return Returns a list of multiple sequences of dynamic networks.
#' @export
#'
#' @examples
#' set.seed(1)
#' SBM_list <- sim_SBM_list(num_seq=1, n=50, rho=0)
sim_SBM_list <- function(num_seq=1, n=50, rho=0){

  output <- list()
  for(rep_iter in 1:num_seq){

    K = 3
    v =  c(26, 51, 76)
    num_time = 100
    df = vector(mode = "list", length = num_time)

    for(t in 1:num_time){

      if( t==1 || t==v[2] ){ # t=1 or t=51

        P =  matrix(0.3,n,n)
        P[1:floor(n/K), 1:floor(n/K)] = 0.5
        P[(1+floor(n/K)):(2*floor(n/K)),(1+floor(n/K)):(2*floor(n/K)) ] = 0.5
        P[(1+2*floor(n/K)):n,(1+2*floor(n/K)):n ] = 0.5
        diag(P) = 0

        A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)

      }

      if(t == v[1] || t == v[3]){ # t=26 or t=76

        Q =  matrix(0.2,n,n)
        Q[1:floor(n/K), 1:floor(n/K)] = 0.45
        Q[(1+floor(n/K)):(2*floor(n/K)),(1+floor(n/K)):(2*floor(n/K)) ] = 0.45
        Q[(1+2*floor(n/K)):n,(1+2*floor(n/K)):n ] = 0.45
        diag(Q) = 0

        A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),Q),n,n)

      }

      if( (t > 1 && t < v[1])  ||  (t > v[2] && t < v[3]) ){ # t=2 to t=25 or t=52 to t=75

        aux1 = (1-P)*rho + P # (1-E(t+1))*rho + E(t+1) if A(t) = 1
        aux2 = P*(1-rho) # (1-E(t+1))*rho + E(t+1) if A(t) = 0

        aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
        aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
        A =  aux1*A + aux2*(1-A)

      }

      if( (t > v[1] && t < v[2]) || ((t > v[3] && t <= num_time)) ){ # t=27 to t=50 or t=77 to t=100

        aux1 = (1-Q)*rho + Q
        aux2 = Q*(1-rho)

        aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
        aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
        A =  aux1*A + aux2*(1-A)

      }

      diag(A) <- 0
      df[[t]] = A

    }

    output[[rep_iter]] <- df

  }

  return(output)
}









#' Simulation of dynamic networks from the Separable Temporal Exponential-family Random Graph Model
#' @description This function simulates multiple sequences of dynamic networks, from the Separable Temporal Exponential-family Random Graph Model (STERGM).
#' @param num_seq The number of sequences of dynamic networks.
#' @param n The number of node.
#' @param network_stats The network statistics for both formation and dissolution models. See search.ergmTerms() for a comprehensive list from the library(ergm).
#' @param coefs_pos The coefficients for the formation model.
#' @param coefs_neg The coefficients for the dissolution model.
#' @param y1_stats The number of edges for the first network.
#' @param node_attr The nodal attributes.
#'
#' @return Returns a list of multiple sequences of dynamic networks
#' @export
#'
#' @examples
#' coefs_pos <- matrix(c(-1, -1, -1, -1,
#'                       -2,  1, -2,  1), nrow=2, ncol=4, byrow=TRUE)
#' coefs_neg <- matrix(c(-1, -1, -1, -1,
#'                       -2, -1, -2, -1), nrow=2, ncol=4, byrow=TRUE)
#'
#' set.seed(1)
#' STERGM_list <- sim_STERGM_list(num_seq=1, n=50, network_stats=c("edges", "mutual"),
#'                                coefs_pos, coefs_neg, y1_stats=250, node_attr=NA)
sim_STERGM_list <- function(num_seq=1, n=50, network_stats,
                             coefs_pos, coefs_neg, y1_stats, node_attr=NA){
  num_nodes <- n
  num_timepts <- 100
  num_changepts <- 3
  form_model <- diss_model <- as.formula(paste("~", paste(network_stats, collapse = "+")))

  output <- list()
  for(rep_iter in 1:num_seq){

    # generate initial network
    g0<-network.initialize(num_nodes, directed=T) # empty
    g1<-san(g0~edges,target.stats=y1_stats, verbose=TRUE)
    if(any(!is.na(node_attr))){
      network::set.vertex.attribute(g1, "node_attr", node_attr)
    }
    ginit <- g1

    time_stable <- num_timepts / (num_changepts+1)
    cur_end <- 0
    res_adj_list <- vector(mode = 'list', length = num_timepts)
    for(i in 1:(num_changepts+1)){
      cat('[INFO] Simulate from ', cur_end+1, ' to ', cur_end + time_stable, '\n')
      stergm.sim <- simulate(
        ginit,
        formation=form_model,
        dissolution=diss_model,
        coef.form=coefs_pos[,i],
        coef.diss=coefs_neg[,i],
        nsim = 1,
        time.slices = time_stable,
        time.start = cur_end # needs to be update every time a new dynamic is generated
      )
      # newstart_nw%n%'net.obs.period' check https://github.com/statnet/tergm/blob/master/R/simulate.stergm.R for details

      for(t in (1 + cur_end) : (time_stable + cur_end)){
        tmpnet <- network.extract(stergm.sim, at = t) %>% network()
        res_adj_list[[t]] <- tmpnet[, ]
      }
      cur_end <- cur_end + time_stable
      ginit <- network.extract(stergm.sim, at = cur_end) %>% network()

    }

    output[[rep_iter]] <- res_adj_list
  }

  return(output)
}








#' Simulation of dynamic networks from the Random Dot Product Graph Model
#' @description This function simulates multiple sequences of dynamic networks, from the Random Dot Product Graph Model (RDPGM).
#' @param num_seq The number of sequences of dynamic networks.
#' @param n The number of node.
#' @param rho The correlation coefficient between consecutive time points for latent position.
#' @param d The latent dimension.
#'
#' @return Returns a list of multiple sequences of dynamic networks.
#' @export
#'
#' @examples
#' set.seed(1)
#' RDPG_list <- sim_RDPG_list(num_seq = 1, n = 50, rho = 0.9, d = 10)
sim_RDPG_list <- function(num_seq = 1, n = 50, rho = 0.9, d = 10) {

  change_points <- c(26, 51, 76)
  lower_bound <- c(0, 1, 0, 1) / 16
  upper_bound <- c(1, 2, 1, 2) / 16
  num_time <- 100
  sigma <- 1
  output <- list()

  generate_latent <- function(n, d) {
    latent <- matrix(rnorm(n * d, mean = 1, sd = 1), n, d)
    return(latent)
  }

  segments <- c(1, change_points, num_time + 1)
  num_segments <- length(segments) - 1

  for (rep_iter in 1:num_seq) {

    dynamic_adj <- list()

    Xt <- generate_latent(n, d)
    Yt <- generate_latent(n, d)

    for (s in 1:num_segments) {
      t_start <- segments[s]
      t_end <- segments[s + 1] - 1

      W <- matrix(runif(d * d, lower_bound[s], upper_bound[s]), d, d)

      for (t in t_start:t_end) {

        if (t > 1) {
          noise_X <- matrix(rnorm(n * d, 0, sigma), n, d)
          Xt <- rho * Xt + (1 - rho) * noise_X

          noise_Y <- matrix(rnorm(n * d, 0, sigma), n, d)
          Yt <- rho * Yt + (1 - rho) * noise_Y
        }

        Xt <- Xt / (sqrt(rowSums(Xt^2)) + 1e-5)
        Yt <- Yt / (sqrt(rowSums(Yt^2)) + 1e-5)

        P <- Xt %*% W %*% t(Yt)
        P[P < 0] <- 1e-5
        P[P > 1] <- 1 - 1e-5

        A <- matrix(rbinom(n * n, size = 1, prob = as.vector(P)), n, n)
        diag(A) <- 0
        dynamic_adj[[t]] <- A
      }
    }

    output[[rep_iter]] <- dynamic_adj
  }

  return(output)
}











#' Generation of input data from dynamic networks
#' @description This function generates the change statistics with formation and dissolution networks, and save them as input data.
#' @param y_data A sequence of dynamic networks.
#' @param directed Whether the networks are directed or not.
#' @param network_stats The network statistics for both formation and dissolution models. See search.ergmTerms() for a comprehensive list from the library(ergm).
#' @param node_attr The nodal attributes.
#'
#' @return Returns a list of change statistics and dynamic networks for both formation and dissolution models.
#' @export
#'
#' @examples
#' set.seed(1)
#' SBM_list <- sim_SBM_list(num_seq=1, n=50, rho=0)
#' input_data <- save_H_y_list(SBM_list[[1]], directed=TRUE, network_stats=c("edges", "mutual"))
save_H_y_list <- function(y_data, directed, network_stats, node_attr=NA){

  output <- list()
  num_time <- length(y_data)
  n <- dim(y_data[[1]])[1]

  H_pos_list <- H_neg_list <- y_pos_list <- y_neg_list <- list()

  f_pos <- as.formula(paste("y_pos", paste(network_stats, collapse = "+"), sep = "~"))
  f_neg <- as.formula(paste("y_neg", paste(network_stats, collapse = "+"), sep = "~"))

  for(t in 2:num_time){
    y_before <- y_data[[t-1]] # matrix
    y_after <- y_data[[t]]    # matrix

    # network object with directed/undirected
    y_pos = network(matrix(as.integer(y_before | y_after), nrow=n, ncol=n, byrow=F), directed = directed)
    y_neg = network(matrix(as.integer(y_before & y_after), nrow=n, ncol=n, byrow=F), directed = directed)

    # assign nodal attributes to network object
    if(any(!is.na(node_attr))){
      network::set.vertex.attribute(y_pos, "node_attr", node_attr)
      network::set.vertex.attribute(y_neg, "node_attr", node_attr)
    }

    # calculate change statistics
    y_pos = ergmMPLE(f_pos, output='array')
    y_neg = ergmMPLE(f_neg, output='array')

    H_pos_output <- c()
    for(i in 1:length(network_stats)){

      this_stats <- y_pos$predictor[,,i]

      # adjust the change statistics
      # when y_ij^t-1 is 1, no effect on change statistics of y^+, hence set it to 0
      # this_stats <- ifelse( (y_before==1) & (!is.na(this_stats)), 0, this_stats) # matrix

      # stack matrix row by row into one vector
      H_pos_output <- cbind(H_pos_output, c(t(this_stats)))
    }


    H_neg_output <- c()
    for(i in 1:length(network_stats)){

      this_stats <- y_neg$predictor[,,i]

      # adjust the change statistics
      # when y_ij^t-1 is 0, no effect on change statistics of y^-, hence set it to 0
      # this_stats <- ifelse( (y_before==0) & (!is.na(this_stats)), 0, this_stats) # matrix

      # stack matrix row by row into one vector
      H_neg_output <- cbind(H_neg_output, c(t(this_stats)))
    }

    # store network in same order for dyad ij
    y_pos_output <- cbind(c(t(y_pos$response)))
    y_neg_output <- cbind(c(t(y_neg$response)))

    if(directed){
      H_pos_output[is.na(H_pos_output)] <- 0 # diagonal NA becomes 0
      H_neg_output[is.na(H_neg_output)] <- 0
      y_pos_output[is.na(y_pos_output)] <- 0
      y_neg_output[is.na(y_neg_output)] <- 0
    }else{
      H_pos_output <- H_pos_output[complete.cases(H_pos_output), ] # upper triangles
      H_neg_output <- H_neg_output[complete.cases(H_neg_output), ]
      y_pos_output <- y_pos_output[complete.cases(y_pos_output), ]
      y_neg_output <- y_neg_output[complete.cases(y_neg_output), ]
    }

    H_pos_list[[t-1]] <- H_pos_output
    H_neg_list[[t-1]] <- H_neg_output
    y_pos_list[[t-1]] <- y_pos_output
    y_neg_list[[t-1]] <- y_neg_output

  }

  output[[1]] <- H_pos_list
  output[[2]] <- H_neg_list
  output[[3]] <- y_pos_list
  output[[4]] <- y_neg_list
  return(output)
}








Evaluation_ADMM_list <- function(H_pos_list, H_neg_list, y_pos_list, y_neg_list,
                                 theta_mat, tau, p1, p2, p, n, true_CP){

  log_lik <- cal_log_likelihood(H_pos_list, H_neg_list, y_pos_list, y_neg_list, theta_mat, tau, p1, p2)

  theta_change <- numeric(tau-1)
  for(i in 1:(tau-1)){theta_change[i] <- norm(theta_mat[i+1,]-theta_mat[i,],"2")}
  t_change <- (theta_change - median(theta_change))/sd(theta_change) # use median
  threshold <- mean(t_change) + qnorm(0.9, lower.tail = T) * sd(t_change) # use mean
  theta_change <- t_change

  est_CP <- c()
  for(i in 1:(tau-1)){
    if(theta_change[i] > threshold & i > 5 & i < (tau-5)) {
      est_CP <- c(est_CP, i) # location of change point
    }
  }

  #min-spacing
  end_i <- 2
  while(end_i <= length(est_CP)){
    prev <- est_CP[end_i-1]
    this <- est_CP[end_i]

    if(this - prev > 5){
      end_i <- end_i + 1
    }else{
      selection <- c(prev,this)
      to_remove <- selection[which.min(theta_change[selection])]
      est_CP <- est_CP[-which(est_CP == to_remove)]
    }
  }

  est_CP <- est_CP + 2 # the 1st theta_diff indicates t=3 is change point
  num_CP <- length(est_CP)

  gt_CP_corrected <- c(1, true_CP, tau+2) # tau+2 = T+1
  est_CP_corrected <- c(1, est_CP, tau+2)

  gt_list <- est_list <- list();
  for(i in 2:length(gt_CP_corrected)){
    gt_list[[i-1]] <- gt_CP_corrected[i-1]:(gt_CP_corrected[i]-1)
  }
  for(i in 2:length(est_CP_corrected)){
    est_list[[i-1]] <- est_CP_corrected[i-1]:(est_CP_corrected[i]-1)
  }


  if(num_CP == 0){
    dist_est_gt <- Inf
    dist_gt_est <- -Inf
    covering_metric <- 0
  }else{

    holder <- c()
    for(i in true_CP){
      dist_diff <- c()
      for(j in est_CP){dist_diff <- c(dist_diff, abs(j-i))}
      holder <- c(holder, min(dist_diff))
    }
    dist_est_gt <- max(holder)

    holder <- c()
    for(i in est_CP){
      dist_diff <- c()
      for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
      holder <- c(holder, min(dist_diff))
    }
    dist_gt_est <- max(holder)

    covering_metric <- 0
    for(i in 1:length(gt_list)){
      A <- gt_list[[i]]
      jaccard <- c()
      for(j in 1:length(est_list)){
        A_prime <- est_list[[j]]
        jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
      }
      covering_metric <- covering_metric + length(A)*max(jaccard)
    }
    covering_metric <- covering_metric/(tau+1) # tau+1 = T

  }

  BIC = -2*log_lik + p*log( (tau+1)*choose(n,2)*2 )*(num_CP+1) # N = choose(n,2)*2 for directed network
  #AIC = -2*log_lik + 2*p*(num_CP+1)

  abs_error <- abs(num_CP - length(true_CP))

  output <- list()
  output[[1]] <- est_CP
  output[[2]] <- BIC
  output[[3]] <- abs_error
  output[[4]] <- dist_est_gt
  output[[5]] <- dist_gt_est
  output[[6]] <- covering_metric
  output[[7]] <- threshold

  return(output)
}










Evaluation_ADMM <- function(H_pos_list, H_neg_list, y_pos_list, y_neg_list,
                            theta_mat, tau, p1, p2, p, n, directed, threshold_alpha){
  output <- list()

  log_lik <- cal_log_likelihood(H_pos_list, H_neg_list, y_pos_list, y_neg_list, theta_mat, tau, p1, p2)

  theta_change <- numeric(tau-1)
  for(i in 1:(tau-1)){theta_change[i] <- norm(theta_mat[i+1,]-theta_mat[i,],"2")}
  med <- median(theta_change); std <- sd(theta_change)
  t_change <- (theta_change - med)/std # normalize
  threshold <- mean(t_change) + qnorm( 1-threshold_alpha, lower.tail = T) * sd(t_change)
  theta_change <- t_change

  est_CP <- c()
  for(i in 1:(tau-1)){
    if(theta_change[i] > threshold & i > 10 & i < (tau-10)) {
      est_CP <- c(est_CP, i) # location of change point
    }
  }

  #min-spacing
  end_i <- 2
  while(end_i <= length(est_CP)){
    prev <- est_CP[end_i-1]
    this <- est_CP[end_i]

    if(this - prev > 5){
      end_i <- end_i + 1
    }else{
      selection <- c(prev,this)
      to_remove <- selection[which.min(theta_change[selection])]
      est_CP <- est_CP[-which(est_CP == to_remove)]
    }
  }

  est_CP <- est_CP + 2 # the 1st theta_diff indicates t=3 is change point
  num_CP <- length(est_CP)

  if(directed){
    net_size <- 2*choose(n,2)
  }else{
    net_size <- choose(n,2)
  }

  BIC = -2*log_lik + p*log( (tau+1)*net_size )*(num_CP+1)
  theta_change[1:10] <- 0; theta_change[(length(theta_change)-9):length(theta_change)] <- 0

  output[[1]] <- est_CP
  output[[2]] <- BIC
  output[[3]] <- theta_change
  output[[4]] <- threshold
  output[[5]] <- log_lik

  return(output)
}











#' Change point detection and performance evaluation for dynamic networks with known ground truth
#' @description This function detects the change points for each sequence of dynamic networks, and it evaluates the performance by calculating the metrics after model selection. The ground truth for the change points must be provided.
#' @param data_list The list of multiple sequences of dynamic networks.
#' @param network_stats The network statistics for both formation and dissolution models. See search.ergmTerms() for a comprehensive list from the library(ergm).
#' @param directed Whether the networks are directed or not.
#' @param node_attr The nodal attributes.
#' @param list_of_lambda The list of tuning parameters \eqn{\lambda}.
#' @param ADMM_iteration The learning iteration for ADMM.
#' @param theta_iter The learning iteration for \eqn{\bm{\theta}}.
#' @param z_iter The learning iteration for \eqn{\bm{z}}.
#' @param theta_tol The tolerance for \eqn{\bm{\theta}} stopping criteria.
#' @param ADMM_tol The tolerance for ADMM stopping criteria.
#' @param true_CP The true change points.
#' @param update_alpha If it is TRUE, the alpha is updated with a schedule.
#' @param verbose If it is TRUE, the information at each ADMM iteration is printed.
#'
#' @return Returns the four evaluation metrics for each sequence of dynamic networks. The four evaluation metrics are (1) absolute error \eqn{|\hat{K}-K|}, (2) one-sided Hausdorff distance \eqn{d(\hat{\mathcal{C}}|\mathcal{C})}, (3) one-sided Hausdorff distance \eqn{d(\mathcal{C}|\hat{\mathcal{C}})}, and (4) coverage metric \eqn{C(\mathcal{G},\mathcal{G'})}.
#' @export
#'
#' @examples
#' set.seed(1)
#' SBM_list <- sim_SBM_list(num_seq=5, n=50, rho=0)
#' network_stats=c("edges", "mutual")
#' result <- CPD_STERGM_list(SBM_list, directed=TRUE, network_stats, list_of_lambda=10^c(4:7))
CPD_STERGM_list <- function(data_list, directed, network_stats, node_attr=NA,
                            list_of_lambda=10^c(-2:7), ADMM_iteration=200, theta_iter=20, z_iter=20,
                            theta_tol=1e-3, ADMM_tol=1e-7, true_CP=c(26, 51, 76), update_alpha=TRUE, verbose=FALSE){
  num_nets <- length(data_list)
  temp <- list()
  BIC_table <- matrix(NA, nrow=num_nets, ncol=length(list_of_lambda))
  output_table <- matrix(NA, nrow=num_nets, ncol=4) # 4 metrics
  p1 <- p2 <- length(network_stats); p <- p1+p2

  for(rep_iter in 1:num_nets){

    temp[[rep_iter]] <- list()
    y_data <- data_list[[rep_iter]]
    input_data <- save_H_y_list(y_data, directed, network_stats, node_attr)
    H_pos_list <- input_data[[1]]
    H_neg_list <- input_data[[2]]
    y_pos_list <- input_data[[3]]
    y_neg_list <- input_data[[4]]

    n <- dim(y_data[[1]])[1]
    tau <- length(y_data)-1
    d_vec <- numeric(tau-1) # fixed
    for(i in 1:(tau-1)){ d_vec[i] = sqrt( tau/ (i*(tau-i)) ) }
    X_mat <- matrix(0, nrow = tau, ncol = tau-1) # fixed
    for(i in 1:tau){ for(j in 1:(tau-1)){ if(i > j) X_mat[i,j] <- d_vec[j] } }
    theta_mat <- z_mat <- u_mat <- matrix(0, nrow=tau, ncol=p)

    for(lambda_index in 1:length(list_of_lambda)){

      lambda <- list_of_lambda[lambda_index]
      output <- CPD_STERGM_cpp(ADMM_iteration, theta_iter, z_iter, H_pos_list, H_neg_list, y_pos_list, y_neg_list,
                               theta_mat, z_mat, u_mat, X_mat, d_vec, alpha=10, tau, p1, p2, lambda, theta_tol, ADMM_tol,
                               update_alpha, verbose)

      if(output[[4]] == "Save"){
        result <- Evaluation_ADMM_list(H_pos_list, H_neg_list, y_pos_list, y_neg_list,
                                        output[[1]], tau, p1, p2, p, n, true_CP)

        temp[[rep_iter]][[lambda_index]] <- list()
        temp[[rep_iter]][[lambda_index]][[1]] <- output[[1]] # theta_mat
        temp[[rep_iter]][[lambda_index]][[2]] <- result # list (est_CP, BIC, four metrics, data-driven threshold)
        BIC_table[rep_iter, lambda_index] <- result[[2]] # BIC
      }

    }

  } # End of all learning

  for(rep_iter in 1:num_nets){
    index_min_BIC <- which.min(BIC_table[rep_iter,])
    output_table[rep_iter,1] <- temp[[rep_iter]][[index_min_BIC]][[2]][[3]]
    output_table[rep_iter,2] <- temp[[rep_iter]][[index_min_BIC]][[2]][[4]]
    output_table[rep_iter,3] <- temp[[rep_iter]][[index_min_BIC]][[2]][[5]]
    output_table[rep_iter,4] <- temp[[rep_iter]][[index_min_BIC]][[2]][[6]]
  }

  colnames(output_table) <- c("abs_err", "C.hat_C", "C_C.hat", "Coverage")
  return(output_table)
}











#' Change point detection for a sequence of dynamic networks
#' @description This function detects the change points for a sequence of dynamic networks
#' @param y_data The sequence of dynamic networks
#' @param network_stats The network statistics for both formation and dissolution models. See search.ergmTerms() for a comprehensive list from the library(ergm).
#' @param directed Whether the networks are directed or not.
#' @param node_attr The nodal attributes.
#' @param list_of_lambda The list of tuning parameters \eqn{\lambda}.
#' @param ADMM_iteration The learning iteration for ADMM.
#' @param theta_iter The learning iteration for \eqn{\bm{\theta}}.
#' @param z_iter The learning iteration for \eqn{\bm{z}}.
#' @param theta_tol The tolerance for \eqn{\bm{\theta}} stopping criteria.
#' @param ADMM_tol The tolerance for ADMM stopping criteria.
#' @param threshold_alpha The alpha level for the data-driven threshold to declare the change points.
#' @param update_alpha If it is TRUE, the alpha is updated with a schedule.
#' @param verbose If it is TRUE, the information at each ADMM iteration is printed.
#'
#' @return Returns a list of results including detected change points, selected lambda, BIC value, log-likelihood, sequential parameter change after standardization \eqn{\Delta \hat{\bm{\zeta}}}, data-driven threshold to declare the change points, and the learned parameters \eqn{\bm{\theta}}.
#' @export
#'
#' @examples
#' library(ecp)
#' data(DJIA)
#' market <- DJIA$market
#' date_vec <- DJIA$dates[1:1138]
#' rownames(market) <- date_vec
#'
#' start <- which(date_vec == '2007-01-01')
#' end <- which(date_vec == '2010-01-04')
#' date_range <- start:end
#' mydata <- array(NA,c(29, 29, length(date_range)))
#' df <- list()
#'
#' for(i in 1:length(date_range)){
#'   temp <- market[date_range[i]:(date_range[i]+3),]
#'   temp <- ifelse(cor(temp)< 0, 1, 0)
#'   diag(temp) <- 0
#'   df[[i]] <- temp
#' }
#'
#' result <- CPD_STERGM(df, directed=FALSE, network_stats=c("edges", "triangles"),
#'                      list_of_lambda=10^c(2:5))
#'
#' seq_date <- rownames(market[start:end,])
#' seq_date[result[[1]]]
CPD_STERGM <- function(y_data, directed, network_stats, node_attr=NA,
                       list_of_lambda=10^c(-2:7), ADMM_iteration=200, theta_iter=20, z_iter=20,
                       theta_tol=1e-3, ADMM_tol=1e-7, threshold_alpha=0.1, update_alpha=TRUE, verbose=FALSE){
  temp <- list()
  BIC_vec <- rep(NA, length(list_of_lambda))
  p1 <- p2 <- length(network_stats); p <- p1+p2

  input_data <- save_H_y_list(y_data, directed, network_stats, node_attr)
  H_pos_list <- input_data[[1]]
  H_neg_list <- input_data[[2]]
  y_pos_list <- input_data[[3]]
  y_neg_list <- input_data[[4]]

  n <- dim(y_data[[1]])[1]
  tau <- length(y_data)-1
  d_vec <- numeric(tau-1) # fixed
  for(i in 1:(tau-1)){ d_vec[i] = sqrt( tau/ (i*(tau-i)) ) }
  X_mat <- matrix(0, nrow = tau, ncol = tau-1) # fixed
  for(i in 1:tau){ for(j in 1:(tau-1)){ if(i > j) X_mat[i,j] <- d_vec[j] } }
  theta_mat <- z_mat <- u_mat <- matrix(0, nrow=tau, ncol=p)

  for(lambda_index in 1:length(list_of_lambda)){

    lambda <- list_of_lambda[lambda_index]
    output <- CPD_STERGM_cpp(ADMM_iteration, theta_iter, z_iter, H_pos_list, H_neg_list, y_pos_list, y_neg_list,
                             theta_mat, z_mat, u_mat, X_mat, d_vec, alpha=10, tau, p1, p2, lambda, theta_tol, ADMM_tol,
                             update_alpha, verbose)

    if(output[[4]] == "Save"){
      result <- Evaluation_ADMM(H_pos_list, H_neg_list, y_pos_list, y_neg_list,
                                output[[1]], tau, p1, p2, p, n, directed, threshold_alpha)
      temp[[lambda_index]] <- list()
      temp[[lambda_index]][[1]] <- output[[1]] # theta_mat
      temp[[lambda_index]][[2]] <- result # list (est_CP, BIC, theta_change, data-driven threshold, log_lik)
      BIC_vec[lambda_index] <- result[[2]] # BIC
    }

  }

  output_list <- list()
  if( sum(is.na(BIC_vec)) < length(list_of_lambda) ){
    index_min_BIC            <- which.min(BIC_vec)
    output_list$est_CP       <- temp[[index_min_BIC]][[2]][[1]] # est_CP
    output_list$lambda       <- list_of_lambda[index_min_BIC]   # selected lambda
    output_list$BIC          <- temp[[index_min_BIC]][[2]][[2]] # BIC value
    output_list$log_lik      <- temp[[index_min_BIC]][[2]][[5]] # log_lik
    output_list$theta_change <- temp[[index_min_BIC]][[2]][[3]] # theta_change
    output_list$threshold    <- temp[[index_min_BIC]][[2]][[4]] # threshold
    output_list$theta_mat    <- temp[[index_min_BIC]][[1]]      # theta_mat
  }else{
    message("The algorithm does not converage for any choice of lambda (tuning parameter). An empty list is returned.")
  }

  return(output_list)
}




