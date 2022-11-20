## usethis namespace: start
#' @useDynLib CPDstergm, .registration = TRUE
## usethis namespace: end
NULL


## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL


#' @import network ndtv dplyr tergm network networkDynamic
#' @importFrom stats as.formula median qnorm rbinom sd simulate
#' @importFrom ergm ergmMPLE san
NULL











#' Simulate networks from SBM
#'
#' @param num_nets number of networks
#' @param n number of node
#' @param rho correlation coefficient
#'
#' @return list of networks
#' @export
#'
#' @examples
#' # set.seed(1)
#' # SBM_list <- sim_SBM_batch(num_nets=10, n=50, rho=0)
sim_SBM_batch <- function(num_nets=10, n=50, rho=0){

  output <- list()
  for(rep_iter in 1:num_nets){
    rho_a = rho
    K = 3
    v =  c(26, 51, 76)
    num_time = 100
    df = vector(mode = "list", length = num_time)

    for(t in 1:num_time){

      if( t==1 || t==v[2] ){ # t=1 or t=51

        #print('1 or 51')
        P =  matrix(0.3,n,n)
        P[1:floor(n/K), 1:floor(n/K)] = 0.5
        P[(1+floor(n/K)):(2*floor(n/K)),(1+floor(n/K)):(2*floor(n/K)) ] = 0.5
        P[(1+2*floor(n/K)):n,(1+2*floor(n/K)):n ] = 0.5
        diag(P) = 0

        A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)

      }

      if(t == v[1] || t == v[3]){ # t=26 or t=76

        #print('26 or 76')
        Q =  matrix(0.2,n,n)
        Q[1:floor(n/K), 1:floor(n/K)] = 0.45
        Q[(1+floor(n/K)):(2*floor(n/K)),(1+floor(n/K)):(2*floor(n/K)) ] = 0.45
        Q[(1+2*floor(n/K)):n,(1+2*floor(n/K)):n ] = 0.45
        diag(Q) = 0

        A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),Q),n,n)

      }

      if( (t > 1 && t < v[1])  ||  (t > v[2] && t < v[3]) ){ # t=2 to t=25 or t=52 to t=75
        #print('t=2 to t=25 or t=52 to t=75')
        aux1 = (1-P)*rho_a + P # (1-E(t+1))*rho + E(t+1) if A(t) = 1
        aux2 = P*(1-rho_a) # (1-E(t+1))*rho + E(t+1) if A(t) = 0

        aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
        aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
        A =  aux1*A + aux2*(1-A)

      }

      if( (t > v[1] && t < v[2]) || ((t > v[3] && t <= num_time)) ){ # t=27 to t=50 or t=77 to t=100
        #print('t=27 to t=50 or t=77 to t=100')
        aux1 = (1-Q)*rho_a + Q
        aux2 = Q*(1-rho_a)

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


#' Simulate networks from STERGM
#'
#' @param num_nets number of networks
#' @param n number of node
#' @param network_stats network statistics
#' @param coefs_pos coefficients for the formation model
#' @param coefs_neg coefficients for the dissolution model
#' @param y1_stats edges for the first network
#' @param node_attr nodal attribute
#'
#' @return list of networks
#' @export
#'
#' @examples
#' # coefs_pos <- matrix(c(-1, -1, -1, -1, -2,  1, -2,  1), nrow=2, ncol=4, byrow = T)
#' # coefs_neg <- matrix(c( -1,  -1, -1,  -1, -2,  -1, -2,  -1), nrow=2, ncol=4, byrow = T)
#'
#' # y1_stats <- 250
#' # network_stats=c("edges", "mutual")
#'
#' # set.seed(1)
#' # STERGM_list <- sim_STERGM_batch(num_nets=10, n=50, network_stats,
#' #                                 coefs_pos, coefs_neg, y1_stats, node_attr=NA)
sim_STERGM_batch <- function(num_nets=10, n=50, network_stats,
                             coefs_pos, coefs_neg, y1_stats, node_attr=NA){
  num_nodes <- n
  num_timepts <- 100
  num_changepts <- 3
  form_model <- diss_model <- as.formula(paste("~", paste(network_stats, collapse = "+")))

  output <- list()
  for(rep_iter in 1:num_nets){

    # generate initial network
    g0<-network.initialize(num_nodes, directed=T) # empty
    g1<-san(g0~edges,target.stats=y1_stats, verbose=TRUE)
    if(any(!is.na(node_attr))){
      network::set.vertex.attribute(g1, "Gender", node_attr)
    }
    ginit <- g1

    time_stable <- num_timepts / (num_changepts+1)
    cur_end <- 0
    res_adj_list <- vector(mode = 'list', length = num_timepts)
    for(i in 1:(num_changepts+1)){
      cat('[INFO] Simulate from ', cur_end+1, ' to ', cur_end + time_stable, '\n')
      stergm.sim <- simulate(
        # if the input graph already has 'net.obs.period', to continue the simulation,
        # we need to specify the time.start to be the end of 'net.obs.period'
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


#' Construct input data to the algorithm
#'
#' @param y_data network data
#' @param network_stats network statistics
#' @param node_attr nodal attribute
#'
#' @return list of input data
save_H_y_list <- function(y_data, network_stats, node_attr=NA){

  output <- list()
  num_time <- length(y_data)
  n <- dim(y_data[[1]])[1]

  #library(ergm)
  #library(network)

  H_pos_list <- H_neg_list <- y_pos_list <- y_neg_list <- list()

  f_pos <- as.formula(paste("y_pos", paste(network_stats, collapse = "+"), sep = "~"))
  f_neg <- as.formula(paste("y_neg", paste(network_stats, collapse = "+"), sep = "~"))

  for(t in 2:num_time){
    y_before <- y_data[[t-1]]
    y_after <- y_data[[t]]

    y_pos = network(matrix(as.integer(y_before | y_after), nrow=n, ncol=n, byrow=F))
    y_neg = network(matrix(as.integer(y_before & y_after), nrow=n, ncol=n, byrow=F))

    if(any(!is.na(node_attr))){
      network::set.vertex.attribute(y_pos, "Gender", node_attr)
      network::set.vertex.attribute(y_neg, "Gender", node_attr)
    }

    y_pos = ergmMPLE(f_pos, output='array')
    y_neg = ergmMPLE(f_neg, output='array')

    H_pos_output <- c()
    for(i in 1:length(network_stats)){
      # stack matrix row by row into one vector
      H_pos_output <- cbind(H_pos_output, c(t(y_pos$predictor[,,i])))
    }

    H_neg_output <- c()
    for(i in 1:length(network_stats)){
      # stack matrix row by row into one vector
      H_neg_output <- cbind(H_neg_output, c(t(y_neg$predictor[,,i])))
    }

    y_pos_output <- cbind(c(t(y_pos$response)))
    y_neg_output <- cbind(c(t(y_neg$response)))

    H_pos_output[is.na(H_pos_output)] <- 0 # diagonal NA becomes 0
    H_neg_output[is.na(H_neg_output)] <- 0
    y_pos_output[is.na(y_pos_output)] <- 0
    y_neg_output[is.na(y_neg_output)] <- 0

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


#' Evaluate the algorithm performance
#'
#' @param H_pos_list change statistics in the formation model
#' @param H_neg_list change statistics in the dissolution model
#' @param y_pos_list formation network
#' @param y_neg_list dissolution network
#' @param theta_mat theta parameter
#' @param tau number of row
#' @param p1 number of network statistics in the formation model
#' @param p2 number of network statistics in the dissolution model
#' @param p total number of network statistics
#' @param n number of node
#' @param true_CP true change points
#'
#' @return result of metrics
Evaluation_ADMM <- function(H_pos_list, H_neg_list, y_pos_list, y_neg_list,
                            theta_mat, tau, p1, p2, p, n, true_CP){
  log_lik <- 0
  for(iter in 1:tau){
    theta_g_pos <- H_pos_list[[iter]] %*% theta_mat[iter,1:p1]
    theta_g_neg <- H_neg_list[[iter]] %*% theta_mat[iter,(p1+1):p]
    log_lik <- log_lik + sum(y_pos_list[[iter]] * theta_g_pos - log(1 + exp(theta_g_pos))) # element-wise
    log_lik <- log_lik + sum(y_neg_list[[iter]] * theta_g_neg - log(1 + exp(theta_g_neg)))
  }

  theta_change <- numeric(tau-1)
  for(i in 1:(tau-1)){theta_change[i] <- norm(theta_mat[i+1,]-theta_mat[i,],"2")}
  t_change <- (theta_change - median(theta_change))/sd(theta_change) # use median
  threshold <- mean(t_change) + qnorm(0.9, lower.tail = T) * sd(t_change) # use mean
  theta_change <- t_change

  #plot(1:(tau-1), theta_change, type='b');abline(h = threshold)

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


#' Perform change point detection
#'
#' @param data_list list of networks
#' @param network_stats network statistics
#' @param node_attr nodal attribute
#' @param list_of_lambda list of tunning parameters
#' @param ADMM_iteration ADMM learning iteration
#' @param theta_iter theta learning iteration
#' @param z_iter z learning iteration
#' @param theta_tol tolerance for theta stopping criteria
#' @param ADMM_tol tolerance for ADMM stopping criteria
#' @param true_CP true change points
#'
#' @return result of metrics
#' @export
#'
#' @examples
#' # network_stats=c("edges", "mutual")
#' # sim1_result1 <- CPD_STERGM_batch(SBM_list, network_stats)
#' # colMeans(sim1_result1)
CPD_STERGM_batch <- function(data_list, network_stats, node_attr=NA,
                             list_of_lambda=10^c(-2:7), ADMM_iteration=200, theta_iter=20, z_iter=20,
                             theta_tol=1e-3, ADMM_tol=1e-7, true_CP=c(26, 51, 76)){
  num_nets <- length(data_list)
  temp <- list()
  BIC_table <- matrix(NA, nrow=num_nets, ncol=length(list_of_lambda))
  output_table <- matrix(NA, nrow=num_nets, ncol=4) # 4 metrics
  p1 <- p2 <- length(network_stats); p <- p1+p2

  for(rep_iter in 1:num_nets){

    temp[[rep_iter]] <- list()
    y_data <- data_list[[rep_iter]]
    input_data <- save_H_y_list(y_data, network_stats, node_attr)
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
      output <- CPD_STERGM(ADMM_iteration, theta_iter, z_iter, H_pos_list, H_neg_list, y_pos_list, y_neg_list,
                           theta_mat, z_mat, u_mat, X_mat, d_vec, alpha=10, tau, p1, p2, lambda, theta_tol, ADMM_tol)

      if(output[[4]] == "Save"){
        result <- Evaluation_ADMM(H_pos_list, H_neg_list, y_pos_list, y_neg_list,
                                  output[[1]], tau, p1, p2, p, n, true_CP)

        temp[[rep_iter]][[lambda_index]] <- list()
        temp[[rep_iter]][[lambda_index]][[1]] <- output[[1]] # theta_mat
        temp[[rep_iter]][[lambda_index]][[2]] <- result # list object
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

  return(output_table)
}



