#include <RcppArmadillo.h>
#include <math.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;

///////////
// utils //
///////////

arma::vec sigmoid(arma::vec x){
  return 1/(1 + exp(-clamp(x, -50, 50))); // element-wise
}


arma::mat vec_to_mat(arma::vec x, int num_row, int num_col){
  // num_row is tau, num_col is p
  arma::mat output;

  output = x;
  output.reshape(num_col, num_row); // fill matrix column-wise, so (p by tau) first then transpose to (tau by p)

  return output.t();
}


List beta_gamma_from_z(arma::vec d_vec, arma::mat z_mat, int p){

  List output;
  arma::mat beta_mat, gamma;

  gamma = z_mat.submat(0,0,0,p-1);
  beta_mat = diff(z_mat, 1, 0) / repelem(d_vec, 1, p); // diff(z_mat, k=1, dim=0) difference between rows
  // repeat 1*nrows but p*ncols for element-wise division

  output.push_back(gamma);
  output.push_back(beta_mat);

  return output;
}


arma::mat z_from_beta_gamma(List beta_gamma_list, int tau, arma::mat X_mat){

  arma::mat output;
  // (tau by p) + (tau by p)
  output = repelem(as<arma::mat>(beta_gamma_list[0]), tau, 1) + X_mat*as<arma::mat>(beta_gamma_list[1]);

  return output;
}

// [[Rcpp::export]]
double cal_log_likelihood(List H_pos_list, List H_neg_list, List y_pos_list, List y_neg_list, arma::mat theta_mat, int tau, int p1, int p2){

  double log_likelihood = 0;
  arma::mat theta_g_pos, theta_g_neg;
  int p = p1+p2;

  for(int iter = 0; iter < tau; ++iter){

    theta_g_pos = as<arma::mat>(H_pos_list[iter]) * theta_mat.submat(iter,0,iter,p1-1).t(); // (E by p1) * (p1 by 1) = (E by 1)
    theta_g_neg = as<arma::mat>(H_neg_list[iter]) * theta_mat.submat(iter,p1,iter,p-1).t(); // (E by p2) * (p2 by 1) = (E by 1)
    // * is matrix multiplication

    log_likelihood += sum( as<arma::vec>(y_pos_list[iter]) % theta_g_pos - log(1 + exp(theta_g_pos)) );
    log_likelihood += sum( as<arma::vec>(y_neg_list[iter]) % theta_g_neg - log(1 + exp(theta_g_neg)) );
    // % is element-wise multiplication

  }

  return log_likelihood;
}


////////////////////
// theta Learning //
////////////////////


arma::vec cal_Gradient(List H_pos_list, List H_neg_list, List y_pos_list, List y_neg_list,
                       arma::mat theta_mat, arma::mat z_mat, arma::mat u_mat, double alpha, int tau, int p1, int p2){

  // tau is num_time-1
  // H_pos_list[iter] is (E by p1) and H_neg_list[iter] is (E by p2)
  // theta_mat is (tau by p)
  int p = p1+p2;
  arma::vec mu_pos, mu_neg, y_mu_pos, y_mu_neg; // (E by 1)
  arma::vec gradient; // (tau p by 1)

  gradient.zeros(tau*p); // gradient is (tau p by 1)

  for(int iter = 0; iter < tau; ++iter){

    mu_pos = sigmoid(as<arma::mat>(H_pos_list[iter]) * theta_mat.submat(iter,0,iter,p1-1).t()); // (E by p1) * (p1 by 1) = (E by 1)
    mu_neg = sigmoid(as<arma::mat>(H_neg_list[iter]) * theta_mat.submat(iter,p1,iter,p-1).t()); // (E by p2) * (p2 by 1) = (E by 1)
    // * is matrix multiplication

    y_mu_pos = as<arma::vec>(y_pos_list[iter]) - mu_pos; // (E by 1) - (E by 1) = (E by 1)
    y_mu_neg = as<arma::vec>(y_neg_list[iter]) - mu_neg; // (E by 1) - (E by 1) = (E by 1)

    gradient.subvec(iter*p, iter*p+p1-1) = -as<arma::mat>(H_pos_list[iter]).t() * y_mu_pos;       // (p1 by E) * (E by 1)  = (p1 by 1)
    gradient.subvec(iter*p+p1, iter*p+p1+p2-1) = -as<arma::mat>(H_neg_list[iter]).t() * y_mu_neg; // (p2 by E) * (E by 1)  = (p2 by 1)
    // copy to corresponding positions for each p1 and p2
    // there is a NEGATIVE sign here
  }

  gradient += alpha * (theta_mat.t().as_col() - z_mat.t().as_col() + u_mat.t().as_col());

  return gradient;

}


arma::mat cal_Hessian(List H_pos_list, List H_neg_list, List y_pos_list, List y_neg_list, arma::mat theta_mat,
                      double alpha, int tau, int p1, int p2){

  // tau is num_time-1
  // H_pos_list[iter] is (E by p1) and H_neg_list[iter] is (E by p2)
  // theta_mat is (tau by p)
  int p = p1+p2;
  arma::mat hessian;
  arma::vec mu_pos, mu_neg;
  arma::vec temp_pos, temp_neg;

  hessian.zeros(tau*p,tau*p);

  for(int iter = 0; iter < tau; ++iter){

    mu_pos = sigmoid(as<arma::mat>(H_pos_list[iter]) * theta_mat.submat(iter,0,iter,p1-1).t()); // (E by p1) * (p1 by 1) = (E by 1)
    mu_neg = sigmoid(as<arma::mat>(H_neg_list[iter]) * theta_mat.submat(iter,p1,iter,p-1).t()); // (E by p2) * (p2 by 1) = (E by 1)

    temp_pos = mu_pos % (1-mu_pos); // (E by 1)
    temp_neg = mu_neg % (1-mu_neg); // (E by 1)
    // % is element-wise multiplication

    hessian.submat(iter*p, iter*p, iter*p+p1-1, iter*p+p1-1) = as<arma::mat>(H_pos_list[iter]).t() * diagmat(temp_pos) * as<arma::mat>(H_pos_list[iter]);
    hessian.submat(iter*p+p1, iter*p+p1, iter*p+p1+p2-1, iter*p+p1+p2-1) = as<arma::mat>(H_neg_list[iter]).t() * diagmat(temp_neg) * as<arma::mat>(H_neg_list[iter]);
    // (p1 by E) * (E by E) * (E by p1) = (p1 by p1)
    // (p2 by E) * (E by E) * (E by p2) = (p2 by p2)

  }

  hessian += alpha * arma::eye<arma::mat>(tau*p,tau*p);
  return hessian;

}


arma::mat theta_learning(int learning_iter, List H_pos_list, List H_neg_list, List y_pos_list, List y_neg_list,
                         arma::mat theta_mat, arma::mat z_mat, arma::mat u_mat, double alpha, int tau, int p1, int p2, double theta_tol, bool verbose){

  // tau is num_time-1
  // H_pos_list[iter] is (E by p1) and H_neg_list[iter] is (E by p2)
  // theta_mat is (tau by p)
  int p = p1+p2;
  arma::vec Grad;
  arma::mat Hess;
  arma::vec theta_vec, theta_vec_old;
  double resnorm;

  for(int iter = 0; iter < learning_iter; ++iter){

    theta_vec = theta_mat.t().as_col(); // (tau by p) mat to (tau*p by 1) vec  (transpose and stack columns)
    theta_vec_old = theta_vec;

    Grad = cal_Gradient(H_pos_list, H_neg_list, y_pos_list, y_neg_list, theta_mat, z_mat, u_mat, alpha, tau, p1, p2);
    Hess = cal_Hessian(H_pos_list, H_neg_list, y_pos_list, y_neg_list, theta_mat, alpha, tau, p1, p2);
    theta_vec -= arma::solve(Hess,arma::eye<arma::mat>(tau*p,tau*p)) * Grad;
    theta_mat = vec_to_mat(theta_vec, tau, p);

    // residual = theta_vec - theta_vec_old;
    resnorm = norm(theta_vec - theta_vec_old, 2);

    if( resnorm < theta_tol ){
      if(verbose){Rcout << "    ---theta converged at iteration " << iter+1 << "\n";}
      break;
    }

  }

  return theta_mat;

}

////////////////
// z Learning //
////////////////


arma::mat z_learning(int learning_iter, arma::mat theta_mat, arma::mat z_mat, arma::mat u_mat,
                     arma::mat X_mat, arma::vec d_vec, double alpha, int tau, int p1, int p2, double lambda, bool verbose){

  int p = p1+p2, counter;
  List beta_gamma_list;
  arma::mat z_updated, gamma, beta_mat, beta_i, condition;
  arma::mat si, temp;
  arma::mat X_without_i, beta_without_i;
  arma::mat zero_row_for_beta, zero_col_for_X;
  List output;
  double si_norm;

  zero_row_for_beta.zeros(1,p);
  zero_col_for_X.zeros(tau,1);

  beta_gamma_list = beta_gamma_from_z(d_vec, z_mat, p);
  gamma = as<arma::mat>(beta_gamma_list[0]);
  beta_mat = as<arma::mat>(beta_gamma_list[1]);

  for(int iter = 0; iter < learning_iter; ++iter){

    // Learning once
    for(int i = 0; i < (tau-1); ++i){
      beta_without_i = beta_mat;
      X_without_i = X_mat;
      beta_without_i.submat(i,0,i,p-1) = zero_row_for_beta;
      X_without_i.submat(0,i,tau-1,i) = zero_col_for_X;

      si = alpha * X_mat.col(i).t() * (theta_mat + u_mat - repelem(gamma, tau, 1) - X_without_i*beta_without_i);
      si_norm = norm(si, 2);

      if(si_norm <= lambda){
        beta_mat.submat(i,0,i,p-1) = zero_row_for_beta;
      }else{
        temp = 1/( alpha*as_scalar(X_mat.col(i).t()*X_mat.col(i)) ) * (1 - lambda/si_norm) * si;
        beta_mat.submat(i,0,i,p-1) = temp;
      }
    }// End of beta learning once

    gamma = arma::mean(theta_mat + u_mat - X_mat*beta_mat,0); // End of gamma learning once

    // Monitor convergence by KKT conditions
    counter = 0;
    for(int i = 0; i < (tau-1); ++i){
      beta_i = beta_mat.submat(i,0,i,p-1);

      if(norm(beta_i, 2) == 0){
        condition = -alpha * X_mat.col(i).t() * (theta_mat + u_mat - repelem(gamma, tau, 1) - X_mat*beta_mat); // negative sign
        if(norm(condition, 2) < lambda){ counter += 1; }
      }else{
        condition = -alpha * X_mat.col(i).t() * (theta_mat + u_mat - repelem(gamma, tau, 1) - X_mat*beta_mat); // negative sign
        condition += lambda * beta_i/norm(beta_i, 2);
        if(norm(condition, 2) == 0){ counter += 1; }
      }

    }// End of checking

    if( counter == (tau-1) ){
      if(verbose){Rcout << "    ---beta converged at iteration " << iter+1 << "\n"; }
      break;
    }

  }// End of beta learning

  output.push_back(gamma);
  output.push_back(beta_mat);
  z_updated = z_from_beta_gamma(output, tau, X_mat);
  return z_updated;
}


////////////////
// CPD STERGM //
////////////////


// [[Rcpp::export]]
List CPD_STERGM_cpp(int ADMM_iter, int theta_iter, int z_iter, List H_pos_list, List H_neg_list, List y_pos_list, List y_neg_list,
                    arma::mat theta_mat, arma::mat z_mat, arma::mat u_mat, arma::mat X_mat, arma::vec d_vec,
                    double alpha, int tau, int p1, int p2, double lambda, double theta_tol, double ADMM_tol, bool verbose){

  Rcpp::Rcout.precision(15);
  int p = p1+p2;
  List output;
  arma::mat z_mat_old;
  arma::mat primal_residual, dual_residual;
  double primal_resnorm, dual_resnorm, converged;
  String s("DoNotSave");

  double log_lik, old_log_lik;
  old_log_lik = arma::datum::inf;

  for(int iter = 0; iter < ADMM_iter; ++iter){
    if(verbose){Rcout << "ADMM iter = " << iter+1 << "\n";}

    theta_mat = theta_learning(theta_iter, H_pos_list, H_neg_list, y_pos_list, y_neg_list, theta_mat, z_mat, u_mat, alpha, tau, p1, p2, theta_tol, verbose);
    if(verbose){Rcout << "    theta updated" << "\n";}

    z_mat_old = z_mat;
    z_mat = z_learning(z_iter, theta_mat, z_mat, u_mat, X_mat, d_vec, alpha, tau, p1, p2, lambda, verbose);
    if(verbose){Rcout << "    z updated" << "\n";}

    dual_residual = z_mat - z_mat_old;
    primal_residual = theta_mat - z_mat;

    u_mat += primal_residual; // End of u_mat learning
    if(verbose){Rcout << "    u updated" << "\n";}

    dual_residual.reshape(tau*p, 1);
    primal_residual.reshape(tau*p, 1);
    dual_resnorm = as_scalar(sqrt(arma::mean(pow(dual_residual,2),0)));
    primal_resnorm = as_scalar(sqrt(arma::mean(pow(primal_residual,2),0)));
    log_lik = cal_log_likelihood(H_pos_list, H_neg_list, y_pos_list, y_neg_list, theta_mat, tau, p1, p2);
    converged = max(dual_resnorm, primal_resnorm);

    if(verbose){
      Rcout << "        primal_resnorm = " << primal_resnorm << "\n";
      Rcout << "        dual_resnorm = " << dual_resnorm << "\n";
      Rcout << "        converged = " << converged << "\n";
      Rcout << "        log_lik = " << log_lik << "\n";
      Rcout << "        tol = " << abs( (log_lik-old_log_lik)/old_log_lik ) << "\n";
    }


    if( abs( (log_lik-old_log_lik)/old_log_lik ) < ADMM_tol ){s="Save";break;}


    if(primal_resnorm > dual_resnorm * 10){
      alpha *= 2;
      u_mat /= 2;
      if(verbose){Rcout << "        ---alpha is increased to " << alpha << "\n";}
    }else if(dual_resnorm > primal_resnorm * 10){
      alpha /= 2;
      u_mat *= 2;
      if(verbose){Rcout << "        ---alpha is decreased to " << alpha << "\n";}
    }

    old_log_lik = log_lik;
    if(converged > 10){break;} // not going to converge, so early stop

  }

  output.push_back(theta_mat);
  output.push_back(z_mat);
  output.push_back(u_mat);
  output.push_back(s);
  output.push_back(log_lik);
  return output;
}



////////////////////
// standard error //
////////////////////

arma::vec cal_Gradient_SE(List H_pos_list, List H_neg_list, List y_pos_list, List y_neg_list,
                          arma::mat theta_mat, arma::mat z_mat, arma::mat u_mat, double alpha, int tau, int p1, int p2){

  // tau is num_time-1
  // H_pos_list[iter] is (E by p1) and H_neg_list[iter] is (E by p2)
  // theta_mat is (tau by p)
  int p = p1+p2;
  arma::vec mu_pos, mu_neg, y_mu_pos, y_mu_neg; // (E by 1)
  arma::vec gradient; // (tau p by 1)

  gradient.zeros(tau*p); // gradient is (tau p by 1)

  for(int iter = 0; iter < tau; ++iter){

    mu_pos = sigmoid(as<arma::mat>(H_pos_list[iter]) * theta_mat.submat(iter,0,iter,p1-1).t()); // (E by p1) * (p1 by 1) = (E by 1)
    mu_neg = sigmoid(as<arma::mat>(H_neg_list[iter]) * theta_mat.submat(iter,p1,iter,p-1).t()); // (E by p2) * (p2 by 1) = (E by 1)
    // * is matrix multiplication

    y_mu_pos = as<arma::vec>(y_pos_list[iter]) - mu_pos; // (E by 1) - (E by 1) = (E by 1)
    y_mu_neg = as<arma::vec>(y_neg_list[iter]) - mu_neg; // (E by 1) - (E by 1) = (E by 1)

    gradient.subvec(iter*p, iter*p+p1-1) = -as<arma::mat>(H_pos_list[iter]).t() * y_mu_pos;       // (p1 by E) * (E by 1)  = (p1 by 1)
    gradient.subvec(iter*p+p1, iter*p+p1+p2-1) = -as<arma::mat>(H_neg_list[iter]).t() * y_mu_neg; // (p2 by E) * (E by 1)  = (p2 by 1)
    // copy to corresponding positions for each p1 and p2
    // there is a NEGATIVE sign here
  }

  return gradient;

}


arma::mat cal_Hessian_SE(List H_pos_list, List H_neg_list, List y_pos_list, List y_neg_list, arma::mat theta_mat,
                         double alpha, int tau, int p1, int p2){

  // tau is num_time-1
  // H_pos_list[iter] is (E by p1) and H_neg_list[iter] is (E by p2)
  // theta_mat is (tau by p)
  int p = p1+p2;
  arma::mat hessian;
  arma::vec mu_pos, mu_neg;
  arma::vec temp_pos, temp_neg;

  hessian.zeros(tau*p,tau*p);

  for(int iter = 0; iter < tau; ++iter){

    mu_pos = sigmoid(as<arma::mat>(H_pos_list[iter]) * theta_mat.submat(iter,0,iter,p1-1).t()); // (E by p1) * (p1 by 1) = (E by 1)
    mu_neg = sigmoid(as<arma::mat>(H_neg_list[iter]) * theta_mat.submat(iter,p1,iter,p-1).t()); // (E by p2) * (p2 by 1) = (E by 1)

    temp_pos = mu_pos % (1-mu_pos); // (E by 1)
    temp_neg = mu_neg % (1-mu_neg); // (E by 1)
    // % is element-wise multiplication

    hessian.submat(iter*p, iter*p, iter*p+p1-1, iter*p+p1-1) = as<arma::mat>(H_pos_list[iter]).t() * diagmat(temp_pos) * as<arma::mat>(H_pos_list[iter]);
    hessian.submat(iter*p+p1, iter*p+p1, iter*p+p1+p2-1, iter*p+p1+p2-1) = as<arma::mat>(H_neg_list[iter]).t() * diagmat(temp_neg) * as<arma::mat>(H_neg_list[iter]);
    // (p1 by E) * (E by E) * (E by p1) = (p1 by p1)
    // (p2 by E) * (E by E) * (E by p2) = (p2 by p2)

  }

  return hessian;

}


