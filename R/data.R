#' @title The MIT Cellphone Data
#'
#' @description The Massachusetts Institute of Technology (MIT) cellphone data consists of the human interactions via cellphone activity, among 96 participants for a duration of 232 days. The data were taken from 15-Sept-2004 to 4-May-2005 inclusive. For participant \eqn{i} and participant \eqn{j}, a connected edge \eqn{\bm{y}_{ij}^t = 1} indicates that they had made at least one phone call on day \eqn{t}, and \eqn{\bm{y}_{ij}^t = 0} indicates that they had made no phone call on day \eqn{t}.
#'
#' @format A list with 232 matrices of dimension of 96 by 96.
#'
#' @references Eagle, Nathan and Pentland, Alex (Sandy). "Reality mining: sensing complex social systems." Personal and ubiquitous computing 10.4 (2006): 255-268.
#'
#' @examples
#' data("MITphone")
"MITphone"
