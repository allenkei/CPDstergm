% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_list.R
\name{sim_RDPG_list}
\alias{sim_RDPG_list}
\title{Simulation of dynamic networks from the Random Dot Product Graph Model}
\usage{
sim_RDPG_list(num_seq = 1, n = 50, rho = 0.9, d = 10)
}
\arguments{
\item{num_seq}{The number of sequences of dynamic networks.}

\item{n}{The number of node.}

\item{rho}{The correlation coefficient between consecutive time points for latent position.}

\item{d}{The latent dimension.}
}
\value{
Returns a list of multiple sequences of dynamic networks.
}
\description{
This function simulates multiple sequences of dynamic networks, from the Random Dot Product Graph Model (RDPGM).
}
\examples{
set.seed(1)
RDPG_list <- sim_RDPG_list(num_seq = 1, n = 50, rho = 0.9, d = 10)
}
