% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/new_MWBS.R
\name{new_MWBS}
\alias{new_MWBS}
\title{new_MWBS function}
\usage{
new_MWBS(y, z, s, e, flag = 0, S = NULL, Dval = NULL, pos = 1, alpha, beta, h)
}
\arguments{
\item{y}{a matrix containing the data. Row correspond to different time points and columns to different variables}

\item{z}{a copy of the matrix y, it can be y itself}

\item{flag}{}
}
\value{
S

Dval

pos
}
\description{
Sub-routine compute change points based on the  MNP method
}
