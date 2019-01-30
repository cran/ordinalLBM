\name{olbm_dat}
\docType{data}
\alias{olbm_dat}
\title{OLBM simulated data}
\description{
It is a list containing i) an ordinal toy data matrix simulated acccording to OLBM and ii) the row/column cluster assignments. To see how the data are simulated, you can type "?simu.olbm" in the R console and look at "Examples". 
}
\usage{data(olbm_dat)}

\format{ 
A list containing three items.
	\describe{
		\item{Y}{: an ordinal data matrix simulated according to OLBM. }
		\item{Rclus}{: the actual row cluster assignments.}
 		\item{Cclust}{: the actual column cluster assignments.}
  }
}