\name{LQE}
\docType{methods}
\alias{LQE}
\alias{LQE-methods}
\alias{LQE,ANY,missing-method}
\alias{LQE,DVH.list,numeric-method}
\alias{LQE,DVH,numeric-method}
\alias{LQE,numeric,numeric-method}
\title{ Linear Quadratic Extrapolated (LQE) Dose Conversion }
\description{
This function and its associated methods calculate LQE-weighted dose(s) for a given \code{DVH} or \code{DVH.list} object, allowing iso-effect conversion between one or more DVH(s) employing different dose fractionation.
}
\usage{
\S4method{LQE}{DVH,numeric}(x, aB, fractions = NULL, dose.units = c("cGy", "Gy"))

\S4method{LQE}{DVH.list,numeric}(x, aB, fractions = NULL, dose.units = NULL)	

\S4method{LQE}{numeric,numeric}(x, aB, fractions = NULL)	
}
\arguments{
  \item{x}{
	Can represent either numeric value(s) or otherwise an object of class \code{DVH} or \code{DVH.list} for which to calculate LQE value(s).
  }
  \item{aB}{
	Numeric value specifying tissue-specific parameter "alpha/beta ratio". Multiple values may be specified for LQE conversion of \code{x}, where \code{x} is a \code{DVH.list} object containing multiple \code{DVH} objects.
  }
  \item{fractions}{
	Numeric value specifying the number of fractions in which doses are deposited. Multiple values may be specified for LQE conversion of \code{x}, where \code{x} is a \code{DVH.list} object containing multiple \code{DVH} objects. Two values should be specified when \code{x} contains numeric value(s), with the first element representing the fractionation for the input dose and the second element used to compute iso-effective output dose.
  }
  \item{dose.units}{
	Value specifying dose units (must be one of \code{"cGy"} or \code{"Gy"}). Multiple values may be specified for LQE conversion of \code{x}, where \code{x} is a \code{DVH.list} object containing multiple \code{DVH} objects. This parameter is not used for LQE conversion when \code{x} contains numeric value(s).
  }
}
\value{
Returns a \code{DVH} or \code{DVH.list} object, corresponding to original input type for parameter \code{x}.
}
\references{
Barendsen, G.W. (1982)
Dose fractionation, dose-rate and iso-effect relationships for normal-tissue response.
\emph{Int J Radiat Oncol Biol Phys} \bold{8}(11):1981-1997.
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:DVH-class]{DVH}}, \code{\link[RadOnc:DVH.list-class]{DVH.list}}
}
\examples{
# 20x 300cGy == 30x 244.6752cGy (isoeffective doses)
LQE(c(4500, 5500, 6000), aB=3, fractions=c(20, 30))
}
\keyword{methods}
\keyword{manip}
