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
\S4method{LQE}{DVH,numeric}(x, aB, fractions = NULL, N = NULL, dose.units = c("cGy", "Gy"))

\S4method{LQE}{DVH.list,numeric}(x, aB, fractions = NULL, N = NULL, dose.units = NULL)	

\S4method{LQE}{numeric,numeric}(x, aB, fractions = NULL, N = NULL, dose.units = c("cGy", "Gy"))	
}
\arguments{
  \item{x}{
	Can represent either numeric value(s) or otherwise an object of class \code{DVH} or \code{DVH.list} for which to calculate isoeffective dose(s) based on LQE. Note that input doses should be specified in units of \code{dose.units} (default is \code{"cGy"}).
  }
  \item{aB}{
	Numeric value specifying tissue-specific parameter "alpha/beta ratio". Multiple values may be specified for LQE conversion of \code{x}, where \code{x} is a \code{DVH.list} object containing multiple \code{DVH} objects.
  }
  \item{fractions}{
	Numeric value specifying the fraction sizes in which doses are deposited. Multiple values may be specified for LQE conversion of \code{x}, where \code{x} is a \code{DVH.list} object containing multiple \code{DVH} objects. Two values should be specified when \code{x} contains numeric value(s), with the first element representing the fractionation for the input dose and the second element used to compute iso-effective output dose. Only a single value should be used to calculate iso-effective doses for a single \code{DVH} input.
  }
  \item{N}{
	Numeric value specifying the number of fractions in which doses are deposited. Multiple values may be specified for LQE conversion of \code{x}, where \code{x} is a \code{DVH.list} object containing multiple \code{DVH} objects. Two values should be specified when \code{x} contains numeric value(s), with the first element representing the fractionation for the input dose and the second element used to compute iso-effective output dose. Only a single value should be used to calculate iso-effective doses for a single \code{DVH} input. Note that parameter \code{N} is ignored unless \code{fractions=NULL}.
  }
  \item{dose.units}{
	Value specifying output dose units (must be one of \code{"cGy"} or \code{"Gy"}). Note that this value also determines input dose units when \code{x} is numeric. Multiple values may be specified for LQE conversion of \code{x}, where \code{x} is a \code{DVH.list} object containing multiple \code{DVH} objects.
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
# 5x 300cGy == 9x 200cGy (1800cGy isoeffective dose)
LQE(1500, aB=3, fractions=c(300, 200), dose.units="cGy")
LQE(15, aB=3, N=c(5, 9), dose.units="Gy")
LQE(18, aB=3, fractions=c(2, 3), dose.units="Gy")

# Multiple simultaneous numerical queries
LQE(c(4500, 5400), aB=3, fractions=c(300, 200), dose.units="cGy")

# DVH processing
LQE(janedoe[["LIVER"]],aB=3,fractions=200,dose.units="cGy")

# DVH list processing (default dose.units is "cGy")
plot(janedoe)
plot(LQE(janedoe,aB=3,fractions=200),col="red",new=FALSE)
}
\keyword{methods}
\keyword{manip}
