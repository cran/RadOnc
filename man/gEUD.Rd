\name{gEUD}
\docType{methods}
\alias{gEUD}
\alias{gEUD-methods}
\alias{gEUD,ANY,missing-method}
\alias{gEUD,DVH.list,numeric-method}
\alias{gEUD,DVH,numeric-method}
\title{ Generalized Equivalent Uniform Dose (gEUD) Calculation}
\description{
This function and its associated methods calculate gEUD value(s) for a given \code{DVH} or \code{DVH.list} object.
}
\usage{
\S4method{gEUD}{DVH,numeric}(x, a, dose.units = c("cGy", "Gy"))

\S4method{gEUD}{DVH.list,numeric}(x, a, dose.units = c("cGy", "Gy"))	
}
\arguments{
  \item{x}{
	Can represent either an object of class \code{DVH} or \code{DVH.list} for which to calculate gEUD value(s).
  }
  \item{a}{
	Numeric value specifying tissue-specific parameter, should be negative for target structures (e.g. tumor) and positive for organs at risk. For \code{a=1}, the EUD is equivalent to mean dose, while for \code{a=Inf} and \code{a=-Inf}, the EUD is equivalent to maximum and minimum doses, respectively.
  }
  \item{dose.units}{
	Value specifying dose units (must be one of \code{"cGy"} or \code{"Gy"}).
  }
}
\value{
Returns a numeric vector containing the gEUD values (in units of \code{dose.units}) for all input objects.
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\references{
Thieke, C., Bortfeld, T., Niemierko, A., and Nill, S. (2003)
From physical dose constraints to equivalent uniform dose constraints in inverse radiotherapy planning.
\emph{Med Phys} \bold{30}(9), 2332-2339.
}
\seealso{
\code{\link[RadOnc:DVH-class]{DVH}}, \code{\link[RadOnc:DVH.list-class]{DVH.list}}
}
\examples{
data("RadOnc")
gEUD(janedoe, 1) == unlist(janedoe$"Dmean")
gEUD(janedoe, Inf) == unlist(janedoe$"Dmax")
gEUD(janedoe, -Inf) == unlist(janedoe$"Dmin")
}
\keyword{methods}
\keyword{manip}
