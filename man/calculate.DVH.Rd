\name{calculate.DVH}
\docType{methods}
\alias{calculate.DVH}
\alias{calculate.DVH-methods}
\alias{calculate.DVH,RTdata,missing-method}
\alias{calculate.DVH,RTdata,array-method}
\alias{calculate.DVH,structure3D,array-method}
\alias{calculate.DVH,structure.list,array-method}
\alias{calculate.DVH,ANY,missing-method}
\alias{calculate.DVH,ANY,array-method}
\alias{calculate.DVH,ANY,ANY-method}
\title{
Calculate DVH data for one or more structure(s) using a dose grid
}
\description{
Function to calculate DVH data for one or more structure(s) using a dose grid.
}
\usage{
\S4method{calculate.DVH}{RTdata,missing}(x, dose, resolution.xyz=c(0.2,0.2,NA), 
	resolution.dose=0.01, method=c("ATC", "surface"), dose.units=c("cGy","Gy"))

\S4method{calculate.DVH}{RTdata,array}(x, dose, resolution.xyz=c(0.2,0.2,NA), 
	resolution.dose=0.01, method=c("ATC", "surface"), dose.units=c("cGy","Gy"))

\S4method{calculate.DVH}{structure3D,array}(x, dose, resolution.xyz=c(0.2,0.2,NA), 
	resolution.dose=0.01, method=c("ATC", "surface"), dose.units=c("cGy","Gy"))

\S4method{calculate.DVH}{structure.list,array}(x, dose, resolution.xyz=c(0.2,0.2,NA), 
	resolution.dose=0.01, method=c("ATC", "surface"), dose.units=c("cGy","Gy"))
}
\arguments{
  \item{x}{
	One of either a \code{"structure3D"} object, a \code{"structure.list"}, or a \code{"RTdata"} object containing one or more structure(s).
}
  \item{dose}{
Object of class \code{"array"} containing dose grid data. This parameter is not required if \code{x} contains a \code{"RTdata"} object including dose grid information.
}
  \item{resolution.xyz}{
A numeric list containing three positive values representing \emph{x}, \emph{y}, and \emph{z} resolution to be used for DVH calculation (specified in units of mm).  The smaller each value, the longer the overall calculation time (default is \code{c(0.2, 0.2, NA)}). \emph{z} resolution defaults to axial slice thickness when \code{method="ATC"}.
}
  \item{resolution.dose}{
A positive numeric value specifying the dose calculation interval (bin width) for DVH calculation (values interpreted in Gy).
}
  \item{method}{
Character value specifying the DVH calculation algorithm to employ. Must be one of \code{"ATC"} or \code{"surface"}. See below for further details about the dose calculation algorithms.
}
  \item{dose.units}{
Value specifying units of dose. Must be one of \code{"cGy"} (default) or \code{"Gy"}.}
}
\value{
Returns a single calculated \code{"DVH"} object or a \code{"DVH.list"} object consisting of multiple DVHs where relevant.
}
\references{
Straube, W., Matthews, J., Bosch, W., and Purdy, J.A. (2005)
DVH Analysis: Consequences for Quality Assurance of Multi-Institutional Clinical Trials.
\emph{Med Phys} \bold{32}(6):2021.
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:DVH-class]{DVH}}, \code{\link[RadOnc:DVH.list-class]{DVH.list}}, \code{\link{new}}
}
\examples{
# Read two DVH example files
# need to add DICOM data to RadOnc.RData so have a dose grid in there!!
# add further details section to discuss algorithm!!  should also have citation for ATC stuff!!
}
\keyword{ methods }
\keyword{ manip }