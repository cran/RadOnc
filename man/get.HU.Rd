\name{get.HU}
\docType{methods}
\alias{get.HU}
\alias{get.HU-methods}
\alias{get.HU,list,array-method}
\alias{get.HU,RTdata,missing-method}
\alias{get.HU,structure.list,array-method}
\alias{get.HU,structure3D,array-method}
\alias{get.HU,ANY,missing-method}
\alias{get.HU,ANY,array-method}
\alias{get.HU,ANY,ANY-method}
\title{ Extract Hounsfield Unit (HU) data for one or more structure(s) from a CT image }
\description{
 Function to extract Hounsfield Unit (HU) data for one or more structure(s) from a CT image.
}
\usage{
\S4method{get.HU}{RTdata,missing}(x, CT, resolution.xyz=NA,	resolution.HU=1, method=NULL)

\S4method{get.HU}{structure3D,array}(x, CT, resolution.xyz=NA, resolution.HU=1, method=NULL)

\S4method{get.HU}{structure.list,array}(x, CT, resolution.xyz=NA, resolution.HU=1, method=NULL)
}
\arguments{
  \item{x}{
	One of either a \code{"structure3D"} object, a \code{"structure.list"}, or a \code{"RTdata"} object containing one or more structure(s).
}
  \item{CT}{
	Object of class \code{"array"} containing CT imaging data. This parameter is not required if \code{x} contains a \code{"RTdata"} object including CT imaging data.
}
  \item{resolution.xyz}{
A numeric list containing three positive values representing \emph{x}, \emph{y}, and \emph{z} resolution to be used for Hounsfield Unit (HU) interpolation from CT imaging data (specified in units of mm).  The smaller each value, the longer the overall calculation time (default is \code{c(NA, NA, NA)} which is interpreted as the native CT image resolution). \emph{z} resolution defaults to axial slice thickness when \code{method="axial"}.
}
  \item{resolution.HU}{
A positive numeric value specifying the HU calculation interval (bin width) for histogram calculation.

}
  \item{method}{
Character value specifying the Hounsfield Unit (HU) calculation algorithm to employ. Must be \code{NULL} or \code{"axial"}. }
}
\value{
Returns a single calculated \code{"histogram"} object or a \code{"list"} of \code{"histogram"} objects consisting of Hounsfield Unit values contained within the input \code{"structure3D"} or \code{"structure.list"}, respectively.
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:RTdata-class]{RTdata}}, \code{\link{hist}}
}
\keyword{ methods }
\keyword{ manip }
