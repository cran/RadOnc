\name{DVH-class}
\Rdversion{1.1}
\docType{class}
\alias{DVH}
\alias{DVH-class}
\alias{[,DVH-method}
\alias{$,DVH-method}
\alias{$<-,DVH-method}
\alias{c,DVH-method}
\alias{initialize,DVH-method}
\alias{lines,DVH-method}
\alias{max,DVH-method}
\alias{mean,DVH-method}
\alias{min,DVH-method}
\alias{names,DVH-method}
\alias{names<-,DVH-method}
\alias{plot,DVH-method}
\alias{points,DVH-method}
\alias{print,DVH-method}
\alias{range,DVH-method}
\alias{show,DVH-method}

\title{Class \code{"DVH"}}
\description{
A data structure containing Dose-Volume Histogram (\code{DVH}) data and associated parameters for a single structure
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DVH", structure.name, structure.volume, type, dose.max, dose.min, dose.mean, dose.median, dose.mode, dose.STD, conf.index, equiv.sphere, gradient, dose.rx, dose.fx, doses, dose.type, volumes, volume.type, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{structure.name}:}{Name of the structure (e.g. "Stomach") }
    \item{\code{structure.volume}:}{Volume of the structure (in cubic centimeters) }
    \item{\code{type}:}{A character string specifying the DVH type (must be one of \code{"cumulative"} or \code{"differential"})}
    \item{\code{dose.max}:}{Maximum (point) dose contained within the structure}
    \item{\code{dose.min}:}{Minimum dose contained within the structure}
    \item{\code{dose.mean}:}{Mean dose to the structure}
    \item{\code{dose.median}:}{Median dose to the structure}
    \item{\code{dose.mode}:}{Modal dose to the structure}
    \item{\code{dose.STD}:}{Standard deviation of dose to the structure }
    \item{\code{conf.index}:}{Conformality index}
    \item{\code{equiv.sphere}:}{Equivalent sphere (diameter in centimeters)}
    \item{\code{gradient}:}{Dose gradient (in centimeters)}
    \item{\code{dose.rx}:}{Prescription dose (in cGy)}
    \item{\code{dose.fx}:}{Fractional dose (in cGy)}
    \item{\code{doses}:}{List of doses corresponding to \code{"volumes"} }
    \item{\code{dose.type}:}{A character string specifying the dose type (must be one of \code{"relative"} or \code{"absolute"})}
    \item{\code{volumes}:}{List of volumes corresponding to \code{"doses"} }
    \item{\code{volume.type}:}{A character string specifying the volume type (must be one of \code{"relative"} or \code{"absolute"})}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "DVH")}: Extract dose or volume parameter(s) from \code{DVH} object }
    \item{$}{\code{signature(x = "DVH")}: ... }
    \item{$<-}{\code{signature(x = "DVH")}: ... }
    \item{c}{\code{signature(x = "DVH")}: Combine multiple \code{DVH} objects into a single list }
    \item{lines}{\code{signature(x = "DVH")}: ... }
    \item{max}{\code{signature(x = "DVH")}: Extract maximum dose from \code{DVH} object }
    \item{mean}{\code{signature(x = "DVH")}: Extract mean dose from \code{DVH} object }
    \item{min}{\code{signature(x = "DVH")}: Extract minimum dose from \code{DVH} object }
    \item{names}{\code{signature(x = "DVH")}: Extract structure name from \code{DVH} object }
    \item{names<-}{\code{signature(x = "DVH")}: Assign structure name to \code{DVH} object }
    \item{plot}{\code{signature(x = "DVH")}: ... }
    \item{points}{\code{signature(x = "DVH")}: ... }
    \item{print}{\code{signature(x = "DVH")}: Display summary of \code{DVH} object }
    \item{range}{\code{signature(x = "DVH")}: Extract dose range from \code{DVH} object }
    \item{show}{\code{signature(x = "DVH")}: Display summary of \code{DVH} object }
	 }
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:DVH.list-class]{DVH.list}}, \code{\link{read.DVH}}, \code{\link[RadOnc:plot.DVH]{plot}}
}
\examples{
# Description of structure/slots in class
showClass("DVH")
}
\keyword{classes}
