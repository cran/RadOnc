\name{zDVH-class}
\Rdversion{1.1}
\docType{class}
\alias{zDVH}
\alias{zDVH-class}
\alias{[,zDVH-method}
\alias{$,zDVH-method}
\alias{$<-,zDVH-method}
\alias{c,zDVH-method}
\alias{initialize,zDVH-method}
\alias{max,zDVH-method}
\alias{mean,zDVH-method}
\alias{min,zDVH-method}
\alias{names,zDVH-method}
\alias{names<-,zDVH-method}
\alias{plot,zDVH-method}
\alias{print,zDVH-method}
\alias{show,zDVH-method}
\title{Class \code{"zDVH"}}
\description{
A data structure containing axially-segmented Dose-Volume Histogram (\code{zDVH}) data and associated parameters for a single structure
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("zDVH", patient, structure.name, structure.volume, type, dose.max, dose.min, dose.mean, dose.median, dose.mode, dose.STD, conf.index, equiv.sphere, gradient, dose.rx, dose.fx, doses, dose.type, dose.units, volumes, volume.type, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{patient}:}{Name of the patient (e.g. "Jane Doe") }
    \item{\code{ID}:}{Additional patient identifier or medical record number (e.g. "123456789") }
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
    \item{\code{dose.rx}:}{Prescription dose (in units specified by \code{dose.units})}
    \item{\code{dose.fx}:}{Fractional dose (in units specified by \code{dose.units})}
    \item{\code{rx.isodose}:}{Isodose line (\%) receiving prescription dose (default is \code{100\%})}
    \item{\code{doses}:}{List of doses corresponding to \code{"volumes"} }
    \item{\code{dose.type}:}{A character string specifying the dose type (must be one of \code{"relative"} or \code{"absolute"})}
    \item{\code{dose.units}:}{A character string specifying the dose units (must be one of \code{"cGy"} or \code{"Gy"})}
    \item{\code{volumes}:}{Matrix containing \emph{N}x\emph{M} dimensions, with \emph{N}\code{=length(doses)} and \emph{M} representing the number of axial segments for which of each a DVH was obtained }
    \item{\code{volume.type}:}{A character string specifying the volume type (must be one of \code{"relative"} or \code{"absolute"})}
  }
}
\section{Methods}{
  \describe{
    \item{[}{Extract dose or volume parameter(s) from \code{zDVH} object. Only one parameter may be specified at a time. Parameter specification should be of the form \code{["<A><B><C>"]}: \code{<A>} is equivalent to \code{"V"} or \code{"D"}, representing a volume or dose, respectively; \code{<B>} usually denotes a numerical value specifying the dose or volume; and \code{<C>} represents the dose or volume units (\code{"cGy"}, \code{"Gy"}, \code{"\%"}, or \code{"cc"}). An example would be \code{["V20Gy"]} which represents the volume of the structure receiving at least 20Gy dose. Specialized dosimetric keywords may also be used:  \code{"Dmax"} (maximum dose), \code{"Dmin"} (minimum dose), \code{"Dmean"} (mean dose), \code{"Dintegral"} (estimated integral dose), \code{"DRx"} (prescription dose), and \code{"volume"} (total structure volume).  If an improper parameter is specified however, \code{NA} results will be returned. See package documentation (vignette) for more details. }
    \item{$}{Extract a given parameter from a \code{zDVH} object}
    \item{$<-}{Assign a value to a given parameter within a \code{zDVH} object}
    \item{c}{Combine multiple \code{zDVH} objects into a single list }
    \item{max}{Extract maximum dose from \code{zDVH} object }
    \item{mean}{Extract mean dose from \code{zDVH} object }
    \item{min}{Extract minimum dose from \code{zDVH} object }
    \item{names}{Extract structure name from \code{zDVH} object }
    \item{names<-}{Assign structure name to \code{zDVH} object }
    \item{plot}{Plot \code{zDVH} object}
    \item{print}{Display summary of \code{zDVH} object }
    \item{range}{Extract dose range from \code{zDVH} object }
    \item{show}{Display summary of \code{zDVH} object }
    \item{sum}{Compute the total (summed) \code{DVH} from two or more \code{zDVH} and/or \code{DVH} objects. \emph{Note that structures are assumed to be non-overlapping; any overlaps in structure volumes may generate inaccurate dose summation.} }
	 }
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:DVH-class]{DVH}}, \code{\link[RadOnc:DVH.list-class]{DVH.list}}, \code{\link{read.DVH}}, \code{\link[RadOnc:plot.DVH]{plot}}
}
\examples{
# Description of structure/slots in class
showClass("zDVH")
}
\keyword{classes}
