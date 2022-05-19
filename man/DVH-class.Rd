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
\alias{sum,DVH-method}

\title{Class \code{"DVH"}}
\description{
A data structure containing Dose-Volume Histogram (\code{DVH}) data and associated parameters for a single structure
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DVH", patient, structure.name, structure.volume, type, dose.max, dose.min, dose.mean, dose.median, dose.mode, dose.STD, conf.index, equiv.sphere, gradient, plan.sum, dose.rx, dose.fx, doses, dose.type, dose.units, volumes, volume.type, ...)}.
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
    \item{\code{plan.sum}:}{Logical specifying whether or not data represents a plan sum (default is \code{FALSE})}
    \item{\code{dose.rx}:}{Prescription dose (in units specified by \code{dose.units})}
    \item{\code{dose.fx}:}{Fractional dose (in units specified by \code{dose.units})}
    \item{\code{rx.isodose}:}{Isodose line (\%) receiving prescription dose (default is \code{100\%})}
    \item{\code{doses}:}{List of doses corresponding to \code{"volumes"} }
    \item{\code{dose.type}:}{A character string specifying the dose type (must be one of \code{"relative"} or \code{"absolute"})}
    \item{\code{dose.units}:}{A character string specifying the dose units (must be one of \code{"cGy"} or \code{"Gy"})}
    \item{\code{volumes}:}{List of volumes corresponding to \code{"doses"} }
    \item{\code{volume.type}:}{A character string specifying the volume type (must be one of \code{"relative"} or \code{"absolute"})}
  }
}
\section{Methods}{
  \describe{
    \item{[}{Extract dose or volume parameter(s) from \code{DVH} object. Only one parameter may be specified at a time. Parameter specification should be of the form \code{["<A><B><C>"]}: \code{<A>} is equivalent to \code{"V"} or \code{"D"}, representing a volume or dose, respectively; \code{<B>} usually denotes a numerical value specifying the dose or volume; and \code{<C>} represents the dose or volume units (\code{"cGy"}, \code{"Gy"}, \code{"\%"}, or \code{"cc"}). An example would be \code{["V20Gy"]} which represents the volume of the structure receiving at least 20Gy dose. Dose ranges may also be specified, for instance \code{["V10-20Gy"]} or \code{["V<20Gy"]}. Specialized dosimetric keywords may also be used:  \code{"Dmax"} (maximum dose), \code{"Dmin"} (minimum dose), \code{"Dmean"} (mean dose), \code{"Dmedian"} (median dose), \code{"Dintegral"} (estimated integral dose), \code{"DRx"} (prescription dose), and \code{"volume"} (total structure volume).  If an improper parameter is specified however, \code{NA} results will be returned. See package documentation (vignette) for more details. }
    \item{$}{Extract a given parameter from a \code{DVH} object}
    \item{$<-}{Assign a value to a given parameter within a \code{DVH} object}
    \item{c}{Combine multiple \code{DVH} objects into a single list }
    \item{lines}{Plot \code{DVH} object (see \code{link{lines}})}
    \item{max}{Extract maximum dose from \code{DVH} object }
    \item{mean}{Extract mean dose from \code{DVH} object }
    \item{min}{Extract minimum dose from \code{DVH} object }
    \item{names}{Extract structure name from \code{DVH} object }
    \item{names<-}{Assign structure name to \code{DVH} object }
    \item{plot}{Plot \code{DVH} object}
    \item{points}{Plot \code{DVH} object (see \code{link{points}}) }
    \item{print}{\code{signature(x = "DVH")}: Display summary of \code{DVH} object }
    \item{range}{\code{signature(x = "DVH")}: Extract dose range from \code{DVH} object }
    \item{show}{\code{signature(x = "DVH")}: Display summary of \code{DVH} object }
    \item{sum}{\code{signature(x = "DVH")}: Compute the total (summed) \code{DVH} from two or more \code{DVH} objects. \emph{Note that structures are assumed to be non-overlapping; any overlaps in structure volumes may generate inaccurate dose summation.} }
	 }
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:DVH.list-class]{DVH.list}}, \code{\link{read.DVH}}
}
\examples{
# Description of structure/slots in class
showClass("DVH")
}
\keyword{classes}
