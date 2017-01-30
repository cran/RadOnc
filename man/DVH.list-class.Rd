\name{DVH.list-class}
\Rdversion{1.1}
\docType{class}
\alias{DVH.list}
\alias{DVH.list-class}
\alias{[,DVH.list-method}
\alias{[[,DVH.list-method}
\alias{[[<-,DVH.list-method}
\alias{$,DVH.list-method}
\alias{as.list,DVH.list-method}
\alias{c,DVH.list-method}
\alias{initialize,DVH.list-method}
\alias{lapply,DVH.list-method}
\alias{length,DVH.list-method}
\alias{mad,DVH.list-method}
\alias{max,DVH.list-method}
\alias{mean,DVH.list-method}
\alias{median,DVH.list-method}
\alias{min,DVH.list-method}
\alias{names,DVH.list-method}
\alias{names<-,DVH.list-method}
\alias{plot,DVH.list-method}
\alias{print,DVH.list-method}
\alias{quantile,DVH.list-method}
\alias{range,DVH.list-method}
\alias{rev,DVH.list-method}
\alias{sd,DVH.list-method}
\alias{show,DVH.list-method}
\alias{sum,DVH.list-method}
\alias{t.test,DVH.list-method}
\alias{var,DVH.list-method}
\alias{wilcox.test,DVH.list-method}

\title{Class \code{"DVH.list"}}
\description{
A data structure containing one or more Dose-Volume Histogram (\code{\link[RadOnc:DVH-class]{DVH}}) objects
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DVH.list", structures, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{structures}:}{List of \code{DVH} objects }
  }
}
\section{Methods}{
  \describe{
    \item{[}{Extract subset of DVH list based on pattern matching with structure names (regular expressions may be specified if desired, see \code{\link{regex}} for more details). Note that case-insensitive matching is supported using the \code{"(?i)"} prefix, such that \code{DVHs["(?i)LIVER"]} will match any combination of upper and lowercase letters spelling "liver". }
    \item{[[}{Extract single \code{DVH} object from DVH list }
    \item{[[<-}{Replace single \code{DVH} object in DVH list }
    \item{[<-}{Replace multiple elements of DVH list }
    \item{$}{Extract dose or volume parameter(s) from \code{DVH} objects within DVH list. Note that multiple parameters may be specified by use of the comma (e.g. \code{\$"V20Gy,D5\%"}). See \code{\link[RadOnc:DVH-class]{DVH}} documentation (\code{[]} usage) for additional details. Note that this functionality can also be used to extract patient name(s) and/or identifier(s) (e.g. \code{\$patients} or \code{$ID}).}
    \item{as.list}{Convert a \code{DVH.list} object to a list containing individual \code{DVH} objects. Note that the reverse conversion can be performed using the \code{\link{as}} command and specifying \code{class="DVH.list"}.}
    \item{c}{Combine two or more DVH lists and/or \code{\link[RadOnc:DVH-class]{DVH}} objects }
    \item{lapply}{Apply function to a list of \code{\link[RadOnc:DVH-class]{DVH}} objects }
    \item{length}{Extract number of \code{\link[RadOnc:DVH-class]{DVH}} objects in DVH list }
    \item{mad}{Compute the median absolute deviation (i.e. the median of the absolute deviations from the median) for all doses extracted from \code{DVH} object list }
    \item{max}{Compute the maximum dose contained within the DVH list }
    \item{mean}{Compute the mean \code{DVH} from all \code{DVH} objects within the list }
    \item{median}{Compute the median \code{DVH} from all \code{DVH} objects within the list }
    \item{min}{Compute the minimum dose contained within the DVH list }
    \item{names}{Extract structure names for \code{DVH} objects in DVH list }
    \item{names<-}{Assign structure name(s) to one or more \code{DVH} objects in DVH list }
    \item{print}{Display summary of DVH list }
    \item{range}{Compute the range (minimum and maximum doses) contained within the DVH list }
    \item{rev}{Return a DVH list whose \code{DVH} objects are in reverse order}
    \item{sd}{Compute the standard deviation for all doses extracted from \code{DVH} object list }
    \item{show}{Display summary of DVH list }
    \item{sum}{Compute the total (summed) \code{DVH} from all \code{DVH} objects within the list. \emph{Note that structures are assumed to be non-overlapping; any overlaps in structure volumes may generate inaccurate dose summation.} }    
    \item{t.test}{Apply Student's t-Test to compare two DVH lists  (see also \code{\link{t.test}})}
    \item{var}{Compute the variance for all doses extracted from \code{DVH} object list }
    \item{wilcox.test}{Apply Wilcoxon Rank Sum and Signed Rank Tests to compare two DVH lists (see also \code{\link{wilcox.test}}) }
	 }
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:DVH-class]{DVH}}, \code{\link[RadOnc:plot.DVH.list]{plot}}
}
\examples{
# Description of structure/slots in class
showClass("DVH.list")
data(list="RadOnc", package="RadOnc")
print(johndoe)
plot(johndoe[c("LIVER", "PTV")], plot.type="i", col=c("red","blue"), lty=1:2, lwd=1:2)
}
\keyword{classes}
