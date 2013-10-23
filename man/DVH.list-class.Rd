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
    \item{[}{\code{signature(x = "DVH.list")}: Extract subset of DVH list }
    \item{[[}{\code{signature(x = "DVH.list")}: Extract single \code{DVH} object from DVH list }
    \item{[[<-}{\code{signature(x = "DVH.list")}: Replace single \code{DVH} object in DVH list }
    \item{[<-}{\code{signature(x = "DVH.list")}: Replace multiple elements of DVH list }
    \item{$}{\code{signature(x = "DVH.list")}: Extract dose or volume parameter(s) from \code{DVH} objects within DVH list. Note that multiple parameters may be specified by use of the comma (e.g. \code{\$"V20Gy,D5\%"}) }
    \item{c}{\code{signature(x = "DVH.list")}: Combine two or more DVH lists and/or \code{\link[RadOnc:DVH-class]{DVH}} objects }
    \item{lapply}{\code{signature(X = "DVH.list")}: Apply function to a list of \code{\link[RadOnc:DVH-class]{DVH}} objects }
    \item{length}{\code{signature(x = "DVH.list")}: Extract number of \code{\link[RadOnc:DVH-class]{DVH}} objects in DVH list }
    \item{mad}{\code{signature(x = "DVH.list")}: Compute the median absolute deviation (i.e. the median of the absolute deviations from the median) for all doses extracted from \code{DVH} object list }
    \item{max}{\code{signature(x = "DVH.list")}: Compute the maximum dose contained within the DVH list }
    \item{mean}{\code{signature(x = "DVH.list")}: Compute the mean \code{DVH} from all \code{DVH} objects within the list }
    \item{median}{\code{signature(x = "DVH.list")}: Compute the median \code{DVH} from all \code{DVH} objects within the list }
    \item{min}{\code{signature(x = "DVH.list")}: Compute the minimum dose contained within the DVH list }
    \item{names}{\code{signature(x = "DVH.list")}: Extract structure names for \code{DVH} objects in DVH list }
    \item{names<-}{\code{signature(x = "DVH.list")}: Assign structure name(s) to one or more \code{DVH} objects in DVH list }
    \item{print}{\code{signature(x = "DVH.list")}: Display summary of DVH list }
    \item{range}{\code{signature(x = "DVH.list")}: Compute the range (minimum and maximum doses) contained within the DVH list }
    \item{rev}{\code{signature(x = "DVH.list")}: Return a DVH list whose \code{DVH} objects are in reverse order}
    \item{sd}{\code{signature(x = "DVH.list")}: Compute the standard deviation for all doses extracted from \code{DVH} object list }
    \item{show}{\code{signature(x = "DVH.list")}: Display summary of DVH list }
    \item{t.test}{\code{signature(x = "DVH.list")}: Apply Student's t-Test to compare two DVH lists }
    \item{var}{\code{signature(x = "DVH.list")}: Compute the variance for all doses extracted from \code{DVH} object list }
    \item{wilcox.test}{\code{signature(x = "DVH.list")}: Apply Wilcoxon Rank Sum and Signed Rank Tests to compare two DVH lists }
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
