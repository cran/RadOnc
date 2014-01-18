\name{RTdata-class}
\Rdversion{1.1}
\docType{class}
\alias{RTdata-class}
\alias{$,RTdata-method}
\alias{$<-,RTdata-method}
\alias{names,RTdata-method}
\alias{names<-,RTdata-method}
\alias{print,RTdata-method}
\alias{show,RTdata-method}

\title{Class \code{"RTdata"}}
\description{
A data structure containing a corresponding CT image, dose grid, and structure set
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RTdata", name, CT, dose, structures, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Name of the dataset (e.g. "Patient XX") }
    \item{\code{CT}:}{Object of class \code{"array"} containing computed tomography (CT) image data in Houndsfield units }
    \item{\code{dose}:}{Object of class \code{"array"} containing dose grid data (pre-calculated using a separate treatment planning system)}
    \item{\code{structures}:}{Object of class \code{"structure.list"} containing the corresponding structure set, indexed to the CT and dose grid coordinate systems }
  }
}
\section{Methods}{
  \describe{
    \item{$}{Extract a given parameter from a \code{RTdata} object}
    \item{$<-}{Assign a value to a given parameter within a \code{RTdata} object}
    \item{names}{Extract dataset name from an \code{RTdata} object }
    \item{names<-}{Assign a name to an \code{RTdata} object}
    \item{print}{Display summary of \code{RTdata} object}
    \item{show}{Display summary of \code{RTdata} object}
	 }
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:structure.list-class]{structure.list}}
}
\examples{
showClass("RTdata")
}
\keyword{classes}
