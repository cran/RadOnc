\name{structure3D-class}
\Rdversion{1.1}
\docType{class}
\alias{structure3D-class}
\alias{$,structure3D-method}
\alias{$<-,structure3D-method}
\alias{c,structure3D-method}
\alias{dim,structure3D-method}
\alias{names,structure3D-method}
\alias{names<-,structure3D-method}
\alias{plot,structure3D,ANY-method}
\alias{plot,structure3D,missing-method}
\alias{print,structure3D-method}
\alias{range,structure3D-method}
\alias{show,structure3D-method}

\title{Class \code{"structure3D"}}
\description{
A data structure containing 3D volumetric (\code{structure3D}) data and associated parameters for a single structure object
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("structure3D", name, volume, volume.units, coordinate.units, vertices, origin, triangles, DVH)}.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Name of the structure (e.g. "Stomach") }
    \item{\code{volume}:}{Volume of the structure (in \code{volume.units}) }
    \item{\code{volume.units}:}{A character string specifying the units of volume (must be \code{"cc"})}
    \item{\code{coordinate.units}:}{A character string specifying the units of the coordinate system (must be one of \code{"cm"} or \code{"mm"})}
    \item{\code{vertices}:}{Specifies all points defining the given structure, with data contained in a \emph{Nx3} matrix where each column represents \emph{x}, \emph{y}, and \emph{z} coordinates, respectively, for each of \emph{N} points}
    \item{\code{origin}:}{A numeric vector containing the \emph{x}, \emph{y}, and \emph{z} coordinates corresponding to the center of the structure}
    \item{\code{triangles}:}{Specifies a complete combination of points which define the triangular mesh surface of the structure; data is contained in a \emph{3xN} matrix where each row represents one of the three vertices of a each triangle within the mesh, with \emph{N} columns corresponding to the number of triangles; note that the value of each element in \code{triangles} references a point in \code{vertices}, such that all values in \code{triangles} must be between \code{1} and \code{dim(vertices)[1]}}
    \item{\code{closed.polys}:}{Specifies a complete set of axially-defined closed polygons which, together, comprise the surface of the structure; data is contained in a list of \emph{Nx3} matrices where each row represents the \emph{x}, \emph{y}, and \emph{z} coordinates of a point in a single polygon}
    \item{\code{DVH}:}{Object of class \code{\link[RadOnc:DVH-class]{DVH}}}
  }
}
\section{Methods}{
  \describe{
    \item{$}{Extract a given parameter from a \code{structure3D} object}
    \item{$<-}{Assign a value to a given parameter within a \code{structure3D} object}
    \item{c}{Combine multiple \code{structure3D} objects into a single \code{structure.list} }
    \item{dim}{Extract dimensions (number of vertices and axial slices) from \code{structure3D} object }
    \item{names}{Extract structure name from \code{structure3D} object  }
    \item{names<-}{Assign structure name to a \code{structure3D} object }
    \item{plot}{Plot \code{structure3D} object}
    \item{print}{ Display summary of \code{structure3D} object }
    \item{range}{ Extract coordinate range from \code{structure3D} object }
    \item{show}{ Display summary of \code{structure3D} object }
	 }
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:DVH-class]{DVH}}
}
\examples{
# Description of structure/slots in class
showClass("structure3D")

# data(list="RadOnc", package="RadOnc")
# plot(mandible)
# points3d(range(mandible), col="red")
}
\keyword{classes}
