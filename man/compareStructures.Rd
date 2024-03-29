\name{compareStructures}
\alias{compareStructures}
\alias{compareStructures.axial}
\alias{compareStructures.EMD}
\alias{compareStructures.hausdorff}
\alias{compareStructures.grid}
\alias{compareStructures.surface}
\title{
Compare Axially-Defined Three-Dimensional Structures
}
\description{
Three-dimensional structural comparison of a \code{structure.list} composed of two or more \code{structure3D} objects.
}
\usage{
compareStructures(structures, method = NULL,
	hausdorff.method = NULL, verbose = TRUE, plot = TRUE, pixels = 100)
}
\arguments{
  \item{structures}{
	An object of class \code{\link[RadOnc:structure.list-class]{structure.list}} containing two or more \code{\link[RadOnc:structure3D-class]{structure3D}} objects.
}
  \item{method}{
One of \code{"axial"}, \code{"surface"}, \code{"hausdorff"}, \code{"DSC"}, or \code{"EMD"} specifying the comparison between structures. For the \code{"axial"} method (default), the volumetric region of comparison is divided into discrete evenly-spaced voxels and the degree of overlap is measured across each voxel in axial cross-sections. For the \code{"surface"} method, overlap is calculated for each point on each structure's surface. For the \code{"hausdorff"} method, Hausdorff distances are calculated for each pairwise combination of structures as the distance between two point clouds. The \code{"EMD"} method is currently not currently supported (in development), but will return a modified earth mover's distance between two structures, accounting for radiation dose. }
  \item{hausdorff.method}{
One of \code{"mean"}, \code{"median"}, \code{"95"}, or \code{"absolute"}, specifying the method to use for Hausdorff distance measurement when \code{method = "hausdorff"}. The \code{"absolute"} Hausdorff distance yields the maximum distance required to connect any point from one point cloud to its closest neighbor in the other. This metric is highly subject to outliers, thus an aggregate metric is implemented by selecting the average (\code{"mean"}, equivalent to the mean surface distance) or median (\code{"median"}) distance required to connect all points in one point cloud to their closest neighboring points in the other. When value is \code{"95"}, the 95-percentile Hausdorff distance is calculated. Note that the Hausdorff distance between two completely superimposable point clouds is zero. This parameter is ignored for other values of \code{method}.}
  \item{verbose}{
Logical value (default is \code{TRUE}) indicating whether or not to print status updates to the active command line interface or workspace.
}
  \item{plot}{
Logical value (default is \code{TRUE}) indicating whether or not to display graphical comparison of structures. This parameter is only relevant when \code{method = "axial"}, as it is otherwise ignored. }
  \item{pixels}{
Integer value (default is \code{100}) specifying the number of pixels in the x- and y-directions used for comparative calculations when \code{method = "axial"}. }
}
\value{
When \code{method = "hausdorff"}, the function returns a \emph{N}x\emph{N} matrix of pairwise Hausdorff distances, where \emph{N} is the number of structures in \code{structures}. When \code{method = "axial"} or \code{method = "surface"}, the function returns a \emph{M}x(\emph{N}+3) dimensional matrix representing a list of \emph{M} points in three-dimensional space and a table of values denoting whether each point is internal or external to each structure (note that the first three columns in the matrix represent the \emph{xyz} coordinates for each point and values of \code{1} and \code{0} represent internal and external points, respectively).
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:structure.list-class]{structure.list}}, \code{\link[RadOnc:structure3D-class]{structure3D}}
}
\examples{
data("RadOnc")
compareStructures(teeth, method="DSC")
# teeth.compare <- compareStructures(teeth, method="axial", plot=TRUE)
# compareStructures(teeth, method="hausdorff", hausdorff.method="mean")
}
\keyword{ hplot }
\keyword{ manip }
