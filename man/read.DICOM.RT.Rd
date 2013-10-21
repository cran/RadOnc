\name{read.DICOM.RT}
\alias{read.DICOM.RT}
\title{
Read DICOM-RT data from an input directory
}
\description{
Function to extract three-dimensional structural data from input file(s) DICOM-RT formats
}
\usage{
read.DICOM.RT(path, exclude = NULL, recursive = TRUE, verbose = TRUE, limit = NULL, ...)
}
\arguments{
  \item{path}{
Location of directory containing relevant DICOM-RT information
}
  \item{exclude}{
Character string specifying exclusion criteria used to skip processing of specific file(s) within the DICOM-RT directory.
}
  \item{recursive}{
Logical (default is \code{TRUE}) specifying whether to recursively extract nested files located within a parent directory
}
  \item{verbose}{
Logical value (default is \code{TRUE}) indicating whether or not to output verbose information and status in text
}
  \item{limit}{
Numerical value specifying size threshold over which to skip processing of specific structures (number corresponds to the number of points used to define each structure surface)
}
  \item{\dots}{
Additional parameters passed to \code{\link{readDICOM}} function
}
}
\value{
Returns a single \code{\link[RadOnc:structure.list-class]{structure.list}} object containing all structures from input DICOM-RT file(s)
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link{readDICOM}}, \code{\link[RadOnc:structure.list-class]{structure.list}}, \code{\link[RadOnc:structure3D-class]{structure3D}}
}
\keyword{ file }