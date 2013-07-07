\name{read.DVH}
\alias{read.DVH}
\title{
Read DVH data from input file
}
\description{
Function to extract DVH data from input file(s) in Aria (v8-11) or DICOM-RT formats
}
\usage{
read.DVH(file, type = c(NA, "aria10", "aria11", "aria8", "dicom"), verbose = TRUE)
}
\arguments{
  \item{file}{
Location of file containing relevant DVH information
}
  \item{type}{
Must be one of \code{"aria8"}, \code{"aria10"}, \code{"aria11"} or \code{"dicom"} (default is \code{NA}).
}
  \item{verbose}{
Logical value indicating whether or not to output verbose information and status in text
}
}
\value{
Returns a single \code{DVH.list} object containing all DVH information from input file
}
\author{
Reid F. Thompson (\email{reid.thompson@gmail.com})
}
\seealso{
\code{\link[RadOnc:DVH-class]{DVH}}, \code{\link[RadOnc:DVH.list-class]{DVH.list}}, \code{\link{new}}
}
\examples{
# Read two DVH example files
file1 <- system.file("extdata/John_Doe.dvh", package="RadOnc")
johndoe <- read.DVH(file=file1, type="aria10", verbose=TRUE)
file2 <- system.file("extdata/Jane_Doe.dvh", package="RadOnc")
janedoe <- read.DVH(file=file2, type="aria10", verbose=TRUE)
}
\keyword{ file }
