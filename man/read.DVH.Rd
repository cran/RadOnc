\name{read.DVH}
\alias{read.DVH}
\title{
Read DVH data from input file
}
\description{
Function to extract DVH data from input file(s). Supported file formats include Aria (v8-13), DICOM-RT, CadPlan, Monaco, RayStation, and TomoTherapy.
}
\usage{
read.DVH(file, type=NA, verbose=TRUE, collapse=TRUE, modality=NULL)
}
\arguments{
  \item{file}{
Location of file containing relevant DVH information. One or more files may be specified. If multiple files are specified, all will be imported simultaneously.
}
  \item{type}{
Character vector specifying the DVH file format corresonding to each element in \code{file}. Value(s) must be one of \code{"aria8"}, \code{"aria10"}, \code{"aria11"}, \code{"aria13"}, \code{"aria15"}, \code{"dicom"}, \code{"cadplan"}, \code{"monaco"}, \code{"raystation"}, or \code{"tomo"} (default is \code{NA}, which will attempt to guess the appropriate file type for each DVH file). Note that multiple different types can be processed so long as the length of \code{type} corresponds to the number of files specified, or otherwise \code{type=NA}.
}
  \item{verbose}{
Single logical value indicating whether or not to output verbose information and status in text
}
  \item{collapse}{
Single logical value indicating whether or not to collapse input from multiple files into a single \code{DVH.list} with all structures (default) or to output a nested list of individual \code{DVH.list} objects with each one corresponding to a single input file
}
  \item{modality}{
When \code{type="dicom"}, this parameter is used to specify the scan type expected in importing DICOM-RT data (must be one one of \code{"CT"} [default] or \code{"MR"}).
}
}
\value{
Returns a single \code{DVH.list} object containing all DVH information from a single input file. If multiple files are specified, a list of \code{DVH.list} objects will be returned.
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
combined <- read.DVH(file=c(file1, file2), type="aria10", collapse=TRUE)
}
\keyword{ file }
