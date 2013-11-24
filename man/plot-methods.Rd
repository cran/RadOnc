\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY,ANY-method}
\alias{plot,DVH.list,ANY-method}
\alias{plot,DVH.list,missing-method}
\alias{plot,DVH,ANY-method}
\alias{plot,DVH,DVH-method}
\alias{plot,DVH,missing-method}
\alias{plot.DVH}
\alias{plot.DVH.list}
\title{ Additional Methods for Function \code{plot} in Package \pkg{graphics} }
\description{
 These additional methods extend function \code{plot} in package \pkg{graphics}, enabling DVH visualization in package \pkg{RadOnc}
}
\usage{
\method{plot}{DVH}(x, \dots, plot.type = c("individual", "grouped", "ttest", "wilcox", "bars"))

\method{plot}{DVH.list}(x, \dots, plot.type = c("individual", "grouped", "ttest", "wilcox", "bars"))
}
\arguments{
  \item{x}{
	Can represent either an object of class \code{DVH} or \code{DVH.list}
}
  \item{\dots}{
	May contain any number of objects of class \code{DVH} or \code{DVH.list}. Note that at least one object of either class is required if \code{plot.type = "ttest"} or \code{"wilcox"}.
}
  \item{plot.type}{
	A character string specifying the type of plot to generate, must be one of \code{"individual"} (default), \code{"grouped"}, \code{"ttest"}, \code{"wilcox"} or \code{"bars"}	
}
}
\section{Graphical Parameters}{
	\describe{
		\item{\code{main}}{
			An overall title for the plot (see \code{\link{title}}, default is \code{""}).
		}
		\item{\code{col}}{
			The color to be used for drawing lines. This parameter can also represent a list of different colors, with length corresponding to the number of \code{DVH} or \code{DVH.list} objects (if \code{plot.type = "individual"}, number of colors should correspond to number of individual \code{DVH} objects).
			
			When \code{plot.type = "bars"}, \code{col} instead represents a list of colors used to define the shading scheme applied to the entire dose range. Default color scheme is a standard rainbow, with colors ranging from \code{"blue"} to \code{"red"}. Actual default is specified as \code{rev(rainbow(n=10, start=0, end=2/3))}. See \code{\link{colors}} for more options and information.
		}
		\item{\code{lty}}{
			The line type. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings \code{"blank"}, \code{"solid"}, \code{"dashed"}, \code{"dotted"}, \code{"dotdash"}, \code{"longdash"}, or \code{"twodash"}, where \code{"blank"} uses 'invisible lines' (i.e., does not draw them).

			Alternatively, a string of up to 8 characters (from c(1:9, "A":"F")) may be given, giving the length of line segments which are alternatively drawn and skipped (see \code{\link{par}} for more details). Note that this paramater can also represent a list of different line types, with length corresponding to the number of \code{DVH} or \code{DVH.list} objects (if \code{plot.type = "individual"}, number of line types should correspond to number of individual \code{DVH} objects).
		}
		\item{\code{lwd}}{
			The line width, a emph{positive} number, defaulting to \code{1}. The interpretation is device-specific, and some devices do not implement line widths less than one (see the help on the device for details of the interpretation). Note that this paramater can also represent a list of different line widths, with length of the list corresponding to the number of \code{DVH} or \code{DVH.list} objects (if \code{plot.type = "individual"}, number of line widths should correspond to number of individual \code{DVH} objects).
		}
		\item{\code{line.transparency}}{
			Factor modifying the transparency of line drawings (value must be between \code{[0,1]}), with \code{fill.transparency = 0} specifying complete transparency and \code{1} (default) specifying complete opacity. Note that this parameter can represent a list of different transparency factors, with length corresponding to the number of \code{DVH} or \code{DVH.list} objects  (if \code{plot.type = "individual"}, number of transparency factors should correspond to number of individual \code{DVH} objects).
		}
		\item{\code{new}}{
			Logical, defaulting to \code{TRUE}. If set to \code{FALSE}, the next high-level plotting command (actually \code{\link{plot.new}}) should \emph{not clean} the frame before drawing \emph{as if it were on a \strong{new} device}. It is an error (ignored with a warning) to try to use \code{new = FALSE} on a device that does not currently contain a high-level plot.
		}
		\item{\code{fill}}{
			Logical, defaulting to \code{TRUE}. Determines whether or not to draw shading around groups of DVH lines. Note that this parameter is relevant when \code{plot.type = "grouped"}, \code{"ttest"} or \code{"wilcox"} (it is ignored when \code{plot.type = "individual"} or \code{"bars"}).
		}	
		\item{\code{density}}{
			The density of shading lines, in lines per inch. The default value of \code{NULL} means that no shading lines are drawn. A zero value of \code{density} means no shading nor filling, whereas negative values and \code{NA} suppress shading lines but allow solid color filling. Note that this paramater can also represent a list of different densities, with length of the list corresponding to the number of \code{DVH} or \code{DVH.list} objects.
		}
		\item{\code{angle}}{
			The slope of shading lines, given as an angle in counter-clockwise degrees (default is \code{45}). This parameter is ignored if \code{density = NULL} or \code{NA} or wherever \code{fill = FALSE} or is ignored (e.g. when \code{plot.type = "individual"} or \code{"bars"}). Also note that this paramater can represent a list of different angles, with length of the list corresponding to the number of \code{DVH} or \code{DVH.list} objects.
		}		
		\item{\code{fill.lty}}{
			The line type used for shading (per \code{density} and \code{angle} parameters). See \code{lty} above for further specifications of line types. Note that this paramater can represent a list of different line types, with length corresponding to the number of \code{DVH} or \code{DVH.list} objects.
		}
		\item{\code{fill.transparency}}{
			Factor modifying the transparency of filling/shading (value must be between \code{[0,1]}), with \code{fill.transparency = 0} specifying complete transparency and \code{1} specifying complete opacity. Note that this parameter can represent a list of different transparency factors, with length corresponding to the number of \code{DVH} or \code{DVH.list} objects. Default value is \code{line.transparency/2}.
		}
		\item{\code{legend}}{
			Value specifying the location to draw a figure legend, must be one of \code{"topright"}, \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"right"}, or \code{"center"}. These keywords place the legend on the inside of the plot frame at the location specified. Partial argument matching is used. If \code{legend = NA} (default), the legend is not drawn.
			
			When \code{plot.type = "bars"}, \code{legend} is interpreted instead as a logical value and must be either \code{TRUE} (default) or \code{FALSE}.	If \code{TRUE}, labels are drawn at the bottom of the plot, outside of the actual plotting frame.		
		}
		\item{\code{legend.labels}}{
			A character or \code{\link{expression}} vector specifying the text to appear in the legend, when relevant. Length of \code{legend.labels} must match the numer of \code{DVH} and \code{DVH.list} objects, otherwise legend text will default to the form: \code{"Group 1", "Group 2", ...}
		}
	}
}
\section{Additional Parameters}{
	\describe{
		\item{\code{type}}{
			Value specifying type of dose-volume histogram to plot (must be one of \code{"cumulative"} or \code{"default"}).
		}
		\item{\code{dose}}{
			Value specifying dose scale (must be one of \code{"relative"} or \code{"absolute"}).
		}
		\item{\code{dose.units}}{
			Value specifying dose units (must be one of \code{"cGy"} or \code{"Gy"}).
		}
		\item{\code{volume}}{
			Value specifying volume scale (must be one of \code{"relative"} or \code{"absolute"}).
		}
	}
}
\details{
	More details about each plotting subtype here . . .
}
\examples{
	# plot(x, ..., plot.type="individual")
	OARs <- c("LIVER","STOMACH","DUODENUM")
	cols <- c("red", "darkgreen", "blue")
	plot(johndoe[OARs], col=cols, lty=1:3, legend="topright", legend.labels=OARs, main="OARs")
	
	# plot(x, ..., plot.type="bars")
	plot(janedoe[2:9], plot.type="bars", volume="absolute",dose="relative")
	
	# plot(x, ..., plot.type="grouped")
	plot(c(johndoe["STOMACH"],janedoe["STOMACH"]), #group 1
		c(johndoe[c("CTV","PTV")],janedoe[c("CTV","PTV")]), #group 2
		c(janedoe["LIVER"],johndoe["LIVER"]), #group 3
		c(johndoe["DUODENUM"],janedoe["DUODENUM"]), #group 4
		plot.type="grouped", col=c("orange","green","blue","red"), lwd=2, dose="relative")
		
	# plot(x, ..., plot.type="ttest")
	# OARs <- c("LIVER", "STOMACH", "SMALL_BOWEL")
	# plot(c(johndoe[c("CTV","PTV")],janedoe[c("CTV","PTV")]), #group 1
	#	c(janedoe[OARs],johndoe[OARs]), #group 2
 	#	plot.type="t", col=c("red", "blue"), lty=2, fill.lty=1, main="Target v. OAR t-Test")
 		
	# plot(x, ..., plot.type="wilcox")
	# plot(c(johndoe[c("CTV","PTV")],janedoe[c("CTV","PTV")]), #group 1
	#	c(janedoe[OARs],johndoe[OARs]), #group 2
 	#	plot.type="w", col=c("red", "blue"), lty=2, fill.lty=1, main="Target v. OAR")
}
\seealso{
\code{\link{plot}}
}
\keyword{methods}
\keyword{aplot}
