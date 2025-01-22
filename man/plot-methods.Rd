\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY,ANY-method}
\alias{plot,DVH.list,ANY-method}
\alias{plot,DVH.list,missing-method}
\alias{plot,DVH,ANY-method}
\alias{plot,DVH,DVH-method}
\alias{plot,DVH,missing-method}
\alias{plot,zDVH,ANY-method}
\alias{plot,zDVH,DVH-method}
\alias{plot,DVH,zDVH-method}
\alias{plot,zDVH,zDVH-method}
\alias{plot.DVH}
\alias{plot.DVH.list}
\alias{plot.zDVH}
\alias{plot,zDVH,missing-method}
\title{ Additional Methods for Function \code{plot} in Package \pkg{graphics} }
\description{
 These additional methods extend function \code{plot} in package \pkg{graphics}, enabling DVH visualization in package \pkg{RadOnc}
}
\usage{
\method{plot}{DVH}(x, \dots, plot.type = NULL)

\method{plot}{DVH.list}(x, \dots, plot.type = NULL)
}
\arguments{
  \item{x}{
	Can represent either an object of class \code{DVH}, \code{zDVH}, or \code{DVH.list}
}
  \item{\dots}{
	May contain any number of objects of class \code{DVH} or \code{DVH.list}. Note that at least one object of either class is required if \code{plot.type = "ttest"} or \code{"wilcox"}. Specification of more than one \code{zDVH} will reduce all input \code{zDVH} objects to \code{DVH} objects for standard plotting.
}
  \item{plot.type}{
	A character string specifying the type of plot to generate, must be one of \code{"individual"} (default), \code{"grouped"}, \code{"ttest"}, \code{"wilcox"}, \code{"bars"}, or \code{"correlation"}. Parameter \code{plot.type} will be ignored if \code{x} is a \code{zDVH} object.
}
}
\section{Graphical Parameters}{
	\describe{
		\item{\code{alpha}}{
			When \code{plot.type="ttest"}, \code{"wilcox"}, or \code{"correlation"}, this specifies the numerical value (threshold) used to display significant p-values.
		}		
		\item{\code{angle}}{
			The slope of shading lines, given as an angle in counter-clockwise degrees (default is \code{45}). This parameter is ignored if \code{density=NULL} or \code{NA} or wherever \code{fill=FALSE} or is ignored (e.g. when \code{plot.type="individual"} or \code{"bars"}). Also note that this paramater can represent a list of different angles, with length of the list corresponding to the number of \code{DVH} or \code{DVH.list} objects.
		}		
		\item{\code{back}}{
			A character value specifying the back-facing surface fill method when \code{x} is a \code{zDVH} object (the parameter is otherwise ignored). Value must be one of \code{"filled"}, \code{"lines"}, \code{"points"}, or \code{"culled"}. Default value is assigned from input parameter \code{front}.
		}		
		\item{\code{center}}{
			A character value specifying where to plot the center of a group when \code{plot.type="grouped"} (the parameter is otherwise ignored). Value must be one of \code{"mean"} or \code{"median"}.
			\emph{Note:  As of package v.1.0.3, this parameter is considered defunct and will be disregarded accordingly}.
		}		
		\item{\code{col}}{
			The color to be used for drawing lines. This parameter can also represent a list of different colors, with length corresponding to the number of \code{DVH}, \code{zDVH}, or \code{DVH.list} objects (if \code{plot.type="individual"}, number of colors should correspond to number of individual \code{DVH} objects).
			
			When \code{plot.type="bars"}, \code{col} instead represents a list of colors used to define the shading scheme applied to the entire dose range. Default color scheme is a standard rainbow, with colors ranging from \code{"blue"} to \code{"red"}. Actual default is specified as \code{rev(rainbow(n=10, start=0, end=2/3))}. See \code{\link{colors}} for more options and information.
		}
		\item{\code{density}}{
			The density of shading lines, in lines per inch. The default value of \code{NULL} means that no shading lines are drawn. A zero value of \code{density} means no shading nor filling, whereas negative values and \code{NA} suppress shading lines but allow solid color filling. Note that this paramater can also represent a list of different densities, with length of the list corresponding to the number of \code{DVH} or \code{DVH.list} objects.
		}
		\item{\code{fill}}{
			Logical, defaulting to \code{TRUE}. Determines whether or not to draw shading around groups of DVH lines. Note that this parameter is relevant when \code{plot.type="grouped"}, \code{"ttest"} or \code{"wilcox"} (it is ignored when \code{plot.type="individual"} or \code{"bars"}).
		}	
		\item{\code{fill.lty}}{
			The line type used for shading (per \code{density} and \code{angle} parameters). See \code{lty} above for further specifications of line types. Note that this paramater can represent a list of different line types, with length corresponding to the number of \code{DVH} or \code{DVH.list} objects.
		}
		\item{\code{fill.transparency}}{
			Factor modifying the transparency of filling/shading (value must be between \code{[0,1]}), with \code{fill.transparency=0} specifying complete transparency and \code{1} specifying complete opacity. Note that this parameter can represent a list of different transparency factors, with length corresponding to the number of \code{DVH} or \code{DVH.list} objects. Default value is \code{line.transparency/2}.
		}
		\item{\code{front}}{
			A character value specifying the front-facing surface fill method when \code{x} is a \code{zDVH} object (the parameter is otherwise ignored). Value must be one of \code{"filled"} (default), \code{"lines"}, \code{"points"}, or \code{"culled"}.
		}		
		\item{\code{grid}}{
			Logical, defaulting to \code{FALSE}. Determines whether or not to draw gridlines on 2D plot. Note that this parameter is relevant when \code{plot.type="individual"}, \code{"grouped"}, \code{"ttest"}, or \code{"wilcox"} (it is ignored when \code{plot.type="bars"} or when plotting \code{zDVH} objects).
		}		
		\item{\code{highlight}}{
			The color to be used for shading p-value data when \code{plot.type="ttest"} or \code{"wilcox"}. See \code{\link{colors}} for more options and information.
		}		
		\item{\code{legend}}{
			Value specifying the location to draw a figure legend, must be one of \code{"topright"}, \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"right"}, or \code{"center"}. These keywords place the legend on the inside of the plot frame at the location specified. Partial argument matching is used. If \code{legend=NA} (default), the legend is not drawn.
			
			When \code{plot.type="bars"}, \code{legend} is interpreted instead as a logical value and must be either \code{TRUE} (default) or \code{FALSE}.	If \code{TRUE}, labels are drawn at the bottom of the plot, outside of the actual plotting frame.	
			
			Note that \code{legend} parameter is ignored when \code{plot.type="wilcox"} and \code{panel.lower="difference"}.	
		}
		\item{\code{legend.labels}}{
			A character or \code{\link{expression}} vector specifying the text to appear in the legend, when relevant. Length of \code{legend.labels} must match the numer of \code{DVH} and \code{DVH.list} objects, otherwise legend text will default to the form: \code{"Group 1", "Group 2", ...}
		}
		\item{\code{line.transparency}}{
			Factor modifying the transparency of line drawings (value must be between \code{[0,1]}), with \code{fill.transparency=0} specifying complete transparency and \code{1} (default) specifying complete opacity. Note that this parameter can represent a list of different transparency factors, with length corresponding to the number of \code{DVH} or \code{DVH.list} objects  (if \code{plot.type="individual"}, number of transparency factors should correspond to number of individual \code{DVH} objects).
		}
		\item{\code{lty}}{
			The line type. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings \code{"blank"}, \code{"solid"}, \code{"dashed"}, \code{"dotted"}, \code{"dotdash"}, \code{"longdash"}, or \code{"twodash"}, where \code{"blank"} uses 'invisible lines' (i.e., does not draw them).

			Alternatively, a string of up to 8 characters (from c(1:9, "A":"F")) may be given, giving the length of line segments which are alternatively drawn and skipped (see \code{\link{par}} for more details). Note that this paramater can also represent a list of different line types, with length corresponding to the number of \code{DVH} or \code{DVH.list} objects (if \code{plot.type="individual"}, number of line types should correspond to number of individual \code{DVH} objects).
		}
		\item{\code{lwd}}{
			The line width, a \emph{positive} number, defaulting to \code{1}. The interpretation is device-specific, and some devices do not implement line widths less than one (see the help on the device for details of the interpretation). Note that this paramater can also represent a list of different line widths, with length of the list corresponding to the number of \code{DVH} or \code{DVH.list} objects (if \code{plot.type="individual"}, number of line widths should correspond to number of individual \code{DVH} objects).
		}
		\item{\code{main}}{
			An overall title for the plot (see \code{\link{title}}, default is \code{""}).
		}
		\item{\code{multiplier}}{
			A single numerical value (default is \code{1}) specifying the degree to which group display width should be shrunken or magnified. This parameter is only processed when \code{plot.type="grouped"} or \code{"ttest"} and \code{width="mad"}, \code{"sd"}, or \code{"var"}.
		}
		\item{\code{new}}{
			Logical, defaulting to \code{TRUE}. If set to \code{FALSE}, the next high-level plotting command (actually \code{\link{plot.new}}) should \emph{not clean} the frame before drawing \emph{as if it were on a \strong{new} device}. It is an error (ignored with a warning) to try to use \code{new=FALSE} on a device that does not currently contain a high-level plot.
		}
		\item{\code{panel.lower}}{
			Value specifying the lower panel plot content when \code{plot.type="wilcox"} (\code{panel.lower} must be one of \code{"grouped"} or \code{"wilcox"} [default]). Partial argument matching is used.			
			When \code{panel.lower="wilcox"}, the median groupwise differences will be displayed with associated confidence intervals obtained from \code{wilcox.test()}. When \code{panel.lower="grouped"}, the individual groups are displayed with surrounding user-specified intervals (see \code{width}, \code{quantile}, and \code{multiplier} parameters).
		}
		\item{\code{quantile}}{
			A numeric vector containing two values in \emph{[0,1]} (default is the interquartile range: \code{(0.25, 0.75)}), specifying the lower and upper probabilities (respectively) to be used when \code{width="quantile"} and \code{plot.type="grouped"} or \code{"ttest"}. See \code{\link{quantile}} for more details.
		}
		\item{\code{width}}{
			A character vector specifying the width of shading to use, when relevant (for a group of DVHs). Value must be one of \code{"range"} (default), \code{"mad"}, \code{"IQR"}, \code{"quantile"}, \code{"var"}, or \code{"sd"}, specifying the absolute range, the mean-absolute-deviation, the interquartile (25-75\%) range, the arbitrary inter-quantile range (see \code{quantile} parameter), the group variance, or standard deviation, respectively. This parameter is ignored unless \code{plot.type="grouped"}. If \code{plot.type="ttest"}, the parameter can be optionally specified (if \code{width=NULL} (default), shading will default to the confidence interval).
		}
		\item{\code{xlim}}{
			A numeric vector of length 2 specifying the minimum and maximum x coordinates for plotting (note that this parameter is ignored when \code{plot.type="bars"}, \code{new=FALSE}, or plotting object(s) of class \code{zDVH}).
		}
		\item{\code{ylim}}{
			A numeric vector of length 2 specifying the minimum and maximum y coordinates for plotting (note that this parameter is ignored when \code{plot.type="bars"}, \code{new=FALSE}, or plotting object(s) of class \code{zDVH}).
		}
	}
}
\section{Additional Parameters}{
	\describe{
		\item{\code{alternative}}{
			When \code{plot.type="ttest"}, \code{"wilcox"}, or \code{"correlation"}, \code{alternative} indicates the alternative hypothesis and must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}. You can specify just the initial letter. \code{"greater"} corresponds to positive association, \code{"less"} to negative association. See \code{\link{cor.test}} for more details.
		}
		\item{\code{dose}}{
			Value specifying dose scale (must be one of \code{"relative"} or \code{"absolute"}).
		}
		\item{\code{dose.units}}{
			Value specifying dose units (must be one of \code{"cGy"} or \code{"Gy"}).
		}
		\item{\code{exact}}{
			When \code{plot.type="wilcox"} or \code{plot.type="correlation"} and \code{method="kendall"} or \code{"spearman"}, \code{exact} specifies a logical indicating whether an exact p-value should be computed.  See \code{\link{cor.test}} or \code{\link{wilcox.test}} for more details.
		}
		\item{\code{method}}{
			When \code{plot.type="correlation"}, \code{method} specifies a character string indicating which correlation coefficient is to be computed. Value must be one of \code{"pearson"} (default), \code{"kendall"}, or \code{"spearman"}.  See \code{\link{cor.test}} for more details. 
		}
		\item{\code{mu}}{
			When \code{plot.type="ttest"} or \code{"wilcox"}, specifies an optional parameter used to form the null hypothesis.  See \code{\link{t.test}} or \code{\link{wilcox.test}} for more details.
		}
		\item{\code{paired}}{
			A logical value (default is \code{FALSE}) specifying whether or not to perform paired groupwise comparisons when \code{plot.type="ttest"} or \code{"wilcox"}.
		}
		\item{\code{type}}{
			Value specifying type of dose-volume histogram to plot (must be one of \code{"cumulative"} or \code{"differential"}).
		}
		\item{\code{volume}}{
			Value specifying volume scale (must be one of \code{"relative"} or \code{"absolute"}).
		}
	}
}
\details{
	See package vignette and examples (below) for more details about each plotting subtype.
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
 	#	plot.type="w", col=c("red", "blue"), lty=2, fill.lty=1, main="Target v. OAR",
 	#	panel.lower="grouped")
 	
 	# plot(x, ..., plot.type="correlation")
	# plot(janedoe, sample(1:100, 10), plot.type="correlation", method="kendall")
	
}
\seealso{
\code{\link{plot}}
}
\keyword{methods}
\keyword{aplot}
