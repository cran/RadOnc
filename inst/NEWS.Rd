% -*- coding: utf-8 -*-
\name{NEWS}
\title{ RadOnc News }
\encoding{UTF-8}

\section{CHANGES IN VERSION 1.0.3}{
  \subsection{NEW FEATURES}{
    \itemize{
    	\item New methods for \code{sum()} implemented to handle DVH data (assumes *non-overlapping* structures); it can be used to calculate total dose from multiple independent structures (e.g. paired/bilateral lungs or kidneys).
    	\item New function \code{gEUD()} to calculate generalized equivalent uniform dose values from \code{DVH} and \code{DVH.list} objects.
    	\item New function \code{LQE()} to perform linear quadratic extrapolation for iso-effective dose conversion among \code{DVH} and \code{DVH.list} objects with different fractionation.
    	\item Extended \code{read.DVH()} function to import CadPlan files (for achival/historical data).
    	\item Extended \code{read.DVH()} function to import multiple files simultaneously, if desired.
    	\item Expanded functionality of DVH parameter access; now a user may specify output type in parentheses (e.g. \code{DVH["V20Gy(cc)"]}).
   }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Updated \code{plot(..., plot.type="grouped")} to ensure that \code{multiplier} contains a positive value (previously inappropriate negative value would result in a plotting error).
      \item Fixed issue with \code{mean()}, \code{median()}, \code{mad()}, \code{sd()}, and \code{var()} calculations on DVH lists where differing dose units were not handled appropriately.
      \item Added method to \code{mean()} for proper handling of \code{DVH} objects (code was previously omitted, leading to an error on call to the function).
      \item Updated NAMESPACE file to export method definitions for \code{dim()}.
      \item Fixed an internal function which was causing errors upon processing of \code{DVH.list} objects containing one or more empty DVHs.
      \item Fixed fatal error in \code{read.DVH()} caused by attempted read of an improperly-formatted file.
      \item Fixed error in \code{read.DICOM.RT()} affecting structure set import when DICOM header "FrameOfReferenceUID" not specified.
    }
  }

  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{
      \item The \code{center} parameter in the DVH-specific \code{plot(..., plot.type="grouped")} function is no longer supported. Instead, the parameter is now implicitly determined from the \code{width} parameter (e.g. \code{center="mean"} when \code{width="sd"}).
    }
  }

  \subsection{OTHER UPDATES}{
    \itemize{
      \item Added version requirement for \code{oro.dicom} package (must be >= 0.4.1) in order to appropriately import high-resolution structure sets.
      \item DVH list method for \code{wilcox.test()} now updated to include \code{estimate} parameter in returned values.
      \item Confidence interval display updated in DVH list methods \code{plot(..., plot.type="ttest")} and \code{plot(..., plot.type="wilcox")}.
    }
   }
  
  \subsection{ACKNOWLEDGEMENTS}{
    \itemize{
      \item Thanks to Daniel Wollschlaeger for his suggestion to include CadPlan support as well as his support in implementation.
      \item Thanks to Brandon Whitcher for his phenomenal and generous support regarding the \code{oro.dicom} pacakge.
   }
  }
}

\section{CHANGES IN VERSION 1.0.2}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Expanded class \code{DVH}, \code{read.DVH()}, and other supporting functions to allow for dose specification (and analysis) in either units of cGy or Gy.
      \item Added ability to calculate integral dose for \code{DVH} objects, with usage and implementation described further in an updated package vignette.
   }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed issue where structure import from DICOM-RT files failed (due to variable ordering of "ReferencedROINumber" DICOM header).
      \item Fixed issue where attempted structure import from an empty DICOM-RT file caused fatal error.
      \item Fixed issue where attempted structure import from an empty structure set file caused fatal error.
      \item Fixed issue where one or more empty structures in structure list input caused \code{compareStructures()} to fail.
      \item Fixed issue within \code{compareStructures()} where fault in internal color specification caused fatal error.
      \item Fixed issue when plotting using \code{compareStructures()} and \code{method="grid"} where improper colors would be displayed for partially overlapping structures.
      \item Fixed issue where axial comparisons were plotted in improper order when using \code{compareStructures()} and \code{method="grid"}.
    }
  }
    
  \subsection{OTHER UPDATES}{
    \itemize{
      \item Added implementation of \code{dim} function to assess size of \code{structure3D} objects.
      \item Updated implementation of \code{read.DVH} function to allow for import of relative dose parameters (e.g. max dose as percent of prescription dose).
      \item \code{new("DVH")} now creates a structure with volume equivalent to zero (instead of prior NULL value default).
      \item \code{johndoe} and \code{janedoe} datasets were updated to conform to newest \code{DVH} specifications.
    }
   }

  \subsection{ACKNOWLEDGEMENTS}{
    \itemize{
      \item Thanks to Daniel Wollschlaeger for his critical feedback and help in updating this version of the package.
   }
  }
}

\section{CHANGES IN VERSION 1.0.1}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New function \code{read.DICOM.RT()} to import 3-dimensional structural information from one or more DICOM-RT files.
      \item New class \code{structure3D} to store 3-dimensional information encoding a structure.
      \item New functions to support interaction with \code{structure3D} class.
      \item New class \code{structure.list} to store a list of \code{structure3D} objects.
      \item New functions to support interaction with \code{structure.list} class.
      \item New function \code{compareStructures()} to assess similarities and differences among two or more\code{structure3D} objects within a \code{structure.list}.
   }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item For \code{plot(..., method="grouped", width="IQR", center="mean")}, shading offset was previously incorrect, resulting in improper plotting.
      \item For \code{plot(..., method="grouped", width="quantile", center="mean")}, shading offset was previously incorrect, resulting in improper plotting.
      \item Fixed error causing overwrite of slot \code{structure.volume}, preventing proper creation of \code{DVH} objects.
      \item Fixed \code{names<-} function for \code{DVH.list} class which had erroneously prevented lookup of newly assigned names
    }
  }
  
  \subsection{OTHER UPDATES}{
    \itemize{
      \item RadOnc.RData file updated to include example 3-dimensional structural data (\code{cord}, \code{mandible}, \code{teeth}).
      \item Implemented coercion from class DVH to DVH.list, \code{as(..., Class="DVH.list")}. 
   }
  }  
}


\section{CHANGES IN VERSION 1.0.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item This is the first release of the \pkg{RadOnc} package.
   }
  }

  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{
      \item This is the first release of the \pkg{RadOnc} package.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item This is the first release of the \pkg{RadOnc} package.
    }
  }
}
