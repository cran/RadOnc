% -*- coding: utf-8 -*-
\name{NEWS}
\title{ RadOnc News }
\encoding{UTF-8}

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
