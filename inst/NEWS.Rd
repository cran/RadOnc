% -*- coding: utf-8 -*-
\name{NEWS}
\title{ RadOnc News }
\encoding{UTF-8}

\section{CHANGES IN VERSION 1.0.1}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New function \code{read.DICOM.RT()} to import 3-dimensional structural information from one or more DICOM-RT files.
      \item New class \code{structure3D} to store 3-dimensional information encoding a structure.
      \item New functions to support interaction with \code{structure3D} class.
      \item New class \code{structure.list} to store a list of \code{structure3D} objects.
      \item New functions to support interaction with \code{structure.list} class.
      \item New function \code{compareContours()} to assess similarities and differences among two or more\code{structure3D} objects within a \code{structure.list}.
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
      \item RadOnc.RData file updated to include example 3-dimensional structural data (\code{cord}, \code{mandible}).
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
