% -*- coding: utf-8 -*-
\name{NEWS}
\title{ RadOnc News }
\encoding{UTF-8}

\section{CHANGES IN VERSION 1.1.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item Extended \code{read.DVH()} function to import Elekta Monaco files.
      \item Added ability to display gridlines when plotting DVH graphs.
      \item Added ability to extract median dose as a dosimetric parameter.
      \item Added ability to check for arbitrary dose constraints (e.g. \code{V20Gy < 15(\%)}).
      \item New function \code{subset()} to extract a desired subset of a \code{DVH.list} object based upon patient identifiers, structure name, and/or dosimetric constraints.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
    	\item \code{read.DVH()} maintained open file connections, problematic for simultaneous import of hundreds of DVH files.
    	\item \code{oro.dicom} package failed to recognize DICOM files formatted without 128 byte headers (patch applied to \code{oro.dicom} package version 0.5.0). This affected import of DICOM-RT data from TomoTherapy systems and other instances where file format departed from strict DICOM standards.
    	\item \code{read.DVH()} failed to calculate structure volumes appropriately when structure volumes were not explicitly specified.
     }
  }

  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{
      \item Function \code{getStructureList()} has been removed and replaced instead with a wider array of functionality using \code{subset()}.
    }
  }

  \subsection{OTHER UPDATES}{
    \itemize{
      \item \code{read.DICOM.RT()} behavior updated, such that CT data is now scaled (using DICOM parameters \code{"RescaleSlope"} and \code{"RescaleIntercept"}) to represent Hounsfield Unit data.
      \item Internal DVH conversion now generates warning when trying to convert from 'relative' to 'absolute' volume units when structure volume is zero.
      \item DVH class behavior changed to enable relative/absolute dose conversion for stored dose parameters (e.g. \code{dose.max}); changes affect \code{read.DVH()}, \code{mean()}, \code{min()}, \code{max()}, \code{range()}, and other internal DVH specification and conversion functions.
      \item Isodose line specification support added to \code{read.DVH(..., type="cadplan")}.
    }
   }

  \subsection{ACKNOWLEDGEMENTS}{
    \itemize{
      \item Thanks to Daniel Wollschlaeger for identifying and fixing the bug in \code{read.DVH()}.
      \item Thanks to Francis Gibbons for his contributions regarding inclusion of Elekta Monaco-formatted DVH data.
      \item Thanks to Usman Lula for helping to identify an issue with DICOM-RT import from certain TomoTherapy workstations and to Brandon Whitcher for patching the \code{oro.dicom} package accordingly.
    }
  }

 }


\section{CHANGES IN VERSION 1.0.9}{

  \subsection{NEW FEATURES}{
    \itemize{
    	\item New parameter \code{panel.lower} to specify plotting behavior when using \code{plot(..., plot.type="wilcox")}.
      \item Enabled \code{xlim} and \code{ylim} specification when plotting \code{DVH} object(s).
      \item Enabled import of a single collapsed \code{DVH.list} from multiple input files (previously only able to import a nested list of individual \code{DVH.list} objects corresponding to each input file).
      \item Enabled dose range specification on parameter export (e.g. can now calculate \code{"V10-20Gy"} or \code{"V<1500cGy"}) for a given \code{DVH} object.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
    	\item \code{range()} calculation for an empty structure list failed, and calculation for a structure defined by a single point inappropriately returned value of \code{NA}.
    	\item \code{calculate.DVH()} failed when structure contained fewer than three points or was confined to a single plane along the x or y axes.
    	\item \code{read.DICOM.RT()} failed when calculating DVHs from dose grid
    	\item \code{calculate.DVH()} subroutines improperly ignored \code{dose.unit} parameter input
    	\item DVHs containing plan sum data without prescription dose/isodose information were considered as invalid DVH objects.
    	\item DVH plotting using \code{plot(..., plot.type="ttest")} or \code{plot(..., plot.type="wilcox")} improperly displayed relative volume information along an absolute volume y-axis when parameter \code{volume="absolute"} was specified.
    	\item DVH plotting using \code{plot(..., plot.type="wilcox")} displayed inaccurate confidence intervals.
     }
  }

  \subsection{OTHER UPDATES}{
    \itemize{
      \item Updated behavior of \code{wilcox.test(..., paired=TRUE)} to return parameter \code{"estimate"} as the true median of pairwise differences in lieu of the statistical (pseudo)median calculated by wilcoxon signed-rank test.
      \item Enhanced legend when using \code{plot(..., plot.type="grouped")} (now supports filling/shading)
      \item Added \code{demo()} functionality as a package tutorial.
    }
   }

  \subsection{ACKNOWLEDGEMENTS}{
    \itemize{
      \item Thanks to Maura Kirk for identifying the plan sum DVH class validation bug.
      \item Thanks to Eric Ojerholm for identifying the DVH plotting bugs, for suggesting user-modifiable \code{xlim}/\code{ylim} specification as well as dose ranges on volumetric parameter export, and for suggesting DVH plotting enhancements when using \code{plot(..., plot.type="wilcox")}.
      \item Thanks to Abodunde Adegunloye for identifying the error in structure \code{range()} calculation as well as the error in \code{calculate.DVH()} and \code{read.DICOM.RT()}.
    }
  }

 }


\section{CHANGES IN VERSION 1.0.8}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item DVH plotting using \code{plot(..., plot.type="ttest")} has been expanded to incorporate width specification as implemented for \code{plot(..., plot.type="grouped")}.
      \item New function \code{getStructureList()} to extract desired structures from a list of DVH lists (generated by \code{read.DVH()} using multiple input files)
   }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Error in validity check for classes \code{DVH} and \code{zDVH} when \code{rx.isodose} not specified.
      \item Error in \code{read.DICOM.RT()} function preventing proper import of RT plan data when the plan contained more than one dose reference point.
      \item Error in \code{read.DICOM.RT()} function preventing proper import of RT plan data when the plan contained no dose grid.
      \item Error in initialization of new \code{RTdata} objects containing dose grid information (attribute \code{dose.units} improperly specified).
      \item Error in sub-selection of \code{DVH.list} and \code{structure.list} objects where individual names were coercible to numerical values (e.g. "45").
      \item Improper z coordinate specification when using \code{read.DICOM.RT()} on an inversely ordered list of DICOM image slices.
     \item DVH plotting using \code{plot(..., plot.type="ttest")} generated improper shading ranges in some instances. This graphical depiction has been adjusted to generate actual confidence intervals with relative contributions scaled by a group's internal variance.
     \item DVH parameter selection (e.g. \code{DVH["V20Gy"]}) failed when called upon an empty DVH object. 
    }
  }

  \subsection{ACKNOWLEDGEMENTS}{
    \itemize{
      \item Thanks to Primoz Peterlin for identifying the bug in \code{DVH} class validation.
      \item Thanks to Nicola Dinapoli for identifying the bug in \code{read.DICOM.RT()} import.
   }
  }

 }


\section{CHANGES IN VERSION 1.0.7}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item New class \code{zDVH} (along with relevant methods) to store one or more axially-segmented and calculated dose-volume histograms. The class otherwise contains all elements of the \code{DVH} class and is largely interchangeable with other \code{DVH} objects.
      \item New plot functionality for axial dose-volume histogram displayed from one or more \code{zDVH} objects.
      \item New method for \code{compareStructures()} ("DSC") which calculates Dice similarity coefficients between pairs of structures.
   }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Fatal error in \code{convert.DVH(..., dose="relative")} when prescription dose not specified.
    }
  }

  \subsection{OTHER UPDATES}{
    \itemize{
      \item Added attribute \code{dose.units} to dose grid in class \code{RTdata} to enable specification and processing of dose units.
      \item Added \code{stomach} to package data contents in order to provide an example of a \code{zDVH} object.
    }
   }

 }


\section{CHANGES IN VERSION 1.0.6}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item Added \code{calculate.DVH()} function to perform DVH calculation given one or more input structure set(s) and dose grid data. This function now also implicitly extends \code{read.DVH(..., method="dicom")} to enable DVH calculation from DICOM-RT data.
   }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed v.1.0.5 Windows build error due to incorrect file path specification in vignette.
    }
  }

  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{
      \item Prior internal function \code{extract.DVH()} has been removed.
    }
  }

  \subsection{OTHER UPDATES}{
    \itemize{
      \item Added \code{janedoe.RTdata} to package data contents in order to provide an example of a \code{RTdata} object.
      \item New method for \code{as()} function to convert a \code{structure3D} object to a \code{DVH} object.
      \item Package vignette has been updated to incorporate some of the changes in recent releases.
    }
   }

 }

\section{CHANGES IN VERSION 1.0.5}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item Functionality expanded for \code{[]} accessor method to \code{DVH.list} and \code{structure.list} objects to enable pattern matching (e.g. \code{janedoe["KIDNEY$"]} and handling of redundant structure names in \code{DVH.list} (see updated package Vignette for further details). 
      \item New method for \code{as()} function to convert a list of \code{DVH} and/or \code{DVH.list} objects to a single combined \code{DVH.list} object.
      \item New method for \code{as()} function to convert a list of \code{structure3D} and/or \code{structure.list} objects to a single combined \code{structure.list} object.
      \item Expanded \code{read.DICOM.RT()} functionality to include calculation of DVHs from dose grid information. DVH calculations can be enabled by specifying logical parameter \code{DVH=TRUE} (default).
      \item \code{LQE()} now supports fractionation specification in either dose per fraction or total number of fractions (parameter \code{N} was introduced for this purpose).
   }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed structure import using \code{read.DICOM.RT()} for structure sets not mapped to a DVH present within the dose file (previously returned empty structures instead of populating the structures with their respective vertices).
      \item Fixed fatal error in DVH extraction from dose grid data (previously occurred when one or more points were duplicated within a structure).
      \item Input of empty directly previously caused fatal error in \code{read.DICOM.RT()}
      \item \code{t.test()}, \code{wilcox.test()} and methods employing this functionality (e.g. \code{plot(..., plot.type="ttest")}) previously failed if any elements in DVH list were empty.
      \item \code{sum()}, \code{mean()}, \code{median()}, \code{mad()}, \code{var()}, and \code{sd()} previously failed if any elements in DVH list were empty.
      \item Fixed integral dose calculation, which previously returned values that were inadvertently scaled by the dose bin width.
      \item Fixed validity checking of \code{DVH} objects containing a single dose-volume value.
      \item x-axis improperly specified dose units as "cGy" when plotting DVHs in units of "Gy". 
      \item Isoeffective dose calculation using \code{LQE()} was previously inaccurate due to incorrect implementation of the Withers isoeffective dose formula -- this has now been corrected.
      \item Previously improper handling of \code{type} specifier in \code{read.DVH()} where specified type was not provided and list of input files was longer than length 6.
    }
  }

  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{
      \item Removed unused package dependency on \pkg{Rcpp} (future releases may reincorporate this dependency in order to support C++ implementations of various algorithms).
    }
  }

  \subsection{OTHER UPDATES}{
    \itemize{
      \item Changed formatting of printed \code{DVH} objects (structure volume is now displayed to precision of one decimal point and min/max doses displayed to two decimals).
      \item \code{dose.units} parameter introduced for \code{LQE()} function.
      \item \code{DVH} class expanded to include patient identifying information (new parameters \code{patient} and \code{ID})
      \item \code{read.DVH(..., type="dicom")} now links to \code{read.DICOM.RT(..., DVH=TRUE)} and returns a DVH list extracted from the relevant DICOM-RT directory.
    }
   }

  \subsection{ACKNOWLEDGEMENTS}{
    \itemize{
      \item Thanks to Daniel Wollschlaeger for identifying the bug in \code{read.DVH()} type specification and for his suggestion to expand the \code{DVH} class specification to include patient identification.
   }
  }

 }


\section{CHANGES IN VERSION 1.0.4}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item Expanded \code{read.DICOM.RT()} functionality to include import of DVHs from dose file where relevant.
      \item Expanded \code{plot()} functionality to include display of a \code{structure.list} object containing multiple structures.
      \item Function \code{compareStructures(..., method="axial")} can now process structures that contain one or more holes (e.g. a hollow sphere)
      \item New wrapper class \code{RTdata} to store a corresponding CT image, dose grid, and structure set.
      \item New functions to support interaction with \code{RTdata} class.
   }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Updated DVH parameter extraction to return maximal dose in cases where extrapolated dose exceeds actual maximum dose.
      \item DVH methods previously failed without graceful warning or exit in cases where prescription dose was not specified (e.g. plan sum DVH data from Eclipse).
    }
  }

  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{
      \item Method \code{compareStructures(..., method="grid")} has been renamed to \code{compareStructures(..., method="axial")}. Code still supports call using "grid" method, however it now generates a warning to use the preferred "axial" method.
    }
  }

  \subsection{OTHER UPDATES}{
    \itemize{
      \item Added plotting support for \code{structure3D} objects without pre-specified triangulation; previously a call to \code{plot()} would generate an empty plotting window.
      \item Updated \code{gEUD()} calculation to improve performance for large values of tissue-specific parameter 'a'.
      \item Added slot \code{rx.isodose} to class \code{DVH} to enable application of prescription dose to variable isodose lines (relevant in select cases).
    }
   }

  \subsection{ACKNOWLEDGEMENTS}{
    \itemize{
      \item Thanks to Daniel Wollschlaeger for identifying the bug in DVH dose extraction as well as for his suggestion regarding variable prescription isodose lines.
      \item Thanks to Brandon Whitcher for his phenomenal and generous support regarding the \code{oro.dicom} pacakge.
   }
  }

 }

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
      \item Added slot \code{dose.fx} to class \code{DVH} to enable specification of dose fractionation (relevant in select cases).
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
