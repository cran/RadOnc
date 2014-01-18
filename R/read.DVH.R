read.DVH <- function (file, type=c(NA, "aria10", "aria11", "aria8", "dicom", "cadplan"), verbose=TRUE) {
	type <- match.arg(tolower(type), choices=c(NA, "aria10", "aria11", "aria8", "dicom", "cadplan"), several.ok=TRUE)
	if (length(file) < 1) {
		warning("argument 'file' is missing, with no default")
		return()
	}
	
	read.DVH.file <- function (file, type, verbose=TRUE) {
		switch(type, 
			aria10 = return(read.DVH.Aria10(file=file, verbose=verbose)),
			aria8 = return(read.DVH.Aria8(file=file, verbose=verbose)),
			aria11 = return(read.DVH.Aria11(file=file, verbose=verbose)),
			dicom = return(read.DVH.DICOM(path=file, verbose=verbose)),
			cadplan = return(read.DVH.CadPlan(file=file, verbose=verbose)),
			warning("DVH file format not specified for file '", file, "'")
		)
		return()
	}		

	
	if (length(file) == 1) {
		if (length(type) > 1) {
			warning("length of 'file' and 'type' do not match")
			type <- type[1]
		}
		return(read.DVH.file(file, type, verbose))
	}
	if (length(type) < length(file)) {
		if (length(type) > 1) {
			warning("length of 'file' and 'type' do not match")
		}
		type <- rep(type[1], length(file))
	}
	else if (length(type) > length(file)) {
		warning("length of 'file' and 'type' do not match")
		type <- type[1:length(file)]	
	}
	return(mapply(function(x, y, z) {list(read.DVH.file(x, y, z))}, file, type, verbose))
}


read.DVH.Aria10 <- function (file, verbose=TRUE) {
	if (!(fid <- file(file, open="r"))) {
		warning(paste("Could not open file '", file, "'", sep=""))		
		return()
	}
	if (verbose) {
		cat("Reading DVH file ('", file, "')... ", sep="")
	}
	data <- readLines(fid)
	close(fid)
	
    ## IDENTIFY STRUCTURES
    struct.start <- grep("^Structure: ", data, perl=TRUE)
    struct.end <- struct.start + diff(c(struct.start, length(data)+1)) - 1
    if (length(struct.start) < 1) {
		warning(paste("File '", file, "' contained no recognizable DVH structure(s)", sep=""))
		if (verbose) {
			cat("ERROR\n")
		}
		return()
    }
    else if (length(struct.start) == 1) {
   	 	structures <- list(data[struct.start:struct.end])
    }
    else {
	    structures <- mapply(function(start, end) list(data[start:end]), struct.start, struct.end)
	}
	# EXTRACT HEADER INFO
    header <- data[1:(struct.start[1]-1)]
    patient <- sub("^.*: (.*?)\\s*([(].*[)]|$).*", "\\1", header[grep("^Patient Name.*: ", header, ignore.case=TRUE, perl=TRUE)], perl=TRUE)
    ID <- sub("^.*: (.+$)", "\\1", header[grep("^Patient ID.*: ", header, ignore.case=TRUE, perl=TRUE)])
    plan.sum <- (grepl("Summed", header[grep("^Comment.*: ", header, ignore.case=TRUE, perl=TRUE)]))
    date <- sub("^.*: (.+$)", "\\1", header[grep("^Date.*: ", header, ignore.case=TRUE, perl=TRUE)])
	plan <- sub("^.*: (.+$)", "\\1", header[grep("^Plan.*: ", header, ignore.case=TRUE, perl=TRUE)])
	DVH.type <- header[grep("^Type.*: ", header, ignore.case=TRUE, perl=TRUE)]
	if (grepl("Cumulative", DVH.type, ignore.case=TRUE, perl=TRUE)) {
		DVH.type <- "cumulative"
	}
	else {
		DVH.type <- "differential"
	}
	# EXTRACT PRESCRIPTION DOSE AND DOSE UNITS
	dose.rx <- header[grep("^Prescribed dose.*: ", header, ignore.case=TRUE, perl=TRUE)]
    dose.units <- toupper(sub("^.*[[](.*)[]].*", "\\1", dose.rx, perl=TRUE))

	if (dose.units == "GY") {
		dose.units <- "Gy"
	}
	else if (dose.units == "CGY") {
		dose.units <- "cGy"
	}
	dose.rx <- suppressWarnings(as.numeric(sub("^Prescribed dose.*: ", "", dose.rx, ignore.case=TRUE, perl=TRUE)))
    rx.isodose <- as.numeric(sub(".*: ", "", header[grep("^[%] for dose.*: ", header, ignore.case=TRUE, perl=TRUE)]))
    
	if (verbose) {
		cat("[exported on ", date, "]\n", sep="")
		cat("  Patient: ", patient, " (", ID, ")\n", sep="")
		cat("  Plan: ", plan, "\n  Dose: ", if (is.na(dose.rx)) {"NOT SPECIFIED"} else {paste(dose.rx, dose.units, " (at ", rx.isodose, "% isodose line)", sep="")}, "\n", sep="")		
	}

	# EXTRACT DVH DATA FOR EACH STRUCTURE
	DVH.list <- lapply(structures,
		function (data) {
		    name <- sub("^.*: (.+$)", "\\1", data[grep("^Structure.*: ", data, ignore.case=TRUE, perl=TRUE)])
		    if (length(name) < 1) {
				warning("Invalid DVH file format, could not extract structure")
				return(new("DVH"))
		    }
		    volume <- suppressWarnings(as.numeric(sub("^.*: (.+$)", "\\1", data[grep("^Volume.*: ", data, ignore.case=TRUE, perl=TRUE)], perl=TRUE)))
		    if (length(volume) < 1) {
				warning("Invalid DVH file format, could not extract structure volume information")
				return(new("DVH"))
		    }
			getDose <- function(dose) {
				if (grepl("[[]%[]]", dose)) {
					dose <- suppressWarnings(as.numeric(sub(".*: ", "", dose)))*dose.rx/100		
				}
				else {
					dose <- suppressWarnings(as.numeric(sub(".*: ", "", dose)))
				}
				return(dose)				
			}

			dose.min <- max(0, getDose(data[grep("^Min Dose.*: ", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.max <- max(0, getDose(data[grep("^Max Dose.*: ", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.mean <- max(0, getDose(data[grep("^Mean Dose.*: ", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			if (verbose) {
				cat("  ..Importing structure: ", name, "  [volume: ", volume, "cc, dose: ", dose.min, " - ", dose.max, dose.units, "]\n", sep="")
			}

			dose.mode <- max(0, getDose(data[grep("^Modal Dose.*: ", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.median <- max(0, getDose(data[grep("^Median Dose.*: ", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.STD <- max(0, getDose(data[grep("^STD.*: ", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)

		    equiv.sphere <- suppressWarnings(as.numeric(sub("^.*: (.+$)", "\\1", data[grep("^Equiv. Sphere Diam.*: ", data, ignore.case=TRUE, perl=TRUE)])))
		    conf.ind <- suppressWarnings(as.numeric(sub("^.*: (.+$)", "\\1", data[grep("^Conformity Index.*: ", data, ignore.case=TRUE, perl=TRUE)])))
		    gradient <- suppressWarnings(as.numeric(sub("^.*: (.+$)", "\\1", data[grep("^Gradient Measure.*: ", data, ignore.case=TRUE, perl=TRUE)])))

        	header <- grep("Dose [[](%|Gy|cGy)[]].*Volume", data, ignore.case=TRUE, perl=TRUE)
			if (grepl("^\\s*Dose [[](cGy|Gy)[]].*Volume", data[header], ignore.case=TRUE, perl=TRUE)) {
				dose.type <- "absolute"
			}
			else {
				dose.type <- "relative"
			}
			if (grepl(".*Volume [[][%][]]", data[header], ignore.case=TRUE, perl=TRUE)) {
				volume.type <- "relative"
			}
			else {
				volume.type <- "absolute"
			}
			dvh <- read.table(textConnection(data[(header+1):length(data)]), header=FALSE, stringsAsFactors=FALSE)
			data.dose <- dvh[, 1]
			if (plan.sum) {
				data <- dvh[, 2]
			}
			else {
				data <- dvh[, 3]
			}
			if (DVH.type == "differential") {
				data <- data * diff(c(-data.dose[1], data.dose))
				temp.doses <- data.dose - diff(c(-data.dose[1], data.dose))/2
				data.dose <- c(temp.doses, (2*data.dose - temp.doses)[length(temp.doses)])
				data <- diffinv(-data, xi=sum(data))
			}
			return(new("DVH", dose.min=dose.min, dose.max=dose.max, dose.mean=dose.mean, dose.mode=dose.mode, dose.median=dose.median, dose.STD=dose.STD, equiv.sphere=equiv.sphere, conf.index=conf.ind, gradient=gradient, dose.rx=dose.rx, rx.isodose=rx.isodose, structure.name=name, structure.volume=volume, doses=data.dose, volumes=data, type="cumulative", dose.type=dose.type, dose.units=dose.units, volume.type=volume.type))	
		}
	)
	
	# RETURN DVH LIST
	names(DVH.list) <- unlist(lapply(DVH.list, names))
	return(new("DVH.list", DVH.list))
}

read.DVH.Aria11 <- function (file, verbose=TRUE) {
	return(read.DVH.Aria10(file, verbose))
}


read.DVH.Aria8 <- function (file, verbose=TRUE) {
	warning("Aria 8 format not currently supported")
	return()
}

read.DVH.DICOM <- function(path, verbose=TRUE) {
	warning("DICOM-RT format not currently supported")
	return()
	
	if (length(list.files(path)) == 0 && file.exists(path)) {
    	filenames <- path
	}
	else {
      filenames <- list.files(path, full.names=TRUE, recursive=TRUE)
    }
	if (verbose) {
		cat("Reading ", length(filenames), " DICOM files from path: '", path, "' ... ", sep="")
	}
	DICOMs <- readDICOM(path, verbose=FALSE, skipSequence=TRUE)
	if (verbose) {
		cat("FINISHED\n")
	}
	modalities <- as.character(unlist(lapply(DICOMs$hdr, function(x) {x[which(x[,"name"]=="Modality"), "value"]})))
	for (i in as.numeric(which(modalities == "RTDOSE"))) {
		if (verbose) {
			cat("Reading DVH data from file: '", filenames[i], "' ... ", sep="")
		}
		DICOMs$hdr[[i]] <- DICOM.i <- readDICOMFile(filenames[i], skipSequence=FALSE)$hdr
		dvh.start <- which(DICOM.i[,"name"] == "DVHType")
     	dvh.end <- dvh.start + diff(c(dvh.start, dim(DICOM.i)[1]+1)) - 1
	    if (length(dvh.start) < 1) {
			dvhs <- list()
	    }
	    else if (length(dvh.start) == 1) {
	    	dvhs <- list(DICOM.i[dvh.start:dvh.end])
	    }
	    else {
		    dvhs <- mapply(function(start, end) list(DICOM.i[start:end]), dvh.start, dvh.end)
		}

#     group element                           name code   length
#1     0002    0000                    GroupLength   UL        4
#2     0002    0001     FileMetaInformationVersion   OB        2
#3     0002    0002        MediaStorageSOPClassUID   UI       30
#4     0002    0003     MediaStorageSOPInstanceUID   UI       50
#5     0002    0010              TransferSyntaxUID   UI       18
#6     0002    0012         ImplementationClassUID   UI       20
#7     0008    0005           SpecificCharacterSet   CS       10
#8     0008    0012           InstanceCreationDate   DA        8
#9     0008    0013           InstanceCreationTime   TM       14
#10    0008    0016                    SOPClassUID   UI       30
#11    0008    0018                 SOPInstanceUID   UI       50
#12    0008    0020                      StudyDate   DA        8
#13    0008    0030                      StudyTime   TM       14
#14    0008    0050                AccessionNumber   SH        0
#15    0008    0060                       Modality   CS        6
#16    0008    0070                   Manufacturer   LO       22
#17    0008    0090        ReferringPhysiciansName   PN        4
#18    0008    1010                    StationName   SH       12
#19    0008    1030               StudyDescription   LO       24
#20    0008    103E              SeriesDescription   LO       14
#21    0008    1048             PhysiciansOfRecord   PN       30
#22    0008    1090         ManufacturersModelName   LO       12
#23    0010    0010                   PatientsName   PN       12
#24    0010    0020                      PatientID   LO       10
#25    0010    0030              PatientsBirthDate   DA        8
#26    0010    0032              PatientsBirthTime   TM        6
#27    0010    0040                    PatientsSex   CS        2
#28    0010    1000                OtherPatientIDs   LO        0
#29    0010    2160                    EthnicGroup   SH        6
#30    0018    0050                 SliceThickness   DS        0
#31    0018    1000             DeviceSerialNumber   LO       10
#32    0018    1020               SoftwareVersions   LO        8
#33    0020    000D               StudyInstanceUID   UI       56
#34    0020    000E              SeriesInstanceUID   UI       50
#35    0020    0010                        StudyID   SH        2
#36    0020    0011                   SeriesNumber   IS        2
#37    0020    0013                 InstanceNumber   IS        0
#38    0020    0032           ImagePositionPatient   DS       32
#39    0020    0037        ImageOrientationPatient   DS       20
#40    0020    0052            FrameOfReferenceUID   UI       56
#41    0020    1040     PositionReferenceIndicator   LO        0
#42    0028    0002                SamplesperPixel   US        2
#43    0028    0004      PhotometricInterpretation   CS       12
#44    0028    0008                 NumberOfFrames   IS        4
#45    0028    0009          FrameIncrementPointer   AT        4
#46    0028    0010                           Rows   US        2
#47    0028    0011                        Columns   US        2
#48    0028    0030                   PixelSpacing   DS       18
#49    0028    0100                  BitsAllocated   US        2
#50    0028    0101                     BitsStored   US        2
#51    0028    0102                        HighBit   US        2
#52    0028    0103            PixelRepresentation   US        2
#53    3004    0002                      DoseUnits   CS        2
#54    3004    0004                       DoseType   CS        8
#55    3004    000A              DoseSummationType   CS        4
#56    3004    000C          GridFrameOffsetVector   DS      440
#57    3004    000E                DoseGridScaling   DS       12
#58    3004    0014  TissueHeterogeneityCorrection   CS       18
#59    3004    0050                    DVHSequence   SQ  4658916
#60    FFFE    E000                           Item   UN    25704
#61    3004    0001                        DVHType   CS       10
#62    3004    0002                      DoseUnits   CS        2
#63    3004    0004                       DoseType   CS        8
#64    3004    0052                 DVHDoseScaling   DS        2
#65    3004    0054                 DVHVolumeUnits   CS        4
#66    3004    0056                DVHNumberOfBins   IS        4
#67    3004    0058                        DVHData   DS    25504
#68    3004    0060       DVHReferencedROISequence   SQ       34
#69    FFFE    E000                           Item   UN       26
#70    3004    0062         DVHROIContributionType   CS        8
#71    3006    0084            ReferencedROINumber   IS        2
#72    3004    0070                 DVHMinimumDose   DS       16
#73    3004    0072                 DVHMaximumDose   DS       16
#74    3004    0074                    DVHMeanDose   DS       16
#75    FFFE    E000                           Item   UN    98978
#76    3004    0001                        DVHType   CS       10
#77    3004    0002                      DoseUnits   CS        2
#78    3004    0004                       DoseType   CS        8
#79    3004    0052                 DVHDoseScaling   DS        2
#80    3004    0054                 DVHVolumeUnits   CS        4
#81    3004    0056                DVHNumberOfBins   IS        4
#82    3004    0058                        DVHData   DS    98778
#83    3004    0060       DVHReferencedROISequence   SQ       34
#84    FFFE    E000                           Item   UN       26
#85    3004    0062         DVHROIContributionType   CS        8
#86    3006    0084            ReferencedROINumber   IS        2
#87    3004    0070                 DVHMinimumDose   DS       16
#88    3004    0072                 DVHMaximumDose   DS       16
#89    3004    0074                    DVHMeanDose   DS       16
#90    FFFE    E000                           Item   UN   129980
#91    3004    0001                        DVHType   CS       10
#92    3004    0002                      DoseUnits   CS        2
#93    3004    0004                       DoseType   CS        8
#94    3004    0052                 DVHDoseScaling   DS        2
#95    3004    0054                 DVHVolumeUnits   CS        4
#96    3004    0056                DVHNumberOfBins   IS        4
#97    3004    0058                        DVHData   DS   129780
#98    3004    0060       DVHReferencedROISequence   SQ       34
#99    FFFE    E000                           Item   UN       26
#100   3004    0062         DVHROIContributionType   CS        8
#101   3006    0084            ReferencedROINumber   IS        2
#102   3004    0070                 DVHMinimumDose   DS       16
#103   3004    0072                 DVHMaximumDose   DS       16
#104   3004    0074                    DVHMeanDose   DS       16
#105   FFFE    E000                           Item   UN    31340
#106   3004    0001                        DVHType   CS       10
#107   3004    0002                      DoseUnits   CS        2
#108   3004    0004                       DoseType   CS        8
#109   3004    0052                 DVHDoseScaling   DS        2
#110   3004    0054                 DVHVolumeUnits   CS        4
#111   3004    0056                DVHNumberOfBins   IS        4
#112   3004    0058                        DVHData   DS    31144
#113   3004    0060       DVHReferencedROISequence   SQ       34
#114   FFFE    E000                           Item   UN       26
#115   3004    0062         DVHROIContributionType   CS        8
#116   3006    0084            ReferencedROINumber   IS        2
#117   3004    0070                 DVHMinimumDose   DS       12
#118   3004    0072                 DVHMaximumDose   DS       16
#119   3004    0074                    DVHMeanDose   DS       16
#120   FFFE    E000                           Item   UN   127398
#121   3004    0001                        DVHType   CS       10
#122   3004    0002                      DoseUnits   CS        2
#123   3004    0004                       DoseType   CS        8
#124   3004    0052                 DVHDoseScaling   DS        2
#125   3004    0054                 DVHVolumeUnits   CS        4
#126   3004    0056                DVHNumberOfBins   IS        4
#127   3004    0058                        DVHData   DS   127198
#128   3004    0060       DVHReferencedROISequence   SQ       34
#129   FFFE    E000                           Item   UN       26
#130   3004    0062         DVHROIContributionType   CS        8
#131   3006    0084            ReferencedROINumber   IS        2
#132   3004    0070                 DVHMinimumDose   DS       16
#133   3004    0072                 DVHMaximumDose   DS       16
#134   3004    0074                    DVHMeanDose   DS       16
#135   FFFE    E000                           Item   UN    33312
#136   3004    0001                        DVHType   CS       10
#137   3004    0002                      DoseUnits   CS        2
#138   3004    0004                       DoseType   CS        8
#139   3004    0052                 DVHDoseScaling   DS        2
#140   3004    0054                 DVHVolumeUnits   CS        4
#141   3004    0056                DVHNumberOfBins   IS        4
#142   3004    0058                        DVHData   DS    33116
#143   3004    0060       DVHReferencedROISequence   SQ       34
#144   FFFE    E000                           Item   UN       26
#145   3004    0062         DVHROIContributionType   CS        8
#146   3006    0084            ReferencedROINumber   IS        2
#147   3004    0070                 DVHMinimumDose   DS       12
#148   3004    0072                 DVHMaximumDose   DS       16
#149   3004    0074                    DVHMeanDose   DS       16
#150   FFFE    E000                           Item   UN   107976
#151   3004    0001                        DVHType   CS       10
#152   3004    0002                      DoseUnits   CS        2
#153   3004    0004                       DoseType   CS        8
#154   3004    0052                 DVHDoseScaling   DS        2
#155   3004    0054                 DVHVolumeUnits   CS        4
#156   3004    0056                DVHNumberOfBins   IS        4
#157   3004    0058                        DVHData   DS   107776
#158   3004    0060       DVHReferencedROISequence   SQ       34
#159   FFFE    E000                           Item   UN       26
#160   3004    0062         DVHROIContributionType   CS        8
#161   3006    0084            ReferencedROINumber   IS        2
#162   3004    0070                 DVHMinimumDose   DS       16
#163   3004    0072                 DVHMaximumDose   DS       16
#164   3004    0074                    DVHMeanDose   DS       16
#165   FFFE    E000                           Item   UN   128072
#166   3004    0001                        DVHType   CS       10
#167   3004    0002                      DoseUnits   CS        2
#168   3004    0004                       DoseType   CS        8
#169   3004    0052                 DVHDoseScaling   DS        2
#170   3004    0054                 DVHVolumeUnits   CS        4
#171   3004    0056                DVHNumberOfBins   IS        4
#172   3004    0058                        DVHData   DS   127876
#173   3004    0060       DVHReferencedROISequence   SQ       34
#174   FFFE    E000                           Item   UN       26
#175   3004    0062         DVHROIContributionType   CS        8
#176   3006    0084            ReferencedROINumber   IS        2
#177   3004    0070                 DVHMinimumDose   DS       12
#178   3004    0072                 DVHMaximumDose   DS       16
#179   3004    0074                    DVHMeanDose   DS       16
#180   FFFE    E000                           Item   UN   107820
#181   3004    0001                        DVHType   CS       10
#182   3004    0002                      DoseUnits   CS        2
#183   3004    0004                       DoseType   CS        8
#184   3004    0052                 DVHDoseScaling   DS        2
#185   3004    0054                 DVHVolumeUnits   CS        4
#186   3004    0056                DVHNumberOfBins   IS        4
#187   3004    0058                        DVHData   DS   107620
#188   3004    0060       DVHReferencedROISequence   SQ       34
#189   FFFE    E000                           Item   UN       26
#190   3004    0062         DVHROIContributionType   CS        8
#191   3006    0084            ReferencedROINumber   IS        2
#192   3004    0070                 DVHMinimumDose   DS       16
#193   3004    0072                 DVHMaximumDose   DS       16
#194   3004    0074                    DVHMeanDose   DS       16
#195   FFFE    E000                           Item   UN   132610
#196   3004    0001                        DVHType   CS       10
#197   3004    0002                      DoseUnits   CS        2
#198   3004    0004                       DoseType   CS        8
#199   3004    0052                 DVHDoseScaling   DS        2
#200   3004    0054                 DVHVolumeUnits   CS        4
#201   3004    0056                DVHNumberOfBins   IS        4
#202   3004    0058                        DVHData   DS   132410


		if (verbose) {
			cat("FINISHED\n")
		}		
	}
}

read.DVH.CadPlan <- function(file, verbose=TRUE) {
	if (!(fid <- file(file, open="r"))) {
		warning(paste("Could not open file '", file, "'", sep=""))		
		return()
	}
	if (verbose) {
		cat("Reading DVH file ('", file, "')... ", sep="")
	}
	data <- readLines(fid)
	close(fid)
	
    ## IDENTIFY STRUCTURES
    struct.start <- grep("^Histogram.*:\\s*", data, perl=TRUE)
    struct.end <- struct.start + diff(c(struct.start, length(data)+1)) - 1
    if (length(struct.start) < 1) {
		warning(paste("File '", file, "' contained no recognizable DVH structure(s)", sep=""))
		if (verbose) {
			cat("ERROR\n")
		}
		return()
    }
    else if (length(struct.start) == 1) {
    	structures <- list(data[struct.start:struct.end])
    }
    else {
	    structures <- mapply(function(start, end) list(data[start:end]), struct.start, struct.end)
	}
	# EXTRACT HEADER INFO
    header <- data[1:(struct.start[1]-1)]
    patient <- sub("^.*: (.*)", "\\1", header[grep("^Patient Name.*:\\s*", header, ignore.case=TRUE, perl=TRUE)], perl=TRUE)
    ID <- sub("^.*: (.+$)", "\\1", header[grep("^Patient ID.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
    date <- sub("^.*: (.+$)", "\\1", header[grep("^Date.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
	plan <- sub("^PLAN\\s*(.+$)", "\\1", header[grep("^PLAN\\s+", header, ignore.case=TRUE, perl=TRUE)])
	DVH.type <- header[grep("Dose Volume Histogram", header, ignore.case=TRUE, perl=TRUE)]
	if (grepl("Cumulative", DVH.type, ignore.case=TRUE, perl=TRUE)) {
		DVH.type <- "cumulative"
	}
	else {
		DVH.type <- "differential"
	}

	if (verbose) {
		cat("[exported on ", date, "]\n", sep="")
		cat("  Patient: ", patient, " (", ID, ")\n", sep="")
		cat("  Plan: ", plan, "\n", sep="")		
	}
	
	# EXTRACT DVH DATA FOR EACH STRUCTURE
	DVH.list <- lapply(structures,
		function (data) {
			# EXTRACT PRESCRIPTION DOSE AND DOSE UNITS
			dose.rx <- data[grep("^Prescr[.] dose.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]
		    dose.units <- toupper(sub("^.*[(](.*)[)].*", "\\1", dose.rx, perl=TRUE))
			if (dose.units == "GY") {
				dose.units <- "Gy"
			}
			else if (dose.units == "CGY") {
				dose.units <- "cGy"
			}
			dose.rx <- suppressWarnings(as.numeric(sub("^Prescr[.] dose.*:\\s*", "", dose.rx, ignore.case=TRUE, perl=TRUE)))
			if (is.na(dose.rx)) {
				warning("Prescription dose not specified")
			}

		    name <- sub("^.*:\\s*(.+$)", "\\1", data[grep("^Histogram.*:\\s*", data, ignore.case=TRUE, perl=TRUE)])
		    if (length(name) < 1) {
				warning("Invalid DVH file format, could not extract structure")
				return(new("DVH"))
		    }
		    volume <- suppressWarnings(as.numeric(sub("^.*:\\s*(.+$)", "\\1", data[grep("^Volume.*:\\s*", data, ignore.case=TRUE, perl=TRUE)], perl=TRUE)))
		    if (length(volume) < 1) {
				warning("Invalid DVH file format, could not extract structure volume information")
				return(new("DVH"))
		    }
			getDose <- function(dose) {
				if (grepl("[(]\\s*[%]\\s*[)]", dose)) {
					dose <- suppressWarnings(as.numeric(sub(".*:\\s*", "", dose)))*dose.rx/100		
				}
				else {
					dose <- suppressWarnings(as.numeric(sub(".*:\\s*", "", dose)))
				}
				return(dose)				
			}

			dose.min <- max(0, getDose(data[grep("^Dose minimum.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.max <- max(0, getDose(data[grep("^Dose maximum.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.mean <- max(0, getDose(data[grep("^Dose mean.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			if (verbose) {
				cat("  ..Importing structure: ", name, "  [volume: ", volume, "cc, dose: ", dose.min, " - ", dose.max, dose.units, "]\n", sep="")
			}

			dose.mode <- max(0, getDose(data[grep("^Dose modal.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.median <- max(0, getDose(data[grep("^Dose median.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.STD <- max(0, getDose(data[grep("^Standard dev.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)

        	header <- grep("Dose\\s*[(]\\s*(%|Gy|cGy)\\s*[)].*Volume", data, ignore.case=TRUE, perl=TRUE)
			if (grepl("^\\s*Dose\\s*[(]\\s*(Gy|cGy)\\s*[)].*Volume", data[header], ignore.case=TRUE, perl=TRUE)) {
				dose.type <- "absolute"
			}
			else {
				dose.type <- "relative"
			}
			if (grepl(".*Volume\\s*[(]\\s*[%]\\s*[)]", data[header], ignore.case=TRUE, perl=TRUE)) {
				volume.type <- "relative"
			}
			else {
				volume.type <- "absolute"
			}
			dvh <- read.table(textConnection(data[(header+1):length(data)]), header=FALSE, stringsAsFactors=FALSE)
			data.dose <- dvh[, 1]
			data <- dvh[, 2]
			if (DVH.type == "differential") {
				data <- data * diff(c(-data.dose[1], data.dose))
				temp.doses <- data.dose - diff(c(-data.dose[1], data.dose))/2
				data.dose <- c(temp.doses, (2*data.dose - temp.doses)[length(temp.doses)])
				data <- diffinv(-data, xi=sum(data))
			}
			return(new("DVH", dose.min=dose.min, dose.max=dose.max, dose.mean=dose.mean, dose.mode=dose.mode, dose.median=dose.median, dose.STD=dose.STD, dose.rx=dose.rx, structure.name=name, structure.volume=volume, doses=data.dose, volumes=data, type="cumulative", dose.type=dose.type, dose.units=dose.units, volume.type=volume.type))	
		}
	)
	
	# RETURN DVH LIST
	names(DVH.list) <- unlist(lapply(DVH.list, names))
	return(new("DVH.list", DVH.list))
}