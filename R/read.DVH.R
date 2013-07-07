read.DVH <- function (file, type=c(NA, "aria10", "aria11", "aria8", "dicom"), verbose=TRUE) {
	type <- tolower(type)
	type <- match.arg(type)
	switch(type, 
		aria10 = return(read.DVH.Aria10(file=file, verbose=verbose)),
		aria8 = return(read.DVH.Aria8(file=file, verbose=verbose)),
		aria11 = return(read.DVH.Aria11(file=file, verbose=verbose)),
		dicom = return(read.DVH.DICOM(path=file, verbose=verbose)),
		warning("DVH file format not specified for file '", file, "'")
	)
	return(NA)
}

read.DVH.Aria10 <- function (file, verbose=TRUE) {
	if (!(fid <- file(file, open="r"))) {
		stop(paste("Could not open file '", file, "'", sep=""))		
	}
	if (verbose) {
		cat("Reading DVH file ('", file, "')... ", sep="")
	}
	DVH.list <- c()
	# EXTRACT HEADER INFO
	patient <- sub("Patient Name.*: ([A-Z, ]*) [(].*[)]", "\\1", readLines(fid, n=1), ignore.case=TRUE)
	ID <- sub("Patient ID.*: ", "", readLines(fid, n=1), ignore.case=TRUE)
	plan.sum <- (grepl("Summed", readLines(fid, n=1), ignore.case=TRUE))
	date <- sub("Date.*: ", "", readLines(fid, n=1), ignore.case=TRUE)
	if (sub("Type.*: ", "", readLines(fid, n=1), ignore.case=TRUE) == "Cumulative Dose Volume Histogram") {
		DVH.type <- "cumulative"
		readLines(fid, n=4)
	}
	else {
		DVH.type <- "differential"
		readLines(fid, n=6)
	}	
	plan <- sub("Plan.*: ", "", readLines(fid, n=1), ignore.case=TRUE)
	dose.rx <- suppressWarnings(as.numeric(sub("Prescribed dose.*: ", "", readLines(fid, n=1), ignore.case=TRUE)))
	if (verbose) {
		cat("[exported on ", date, "]\n", sep="")
		cat("  Patient: ", patient, " (", ID, ")\n", sep="")
		cat("  Plan: ", plan, "\n  Dose: ", dose.rx, "cGy\n", sep="")
	}
	readLines(fid, n=2)
	while (length(structure <- readLines(fid, n=1)) > 0) {
		if (!grepl("Structure.*: ", structure, ignore.case=TRUE)) {
			if (structure == "") { break }
			stop("Invalid DVH file format, could not extract structure")
		}
		structure <- sub(".*: ", "", structure)
		readLines(fid, n=3)
		volume <- readLines(fid, n=1)
		sub("cm\xb3", "cc", volume)
		if (!grepl("Volume.*: ", volume, ignore.case=TRUE, perl=TRUE)) {
			stop("Invalid DVH file format, could not extract structure volume information")
		}
		volume <- suppressWarnings(as.numeric(sub(".*: ", "", volume)))
		readLines(fid, n=2)
		dose.min <- readLines(fid, n=1)
		if (!grepl("Min Dose.*: ", dose.min, ignore.case=TRUE)) {
			stop("Invalid DVH file format, could not extract structure minimum dose")
		}
		dose.min <- suppressWarnings(as.numeric(sub(".*: ", "", dose.min)))
		dose.max <- readLines(fid, n=1)
		if (!grepl("Max Dose.*: ", dose.max, ignore.case=TRUE)) {
			stop("Invalid DVH file format, could not extract structure maximum dose")
		}
		dose.max <- suppressWarnings(as.numeric(sub(".*: ", "", dose.max)))
		dose.mean <- readLines(fid, n=1)
		if (!grepl("Mean Dose.*: ", dose.mean, ignore.case=TRUE)) {
			stop("Invalid DVH file format, could not extract structure mean dose")
		}
		dose.mean <- suppressWarnings(as.numeric(sub(".*: ", "", dose.mean)))
		if (verbose) {
			cat("  ..Importing structure: ", structure, "  [volume: ", volume, "cc, dose: ", dose.min, " - ", dose.max, "cGy]\n", sep="")
		}
		dose.mode <- suppressWarnings(as.numeric(sub(".*: ", "", readLines(fid, n=1))))
		dose.median <- suppressWarnings(as.numeric(sub(".*: ", "", readLines(fid, n=1))))
		dose.STD <- suppressWarnings(as.numeric(sub(".*: ", "", readLines(fid, n=1))))
		equiv.sphere <- suppressWarnings(as.numeric(sub(".*: ", "", readLines(fid, n=1))))
		conf.ind <- suppressWarnings(as.numeric(sub(".*: ", "", readLines(fid, n=1))))
		gradient <- suppressWarnings(as.numeric(sub(".*: ", "", readLines(fid, n=1))))
		readLines(fid, n=1)
		header <- readLines(fid, n=1)
		if (!grepl("Dose.*Volume", header, ignore.case=TRUE, perl=TRUE)) {
			stop("Invalid DVH file format, could not isolate dose and volumetric data")
		}
		if (grepl("^[ ]*Dose [[]cGy[]].*Volume", header, ignore.case=TRUE, perl=TRUE)) {
			dose.type <- "absolute"
		}
		else {
			dose.type <- "relative"
		}
		if (grepl(".*Volume [[][%][]]", header, ignore.case=TRUE, perl=TRUE)) {
			volume.type <- "relative"
		}
		else {
			volume.type <- "absolute"
		}
		plan.sum <- (!grepl("Dose.*Dose.*Volume", header, ignore.case=TRUE, perl=TRUE))
		# READ DVH INFORMATION

		data <- c()
		data.dose <- c()

		while (length(line <- readLines(fid, n=1)) > 0) {
			if (line == "") { break }
			data.dose <- c(data.dose, as.numeric(sub("[ ]*([0-9.]*).*", "\\1", line)))
			if (plan.sum) {
				data <- c(data, as.numeric(sub("[ ]*[0-9.]*[ ]*([-0-9.e]*).*", "\\1", line)))
			}
			else {
				data <- c(data, as.numeric(sub("[ ]*[0-9.]*[ ]*[0-9.]*[ ]*([-0-9.e]*).*", "\\1", line)))
			}
		}		
		if (DVH.type == "differential") {
			data <- data * diff(c(-data.dose[1], data.dose))
			temp.doses <- data.dose - diff(c(-data.dose[1], data.dose))/2
			data.dose <- c(temp.doses, (2*data.dose - temp.doses)[length(temp.doses)])
			data <- diffinv(-data, xi=sum(data))
		}
		DVH.list <- c(DVH.list, new("DVH", dose.min=dose.min, dose.max=dose.max, dose.mean=dose.mean, dose.mode=dose.mode, dose.median=dose.median, dose.STD=dose.STD, equiv.sphere=equiv.sphere, conf.index=conf.ind, gradient=gradient, dose.rx=dose.rx, structure.name=structure, structure.volume=volume, doses=data.dose, volumes=data, type="cumulative", dose.type=dose.type, volume.type=volume.type))
	}
	close (fid)
	names(DVH.list) <- unlist(lapply(DVH.list, names))
	return(new("DVH.list", DVH.list))
}

read.DVH.Aria11 <- function (file, verbose=TRUE) {
	return(read.DVH.Aria10(file, verbose))
}


read.DVH.Aria8 <- function (file, verbose=TRUE) {
	warning("Aria 8 format not currently supported")
	return(NA)
}

read.DVH.DICOM <- function(path, verbose=TRUE) {
	warning("DICOM-RT format not currently supported")
	return(NA)
}
