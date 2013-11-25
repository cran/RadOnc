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
	data <- readLines(fid)
	close(fid)
	
    ## IDENTIFY STRUCTURES
    struct.start <- grep("^Structure: ", data, perl=TRUE)
    struct.end <- struct.start + diff(c(struct.start, length(data))) - 1
    structures <- mapply(function(start, end) data[start:end], struct.start, struct.end)
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
	if (verbose) {
		cat("[exported on ", date, "]\n", sep="")
		cat("  Patient: ", patient, " (", ID, ")\n", sep="")
		cat("  Plan: ", plan, "\n  Dose: ", dose.rx, dose.units, "\n", sep="")		
	}

	if (length(structures) < 1) {
		warning("No available structures to import")
		return()
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
			return(new("DVH", dose.min=dose.min, dose.max=dose.max, dose.mean=dose.mean, dose.mode=dose.mode, dose.median=dose.median, dose.STD=dose.STD, equiv.sphere=equiv.sphere, conf.index=conf.ind, gradient=gradient, dose.rx=dose.rx, structure.name=name, structure.volume=volume, doses=data.dose, volumes=data, type="cumulative", dose.type=dose.type, dose.units=dose.units, volume.type=volume.type))	
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
	return(NA)
}

read.DVH.DICOM <- function(path, verbose=TRUE) {
	warning("DICOM-RT format not currently supported")
	return(NA)
}
