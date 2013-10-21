
read.DICOM.RT <- function(path, exclude=NULL, recursive=TRUE, verbose=TRUE, limit=NULL, ...) {
	if (length(list.files(path)) == 0 && file.exists(path)) {
    	filenames <- path
	}
	else {
      filenames <- list.files(path, full.names=TRUE, recursive=recursive)
    }
	if (! is.null(exclude)) {
    	filenames <- grep(exclude, filenames, ignore.case=TRUE, value=TRUE, invert=TRUE)
  	}
	if (verbose) {
		cat("Reading ", length(filenames), " DICOM files from path: '", path, "' ... ", sep="")
	}
	DICOMs <- readDICOM(path, verbose=FALSE, exclude=exclude, recursive=recursive, ...)
	if (verbose) {
		cat("FINISHED\nExtracting CT data ... ", sep="")
	}
	modalities <- as.character(unlist(lapply(DICOMs$hdr, function(x) {x[which(x[,"name"]=="Modality"), "value"]})))
	CT <- as.numeric(which(modalities == "CT"))
	frame.ref.CT <- as.character(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "FrameOfReferenceUID"), "value"])
	## ASSUMES CONSTANT VOXEL SIZE AND SLICE THICKNESS FOR ALL DICOM FILES IN CT!!!!
	voxel.size <- as.numeric(unlist(strsplit(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "PixelSpacing"), "value"], " ")))
	image.position <- as.numeric(unlist(strsplit(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "ImagePositionPatient"), "value"], " ")))
	z.slices <- lapply(DICOMs$hdr[CT], function(x) { as.numeric(unlist(strsplit(x[which(x[,"name"] == "ImagePositionPatient"), "value"], " "))[3]) })
	voxel.size <- c(voxel.size, as.numeric(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "SliceThickness"), "value"]))
	CT <- create3D(list(hdr=DICOMs$hdr[CT], img=DICOMs$img[CT]))
	dimnames(CT) <- list((1:dim(CT)[1]-1)*voxel.size[1]+image.position[1], (1:dim(CT)[2]-1)*voxel.size[2]+image.position[2], image.position[3]-(1:dim(CT)[3]-1)*voxel.size[3]) 
	if (verbose) {
		cat("FINISHED\n")
	}
	first <- TRUE
	for (i in as.numeric(which(modalities == "RTSTRUCT"))) {
		if (verbose) {
			cat("Reading structure set from file: '", filenames[i], "' ... ", sep="")
		}
		DICOMs$hdr[[i]] <- DICOM.i <- readDICOMFile(filenames[i], skipSequence=FALSE)$hdr
			
		structures <- DICOM.i[which(DICOM.i[,"name"] %in% c("ROIName", "ROINumber")),]
		N <- dim(structures)[1]/2
		if (frame.ref.CT != as.character(DICOM.i[which(DICOM.i[,"name"] == "FrameOfReferenceUID"), "value"])) {
			if (verbose) {
				warning("Reference frame mismatch")
				cat("ERROR\n")
			}			
			next
		}
		structureset <- as.character(DICOM.i[which(DICOM.i[,"name"] == "StructureSetName"), "value"])
		if (length(structureset) == 0) {
			structureset <- as.character(DICOM.i[which(DICOM.i[,"name"] == "StructureSetLabel"), "value"])
		}
		if (verbose) {
			cat("(", N, " structures identified) ", sep="")
		}
		structure.IDs <- as.numeric(structures[1:N*2-1, "value"])
		names(structure.IDs) <- structures[1:N*2, "value"]
		colors <- as.numeric(which(DICOM.i[,"name"] == "ROIDisplayColor"))
		col <- c()
		structures <- as.numeric(which(DICOM.i[,"name"] == "ReferencedROINumber")[1:N])
		contours <- as.numeric(which(DICOM.i[,"name"] == "ContourData"))
		if (!first) {
			data.old$set <- c(data.old$set, data$set)
			data.old$name <- c(data.old$name, list(data$name))
			data.old$points <- c(data.old$points, list(data$points))
		}
		else {
			data.old <- list(set=NULL, name=NULL, points=NULL)
			first <- FALSE
		}
		data <- list(set=structureset, name=names(structure.IDs), points=vector("list", N))
		for (j in 1:N) {
			data.j <- strsplit(DICOM.i[intersect(colors[j]:structures[j], contours), "value"], " ")
			if (length(data.j) < 1) {
				warning(paste("Structure '", names(structure.IDs)[j], "' is empty", sep=""))
				data$points[[which(structure.IDs == as.numeric(DICOM.i[structures[j], "value"]))]] <- NA
				col <- c(col, DICOM.i[colors[j], "value"])
				next
			}
			data.j <- lapply(data.j,
				function(x) {
					x <- as.numeric(x)
					if (length(x) < 3) {
						warning(paste("Structure '", names(structure.IDs)[j], "' is missing slices", sep=""))
						return(NA)
					}
					x <- cbind(x[1:(length(x)/3)*3-2], x[1:(length(x)/3)*3-1], x[1:(length(x)/3)*3])
					x[,2] <- sum(range(as.numeric(dimnames(CT)[[2]]))) - x[,2]
					return(x)
				}
			)			
			data$points[[which(structure.IDs == as.numeric(DICOM.i[structures[j], "value"]))]] <- data.j
			col <- c(col, DICOM.i[colors[j], "value"])
		}
		if (verbose) {
			cat("FINISHED\n")
		}
	}
	data$set <- c(data.old$set, data$set)
	data$name <- c(data.old$name, list(data$name))
	data$points <- c(data.old$points, list(data$points))

#	return(list(structures=list(col=col, data=data), DICOM=DICOMs))

	if (verbose) {
		cat("Processing (", length(unlist(data$name)), ") structures:\n", sep="")
	}
	N <- length(data$name)
	struct.list <- new("structure.list")
	if (is.null(limit)) {
		limit <- Inf
	}
	for (i in 1:N) {
		for (j in 1:length(data$name[[i]])) {
			struct.i <- data$points[[i]][[j]]
			if (length(unlist(struct.i, recursive=FALSE)) > limit) {
				if (verbose) {
					cat("  ", data$set[[i]], ": ", data$name[[i]][j], " [", length(struct.i), " axial slice(s), ", length(unlist(struct.i, recursive=FALSE)), " point(s)] ... skipped\n", sep="")
				}
				struct.list <- c(struct.list, new("structure3D", name=paste(data$name[[i]][j], data$set[[i]])))
				next
			}
			else if (length(unlist(struct.i, recursive=FALSE)) > 1) {
				if (verbose) {
					cat("  ", data$set[[i]], ": ", data$name[[i]][j], " [", length(struct.i), " axial slice(s), ", length(unlist(struct.i, recursive=FALSE)), " point(s)] ... ", sep="")
				}
				if (identical(struct.i, NA)) {
					pts.i <- matrix(nrow=0, ncol=3)
				}
				else if (is.null(dim(struct.i))) {
					pts.i <- matrix(NA, nrow=0, ncol=3)
					for (k in 1:length(struct.i)) {
						pts.i <- rbind(pts.i, struct.i[[k]])
					}
				}						
				struct.list <- c(struct.list, new("structure3D", name=paste(data$name[[i]][j], data$set[[i]]), vertices=pts.i, closed.polys=struct.i))
				if (verbose) {
					cat("FINISHED\n")
				}
			}
			else if ((length(unlist(struct.i, recursive=FALSE)) == 1) & (!is.na(struct.i))) {
				if (verbose) {
					cat("  ", data$set[[i]], ": ", data$name[[i]][j], " [", length(struct.i), " axial slice(s), ", length(unlist(struct.i, recursive=FALSE)), " point(s)] ... ", sep="")
				}
				if (identical(struct.i, NA)) {
					pts.i <- matrix(nrow=0, ncol=3)
				}
				else if (is.null(dim(struct.i))) {
					pts.i <- matrix(NA, nrow=0, ncol=3)
					for (k in 1:length(struct.i)) {
						pts.i <- rbind(pts.i, struct.i[[k]])
					}
				}						
				struct.list <- c(struct.list, new("structure3D", name=paste(data$name[[i]][j], data$set[[i]]), vertices=pts.i, closed.polys=struct.i))
				if (verbose) {
					cat("FINISHED\n")
				}
			}
			else {
				if (verbose) {
					cat("  ", data$set[[i]], ": ", data$name[[i]][j], " [EMPTY] ... FINISHED\n", sep="")
				}
				struct.list <- c(struct.list, new("structure3D", name=paste(data$name[[i]][j], data$set[[i]])))				
			}
		}
	}

#	return(list(structures=list(col=col, data=data), DICOM=DICOMs, struct.list=struct.list))
	return(struct.list)
	## FOR OTHER FILES LOAD SPECIFIC DICOM files with skipSequwnce=FALSE and re-store hdr info in DICOM list!

}

