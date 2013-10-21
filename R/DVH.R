setClass("DVH",
	representation(
		structure.name = "character",
		structure.volume = "numeric",
		type = "character",
		dose.max = "numeric",
		dose.min = "numeric",
		dose.mean = "numeric",
		dose.median = "numeric",
		dose.mode = "numeric",
		dose.STD = "numeric",
		conf.index = "numeric",
		equiv.sphere = "numeric",
		gradient = "numeric",
		dose.rx = "numeric",
		dose.fx = "numeric",
		doses = "numeric",
		dose.type = "character",
		volumes = "numeric",
		volume.type = "character"
	),
	prototype(
		structure.name = character(),
		structure.volume = numeric(),
		type = character(),
		dose.max = numeric(),
		dose.min = numeric(),
		dose.mean = numeric(),
		dose.median = numeric(),
		dose.mode = numeric(),
		dose.STD = numeric(),
		conf.index = numeric(),
		equiv.sphere = numeric(),
		gradient = numeric(),
		dose.rx = numeric(),
		dose.fx = numeric(),
		doses = numeric(),
		dose.type = character(),
		volumes = numeric(),
		volume.type = character()
	)
)

setMethod("initialize",
	"DVH",
	function (.Object,
		structure.name = "",
		structure.volume = numeric(),
		type = c("cumulative", "differential"),
		dose.max = NA,
		dose.min = NA,
		dose.mean = numeric(),
		dose.median = numeric(),
		dose.mode = numeric(),
		dose.STD = numeric(),
		conf.index = numeric(),
		equiv.sphere = numeric(),
		gradient = numeric(),
		dose.rx = numeric(),
		dose.fx = numeric(),
		doses = numeric(),
		dose.type = c("absolute", "relative"),
		volumes = numeric(),
		volume.type = c("relative", "absolute"),
		...
	) {
		.Object@structure.name <- as.character(structure.name)
		.Object@structure.volume <- as.numeric(structure.volume)
		.Object@type <- match.arg(type)
		if (length(doses) > 0) {
			if (is.na(dose.max)) dose.max <- range(doses)[2]
			if (is.na(dose.min)) dose.min <- range(doses)[1]			
			.Object@doses <- doses	
		}
		else {
			.Object@doses <- numeric()
		}
		.Object@dose.max <- max(0, dose.max, na.rm=TRUE)		
		.Object@dose.min <- max(0, dose.min, na.rm=TRUE)
		.Object@dose.mean <- max(0, dose.mean, na.rm=TRUE)
		.Object@dose.median <- max(0, dose.median, na.rm=TRUE)
		.Object@dose.mode <- max(0, dose.mode, na.rm=TRUE)
		.Object@dose.STD <- max(0, dose.STD, na.rm=TRUE)
		.Object@conf.index <- max(0, conf.index, na.rm=TRUE)
		.Object@equiv.sphere <- max(0, equiv.sphere, na.rm=TRUE)
		.Object@gradient <- max(0, gradient, na.rm=TRUE)
		.Object@dose.rx <- max(0, dose.rx, na.rm=TRUE)
		.Object@dose.fx <- max(0, dose.fx, na.rm=TRUE)
		.Object@dose.type <- match.arg(dose.type)
		.Object@volume.type <- match.arg(volume.type)
		if (length(volumes) > 0) {
			.Object@volumes <- as.numeric(volumes)
			if ((.Object@structure.volume < max(.Object@volumes, na.rm=TRUE)) & (.Object@volume.type == "absolute")) {
				.Object@structure.volume <- max(.Object@volumes, na.rm=TRUE)
			}
		}
		else {
			.Object@volumes <- numeric()	
		}
		return(.Object)
	}
)


setValidity("DVH",
	function(object) {
		if (length(object@doses) != length(object@volumes)) return(FALSE)
		if (length(object@doses) == 0) return(TRUE)
		if (length(object@doses) < 2) return(FALSE)
		if (any(is.na(object@doses))) return(FALSE)
		if (object@dose.rx <= 0) return(FALSE)
		if (object@dose.min > object@dose.max) return(FALSE)
		if ((object@dose.mean > object@dose.max) | (object@dose.mean < object@dose.min)) return(FALSE)
		if (!identical(order(object@doses, decreasing=FALSE), 1:length(object@doses))) return(FALSE)
		if ((object@dose.type == "relative") & (range(object@doses, na.rm=TRUE)[2] > 250)) return(FALSE)
		if (any(is.na(object@volumes))) return(FALSE)		
		# ENSURE RELATIVE DVH VOLUMES ARE ON SCALE UP TO 100% MAXIMUM (VOLUME SHOULD NEVER BE >100%)			
		if ((object@volume.type == "relative") & (max(object@volumes, na.rm=TRUE) > 100.00000000001)) return(FALSE)	
		# ENSURE STRUCTURE VOLUME IS SUFFICIENTLY LARGE TO ENCOMPASS ALL LISTED DVH INFORMATION	
		if ((object@volume.type == "absolute") & (object@structure.volume < floor(max(object@volumes, na.rm=TRUE)))) return(FALSE)
		# ENSURE CUMULATIVE DOSE HISTOGRAM HAS APPROPRIATE DATA (DOSE RANGE MUST START AT 0)	
		if ((object@type == "cumulative") & (object@doses[1] != 0)) return(FALSE)
		# ENSURE STRUCTURE VOLUME AND DVH VOLUME DATA ARE EQUIVALENT (TOLERANCE=0.1%)
		if (!grepl("(mean|median)[(].*[)]", object@structure.name)) {
			if ((object@type == "differential") & (object@volume.type == "relative") & (abs(sum(object@volumes, na.rm=TRUE) - 100) > 0.1)) return(FALSE)
			if ((object@type == "differential") & (object@volume.type == "absolute") & (abs(sum(object@volumes, na.rm=TRUE) - object@structure.volume) / object@structure.volume > 0.001)) return(FALSE)
		}
		return(TRUE)
	}
)


setMethod("$", "DVH",
	function (x, name) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			return(NULL)	
		}
		else {
			return(slot(x, name))	
		}
	}
)

setMethod("names", "DVH",
	function (x) {
		return(x$structure.name)
	}
)

setMethod("names<-", "DVH",
 	function (x, value) {
 		x$structure.name <- value
 		return(x)
 	}
)

setMethod("$<-", "DVH",
	function (x, name, value) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'DVH'")
		}
		else {
			slot(x, name) <- value
		}
		return(x)
	}
)

setMethod("[", "DVH",
	function (x, i, ...) {
		if (!validObject(x)) {
			stop("not a valid object of 'DVH'")
		}
		if (!missing("i")) {
			if (is.null(i)) return(NULL)
			if (length(i) < 1) return(numeric())
			if (is.logical(i)) {
				if (length(i) != length(x@doses)) {
					stop("(subscript 'i') logical subscript does not match length of 'DVH' doses")
				}
				i <- which(i)
			}
			result <- c()
			result.units <- c()
			i <- toupper(as.character(i))
			type <- sub("(V|D).*", "\\1", i)
			volume <- grepl("VOL", i)
			if (any(volume)) {
				type[volume] <- "VOLUME"
			}
			value <- sub("[VD]([.0-9]+|MAX|MIN|MEAN|RX)[^0-9]*", "\\1", i)
			type2 <- sub("[VD]([.0-9]+|MAX|MIN|MEAN|RX)([%]|GY|CGY|CC)$", "\\2", i)
			for (count in 1:length(i)) {
				switch(type[count],
					VOLUME = {
						result <- c(result, x@structure.volume)
						result.units <- c(result.units, "cc")
					},
					V = {
						value[count] <- suppressWarnings(as.numeric(value[count]))
						if (is.na(value[count])) {
							result <- c(result, NA)
							result.units <- c(result.units, NA)
							next
						}
						switch(x@type,
							cumulative = {
								switch(x@dose.type,
									absolute = {
										switch(type2[count],
											CGY = TRUE,
											"%" = value[count] <- as.numeric(value[count]) * x@dose.rx / 100,
											GY = value[count] <- as.numeric(value[count]) * 100,
											CC = value[count] <- NA,
											value[count] <- NA
										)
									},
									relative = {
										switch(type2[count],
											"%" = TRUE,
											GY = value[count] <- (as.numeric(value[count]) * 100 / x@dose.rx) * 100,
											CGY = value[count] <- (as.numeric(value[count]) / x@dose.rx) * 100,
											CC = value[count] <- NA,
											value[count] <- NA
										)
									},
									value[count] <- NA
								)
								result <- c(result, approx(x@doses, x@volumes, value[count], yright=0)$y)
								if (is.na(value[count])) {
									result.units <- c(result.units, NA)
								}
								else if (x@volume.type == "absolute") {
									result.units <- c(result.units, "cc")
								}
								else {
									result.units <- c(result.units, "%")										
								}
							},
							differential = {
								warning("No method available to extract volume given differential doses")
								result <- c(result, NA)
								result.units <- c(result.units, NA)
							}
						)
					},
					D = {
						switch(x@type,
							cumulative = {
								if (value[count] %in% c("MAX", "MIN", "MEAN", "RX")) {
									switch(value[count],
										MAX = value[count] <- x@dose.max,
										MIN = value[count] <- x@dose.min,
										MEAN = value[count] <- x@dose.mean,
										RX = value[count] <- x@dose.rx
									)
									if (type2[count] == "%") {
										result <- c(result, 100 * as.numeric(value[count]) / x@dose.rx)
									}
									else {
										result <- c(result, as.numeric(value[count]))
									}
									if ((x@dose.type == "absolute") & (type2[count] != "%")) {
										result.units <- c(result.units, "Gy")								
									}
									else {
										result.units <- c(result.units, "%")
									}
									next
								}
								value[count] <- suppressWarnings(as.numeric(value[count]))
								if (is.na(value[count])) {
									result <- c(result, NA)
									result.units <- c(result.units, NA)
									next
								}
								switch(x@volume.type,
									absolute = {
										switch(type2[count],
											CC = TRUE,
											"%" = value[count] <- as.numeric(value[count]) * x@structure.volume / 100,
											GY = value[count] <- NA,
											CGY = value[count] <- NA,
											value[count] <- NA
										)									
									},
									relative = {
										switch(type2[count],
											"%" = TRUE,
											CC = value[count] <- as.numeric(value[count]) * 100 / x@structure.volume,
											GY = value[count] <- NA,
											CGY = value[count] <- NA,
											value[count] <- NA
										)
									},
									value[count] <- NA
								)
								value.count <- as.numeric(value[count])
								if (((value.count >= 100) & (x@volume.type == "relative")) | ((value.count >= x@structure.volume) & (x@volume.type == "absolute"))) {
									result <- c(result, x@dose.min)
								}
								else {
									result <- c(result, approx(x@volumes, x@doses, value.count)$y)			
								}
								if (x@dose.type == "absolute") {
									result.units <- c(result.units, "cGy")
								}
								else {
									result.units <- c(result.units, "%")										
								}
							},
							differential = {
								warning("No method available to extract dose given differential doses")
								result <- c(result, NA)
								result.units <- c(result.units, NA)
							},
							{
								result <- c(result, NA)								
								result.units <- c(result.units, NA)
							}
						)
					},
					{
						result <- c(result, NA)
						result.units <- c(result.units, NA)
					}
				)
			}
		}
		else {
			return()
		}
		names(result) <- result.units
		return(result)
	}
)


setMethod("c", "DVH",
	function (x, ..., recursive = FALSE) {
		return(c(as(x, "DVH.list"), ..., recursive=FALSE))
	}
)

setGeneric("print",
	print
)

setMethod("print", "DVH",
	function (x, ...) {
		if (x@dose.type == "relative") { dose.type <- "%" } else { dose.type <- "cGy" }
		print(paste("Structure: ", x@structure.name, " (", x@structure.volume, " cc), Dose: ", x@dose.min, "-", x@dose.max, dose.type, " (", x@dose.rx, "cGy prescribed), DVH: ", x@type, ", Volume: ", x@volume.type, sep=""))	
	}
)

setMethod("show", "DVH",
	function (object) {
		print(object)
	}
)

convert.DVH <- function(..., type=c("cumulative", "differential"), dose=c("absolute", "relative"), volume=c("relative", "absolute")) {
	type <- match.arg(type)
	dose <- match.arg(dose)
	volume <- match.arg(volume)
	arglist <- c(...)
	arglist <- arglist[unlist(lapply(arglist, function(arg) { ((class(arg)[1] == "DVH") & validObject(arg)) }))]
	N <- length(arglist)
	if (N <=0) {
		return(NULL)
	}
	
	for (i in 1:N) {
		x <- arglist[[i]]
		if (dose != x@dose.type) {
			if (dose == "absolute") {
				x@doses <- x@doses * x@dose.rx / 100
				x@dose.max <- x@dose.max * x@dose.rx / 100
				x@dose.min <- x@dose.min * x@dose.rx / 100
				x@dose.mean <- x@dose.mean * x@dose.rx / 100
				x@dose.median <- x@dose.median * x@dose.rx / 100
				x@dose.mode <- x@dose.mode * x@dose.rx / 100
				x@dose.STD <- x@dose.STD * x@dose.rx / 100
				x@dose.type <- "absolute"
			}
			else {
				x@doses <- 100 * x@doses / x@dose.rx
				x@dose.max <- 100 * x@dose.max / x@dose.rx
				x@dose.min <- 100 * x@dose.min / x@dose.rx
				x@dose.mean <- 100 * x@dose.mean / x@dose.rx
				x@dose.median <- 100 * x@dose.median / x@dose.rx
				x@dose.mode <- 100 * x@dose.mode / x@dose.rx
				x@dose.STD <- 100 * x@dose.STD / x@dose.rx
				x@dose.type <- "relative"
			}
		}
		if (volume != x@volume.type) {
			if (volume == "absolute") {
				x@volumes <- x@volumes * x@structure.volume / 100
				x@volume.type <- "absolute"
			}
			else {
				x@volumes <- 100 * x@volumes / x@structure.volume
				x@volume.type <- "relative"
			}
		}
		if (type != x@type) {
			if (type == "cumulative") {
				temp.doses <- x@doses - diff(c(-x@doses[1], x@doses))/2
				x@doses <- c(temp.doses, (2*x@doses - temp.doses)[length(temp.doses)])
				if (volume == "relative") {
					x@volumes <- diffinv(-x@volumes, xi=100)
				}
				else {
					x@volumes <- diffinv(-x@volumes, xi=x@structure.volume)
				}
				x@type <- "cumulative"
			}
			else {
				x@volumes <- -diff(x@volumes)
				x@doses <- x@doses[1:(length(x@doses)-1)] + diff(x@doses)/2
				x@type <- "differential"
			}
		}
		arglist[[i]] <- x
	}
	if (N == 1) { return(arglist[[1]]) }
	return(arglist)
}
