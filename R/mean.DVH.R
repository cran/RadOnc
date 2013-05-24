setGeneric("mean",
	mean
)

setMethod("mean", "DVH.list",
	function (x, ..., type=c("cumulative", "differential"), dose=c("absolute", "relative"), volume=c("relative", "absolute"), weighted=FALSE) {
		type <- match.arg(type)
		dose <- match.arg(dose)
		volume <- match.arg(volume)
		x <- new("DVH.list", lapply(x, convert.DVH, type=type, dose=dose, volume=volume))
		N <- length(x)
		structure.name <- paste("mean(", paste(names(x), collapse=", ", sep=""), ")", sep="")
		structure.volumes <- as.numeric(lapply(x, slot, "structure.volume"))
		structure.means <- as.numeric(lapply(x, slot, "dose.mean"))
		if (weighted) {
			dose.mean <- sum(structure.means * structure.volumes, na.rm=TRUE)/sum(structure.volumes, na.rm=TRUE)
		}
		else {
			dose.mean <- mean(structure.means, na.rm=TRUE)
		}
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		dose.rx <- max(unlist(lapply(x,slot,"dose.rx")), na.rm=TRUE)
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		if (weighted) {
			volumes.new <- apply(t(t(volume.matrix)*structure.volumes), 1, sum, na.rm=TRUE)/sum(structure.volumes)
		}
		else {
			volumes.new <- apply(volume.matrix, 1, mean, na.rm=TRUE)
		}
		new("DVH", type=type, dose.type=dose, volume.type=volume, structure.name=structure.name, structure.volume=mean(structure.volumes), dose.min=dose.min, dose.rx=dose.rx, dose.max=dose.max, dose.mean=dose.mean, doses=doses.new, volumes=volumes.new)
	}
)

setGeneric("median",
	median
)

setMethod("median", "DVH.list",
	function (x, na.rm=TRUE) {
		N <- length(x)
		if (N < 1) {
			stop("argument 'x' is missing with no default")
		}
		if (length(unique(unlist(lapply(x, slot, "type")))) > 1) {
			stop("'type' parameter must be the same for all DVH data (e.g. cumulative or differential)")
		}
		if (length(unique(unlist(lapply(x, slot, "dose.type")))) > 1) {
			stop("'dose.type' parameter must be the same for all DVH data (e.g. absolute or relative)")
		}
		if (length(unique(unlist(lapply(x, slot, "volume.type")))) > 1) {
			stop("'volume.type' parameter must be the same for all DVH data (e.g. absolute or relative)")
		}
		structure.name <- paste("median(", paste(names(x), collapse=", ", sep=""), ")", sep="")
		structure.volume <- median(as.numeric(lapply(x, slot, "structure.volume")), na.rm=na.rm)
		structure.means <- as.numeric(lapply(x, slot, "dose.mean"))
		dose.mean <- median(structure.means, na.rm=na.rm)
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		dose.rx <- max(unlist(lapply(x,slot,"dose.rx")), na.rm=TRUE)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		volumes.new <- apply(volume.matrix, 1, median, na.rm=na.rm)
		new("DVH", type=x[[1]]$type, dose.type=x[[1]]$dose.type, volume.type=x[[1]]$volume.type, structure.name=structure.name, structure.volume=structure.volume, dose.rx=dose.rx, dose.min=dose.min, dose.max=dose.max, dose.mean=dose.mean, doses=doses.new, volumes=volumes.new)
	}
)


setGeneric("mad", 
	mad
)

setMethod("mad", "DVH.list",
	function (x) {
		N <- length(x)
		if (N < 1) {
			return(NA)
		}
		if (length(unique(unlist(lapply(x, slot, "type")))) > 1) {
			stop("'type' parameter must be the same for all DVH data (e.g. cumulative or differential)")
		}
		if (length(unique(unlist(lapply(x, slot, "dose.type")))) > 1) {
			stop("'dose.type' parameter must be the same for all DVH data (e.g. absolute or relative)")
		}
		if (length(unique(unlist(lapply(x, slot, "volume.type")))) > 1) {
			stop("'volume.type' parameter must be the same for all DVH data (e.g. absolute or relative)")
		}
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		if (N == 1) {
			return(list(dose=doses.new, mad=rep(0, length(doses.new))))
		}
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		return(list(dose=doses.new, mad=apply(volume.matrix, 1, mad)))
	}
)


setGeneric("quantile", 
	quantile
)

setMethod("quantile", "DVH.list",
	function (x, type=7, ...) {
		N <- length(x)
		if (N < 1) {
			return(NA)
		}
		if (length(unique(unlist(lapply(x, slot, "type")))) > 1) {
			stop("'type' parameter must be the same for all DVH data (e.g. cumulative or differential)")
		}
		if (length(unique(unlist(lapply(x, slot, "dose.type")))) > 1) {
			stop("'dose.type' parameter must be the same for all DVH data (e.g. absolute or relative)")
		}
		if (length(unique(unlist(lapply(x, slot, "volume.type")))) > 1) {
			stop("'volume.type' parameter must be the same for all DVH data (e.g. absolute or relative)")
		}
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		return(list(dose=doses.new, quantiles=apply(volume.matrix, 1, quantile, type=type, ...)))
	}
)

setGeneric("var", 
	var
)

setMethod("var", "DVH.list",
	function (x, na.rm=TRUE) {
		N <- length(x)
		if (N < 1) {
			return(NA)
		}
		if (length(unique(unlist(lapply(x, slot, "type")))) > 1) {
			stop("'type' parameter must be the same for all DVH data (e.g. cumulative or differential)")
		}
		if (length(unique(unlist(lapply(x, slot, "dose.type")))) > 1) {
			stop("'dose.type' parameter must be the same for all DVH data (e.g. absolute or relative)")
		}
		if (length(unique(unlist(lapply(x, slot, "volume.type")))) > 1) {
			stop("'volume.type' parameter must be the same for all DVH data (e.g. absolute or relative)")
		}
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		if (N == 1) {
			return(list(dose=doses.new, var=rep(0, length(doses.new))))
		}
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		return(list(dose=doses.new, var=apply(volume.matrix, 1, var, na.rm=na.rm)))
	}
)

setGeneric("sd", 
	sd
)

setMethod("sd", "DVH.list",
	function (x, na.rm=TRUE) {
		N <- length(x)
		if (N < 1) {
			return(NA)
		}
		if (length(unique(unlist(lapply(x, slot, "type")))) > 1) {
			stop("'type' parameter must be the same for all DVH data (e.g. cumulative or differential)")
		}
		if (length(unique(unlist(lapply(x, slot, "dose.type")))) > 1) {
			stop("'dose.type' parameter must be the same for all DVH data (e.g. absolute or relative)")
		}
		if (length(unique(unlist(lapply(x, slot, "volume.type")))) > 1) {
			stop("'volume.type' parameter must be the same for all DVH data (e.g. absolute or relative)")
		}
		size <- ceiling(mean(as.numeric(lapply(x, function(DVH) { length(DVH@doses) })), na.rm=TRUE))
		dose.min <- min(x)
		dose.max <- max(x)
		doses.new <- diffinv(rep((ceiling(dose.max)-0)/(size-1), size-1), xi=0)
		if (N == 1) {
			return(list(dose=doses.new, sd=rep(0, length(doses.new))))
		}
		volume.matrix <- matrix(NA, ncol=N, nrow=size)
		for (i in 1:N) {
			volume.matrix[,i] <- approx(x[[i]]$doses, x[[i]]$volumes, doses.new, rule=2)$y
		}
		return(list(dose=doses.new, sd=apply(volume.matrix, 1, sd, na.rm=na.rm)))
	}
)

setMethod("max", "DVH.list",
	function(x, ..., na.rm=TRUE) {
		x <- c(x, ...)
		return(max(as.numeric(lapply(x, max, na.rm=na.rm)), na.rm=na.rm))
	}
)


setMethod("max", "DVH",
	function (x, na.rm=TRUE) {
		if (is.na(x@dose.max)) {
			return(max(x@doses, na.rm=na.rm))
		}
		else {
			return(x@dose.max)
		}
	}
)

setMethod("min", "DVH.list",
	function(x, ..., na.rm=TRUE) {
		x <- c(x, ...)
		return(min(as.numeric(lapply(x, min, na.rm=na.rm)), na.rm=na.rm))
	}
)


setMethod("min", "DVH",
	function (x, na.rm=TRUE) {
		if (is.na(x@dose.min)) {
			return(min(x@doses, na.rm=na.rm))
		}
		else {
			return(x@dose.min)
		}
	}
)


setMethod("range", "DVH.list",
	function (x, ..., na.rm=TRUE) {
		x <- c(x, ...)
		return(range(unlist(lapply(x, range, na.rm=na.rm)), na.rm=na.rm))
	}
)


setMethod("range", "DVH",
	function (x, ..., na.rm=TRUE) {
		if (is.na(x@dose.min) | is.na(x@dose.max)) {
			return(range(x@doses, na.rm=na.rm))
		}
		else {
			return(c(x@dose.min, x@dose.max))
		}
	}
)
