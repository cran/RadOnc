setGeneric("LQE",
	function (x, aB, ...) {
		standardGeneric("LQE") 
	}
)

setMethod("LQE", c("ANY", "missing"),
	function (x, aB, ...) {
		stop("argument 'aB' is missing, with no default")
	}
)

setMethod("LQE", c("numeric", "numeric"),
	function (x, aB, fractions=NULL) {
		if (is.null(fractions)) {
			stop("argument 'fractions' is missing, with no default")		
		}
		if (length(fractions) != 2) {
			warning("argument 'fractions' must be of length two")
			return(NA)
		}
		fractions <- as.integer(fractions)
		if (any(fractions < 1)) {
			warning("argument 'fractions' must specify two positive integer values")
			return(NA)
		}
		if (length(aB) < 1) {
			warning("argument 'aB' is of zero length")
			return(NA)
		}
		if (aB == 0) {
			warning("argument 'aB' must be non-zero")
			return(NA)
		}
		x <- x * (1 + x / (fractions[1] * aB))
		x <- suppressWarnings(unlist(lapply(x,
				function(dose) {
					dose <- as.numeric(polyroot(c(-dose, 1, 1/(fractions[2] * aB))))
					return(dose[dose >= 0])
				}
			)))
		return(x)
	}
)

setMethod("LQE", c("DVH", "numeric"),
	function (x, aB, fractions=NULL, dose.units=c("cGy", "Gy")) {
		dose.units <- match.arg(dose.units)
		if (is.null(fractions) || (as.integer(fractions) < 1)) {
			warning("argument 'fractions' must be positive integer value")
			return(NA)
		}
		else {
			fractions <- as.integer(fractions)
		}
		if (length(aB) < 1) {
			warning("argument 'aB' is of zero length")
			return(NA)
		}
		else if (length(aB) > 1) {
			warning(paste("length of 'aB' exceeds length of 'x', will use single value for aB=",aB[1], sep=""))
			aB <- aB[1]
		}
		if (aB == 0) {
			warning("argument 'aB' must be non-zero")
			return(NA)
		}
		if (is.empty(x)) {
			warning("argument 'x' is an empty DVH")
			return(NA)
		}
		x <- convert.DVH(x, dose="absolute", dose.units=dose.units)
		if (x$dose.fx == 0) {
			x$dose.fx <- fractions
		}
		else if (x$dose.fx != fractions) {
			LQE.calc <- function (doses) {
				doses <- doses * (1 + doses / (x$dose.fx * aB))
				doses <- suppressWarnings(unlist(lapply(doses,
					function(dose) {
						dose <- as.numeric(polyroot(c(-dose, 1, 1/(fractions * aB))))
						return(dose[dose >= 0])
					}
				)))
			}

			x$doses <- LQE.calc(x$doses)
			x$dose.max <- LQE.calc(x$dose.max)
			x$dose.min <- LQE.calc(x$dose.min)
			x$dose.mean <- LQE.calc(x$dose.mean)
			x$dose.median <- LQE.calc(x$dose.median)
			x$dose.mode <- LQE.calc(x$dose.mode)
			x$dose.STD <- LQE.calc(x$dose.STD)
			x$dose.rx <- LQE.calc(x$dose.rx)
			x$dose.fx <- fractions		
		}
		return(x) 
	}
)

setMethod("LQE", c("DVH.list", "numeric"),
	function (x, aB, fractions=NULL, dose.units=NULL) {
		if (length(aB) != length(x)) {
			if (length(aB) > 1) {
				warning(paste("length of 'x' and 'aB' do not match, will use single value for aB=",aB[1], sep=""))
			}
			aB <- rep(aB[1], length(x))
		}
		if (is.null(fractions)) {
			warning("argument 'fractions' is missing, with no default")
			return(NA)
		}
		if (length(fractions) != length(x)) {
			if (length(fractions) > 1) {
				warning(paste("length of 'x' and 'fractions' do not match, will use single value for fractions=",fractions[1], sep=""))
			}
			fractions <- rep(fractions[1], length(x))
		}
		dose.units <- match.arg(dose.units, choices=c("cGy","Gy"), several.ok=TRUE)
		if (length(dose.units) != length(x)) {
			if (length(dose.units) > 1) {
				warning(paste("length of 'x' and 'dose.units' do not match, will use single value for dose.units=",dose.units[1], sep=""))
			}
			dose.units <- rep(dose.units[1], length(x))
		}
		return(new("DVH.list", mapply(LQE, x, aB=aB, fractions=fractions, dose.units=dose.units)))
	}
)
