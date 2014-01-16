setGeneric("gEUD",
	function (x, a, ...) {
		standardGeneric("gEUD") 
	}
)

setMethod("gEUD", c("ANY", "missing"),
	function (x, a, ...) {
		stop("argument 'a' is missing, with no default")
	}
)

setMethod("gEUD", c("DVH", "numeric"),
	function (x, a, dose.units=c("cGy", "Gy")) {
		dose.units <- match.arg(dose.units)
		if (length(a) < 1) {
			warning("argument 'a' is of zero length")
			return(NA)
		}
		else if (length(a) > 1) {
			warning(paste("length of 'a' exceeds length of 'x', will use single value for a=",a[1], sep=""))
			a <- a[1]
		}
		if (a == 0) {
			warning("argument 'a' must be non-zero")
			return(NA)
		}
		if (is.empty(x)) {
			warning("argument 'x' is an empty DVH")
			return(NA)
		}
		x <- convert.DVH(x, type="differential", dose="absolute", volume="absolute", dose.units=dose.units)
		if (a == 1) {
			return(x$dose.mean)
		}
		else if (a == Inf) {
			return(x$dose.max)
		}
		else if (a == -Inf) {
			return(x$dose.min)
		}
		res <- sum(x$volumes * (x$doses ** a), na.rm=TRUE) / x$structure.volume
		if ((res == Inf) & (a > 0)) {
			return(x$dose.max)
		}
		else if ((res == 0) & (a < 0)) {
			return(x$dose.min)
		}
		return(res ** (1/a)) 
	}
)

setMethod("gEUD", c("DVH.list", "numeric"),
	function (x, a, dose.units=c("cGy", "Gy")) {
		dose.units <- match.arg(dose.units)
		if (length(a) < 1) {
			warning("argument 'a' is of zero length")
			return(rep(NA, length(x)))
		}
		x <- new("DVH.list", lapply(x, convert.DVH, type="differential", dose="absolute", volume="absolute", dose.units=dose.units))
		if (length(a) != length(x)) {
			if (length(a) > 1) {
				warning(paste("length of 'x' and 'a' do not match, will use single value for a=",a[1], sep=""))
			}
			a <- rep(a[1], length(x))
		}
		return(unlist(mapply(
			function (dvh, a) {
				if (a == 1) {
					return(dvh$dose.mean)
				}
				if (a == Inf) {
					return(dvh$dose.max)
				}
				if (a == -Inf) {
					return(dvh$dose.min)
				}
				if (a == 0) {
					return(NA)
				}
				if (is.empty(dvh)) {
					return(NA)
				}
				res <- sum(dvh$volumes * (dvh$doses ** a), na.rm=TRUE) / dvh$structure.volume
				if ((res == Inf) & (a > 0)) {
					return(dvh$dose.max)
				}
				else if ((res == 0) & (a < 0)) {
					return(dvh$dose.min)
				}
				return(res ** (1/a)) 
			}, x, a))
		)
	}
)
