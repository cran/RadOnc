convert.DVH <- function(..., type=c(NA, "cumulative", "differential"), dose=c(NA, "absolute", "relative"), volume=c(NA, "relative", "absolute"), dose.units=c(NA, "cGy", "Gy")) {
	type <- match.arg(type)
	dose <- match.arg(dose)
	volume <- match.arg(volume)
	dose.units <- match.arg(dose.units)
	arglist <- c(...)
	arglist <- arglist[unlist(lapply(arglist, function(arg) { ((class(arg)[1] == "DVH") & validObject(arg)) }))]
	N <- length(arglist)
	if (N <=0) {
		return(NULL)
	}
	
	for (i in 1:N) {
		x <- arglist[[i]]
		if (is.empty(x)) { 
			x@dose.type <- dose
			x@volume.type <- volume
			x@type <- type
			x@dose.units <- dose.units
			arglist[[i]] <- x
			next	
		}
		if ((!is.na(dose)) & (dose != x@dose.type)) {
			if (is.na(x@dose.rx)) {
				warning(paste("Cannot convert DVH (", x@name, ") because prescription dose is not specified", sep=""))
				arglist[[i]] <- x
				next	
			}
			if (dose == "absolute") {
				x@doses <- x@doses * x@dose.rx / x@rx.isodose
				x@dose.type <- "absolute"
			}
			else {
				x@doses <- x@doses * x@rx.isodose / x@dose.rx
				x@dose.type <- "relative"
			}
		}
		if ((!is.na(volume)) & (volume != x@volume.type)) {
			if (volume == "absolute") {
				x@volumes <- x@volumes * x@structure.volume / 100
				x@volume.type <- "absolute"
			}
			else {
				x@volumes <- 100 * x@volumes / x@structure.volume
				x@volume.type <- "relative"
			}
		}
		if ((!is.na(type)) & (type != x@type)) {
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
		if ((!is.na(dose.units)) & (dose.units != x@dose.units)) {
			if (dose.units == "cGy") {
				if (x@dose.type == "absolute") {
					x@doses <- x@doses * 100
				}
				x@dose.rx <- x@dose.rx * 100
				x@dose.max <- x@dose.max * 100
				x@dose.min <- x@dose.min * 100
				x@dose.mean <- x@dose.mean * 100
				x@dose.median <- x@dose.median * 100
				x@dose.mode <- x@dose.mode * 100
				x@dose.STD <- x@dose.STD * 100
			}
			else {
				if (x@dose.type == "absolute") {
					x@doses <- x@doses / 100
				}
				x@dose.rx <- x@dose.rx / 100
				x@dose.max <- x@dose.max / 100
				x@dose.min <- x@dose.min / 100
				x@dose.mean <- x@dose.mean / 100
				x@dose.median <- x@dose.median / 100
				x@dose.mode <- x@dose.mode / 100
				x@dose.STD <- x@dose.STD / 100
			}
			x@dose.units <- dose.units
		}
		arglist[[i]] <- x
	}
	if (N == 1) { return(arglist[[1]]) }
	return(arglist)
}