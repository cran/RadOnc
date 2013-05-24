setClass("DVH.list",
	representation(
		structures = "list"
	),
	prototype(
		structures = list()
	)
)


setMethod("initialize",
	"DVH.list",
	function (.Object,
		structures = list(),
		...
	) {
		if ((length(structures) == 1) & (class(structures) == "DVH")) {
			structures <- list(structures)
		}
		DVHs <- which(unlist(lapply(structures, class)) == "DVH")
		if (length(DVHs) >= 1) {
			.Object@structures <- structures[DVHs]
		}
		else {
			.Object@structures <- list()
		}
		return(.Object)
	}
)

setValidity("DVH.list",
	function(object) {
		if (!is.list(object)) return(FALSE)
		if (length(object) == 0) return(TRUE)
		if (!all(unlist(lapply(structures, class)) == "DVH")) return(FALSE)
		return(TRUE)
	}
)


setGeneric("as.list",
	as.list
)

setMethod("as.list", "DVH.list",
	function(x, ...) {
		return(attr(x,"structures"))
	}
)

setGeneric("lapply",
	lapply
)

setMethod("lapply", "DVH.list",
	function (X, FUN, ...) {
    	X <- as.list(X)
    	.Internal(lapply(X, FUN))
	}
)



setMethod("length", "DVH.list",
	function (x) {
		return(length(attr(x,"structures")))
	}
)


setMethod("[", "DVH.list",
	function (x, i, ...) {
		x <- attr(x,"structures")
		return(new("DVH.list", x[i]))
	}
)

setMethod("$", "DVH.list",
	function (x, name) {
		name <- unlist(strsplit(name, ","))
		return(lapply(x, function (DVH) { DVH[name] }))		
	}
)


setMethod("[[", "DVH.list",
	function (x, i, exact=TRUE) {
		x <- attr(x,"structures")
		return(x[[i]])
	}
)

setMethod("[[<-", "DVH.list",
	function (x, i, value) {
		x <- attr(x,"structures")
		if (class(value) == "DVH") {
			x[[i]] <- value
		}
		else {
			stop("'value' must be an object of class 'DVH'")
		}
		return(new("DVH.list", x))
	}
)

setMethod("c", "DVH.list",
	function (x, ..., recursive = FALSE) {
		return(new("DVH.list", c(as.list(x), as.list(c(... , recursive=FALSE)))))
	}
)

setGeneric("rev",
	rev
)

setMethod("rev", "DVH.list",
	function (x) {
		if (length(x) <= 1) {
			return(x)
		}
		else {
			return(x[length(x):1])
		}
	}
)

setGeneric("print",
	print
)

setMethod("print", "DVH.list",
	function (x, ...) {
		print(paste("List containing ", length(x), " DVH objects (", paste(names(x), collapse=", ", sep=""), ")", sep=""))
	}
)

setMethod("show", "DVH.list",
	function (object) {
		print(object)
	}
)

setMethod("names", "DVH.list",
	function (x) {
		return(as.character(unlist(lapply(x, names))))
	}
)

setMethod("names<-", "DVH.list",
 	function (x, value) {
		if (length(x) != length(value)) {
			stop(paste("'names' attribute [", length(value), "] must be the same length as the DVH list [", length(x), "]", sep=""))
		}
		return(new("DVH.list", mapply(function(DVH, name) {
				DVH$structure.name <- name
				return(DVH)
			},
			x, value
		)))
  	}
)

