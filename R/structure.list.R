setClass("structure.list",
	representation(
		structures = "list"
	),
	prototype(
		structures = list()
	)
)


setMethod("initialize",
	"structure.list",
	function (.Object,
		structures = list(),
		...
	) {
		if ((length(structures) == 1) & (class(structures) == "structure3D")) {
			structures <- list(structures)
		}
		structs <- which(unlist(lapply(structures, class)) == "structure3D")
		if (length(structs) >= 1) {
			.Object@structures <- structures[structs]
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
		if (!all(unlist(lapply(structures, class)) == "structure3D")) return(FALSE)
		return(TRUE)
	}
)

setGeneric("as.list",
	as.list
)

setMethod("as.list", "structure.list",
	function(x, ...) {
		return(attr(x,"structures"))
	}
)

setGeneric("lapply",
	lapply
)

setMethod("lapply", "structure.list",
	function (X, FUN, ...) {
    	X <- as.list(X)
    	.Internal(lapply(X, FUN))
	}
)


setMethod("length", "structure.list",
	function (x) {
		return(length(attr(x,"structures")))
	}
)


setMethod("[", "structure.list",
	function (x, i, ...) {
		x <- attr(x,"structures")
		return(new("structure.list", x[i]))
	}
)

setMethod("$", "structure.list",
	function (x, name) {
		name <- unlist(strsplit(name, ","))
		return(lapply(x, function (struct) { struct[name] }))		
	}
)


setMethod("[[", "structure.list",
	function (x, i, exact=TRUE) {
		x <- attr(x,"structures")
		return(x[[i]])
	}
)


setMethod("[[<-", "structure.list",
	function (x, i, value) {
		x <- attr(x,"structures")
		if (class(value) == "structure3D") {
			x[[i]] <- value
		}
		else {
			stop("'value' must be an object of class 'structure3D'")
		}
		return(new("structure.list", x))
	}
)

setMethod("c", "structure.list",
	function (x, ..., recursive = FALSE) {
		return(new("structure.list", c(as.list(x), as.list(c(... , recursive=FALSE)))))
	}
)

setGeneric("rev",
	rev
)

setMethod("rev", "structure.list",
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

setMethod("print", "structure.list",
	function (x, ...) {
		print(paste("List containing ", length(x), " structure3D objects (", paste(names(x), collapse=", ", sep=""), ")", sep=""))
	}
)


setMethod("show", "structure.list",
	function (object) {
		print(object)
	}
)


setMethod("names", "structure.list",
	function (x) {
		return(as.character(unlist(lapply(x, names))))
	}
)


setMethod("names<-", "structure.list",
 	function (x, value) {
		if (length(x) != length(value)) {
			stop(paste("'names' attribute [", length(value), "] must be the same length as the structure3D list [", length(x), "]", sep=""))
		}
		struct.list <- new("structure.list", mapply(function(struct, name) {
				struct$name <- name
				return(struct)
			},
			x, value
		))
		names(attr(struct.list,"structures")) <- value
		return(struct.list)
  	}
)

setMethod("range", "structure.list",
	function (x, ..., na.rm=TRUE) {
		ranges <- lapply(x, range)
		range <- matrix(rep(c(Inf, -Inf), 3), nrow=2, ncol=3, dimnames=list(c("min", "max"), c("x", "y", "z")))
		for (i in 1:length(ranges)) {
			range[1, ] <- pmin(range[1, ], ranges[[i]][1, ], na.rm=na.rm)
			range[2, ] <- pmax(range[2, ], ranges[[i]][2, ], na.rm=na.rm)
		}
		return(range)
	}
)
