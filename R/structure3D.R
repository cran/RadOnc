setClass("structure3D",
	representation(
		name = "character",
		volume = "numeric",
		volume.units = "character",
		coordinate.units = "character",
		vertices = "matrix",
		origin = "numeric",
		triangles = "matrix",
		closed.polys = "matrix",
		DVH = "DVH"
	),
	prototype(
		name = character(),
		volume = numeric(),
		volume.units = character(),
		coordinate.units = character(),
		vertices = matrix(),
		origin = numeric(),
		triangles = matrix(),
		closed.polys = matrix(),
		DVH = new("DVH")
	)
)

setMethod("initialize",
	"structure3D",
	function(.Object,
		name = "",
		volume = NULL,
		volume.units = c("cc"),
		coordinate.units = c("cm", "mm"),
		vertices = matrix(nrow=0, ncol=3),
		origin = NULL,
		triangles = matrix(nrow=3, ncol=0),
		closed.polys = matrix(nrow=0, ncol=3),
		DVH = new("DVH")
	) {
		.Object@name <- as.character(name)
		if (is.null(volume)) {
			.Object@volume <- 0
			# calculate volume of structure3D
			# .Object@volume <- as.numeric(volume)
		}
		else {
			.Object@volume <- as.numeric(volume)		
		}		
		if (is.null(vertices)) {
			vertices <- matrix(nrow=0, ncol=3)
		}
		if (is.null(closed.polys)) {
			closed.polys <- matrix(nrow=0, ncol=3)
		}
		if (is.null(origin)) {
			if (dim(vertices)[1] <= 1) {
				origin <- as.numeric(vertices)
			}
			else {
				origin <- apply(vertices, 2, mean)
			}
		}
		if (length(origin) != 3) {
			.Object@origin <- c(0, 0, 0)
		}
		else {
			.Object@origin <- origin
		}
		volume.units <- match.arg(volume.units)
		.Object@volume.units <- as.character(volume.units)
		coordinate.units <- match.arg(coordinate.units)
		.Object@coordinate.units <- as.character(coordinate.units)
		.Object@vertices <- as.matrix(vertices)
		.Object@triangles <- as.matrix(triangles)
		.Object@closed.polys <- as.matrix(closed.polys)
		return(.Object)
	}
)

setValidity("structure3D",
	function(object) {
		if (object@volume < 0) return(FALSE)
		if (!is.matrix(object@vertices)) return(FALSE)
		if (!is.matrix(object@triangles)) return(FALSE)
		if (dim(object@vertices)[2] != 3) return(FALSE)
		if (dim(object@triangles)[1] != 3) return(FALSE)
		if (length(object@origin) != 3) return(FALSE)
		if ((dim(object@triangles)[2] > 0) & (dim(object@vertices)[1] == 0)) return(FALSE)
#		if ((dim(object@vertices)[1] > 0) & (dim(object@triangles)[2] == 0)) return(FALSE)
		if ((dim(object@vertices)[1] > 0) & (dim(object@triangles)[2] > 0)) {
			range.triangles <- suppressWarnings(range(object@triangles))
			if (range.triangles[1] < 1) return(FALSE)
			if (range.triangles[2] > dim(object@vertices)[1]) return(FALSE)			
		}
		return(validObject(object@DVH))
	}
)

setMethod("$", "structure3D",
	function (x, name) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			return(NULL)	
		}
		else {
			return(slot(x, name))	
		}
	}
)


setMethod("names", "structure3D",
	function (x) {
		return(x$name)
	}
)


setMethod("names<-", "structure3D",
 	function (x, value) {
 		x$name <- value
 		return(x)
 	}
)


setMethod("$<-", "structure3D",
	function (x, name, value) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'structure3D'")
		}
		else {
			slot(x, name) <- value
		}
		return(x)
	}
)

setMethod("c", "structure3D",
	function (x, ..., recursive = FALSE) {
		return(c(as(x, "structure.list"), ..., recursive=FALSE))
	}
)


setMethod("range", "structure3D",
	function (x, ..., na.rm=TRUE) {
		if (dim(x$vertices)[1] > 1) {
			range <- apply(x$vertices, 2, range, na.rm=na.rm)
		}
		else {
			range <- matrix(NA, nrow=2, ncol=3)
		}
		dimnames(range) <- list(c("min", "max"), c("x", "y", "z"))
		return(range)
	}
)


setMethod("plot", c("structure3D", "missing"),
	function(x, col="gray", alpha=1, ...) {
		open3d()
		triangles3d(x$vertices[x$triangles,1], x$vertices[x$triangles,2], x$vertices[x$triangles,3], col=col, alpha=alpha)
	}
)

setMethod("plot", c("structure3D", "ANY"),
	function(x, y, col="gray", alpha=1, ...) {
		open3d()
		triangles3d(x$vertices[x$triangles,1], x$vertices[x$triangles,2], x$vertices[x$triangles,3], col=col, alpha=alpha)
	}
)

setAs("structure3D", "structure.list", 
	function(from) {
		return(new("structure.list", structures=from))
	}
)

setMethod("print", "structure3D",
	function (x, ...) {
		print(paste("Structure (", names(x), ") defined by ", dim(x$vertices)[1], " points in ", length(x$closed.polys), " axial slices", sep=""))
	}
)


setMethod("show", "structure3D",
	function (object) {
		print(object)
	}
)

setMethod("dim", "structure3D",
	function (x) {
		return(c(dim(attr(x,"vertices"))[1], length(attr(x,"closed.polys"))))
	}
)
