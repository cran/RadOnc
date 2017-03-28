setGeneric("get.HU",
	function (x, CT, ...) {
		standardGeneric("get.HU")
	}
)

setMethod("get.HU", c("RTdata", "missing"),
	function (x, CT, resolution.xyz=NA, resolution.HU=1, method=NULL) {
		return(get.HU(x$structures, x$CT, resolution.xyz, resolution.HU, method))	
	}
)

setMethod("get.HU", c("structure3D", "array"),
	function (x, CT, resolution.xyz=NA, resolution.HU=1, method=NULL) {
		if (any(!is.na(resolution.xyz)) & length(resolution.xyz) != 3) {
			warning("Argument 'resolution.xyz' must be of length 3")
			return()
		}
		method <- match.arg(method, choices=c("axial"))
		switch(method,
			axial = return(calc.HU.axial(x, CT, resolution.xyz, resolution.HU))
		)
	}
)

setMethod("get.HU", c("structure.list", "array"),
	function(x, CT, resolution.xyz=NA, resolution.HU=1, method=NULL) {
		if (any(!is.na(resolution.xyz)) & length(resolution.xyz) != 3) {
			warning("Argument 'resolution.xyz' must be of length 3")
			return()
		}
		method <- match.arg(method, choices=c("axial"))	
		switch(method,
			axial = {
				return(as(lapply(x, function(structure) {
					return(calc.HU.axial(structure, CT, resolution.xyz, resolution.HU))
				}), "list"))
			},
			{
				warning("Inappropriate method selection (", method, ")")
				return()
			}
		)
	}
)

setMethod("get.HU", c("ANY", "missing"),
	function (x, CT, ...) {
		warning("Argument 'CT' is missing with no default")
		return()
	}
)

setMethod("get.HU", c("ANY", "array"),
	function (x, CT, ...) {
		warning("Argument 'x' is not an object of class structure3D, structure.list, or RTdata")
		return()
	}
)

setMethod("get.HU", c("ANY", "ANY"),
	function (x, CT, ...) {
		warning("Improper input(s) 'x' and/or 'CT', please refer to RadOnc package documentation for further information")
		return()
	}
)


calc.HU.axial <- function(x, CT, resolution.xyz=NA, resolution.HU=1) {
	if (length(x$closed.polys) < 1) {
		warning("Structure '", names(x), "' is empty (it contains no pre-defined axial slices)")	
		return()	
	}
	if (any(!is.na(resolution.xyz)) & length(resolution.xyz) != 3) {
		warning("Argument 'resolution.xyz' must be of length 3")
		return()
	}
	CT.xrange <- range(as.numeric(dimnames(CT)[[1]]), na.rm=TRUE)
	CT.yrange <- range(as.numeric(dimnames(CT)[[2]]), na.rm=TRUE)
	CT.zrange <- range(as.numeric(dimnames(CT)[[3]]), na.rm=TRUE)
	range.struct <- range(x)
	if ((range.struct[1,1] < CT.xrange[1]) | (range.struct[2,1] > CT.xrange[2]) | 
		(range.struct[1,2] < CT.yrange[1]) | (range.struct[2,2] > CT.yrange[2]) |
		(range.struct[1,3] < CT.zrange[1]) | (range.struct[2,3] > CT.zrange[2])) {
			warning("Structure '", names(x), "' extends beyond CT coordinates")
			return()
	}
	if (is.na(resolution.xyz[1])) {
		xseq <- as.numeric(dimnames(CT)[[1]])
	}
	else {
		offset.x <- ((range.struct[2,1]-range.struct[1,1]) %% resolution.xyz[1]) / 2
		xseq <- seq(from=range.struct[1,1]+offset.x, to=range.struct[2,1]-offset.x, by=resolution.xyz[1])
	}
	N.x <- length(xseq)
	if (is.na(resolution.xyz[2])) {
		yseq <- as.numeric(dimnames(CT)[[2]])
	}
	else {
		offset.y <- ((range.struct[2,2]-range.struct[1,2]) %% resolution.xyz[2]) / 2
		yseq <- seq(from=range.struct[1,2]+offset.y, to=range.struct[2,2]-offset.y, by=resolution.xyz[2])
	}
	N.y <- length(yseq)
	z.unique <- sort(unique(x$vertices[,3]))
	if (is.na(resolution.xyz[3])) {
		resolution.xyz[3] <- median(abs(diff(z.unique)))
	}
	poly.z <- unlist(lapply(x$closed.polys, function(poly) {return(poly[1,3])}))
	N.z <- length(z.unique)
	HUs <- seq(from=min(CT), to=max(CT), by=resolution.HU)
	voxels <- c()
	for (i in z.unique) {
		poly.i <- which(poly.z == i)
		# TEST WHETHER EACH POINT IN GRID IS CONTAINED WITHIN POLYGON(S)
		pts.xyz <- cbind(rep(xseq, each=N.y), rep(yseq, N.x), i)
		results <- rep(0, N.x*N.y)
		lapply(x$closed.polys[poly.i], function(poly) {
			results <<- results+pointInPoly2D(pts.xyz[,1:2], poly[,1:2])
		})
		pts.xyz <- pts.xyz[which(results %%2 != 0),]
		# CALCULATE (INTERPOLATE) DOSES FOR EACH POINT CONTAINED IN STRUCTURE
		voxels <- c(voxels, approx3D(CT, pts.xyz))
	}
	return(hist(voxels, breaks=HUs, plot=FALSE, right=FALSE))
#	HU <- hist(voxels, breaks=HUs,plot=FALSE,right=FALSE)$counts*prod(resolution.xyz)/1000
#	print(sum(HU))
#	return(HU)
}