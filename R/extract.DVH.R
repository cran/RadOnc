extract.DVH <- function(structure, dose, resolution.xyz=c(0.2,0.2,NA), resolution.dose=0.01, method=c("ATC"), dose.units=c("cGy","Gy")) {
#	reference is Straube and Matthews and Bosche and Prudy 2005 (DVH Analysis: consequences for quality assurance of multi-insitutional clinial trials)
#	placeholder for implementation of multiple DVH calculation algorithms (currrently implemented algorithm makes use of evenly spaced voxel grid and tests every point to see whether or not inside every slice on axial-by-axial basis -- note that this assumes "straight walls" i.e. stepwise surface, no angles/curves between z slices)
#	next step will be to take 3D surface and interpolate dose along z direction as well!!!!
#   also consider different bounding box for each slice?  or use same larger bounding box for overall?  I'll implement both ways and do a timing test to see which one compares better . . . keep current function as-is as baseline check to see performance difference!!!!!!!!!
#	extract.DVH <- function(..., method=c("straube", ...)) {
#	method <- match.arg(method)	
	dose.units <- match.arg(dose.units)	
	if (!class(structure) == "structure3D") {
		warning("Input 'structure' not an object of class structure3D")
		return()
	}
	if (is.null(dose)) {
		warning("Argument 'dose' is missing with no default")
		return()
	}
	if (length(structure$closed.polys) < 1) {
		warning("Structure '", names(structure), "' is empty (it contains no pre-defined axial slices)")	
		return()	
	}
	dose.xrange <- range(as.numeric(dimnames(dose)[[1]]), na.rm=TRUE)
	dose.yrange <- range(as.numeric(dimnames(dose)[[2]]), na.rm=TRUE)
	dose.zrange <- range(as.numeric(dimnames(dose)[[3]]), na.rm=TRUE)
	range.struct <- range(structure)
	if ((range.struct[1,1] < dose.xrange[1]) | (range.struct[2,1] > dose.xrange[2]) | 
		(range.struct[1,2] < dose.yrange[1]) | (range.struct[2,2] > dose.yrange[2]) |
		(range.struct[1,3] < dose.zrange[1]) | (range.struct[2,3] > dose.zrange[2])) {
			warning("Structure '", names(structure), "' extends beyond calculated dose grid")
			return()
	}
	z.unique <- sort(unique(structure$vertices[,3]))
	if (is.na(resolution.xyz[3])) {
		resolution.xyz[3] <- median(abs(diff(z.unique)))
	}
	offset.x <- ((range.struct[2,1]-range.struct[1,1]) %% resolution.xyz[1]) / 2
	x <- seq(from=range.struct[1,1]+offset.x, to=range.struct[2,1]-offset.x, by=resolution.xyz[1])
	N.x <- length(x)
	offset.y <- ((range.struct[2,2]-range.struct[1,2]) %% resolution.xyz[2]) / 2
	y <- seq(from=range.struct[1,2]+offset.y, to=range.struct[2,2]-offset.y, by=resolution.xyz[2])
	N.y <- length(y)
#	offset.z <- ((range.struct[2,3]-range.struct[1,3]) %% resolution.xyz[3]) / 2
#	z <- seq(from=range.struct[1,3]+offset.z, to=range.struct[2,3]-offset.z, by=resolution.xyz[3])
#	N.z <- length(z)
	poly.z <- unlist(lapply(structure$closed.polys, function(poly) {return(poly[1,3])}))
	N.z <- length(z.unique)
	doses <- seq(from=min(dose), to=max(dose), by=resolution.dose)
	voxels <- c()
	for (i in z.unique) {
		poly.i <- which(poly.z == i)
		# TEST WHETHER EACH POINT IN GRID IS CONTAINED WITHIN POLYGON(S)
		pts.xyz <- cbind(rep(x, each=N.y), rep(y, N.x), i)
		results <- rep(0, N.x*N.y)
		lapply(structure$closed.polys[poly.i], function(poly) {
			results <<- results+pointInPoly2D(pts.xyz[,1:2], poly[,1:2])
		})
		pts.xyz <- pts.xyz[which(results %%2 != 0),]
		# CALCULATE (INTERPOLATE) DOSES FOR EACH POINT CONTAINED IN STRUCTURE
		dose.xyz <- approx3D(dose, pts.xyz)
		voxels <- c(voxels, dose.xyz)
	}
	dvh <- hist(voxels, breaks=c(doses, max(dose)),plot=FALSE,right=FALSE)$counts*prod(resolution.xyz)/1000
	dvh.volume <- sum(dvh)
	return(new("DVH", type="differential", dose.type="absolute", volume.type="absolute", structure.volume=dvh.volume, doses=doses, volumes=dvh, dose.max=max(voxels,na.rm=TRUE), dose.min=min(voxels,na.rm=TRUE), dose.mean=sum(dvh*doses)/dvh.volume, dose.units=dose.units, structure.name=names(structure)))
}