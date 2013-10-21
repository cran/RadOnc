compareStructures <- function(structures, method=c("grid", "surface", "hausdorff"), hausdorff.method=c("mean", "median", "absolute"), verbose=TRUE, plot=TRUE, pixels=100) {
	if (class(structures) != "structure.list") {
		stop("Input 'structures' must be of class 'structure.list'")
	}
	if (length(structures) < 2) {
		stop("Need at least 2 structures to perform comparison")
	}
	method <- match.arg(method)
	switch(method,
		grid = contours <- compareStructures.grid(structures, pixels=pixels),
		surface = contours <- compareStructures.surface(structures),
		hausdorff = contours <- compareStructures.hausdorff(structures, method=hausdorff.method, verbose=verbose)
	)
	if ((plot) & (method == "grid")) {
		mar.old <- par()$mar
		par(mar=c(0, 0, 0, 0))
		N.z <- length(unique(contours[,3]))
		layout(matrix(c(1:N.z*2, 1:N.z*2-1), nrow=N.z, ncol=2), widths=c(1, 10), heights=1)
		for (i in 1:N.z) {
			z.i <- unique(contours[,3])[i]
			contours.i <- contours[which(contours[,3] == z.i), ]
			sum.i <- apply(contours.i[, 4:(dim(contours.i)[2])], 1, sum)
			x <- unique(contours.i[, 1])
			y <- unique(contours.i[, 2])
			lvl.i <- matrix(sum.i, nrow=length(x), ncol=length(y))
			lvl.i <- contourLines(x=x, y=y, z=lvl.i)
			plot(range(x), range(y), type="n", xaxt="n", yaxt="n")
			for (j in 1:length(lvl.i)) {
				cl.j <- lvl.i[[j]]
				polygon(cl.j$x, cl.j$y, rep(z.i, length(cl.j$x)), col=rev(heat.colors(length(structures)*2+1))[cl.j$level*2], border=NA)
			}
			plot(1,type="n",xaxt="n",yaxt="n")
			text(1, labels=i)
		}
		par(mar=mar.old)
	}
	return(contours)
}	

compareStructures.surface <- function (structures) {

	pointInPoly2D <- function (x, y, points) {
		if (length(x) != length(y)) {
			warning("'x' and 'y' lengths differ")
		}
		if (length(x) + length(y) < 1) {
			return()
		}
		n <- dim(points)[1]
		c <- v <- rep(0, min(length(x), length(y), na.rm=TRUE))
		j <- n
		for (i in 1:n) {
			c <- abs(c - as.numeric(((y<points[i,2]) != (y<points[j,2])) * ((x<(points[j,1]-points[i,1])*(y-points[i,2])/(points[j,2]-points[i,2])+points[i,1]))))
			v <- v+(y==points[i,2])*(x==points[i,1])
			j <- i
		}
		c[which(is.na(c))] <- FALSE
		c[v>0] <- TRUE
		return(as.logical(c)[1:min(length(x), length(y))])
	}
	
	N <- length(structures)
	z <- as.list(rep(NA, N))
	pts <- matrix(nrow=0, ncol=3, dimnames=list(NULL, c("X", "Y", "Z")))
	for (i in 1:N) {
		if (length(structures[[i]]$vertices) < 1) {
			next
		}
		z[[i]] <- unlist(lapply(structures[[i]]$closed.polys, function(closed.poly) {return(unique(closed.poly[,3]))}))
		pts <- rbind(pts, structures[[i]]$vertices)
	}
	results <- matrix(0, nrow=dim(pts)[1], ncol=N, dimnames=list(NULL, names(structures)))
	for (i in 1:N) {
#		plot3d(structures[[i]]$vertices,col="gray",cex=0.2)
		for (j in unique(z[[i]])) {
			pts.j <- pts[which(pts[, 3]== j), 1:2]
			results.j <- rep(0, dim(pts.j)[1])
			z.j <- which(z[[i]] == j)
			## THIS LOOP ACCOUNTS FOR AXIAL SLICES WITH MULTIPLE SEPARATE CLOSED POLYGONS (e.g. 3 ROOTS FOR SINGLE TOOTH)
			## THIS LOOP DOES NOT(!!!) ACCOUNT FOR DONUTS (E.G. STRUCTURES WITH HOLE IN THEM -- NEED TO FIGURE OUT HOW THOSE ARE STORED FIRST) -- IF STRUCTURE HAS A HOLE, ALL BETS ARE OFF AT THE MOMENT... SOLUTION WILL BE TO DO LOGICAL SUBTRACTION RATHER THAN ADDITION OF RESULTS
			for (k in 1:length(z.j)) {
#				poly.jk <- structures[[i]]$closed.polys[[z.j[k]]][, 1:2]
				results.j <- results.j + as.numeric(pointInPoly2D(pts.j[,1], pts.j[,2], structures[[i]]$closed.polys[[z.j[k]]][, 1:2]))
			}
			results[which(pts[, 3]== j), i] <- results[which(pts[, 3]== j), i]+results.j
		}
	}
#	points3d(pts, col=rainbow(n=3)[apply(results,1,sum)])
#	points3d(pts[which(apply(results,1,sum)==3),],col="black",cex=2)
	return(cbind(pts, results))
}

compareStructures.grid <- function (structures, pixels=100) {

	pointInPoly2D <- function (x, y, points) {
		if (length(x) != length(y)) {
			warning("'x' and 'y' lengths differ")
		}
		if (length(x) + length(y) < 1) {
			return()
		}
		n <- dim(points)[1]
		c <- v <- rep(0, min(length(x), length(y), na.rm=TRUE))
		j <- n
		for (i in 1:n) {
			c <- abs(c - as.numeric(((y<points[i,2]) != (y<points[j,2])) * ((x<(points[j,1]-points[i,1])*(y-points[i,2])/(points[j,2]-points[i,2])+points[i,1]))))
			v <- v+(y==points[i,2])*(x==points[i,1])
			j <- i
		}
		c[which(is.na(c))] <- FALSE
		c[v>0] <- TRUE
		return(as.logical(c)[1:min(length(x), length(y))])
	}
	
	N <- length(structures)
	z <- as.list(rep(NA, N))
	bounds <- range(structures)
	x.coords <- seq(from=bounds[1,1], to=bounds[2,1], length.out=pixels)
	y.coords <- seq(from=bounds[1,2], to=bounds[2,2], length.out=pixels)
	for (i in 1:N) {
		if (length(structures[[i]]$vertices) < 1) {
			next
		}
		z[[i]] <- unlist(lapply(structures[[i]]$closed.polys, function(closed.poly) {return(unique(closed.poly[,3]))}))
	}
	z.coords <- unique(unlist(z))
	pts <- matrix(nrow=length(x.coords)*length(y.coords)*length(z.coords), ncol=3, dimnames=list(NULL, c("X", "Y", "Z")))
	pts <- matrix(c(rep(x.coords, each=length(y.coords)*length(z.coords)), rep(rep(y.coords, each=length(z.coords)), length(x.coords)), rep(z.coords, length(x.coords)*length(y.coords))), nrow=length(x.coords)*length(y.coords)*length(z.coords), ncol=3, dimnames=list(NULL, c("X", "Y", "Z"))) 
	results <- matrix(0, nrow=dim(pts)[1], ncol=N, dimnames=list(NULL, names(structures)))
	for (i in 1:N) {
		for (j in unique(z[[i]])) {
			pts.j <- pts[which(pts[, 3]== j), 1:2]
			results.j <- rep(0, dim(pts.j)[1])
			z.j <- which(z[[i]] == j)
			## THIS LOOP ACCOUNTS FOR AXIAL SLICES WITH MULTIPLE SEPARATE CLOSED POLYGONS (e.g. 3 ROOTS FOR SINGLE TOOTH)
			## THIS LOOP DOES NOT(!!!) ACCOUNT FOR DONUTS (E.G. STRUCTURES WITH HOLE IN THEM -- NEED TO FIGURE OUT HOW THOSE ARE STORED FIRST) -- IF STRUCTURE HAS A HOLE, ALL BETS ARE OFF AT THE MOMENT... SOLUTION WILL BE TO DO LOGICAL SUBTRACTION RATHER THAN ADDITION OF RESULTS
			for (k in 1:length(z.j)) {
#				poly.jk <- structures[[i]]$closed.polys[[z.j[k]]][, 1:2]
				results.j <- results.j + as.numeric(pointInPoly2D(pts.j[,1], pts.j[,2], structures[[i]]$closed.polys[[z.j[k]]][, 1:2]))
			}
			results[which(pts[, 3]== j), i] <- results[which(pts[, 3]== j), i]+results.j
		}
	}
	return(cbind(pts, results))
}


compareStructures.hausdorff <- function (structures, verbose=TRUE, method=c("mean", "median", "absolute")) {

	hausdorff.dist <- function (A, B, method) {
		if (ncol(A) != ncol(B)){
			warning("Dimensionality of A and B must be the same")
			return(NA)
		}
		compute.dist = function (a0, B0){
			C0 <- matrix(rep(a0, each=nrow(B0)),byrow=F, ncol=ncol(B0))
			return(min(apply(C0-B0, 1, function(x) {sqrt(sum(t(x)*x))}), na.rm=TRUE))
		}
	
		if (method == "mean") {
			d1 <- apply(A, 1, compute.dist, B0=B)
			d2 <- apply(B, 1, compute.dist, B0=A)
			return(mean(c(d1, d2), na.rm=TRUE))
		}
		else if (method == "median") {
			d1 <- apply(A, 1, compute.dist, B0=B)
			d2 <- apply(B, 1, compute.dist, B0=A)
			return(median(c(d1, d2), na.rm=TRUE))
		}
		else if (method == "absolute") {
			d1 <- max(apply(A, 1, compute.dist, B0=B))
			d2 <- max(apply(B, 1, compute.dist, B0=A))
			return(max(d1, d2, na.rm=TRUE))
		}
		else {
			warning("Invalid 'method' argument; must be one of 'mean', 'median', or 'absolute'")
			return(NA)
		}
	}

	method <- match.arg(method)
	N <- length(structures)
	results <- matrix(0, nrow=N, ncol=N, dimnames=list(names(structures), names(structures)))
	for (i in 1:N) {
		if (verbose) {
			cat("Analyzing structure ", i, "/", N, " (", structures[[i]]$name, ") ... ", sep="")
		}
		for (j in i:N) {
			results[i, j] <- hausdorff.dist(structures[[i]]$vertices, structures[[j]]$vertices, method=method)
			results[j, i] <- results[i, j]
		}
		if (verbose) {
			cat("FINISHED\n")
		}
	}
	return(results)
}
