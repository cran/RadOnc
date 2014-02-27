plot.DVH <- plot.DVH.list <- function(x, ..., plot.type=c("individual", "grouped", "ttest", "wilcox", "bars")) {
	plot.type <- match.arg(plot.type)
	switch(plot.type,
		individual = {
			plot.DVH.individual(x, ...)
		},
		grouped = {
			plot.DVH.groups(x, ...)
		},
		ttest = {
			plot.DVH.ttest(x, ...)
		},
		wilcox = {
			plot.DVH.wilcox(x, ...)
		},
		bars = {
			plot.DVH.bars(x, ...)
		}
	)
}

setGeneric("plot",
	plot
)

setMethod("plot", c("DVH", "missing"),
	function(x, ...) {
		plot.DVH(new("DVH.list",list(x)), ...)
	}
)

setMethod("plot", c("zDVH", "missing"),
	function(x, ...) {
		plot.zDVH(x, ...)
	}
)

setMethod("plot", c("DVH", "DVH"),
	function(x, y, ...) {
		plot.DVH(new("DVH.list",list(x)), new("DVH.list",list(y)), ...)
	}
)

setMethod("plot", c("zDVH", "DVH"),
	function(x, y, ...) {
		plot.DVH(new("DVH.list",list(x)), new("DVH.list",list(y)), ...)
	}
)

setMethod("plot", c("DVH", "zDVH"),
	function(x, y, ...) {
		plot.DVH(new("DVH.list",list(x)), new("DVH.list",list(y)), ...)
	}
)

setMethod("plot", c("zDVH", "zDVH"),
	function(x, y, ...) {
		plot.DVH(new("DVH.list",list(x)), new("DVH.list",list(y)), ...)
	}
)


setMethod("plot", c("DVH", "ANY"),
	function(x, y, ...) {
		plot.DVH(new("DVH.list",list(x)), y, ...)
	}
)

setMethod("plot", c("zDVH", "ANY"),
	function(x, y, ...) {
		plot.DVH(new("DVH.list",list(x)), y, ...)
	}
)

setMethod("plot", c("DVH.list", "missing"),
	function (x, y=NULL, ...) {
		plot.DVH(x, ...)
	}
) 

setMethod("plot", c("DVH.list", "ANY"),
	function (x, y, ...) {
		plot.DVH(x, y, ...)
	}
) 

plot.DVH.ttest <- function(x, y, ..., paired=FALSE, col="black", lty="solid", lwd=1, alpha=0.05, dose=c("absolute", "relative"), dose.units=c("cGy", "Gy"), volume=c("relative", "absolute"), type=c("cumulative", "differential"), main="", line.transparency=1, fill.transparency=line.transparency/2, angle=45, density=NULL, fill.lty=lty, fill=TRUE, legend=c(NA, "topright", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "right", "center"), legend.labels=NULL, new=TRUE, highlight="lightyellow") {
	dose <- match.arg(dose)
	dose.units <- match.arg(dose.units)
	volume <- match.arg(volume)
	type <- match.arg(type)
	legend <- match.arg(legend)
	x <- new("DVH.list", lapply(x, convert.DVH, type=type, dose=dose, dose.units=dose.units, volume=volume))
	x <- x[!unlist(lapply(x, is.empty))]
	N.x <- length(x)
	y <- new("DVH.list", lapply(y, convert.DVH, type=type, dose=dose, dose.units=dose.units, volume=volume))
	y <- y[!unlist(lapply(y, is.empty))]
	N.y <- length(y)
	if (N.x < 2) {
		stop("not enough 'x' observations")
	}
	if (N.y < 2) {
		stop("not enough 'y' observations")
	}
	if (paired & (N.x != N.y)) {
		stop("'x' and 'y' DVH comparison groups have differing lengths -- cannot compute pairwise statistics")
	}
	if (length(col) != 2) {
		col <- rep(col[1], 2)
	}
	col.x <- col2rgb(col[1])/255
	col.y <- col2rgb(col[2])/255	
	if (length(lty) != 2) {
		lty <- rep(lty[1], 2)
	}
	if (length(lwd) != 2) {
		lwd <- rep(lwd[1], 2)
	}
	if (length(fill.transparency) != 2) {
		fill.transparency <- rep(fill.transparency[1], 2)
	}
	if (length(line.transparency) != 2) {
		line.transparency <- rep(line.transparency[1], 2)
	}
	if (length(fill.lty) != 2) {
		fill.lty <- rep(fill.lty[1], 2)
	}
	if (length(angle) != 2) {
		angle <- rep(angle[1], 2)
	}
	if (length(density) != 2) {
		density <- rep(density[1], 2)
	}	
	data.ttest <- t.test(x, y, paired=paired, conf.level=1-alpha)
	doses <- data.ttest$dose
	if (new) {
		layout(c(2,1),heights=c(1,4))
		par(mar=c(4.1, 4.1, 0.5, 2.1))
		plot(NULL, xlim=range(doses), ylim=if (volume == "relative") {c(0, 100)} else {range(unlist(lapply(c(x, y), slot, "volumes")), na.rm=TRUE)}, xlab=if (dose == "absolute") {paste("Dose (", dose.units, ")",sep="")} else {"Dose (%)"}, ylab=if (volume == "relative") {"Volume (%)"} else {"Volume (cc)"}, main="")
	}
	conf.int <- abs(data.ttest$conf.int1 - (data.ttest$x.mean - data.ttest$y.mean))/2
	conf.int[which(is.na(conf.int))] <- 0
	var.x <- var(x)
	var.x <- approx(var.x$dose, var.x$var, doses, yright=0, yleft=0)$y
	var.y <- var(y)
	var.y <- approx(var.y$dose, var.y$var, doses, yright=0, yleft=0)$y
	sum.var <- var.x + var.y
	var.x[which(sum.var != 0)] <- (var.x / sum.var)[which(sum.var != 0)]
	var.y[which(sum.var != 0)] <- (var.y / sum.var)[which(sum.var != 0)]
	y.max <- par("usr")[3]+par("usr")[4]
	x.upper <- pmin(data.ttest$x.mean+conf.int * var.x, y.max)
	x.lower <- pmax(data.ttest$x.mean-conf.int * var.x, 0)
	y.upper <- pmin(data.ttest$y.mean+conf.int * var.y, y.max)
	y.lower <- pmax(data.ttest$y.mean-conf.int * var.y, 0)
	points(doses, x.upper, type="l",col=rgb(col.x[1],col.x[2],col.x[3],fill.transparency[1]), lty=fill.lty[1])
	points(doses, x.lower, type="l",col=rgb(col.x[1],col.x[2],col.x[3],fill.transparency[1]), lty=fill.lty[1])
	points(doses, y.upper, type="l",col=rgb(col.y[1],col.y[2],col.y[3],fill.transparency[2]), lty=fill.lty[2])
	points(doses, y.lower, type="l",col=rgb(col.y[1],col.y[2],col.y[3],fill.transparency[2]), lty=fill.lty[2])
	use.x <- (!(is.na(x.upper) | is.na(x.lower)))
	use.y <- (!(is.na(y.upper) | is.na(y.lower)))
	if (fill) {
		polygon(c(doses[use.x], rev(doses[use.x])), c(x.upper[use.x], rev(x.lower[use.x])), col=rgb(col.x[1],col.x[2],col.x[3],fill.transparency[1]), border=NA, angle=angle[1], density=density[1], lty=fill.lty[1])
   		polygon(c(doses[use.y], rev(doses[use.y])), c(y.upper[use.y], rev(y.lower[use.y])), col=rgb(col.y[1],col.y[2],col.y[3],fill.transparency[2]), border=NA, angle=angle[2], density=density[2],lty=fill.lty[2])
	}
	points(doses, data.ttest$x.mean, type="l", col=rgb(col.x[1],col.x[2],col.x[3],line.transparency[1]), lty=lty[1], lwd=lwd[1])
	points(doses, data.ttest$y.mean, type="l", col=rgb(col.y[1],col.y[2],col.y[3],line.transparency[2]), lty=lty[2], lwd=lwd[2])
	if (!is.na(legend)) {
		if (length(legend.labels) >= 2) { 
			legend(legend, legend=legend.labels[1:2], lty=lty, lwd=lwd, col=c(rgb(col.x[1],col.x[2],col.x[3],line.transparency[1]), rgb(col.y[1],col.y[2],col.y[3],line.transparency[2])), fill=if (fill) {c(rgb(col.x[1],col.x[2],col.x[3],fill.transparency[1]), rgb(col.y[1],col.y[2],col.y[3],fill.transparency[2]))} else {NULL}, density=density, angle=angle)
		}
		else {
			legend(legend, legend=c("Group 1", "Group 2"), lty=lty, lwd=lwd, col=c(rgb(col.x[1],col.x[2],col.x[3],line.transparency[1]), rgb(col.y[1],col.y[2],col.y[3],line.transparency[2])), fill=if (fill) {c(rgb(col.x[1],col.x[2],col.x[3],fill.transparency[1]), rgb(col.y[1],col.y[2],col.y[3],fill.transparency[2]))} else {NULL}, density=density, angle=angle)
		}
	}
	if (new) {
		par(mar=c(0.1,4.1,2.1,2.1))
		p <- data.ttest$p
		plot(range(doses), c(min(p, na.rm=TRUE)/10,1), type="n", xlab="", ylab="P-value", log="y", xaxt="n", yaxt="n", main=main)
		rect(doses[which(p<alpha)],rep(alpha, length(which(p<alpha))), doses[which(p<alpha)], p[which(p<alpha)], col=highlight, border=highlight)
		abline(h=alpha,lty="dotted",col="gray")
		points(doses,p, type="l")
		ticks <- axTicks(2, log=TRUE)	
		ticks <- ticks[which(as.integer(log10(ticks))==log10(ticks))]
		for (i in 1:length(ticks)) {
			j <- log10(ticks[i])
			axis(2, at=ticks[i], labels=substitute(10^j), las=1)
		}
	}
}


plot.DVH.wilcox <- function(x, y, ..., alternative=c("two.sided", "greater", "less"), mu=0, paired=FALSE, exact=TRUE, correct=TRUE, alpha=0.05, col="black", lty="solid", lwd=1, line.transparency=1, fill.transparency=line.transparency/2, angle=45, density=NULL, fill=TRUE, fill.lty=lty, dose=c("absolute", "relative"), dose.units=c("cGy", "Gy"), volume=c("relative", "absolute"), type=c("cumulative", "differential"), main="", legend=c(NA, "topright", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "right", "center"), legend.labels=NULL, new=TRUE, highlight="lightyellow") {
	dose <- match.arg(dose)
	dose.units <- match.arg(dose.units)
	volume <- match.arg(volume)
	type <- match.arg(type)
	alternative <- match.arg(alternative)
	legend <- match.arg(legend)	
	x <- new("DVH.list", lapply(x, convert.DVH, type=type, dose=dose, dose.units=dose.units, volume=volume))
	x <- x[!unlist(lapply(x, is.empty))]
	N.x <- length(x)
	y <- new("DVH.list", lapply(y, convert.DVH, type=type, dose=dose, dose.units=dose.units, volume=volume))
	y <- y[!unlist(lapply(y, is.empty))]
	N.y <- length(y)
	if (N.x < 2) {
		stop("not enough 'x' observations")
	}
	if (N.y < 2) {
		stop("not enough 'y' observations")
	}
	if (length(col) != 2) {
		col <- rep(col[1], 2)
	}
	col.x <- col2rgb(col[1])/255
	col.y <- col2rgb(col[2])/255	
	if (length(lty) != 2) {
		lty <- rep(lty[1], 2)
	}
	if (length(lwd) != 2) {
		lwd <- rep(lwd[1], 2)
	}
	if (length(fill.transparency) != 2) {
		fill.transparency <- rep(fill.transparency[1], 2)
	}
	if (length(line.transparency) != 2) {
		line.transparency <- rep(line.transparency[1], 2)
	}
	if (length(fill.lty) != 2) {
		fill.lty <- rep(fill.lty[1], 2)
	}
	if (length(angle) != 2) {
		angle <- rep(angle[1], 2)
	}
	if (length(density) != 2) {
		density <- rep(density[1], 2)
	}	
	data.wilcox <- wilcox.test(x, y, paired=paired, conf.level=1-alpha, alternative=alternative, mu=mu, exact=exact, correct=correct)
	doses <- data.wilcox$dose
	if (new) {
		layout(c(2,1),heights=c(1,4))
		par(mar=c(4.1, 4.1, 0.5, 2.1))
		plot(NULL, xlim=range(doses), ylim=if (volume == "relative") {c(0, 100)} else {range(unlist(lapply(c(x, y), slot, "volumes")), na.rm=TRUE)}, xlab=if (dose == "absolute") {paste("Dose (", dose.units, ")",sep="")} else {"Dose (%)"}, ylab=if (volume == "relative") {"Volume (%)"} else {"Volume (cc)"}, main="")
	}
	conf.int <- abs(data.wilcox$conf.int2 - data.wilcox$conf.int1)/2
	conf.int[which(is.na(conf.int))] <- 0
	mad.x <- mad(x)
	mad.x <- approx(mad.x$dose, mad.x$mad, doses, yright=0, yleft=0)$y
	mad.y <- mad(y)
	mad.y <- approx(mad.y$dose, mad.y$mad, doses, yright=0, yleft=0)$y
	sum.mad <- mad.x + mad.y
	mad.x[which(sum.mad != 0)] <- (mad.x / sum.mad)[which(sum.mad != 0)]
	mad.y[which(sum.mad != 0)] <- (mad.y / sum.mad)[which(sum.mad != 0)]
	y.max <- par("usr")[3]+par("usr")[4]
	x.upper <- pmin(data.wilcox$x.med+conf.int * mad.x, y.max)
	x.lower <- pmax(data.wilcox$x.med-conf.int * mad.x, 0)
	y.upper <- pmin(data.wilcox$y.med+conf.int * mad.y, y.max)
	y.lower <- pmax(data.wilcox$y.med-conf.int * mad.y, 0)
	points(doses, x.upper, type="l",col=rgb(col.x[1],col.x[2],col.x[3],fill.transparency[1]), lty=fill.lty[1])
	points(doses, x.lower, type="l",col=rgb(col.x[1],col.x[2],col.x[3],fill.transparency[1]), lty=fill.lty[1])
	points(doses, y.upper, type="l",col=rgb(col.y[1],col.y[2],col.y[3],fill.transparency[2]), lty=fill.lty[2])
	points(doses, y.lower, type="l",col=rgb(col.y[1],col.y[2],col.y[3],fill.transparency[2]), lty=fill.lty[2])
	use.x <- (!(is.na(x.upper) | is.na(x.lower)))
	use.y <- (!(is.na(y.upper) | is.na(y.lower)))
	if (fill) {
		polygon(c(doses[use.x], rev(doses[use.x])), c(x.upper[use.x], rev(x.lower[use.x])), col=rgb(col.x[1],col.x[2],col.x[3],fill.transparency[1]), border=NA, angle=angle[1], density=density[1], lty=fill.lty[1])
   		polygon(c(doses[use.y], rev(doses[use.y])), c(y.upper[use.y], rev(y.lower[use.y])), col=rgb(col.y[1],col.y[2],col.y[3],fill.transparency[2]), border=NA, angle=angle[2], density=density[2],lty=fill.lty[2])
	}
	points(doses, data.wilcox$x.med, type="l", col=rgb(col.x[1],col.x[2],col.x[3],line.transparency[1]), lty=lty[1], lwd=lwd[1])
	points(doses, data.wilcox$y.med, type="l", col=rgb(col.y[1],col.y[2],col.y[3],line.transparency[2]), lty=lty[2], lwd=lwd[2])
	if (!is.na(legend)) {
		if (length(legend.labels) >= 2) { 
			legend(legend, legend=legend.labels[1:2], lty=lty, lwd=lwd, col=c(rgb(col.x[1],col.x[2],col.x[3],line.transparency[1]), rgb(col.y[1],col.y[2],col.y[3],line.transparency[2])), fill=if (fill) {c(rgb(col.x[1],col.x[2],col.x[3],fill.transparency[1]), rgb(col.y[1],col.y[2],col.y[3],fill.transparency[2]))} else {NULL}, density=density, angle=angle)
		}
		else {
			legend(legend, legend=c("Group 1", "Group 2"), lty=lty, lwd=lwd, col=c(rgb(col.x[1],col.x[2],col.x[3],line.transparency[1]), rgb(col.y[1],col.y[2],col.y[3],line.transparency[2])), fill=if (fill) {c(rgb(col.x[1],col.x[2],col.x[3],fill.transparency[1]), rgb(col.y[1],col.y[2],col.y[3],fill.transparency[2]))} else {NULL}, density=density, angle=angle)
		}
	}
	if (new) {
		par(mar=c(0.1,4.1,2.1,2.1))
		p <- data.wilcox$p
		plot(range(doses), c(min(p, na.rm=TRUE)/10,1), type="n", xlab="", ylab="P-value", log="y", xaxt="n", yaxt="n", main=main)
		rect(doses[which(p<alpha)],rep(alpha, length(which(p<alpha))), doses[which(p<alpha)], p[which(p<alpha)], col=highlight, border=highlight)
		abline(h=alpha,lty="dotted",col="gray")
		suppressWarnings(points(doses,p, type="l"))
		ticks <- axTicks(2, log=TRUE)	
		ticks <- ticks[which(as.integer(log10(ticks))==log10(ticks))]
		for (i in 1:length(ticks)) {
			j <- log10(ticks[i])
			axis(2, at=ticks[i], labels=substitute(10^j), las=1)
		}
	}	
}




plot.DVH.bars <- function(x, ..., new=TRUE, legend=TRUE, legend.labels=NULL, dose=c("absolute", "relative"), dose.units=c("cGy", "Gy"), volume=c("relative", "absolute"), type=c("cumulative", "differential"), main="", col=rev(rainbow(n=10, start=0, end=2/3))) {
	dose <- match.arg(dose)
	dose.units <- match.arg(dose.units)
	volume <- match.arg(volume)
	x <- c(x, ...)
	x <- new("DVH.list", lapply(x, convert.DVH, type="differential", dose=dose, dose.units=dose.units, volume=volume))
	N <- length(x)
	if (N <= 0) {
		stop("no data to plot DVH")
	}
	doses <- var(x)$dose
	N.dose <- length(doses)	
	colors <- colorRampPalette(col)(N.dose)
	if (volume == "relative") {
		max.vol <- rep(100, N)
		plot(NULL, xlim=c(0,N+1), ylim=c(0,100), xlab="", ylab="Volume (%)", xaxt="n", main=main)
		for (i in 1:N) {
			x.i <- approx(x[[i]]$doses, x[[i]]$volumes, doses, yleft=0, yright=0)$y
			x.i <- diffinv(100*x.i/sum(x.i))
			rect(i-0.35, x.i[1:N.dose], i+0.35, x.i[2:length(x.i)], border=NA, col=colors)
			rect(i-0.35,0,i+0.35,100)
		}
	}
	else {
		max.vol <- max(unlist(lapply(x, slot, "structure.volume")), na.rm=TRUE)
		plot(NULL, xlim=c(0,N+1), ylim=c(0, max.vol), xlab="", ylab="Volume (cc)", xaxt="n", main="")
		for (i in 1:N) {
			vol.i <- x[[i]]$structure.volume
			x.i <- approx(x[[i]]$doses, x[[i]]$volumes, doses, yleft=0, yright=0)$y
			x.i <- diffinv(vol.i*x.i/sum(x.i))
			rect(i-0.35, x.i[1:N.dose], i+0.35, x.i[2:length(x.i)], border=NA, col=colors)
			rect(i-0.35,0,i+0.35,vol.i)
		}
	}
	left <- par("usr")[1]
	right <- par("usr")[2]
	dist <- right - left
	rect(left + dist*(1:N.dose - 1)/N.dose, par("usr")[4], left+dist*(1:N.dose/N.dose), par("usr")[4]+max.vol/25, border=NA, col=colors, xpd=TRUE)
	if (legend) {
		if (length(legend.labels) != N) {
			legend.labels <- unlist(lapply(x, slot, "structure.name"))
		}
		text(1:N, par("usr")[3], labels=paste(substr(legend.labels, 1, 12), " ", sep=""), srt=45, xpd=TRUE, adj=1)
	}
	ax.labs <- pretty(doses)
	at.labs <- left +(dist)*ax.labs/max(doses, na.rm=TRUE)
	axis(3, at=at.labs, labels=ax.labs, tick=FALSE)
	for (i in 1:length(ax.labs)) {
		rect(c(left,at.labs)[i], par("usr")[4], min(at.labs[i],right), par("usr")[4]+max.vol/25, xpd=TRUE)	
	}
	mtext(if (dose == "relative") {"Dose (%)"} else {paste("Dose (", dose.units, ")",sep="")}, at=(left + right)/2, line=2, xpd=TRUE)
	mtext(main, at=left, line=2, xpd=TRUE, font=2)
}


plot.DVH.individual <- function(x, ..., col="black", lty ="solid", lwd=1, line.transparency=1, new=TRUE, legend=c(NA, "topright", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "right", "center"), legend.labels=NULL, dose=c("absolute", "relative"), dose.units=c("cGy", "Gy"), volume=c("relative", "absolute"), type=c("cumulative", "differential"), main="") {
	dose <- match.arg(dose)
	dose.units <- match.arg(dose.units)
	volume <- match.arg(volume)
	type <- match.arg(type)
	legend <- match.arg(legend)
	x <- c(x, ...)
	x <- new("DVH.list", lapply(x, convert.DVH, type=type, dose=dose, dose.units=dose.units, volume=volume))
	N <- length(x)
	if (N <= 0) {
		stop("no data to plot DVH")
	}
	if (length(col) != N) {
		col <- rep(col[1], N)
	}
	if (length(lty) != N) {
		lty <- rep(lty[1], N)
	}
	if (length(lwd) != N) {
		lwd <- rep(lwd[1], N)
	}
	if (length(line.transparency) != N) {
		line.transparency <- rep(line.transparency[1], N)
	}
	if (new) {	
		plot(NULL, xlim=range(x), ylim=range(unlist(lapply(x, slot, "volumes"))), xlab=if (dose == "absolute") {paste("Dose (", dose.units, ")",sep="")} else {"Dose (%)"}, ylab=if (volume == "relative") {"Volume (%)"} else {"Volume (cc)"}, main=main)
	}
	for (i in 1:N) {
		col.i <- col2rgb(col[i])/255
		points(x[[i]]$doses, x[[i]]$volumes, type="l", lty=lty[i], lwd=lwd[i], col=rgb(col.i[1],col.i[2],col.i[3],line.transparency[i]))
	}
	if (!is.na(legend)) {		 
		legend(legend, legend=if (length(legend.labels) >= N) {legend.labels[1:N]} else {paste("Structure", 1:N)}, lty=lty, lwd=lwd, col=col)
	}
}

plot.DVH.groups <- function(x, ..., col="black", lty ="solid", lwd=1, line.transparency=1, fill.transparency=line.transparency/2, new=TRUE, legend=c(NA, "topright", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "right", "center"), legend.labels=NULL, dose=c("absolute", "relative"), dose.units=c("cGy", "Gy"), volume=c("relative", "absolute"), type=c("cumulative", "differential"), width=c("range", "mad", "IQR", "quantile", "var", "sd"), main="", multiplier=1, quantile=c(0.25, 0.75), fill=TRUE, angle=45, density=NULL, fill.lty=lty) {
	dose <- match.arg(dose)
	dose.units <- match.arg(dose.units)
	volume <- match.arg(volume)
	type <- match.arg(type)
	width <- match.arg(width)
	multiplier <- max(0, multiplier, na.rm=TRUE)
	legend <- match.arg(legend)
	groups <- c(list(x), list(...))
	classes <- unlist(lapply(groups, class))
	which.DVH <- which(classes %in% c("DVH", "DVH.list"))
	classes <- classes[which.DVH]
	groups <- groups[which.DVH]
	N <- length(classes)
	if (length(col) != N) {
		col <- rep(col[1], N)
	}
	if (length(lty) != N) {
		lty <- rep(lty[1], N)
	}
	if (length(lwd) != N) {
		lwd <- rep(lwd[1], N)
	}
	if (length(fill.transparency) != N) {
		fill.transparency <- rep(fill.transparency[1], N)
	}
	if (length(line.transparency) != N) {
		line.transparency <- rep(line.transparency[1], N)
	}
	if (length(fill.lty) != N) {
		fill.lty <- rep(fill.lty[1], N)
	}
	if (length(angle) != N) {
		angle <- rep(angle[1], N)
	}
	if (length(density) != N) {
		density <- rep(density[1], N)
	}	
	range.dose <- NA
	range.volume <- NA
	for (i in 1:N) {
		if (classes[i] == "DVH.list") {
			groups[[i]] <- new("DVH.list", lapply(groups[[i]], convert.DVH, type=type, dose=dose, dose.units=dose.units, volume=volume))
			range.dose <- range(range.dose, range(groups[[i]]), na.rm=TRUE)
			range.volume <- range(range.volume, range(unlist(lapply(groups[[i]], slot, "volumes"))), na.rm=TRUE)
		}
		else {
			groups[[i]] <- convert.DVH(groups[[i]], type=type, dose=dose, dose.units=dose.units, volume=volume)
			range.dose <- range(range.dose, range(groups[[i]]), na.rm=TRUE)
			range.volume <- range(range.volume, range(slot(groups[[i]], "volumes")), na.rm=TRUE)
		}
	}
	if (new) {
		plot(NULL, xlim=range.dose, ylim=range.volume, xlab=if (dose == "absolute") {paste("Dose (", dose.units, ")",sep="")} else {"Dose (%)"}, ylab=if (volume == "relative") {"Volume (%)"} else {"Volume (cc)"}, main=main)
	}
	y.max <- par("usr")[3]+par("usr")[4]
	for (i in 1:N) {
		col.i <- col2rgb(col[i])/255
		if (classes[i] == "DVH.list") {
			if (fill) {
				switch(width,
					range = {
						DVH.center <- mean(groups[[i]], type=type, dose=dose, volume=volume)
						DVH.range <- quantile(groups[[i]], probs=c(0, 1), type=7, na.rm=TRUE)
						polygon(c(DVH.range$dose, rev(DVH.range$dose)), c(pmax(DVH.range$quantiles[1,], 0), rev(pmin(DVH.range$quantiles[2,], y.max))), col=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), border=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), angle=angle[i], density=density[i], lty=fill.lty[i])
					},
					mad = {
						DVH.center <- median(groups[[i]])
						DVH.range <- mad(groups[[i]])
						polygon(c(DVH.range$dose, rev(DVH.range$dose)), c(pmax(DVH.center$volumes-DVH.range$mad*multiplier, 0), rev(pmin(DVH.center$volumes+DVH.range$mad*multiplier, y.max))), col=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), border=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), angle=angle[i], density=density[i], lty=fill.lty[i])
					},
					IQR = {
						DVH.center <- median(groups[[i]])
						DVH.range <- quantile(groups[[i]], probs=c(0.25, 0.75), type=7, na.rm=TRUE)
						polygon(c(DVH.range$dose, rev(DVH.range$dose)), c(pmax(DVH.range$quantiles[1,], 0), rev(pmin(DVH.range$quantiles[2,], y.max))), col=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), border=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), angle=angle[i], density=density[i], lty=fill.lty[i])
					},
					quantile = {
						if (length(quantile) != 2) {
							quantile <- c(0.25, 0.75)
						}
						DVH.center <- quantile(groups[[i]], probs=mean(quantile), type=7, na.rm=TRUE)
						names(DVH.center) <- c("doses", "volumes")
						DVH.range <- quantile(groups[[i]], probs=quantile, type=7, na.rm=TRUE)
						polygon(c(DVH.range$dose, rev(DVH.range$dose)), c(pmax(DVH.range$quantiles[1,], 0), rev(pmin(DVH.range$quantiles[2,], y.max))), col=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), border=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), angle=angle[i], density=density[i], lty=fill.lty[i])
					},
					sd = {
						DVH.center <- mean(groups[[i]], type=type, dose=dose, volume=volume)
						DVH.range <- sd(groups[[i]])
						polygon(c(DVH.range$dose, rev(DVH.range$dose)), c(pmax(DVH.center$volumes-DVH.range$sd*multiplier, 0), rev(pmin(DVH.center$volumes+DVH.range$sd*multiplier, y.max))), col=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), border=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), angle=angle[i], density=density[i], lty=fill.lty[i])
					},
					var = {
						DVH.center <- mean(groups[[i]], type=type, dose=dose, volume=volume)
						DVH.range <- var(groups[[i]])
						polygon(c(DVH.range$dose, rev(DVH.range$dose)), c(pmax(DVH.center$volumes-DVH.range$var*multiplier, 0), rev(pmin(DVH.center$volumes+DVH.range$var*multiplier, y.max))), col=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), border=rgb(col.i[1],col.i[2],col.i[3],fill.transparency[i]), angle=angle[i], density=density[i], lty=fill.lty[i])
					}
				)
			}
			points(DVH.center$doses, DVH.center$volumes, type="l", col=rgb(col.i[1],col.i[2],col.i[3],line.transparency[i]), lty=lty[i], lwd=lwd[i])
		}
		else {
			points(groups[[i]]$doses, groups[[i]]$volumes, type="l", lty=lty[i], lwd=lwd[i], col=rgb(col.i[1],col.i[2],col.i[3],line.transparency[i]))
		}		
	}	
	if (!is.na(legend)) {		 
		legend(legend, legend=if (length(legend.labels) >= N) {legend.labels[1:N]} else {paste("Group", 1:N)}, lty=lty, lwd=lwd, col=col)
	}
}


plot.zDVH <- function(x, ..., col="black", front=NULL, back=front, new=TRUE, dose=NULL, dose.units=NULL, volume=NULL, type=NULL, main="") {
	dose.units <- match.arg(dose.units, choices=c(NA, "cGy", "Gy"))
	type <- match.arg(type, choices=c(NA, "cumulative", "differential"))
	volume <- match.arg(volume, choices=c(NA, "relative", "absolute"))
	dose <- match.arg(dose, choices=c(NA, "absolute", "relative"))
	front <- match.arg(front, choices=c("filled", "lines", "points", "culled"))
	back <- match.arg(back, choices=c("filled", "lines", "points", "culled"))
	x <- convert.DVH(x, type=type, dose=dose, dose.units=dose.units, volume=volume)
	if (length(unique(diff(x$doses))) > 1) {
		persp3d(x$doses[2:length(x$doses)], as.numeric(colnames(x$volumes)), x$volumes[2:length(x$doses),], col=col, border=NA, shade=0.5, xlab=paste("Dose (", x$dose.units, ")", sep=""), ylab="z (mm)", zlab=paste("Volume (", if (x$volume.type == "relative") {"%"} else {"cc"}, ")", sep=""), add=!new)
	}
	persp3d(x$doses, as.numeric(colnames(x$volumes)), x$volumes, col=col, border=NA, shade=0.5, xlab=paste("Dose (", x$dose.units, ")", sep=""), ylab="z (mm)", zlab=paste("Volume (", if (x$volume.type == "relative") {"%"} else {"cc"}, ")", sep=""), add=!new, front=front, back=back)
}