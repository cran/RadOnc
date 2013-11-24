### R code from vignette source 'RadOnc.Rnw'

###################################################
### code chunk number 1: RadOnc.Rnw:49-50
###################################################
library(RadOnc)


###################################################
### code chunk number 2: RadOnc.Rnw:52-55
###################################################
data.path <- system.file("extdata", package="RadOnc")
initial.path <- getwd()
options(width=75)


###################################################
### code chunk number 3: RadOnc.Rnw:85-89
###################################################
temp <- c(readLines(paste(data.path, 'Jane_Doe.dvh', sep='/'),n=50), "...", "...")
for (i in 1:52) {
cat(temp[i], "\n")
}


###################################################
### code chunk number 4: RadOnc.Rnw:92-93
###################################################
johndoe <- read.DVH(file=system.file("extdata/John_Doe.dvh", package="RadOnc"), type="aria10", verbose=TRUE)


###################################################
### code chunk number 5: RadOnc.Rnw:96-97 (eval = FALSE)
###################################################
## read.DVH(file="Jane_Doe.dvh", type="aria10", verbose=TRUE)


###################################################
### code chunk number 6: RadOnc.Rnw:99-100
###################################################
janedoe <- read.DVH(file=system.file("extdata/Jane_Doe.dvh", package="RadOnc"), type="aria10", verbose=TRUE)


###################################################
### code chunk number 7: RadOnc.Rnw:107-109
###################################################
janedoe[1:4]
c(janedoe[c("PTV")], johndoe[c("CTV", "DUODENUM")])


###################################################
### code chunk number 8: RadOnc.Rnw:111-112
###################################################
temp <- janedoe


###################################################
### code chunk number 9: RadOnc.Rnw:114-117
###################################################
johndoe[["CTV"]]
janedoe[[1]] <- johndoe[["CTV"]]
janedoe[1:4]


###################################################
### code chunk number 10: RadOnc.Rnw:119-120
###################################################
janedoe <- temp


###################################################
### code chunk number 11: RadOnc.Rnw:126-130
###################################################
names(janedoe)[1:4] <- c("A1", "B2", "C3", "D4")
names(rev(janedoe[1:4]))
length(johndoe)
lapply(johndoe, function(DVH) { DVH[c("DMIN", "D50%", "DMAX", "V20%")] })


###################################################
### code chunk number 12: RadOnc.Rnw:132-133
###################################################
janedoe <- temp


###################################################
### code chunk number 13: RadOnc.Rnw:140-144
###################################################
johndoe[["DUODENUM"]]["V20Gy"]
johndoe[["DUODENUM"]]["D2.5%"]
johndoe[["DUODENUM"]]["volume"] * 0.025
johndoe[["DUODENUM"]]["D2.3286cc"]


###################################################
### code chunk number 14: RadOnc.Rnw:147-149
###################################################
johndoe[["DUODENUM"]][c("V5%", "V20Gy", "D2.5%", "D2.3286cc", "Dmax")]
johndoe[1:4]$"V20Gy,Dmax"


###################################################
### code chunk number 15: RadOnc.Rnw:152-153
###################################################
johndoe[["DUODENUM"]][c("V5", "VGy", "volume", 2.5, "", "Dmax%")]


###################################################
### code chunk number 16: fig1
###################################################
plot(janedoe[[3]], volume="relative", dose="absolute", type="cumulative")


###################################################
### code chunk number 17: fig2
###################################################
plot(janedoe[1:3], plot.type="i", col=c("red", "green", "blue"), 
legend="topright", legend.labels=names(janedoe[1:3]))


###################################################
### code chunk number 18: fig3
###################################################
plot(c(johndoe["STOMACH"],janedoe["STOMACH"]), #group 1
c(janedoe["LIVER"],johndoe["LIVER"]), #group 2
c(johndoe["DUODENUM"],janedoe["DUODENUM"]), #group 3
plot.type="g", dose="relative", col=c("blue", "red", "green"), 
lwd=2, lty="dashed", fill.lty="solid", fill.transparency=0.3)


###################################################
### code chunk number 19: fig4
###################################################
group1 <- c("CTV", "PTV")
group2 <- c("LIVER", "STOMACH", "SMALL_BOWEL")
plot(c(johndoe[group1],janedoe[group1]), 
c(janedoe[group2],johndoe[group2]),
plot.type="t", main="Target v. OAR t-Test", alpha=0.001, 
col=c("red", "blue"), lty="dashed", fill.lty="solid")


###################################################
### code chunk number 20: fig5
###################################################
plot(janedoe[2:9], plot.type="b", volume="abs", dose="rel")


###################################################
### code chunk number 21: fig6
###################################################
plot(janedoe)
plot(median(janedoe), new=FALSE, col="red", lwd=2)
plot(mean(janedoe), new=FALSE, col="blue", lwd=2, lty="dashed")


###################################################
### code chunk number 22: RadOnc.Rnw:243-246
###################################################
groupA <- janedoe[c("LIVER","LEFT_KIDNEY","RIGHT_KIDNEY","CORD")]
groupB <- janedoe[c("CTV", "PTV")]
t.test(unlist(groupA$"V20Gy"), unlist(groupB$"V20Gy"))


###################################################
### code chunk number 23: RadOnc.Rnw:251-253 (eval = FALSE)
###################################################
## AvB <- t.test(groupA, groupB)
## plot(AvB$dose, AvB$p, type="l", log="y", xlab="Dose (cGy)", ylab="p-value")


###################################################
### code chunk number 24: fig7
###################################################
AvB <- t.test(groupA, groupB)
plot(AvB$dose, AvB$p, type="l", log="y", xlab="Dose (cGy)", ylab="p-value")
abline(v=2000,col="gray", lty="dashed")
points(2000,approx(AvB$dose, AvB$p, 2000)$y, col="red")
text(2000,approx(AvB$dose, AvB$p, 2000)$y, col="red", labels="V20Gy (p=5.347e-05)",pos=4)


###################################################
### code chunk number 25: RadOnc.Rnw:278-279 (eval = FALSE)
###################################################
## data <- read.DICOM.RT(path="<<DICOM directory>>", verbose=TRUE)


###################################################
### code chunk number 26: RadOnc.Rnw:283-284
###################################################
data("RadOnc")


###################################################
### code chunk number 27: RadOnc.Rnw:291-293
###################################################
teeth[1:2]
c(cord, mandible)


###################################################
### code chunk number 28: RadOnc.Rnw:295-296
###################################################
temp <- teeth


###################################################
### code chunk number 29: RadOnc.Rnw:298-301
###################################################
teeth[[1]]
teeth[[1]] <- teeth[["Tooth #3"]]
teeth


###################################################
### code chunk number 30: RadOnc.Rnw:303-304
###################################################
teeth <- temp


###################################################
### code chunk number 31: RadOnc.Rnw:310-314
###################################################
names(teeth) <- c("Larry", "Curly", "Moe")
names(rev(teeth[1:3]))
length(teeth)
lapply(teeth, function(tooth) { range(tooth) })


###################################################
### code chunk number 32: RadOnc.Rnw:316-317
###################################################
teeth <- temp


###################################################
### code chunk number 33: fig8 (eval = FALSE)
###################################################
## plot(mandible)


###################################################
### code chunk number 34: fig9 (eval = FALSE)
###################################################
## plot(cord)


###################################################
### code chunk number 35: RadOnc.Rnw:355-356 (eval = FALSE)
###################################################
## compareStructures(teeth, method="grid", plot=TRUE)


###################################################
### code chunk number 36: fig10
###################################################
compareStructures(teeth, method="grid", plot=TRUE, pixels=40)


###################################################
### code chunk number 37: RadOnc.Rnw:368-369
###################################################
teeth <- teeth[c(1,3)]


###################################################
### code chunk number 38: RadOnc.Rnw:371-372
###################################################
compareStructures(teeth, method="hausdorff", hausdorff.method="mean")


###################################################
### code chunk number 39: RadOnc.Rnw:374-375
###################################################
teeth <- temp


