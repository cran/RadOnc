### R code from vignette source 'RadOnc.Rnw'

###################################################
### code chunk number 1: RadOnc.Rnw:50-51
###################################################
library(RadOnc)


###################################################
### code chunk number 2: RadOnc.Rnw:53-56
###################################################
data.path <- system.file("extdata", package="RadOnc")
initial.path <- getwd()
options(width=75)


###################################################
### code chunk number 3: RadOnc.Rnw:86-90
###################################################
temp <- c(readLines(paste(data.path, 'Jane_Doe.dvh', sep='/'),n=50), "...", "...")
for (i in 1:52) {
cat(temp[i], "\n")
}


###################################################
### code chunk number 4: RadOnc.Rnw:93-94
###################################################
johndoe <- read.DVH(file=system.file("extdata/John_Doe.dvh", package="RadOnc"), type="aria10", verbose=TRUE)


###################################################
### code chunk number 5: RadOnc.Rnw:97-98 (eval = FALSE)
###################################################
## read.DVH(file="Jane_Doe.dvh", type="aria10", verbose=TRUE)


###################################################
### code chunk number 6: RadOnc.Rnw:100-101
###################################################
janedoe <- read.DVH(file=system.file("extdata/Jane_Doe.dvh", package="RadOnc"), type="aria10", verbose=TRUE)


###################################################
### code chunk number 7: RadOnc.Rnw:108-110
###################################################
janedoe[1:4]
c(janedoe[c("PTV")], johndoe[c("CTV", "DUODENUM")])


###################################################
### code chunk number 8: RadOnc.Rnw:112-113
###################################################
temp <- janedoe


###################################################
### code chunk number 9: RadOnc.Rnw:115-118
###################################################
johndoe[["CTV"]]
janedoe[[1]] <- johndoe[["CTV"]]
janedoe[1:4]


###################################################
### code chunk number 10: RadOnc.Rnw:120-121
###################################################
janedoe <- temp


###################################################
### code chunk number 11: RadOnc.Rnw:127-131
###################################################
names(janedoe)[1:4] <- c("A1", "B2", "C3", "D4")
names(rev(janedoe[1:4]))
length(johndoe)
lapply(johndoe, function(DVH) { DVH[c("DMIN", "D50%", "DMAX", "V20%")] })


###################################################
### code chunk number 12: RadOnc.Rnw:133-134
###################################################
janedoe <- temp


###################################################
### code chunk number 13: RadOnc.Rnw:141-145
###################################################
johndoe[["DUODENUM"]]["V20Gy"]
johndoe[["DUODENUM"]]["D2.5%"]
johndoe[["DUODENUM"]]["volume"] * 0.025
johndoe[["DUODENUM"]]["D2.3286cc"]


###################################################
### code chunk number 14: RadOnc.Rnw:148-150
###################################################
johndoe[["DUODENUM"]][c("V5%", "V20Gy", "D2.5%", "D2.3286cc", "Dmax")]
johndoe[1:4]$"V20Gy,Dmax"


###################################################
### code chunk number 15: RadOnc.Rnw:153-154
###################################################
johndoe[["DUODENUM"]][c("V5", "VGy", "volume", 2.5, "", "Dmax%")]


###################################################
### code chunk number 16: RadOnc.Rnw:157-158
###################################################
johndoe[["LIVER"]][c("V10Gy(%)","D25%","D25%(Gy)")]


###################################################
### code chunk number 17: RadOnc.Rnw:161-163
###################################################
johndoe[["LIVER"]][c("Dintegral","Dintegral(>0cGy)")]
johndoe[["LIVER"]][c("Dintegral(<20Gy)","Dintegral(10-20Gy)")]


###################################################
### code chunk number 18: RadOnc.Rnw:171-175
###################################################
gEUD(janedoe[1:3], 6:8)
gEUD(janedoe[1:3], 1) == unlist(janedoe[1:3]$"Dmean")
gEUD(janedoe[1:3], Inf) == unlist(janedoe[1:3]$"Dmax")
gEUD(janedoe[1:3], -Inf) == unlist(janedoe[1:3]$"Dmin")


###################################################
### code chunk number 19: RadOnc.Rnw:182-183
###################################################
LQE(c(4500, 5500, 6000), aB=3, fractions=c(20, 30))


###################################################
### code chunk number 20: fig1
###################################################
plot(janedoe[[3]], volume="relative", dose="absolute", type="cumulative")


###################################################
### code chunk number 21: fig2
###################################################
plot(janedoe[1:3], plot.type="i", col=c("red", "green", "blue"), 
legend="topright", legend.labels=names(janedoe[1:3]))


###################################################
### code chunk number 22: fig3
###################################################
plot(c(johndoe["STOMACH"],janedoe["STOMACH"]), #group 1
c(janedoe["LIVER"],johndoe["LIVER"]), #group 2
c(johndoe["DUODENUM"],janedoe["DUODENUM"]), #group 3
plot.type="g", dose="relative", col=c("blue", "red", "green"), 
lwd=2, lty="dashed", fill.lty="solid", fill.transparency=0.3)


###################################################
### code chunk number 23: fig4
###################################################
group1 <- c("CTV", "PTV")
group2 <- c("LIVER", "STOMACH", "SMALL_BOWEL")
plot(c(johndoe[group1],janedoe[group1]), 
c(janedoe[group2],johndoe[group2]),
plot.type="t", main="Target v. OAR t-Test", alpha=0.001, 
col=c("red", "blue"), lty="dashed", fill.lty="solid")


###################################################
### code chunk number 24: fig5
###################################################
plot(janedoe[2:9], plot.type="b", volume="abs", dose="rel")


###################################################
### code chunk number 25: fig6
###################################################
plot(janedoe)
plot(median(janedoe), new=FALSE, col="red", lwd=2)
plot(mean(janedoe), new=FALSE, col="blue", lwd=2, lty="dashed")


###################################################
### code chunk number 26: RadOnc.Rnw:274-277
###################################################
L.kidney <- janedoe[["LEFT_KIDNEY"]]
R.kidney <- janedoe[["RIGHT_KIDNEY"]]
total.kidney <- sum(janedoe[c("LEFT_KIDNEY", "RIGHT_KIDNEY")])


###################################################
### code chunk number 27: RadOnc.Rnw:279-285 (eval = FALSE)
###################################################
## L.kidney <- janedoe[["LEFT_KIDNEY"]]
## R.kidney <- janedoe[["RIGHT_KIDNEY"]]
## total.kidney <- sum(janedoe[c("LEFT_KIDNEY", "RIGHT_KIDNEY")])
## plot(total.kidney, type="diff", volume="abs")
## plot(L.kidney, new=FALSE, type="diff", volume="abs", col="red")
## plot(R.kidney, new=FALSE, type="diff", volume="abs", col="blue")


###################################################
### code chunk number 28: fig7
###################################################
plot(total.kidney, type="diff", volume="abs")
plot(L.kidney, new=FALSE, type="diff", volume="abs", col="red")
plot(R.kidney, new=FALSE, type="diff", volume="abs", col="blue")


###################################################
### code chunk number 29: RadOnc.Rnw:301-304
###################################################
groupA <- janedoe[c("LIVER","LEFT_KIDNEY","RIGHT_KIDNEY","CORD")]
groupB <- janedoe[c("CTV", "PTV")]
t.test(unlist(groupA$"V20Gy"), unlist(groupB$"V20Gy"))


###################################################
### code chunk number 30: RadOnc.Rnw:309-311 (eval = FALSE)
###################################################
## AvB <- t.test(groupA, groupB)
## plot(AvB$dose, AvB$p, type="l", log="y", xlab="Dose (cGy)", ylab="p-value")


###################################################
### code chunk number 31: fig8
###################################################
AvB <- t.test(groupA, groupB)
plot(AvB$dose, AvB$p, type="l", log="y", xlab="Dose (cGy)", ylab="p-value")
abline(v=2000,col="gray", lty="dashed")
points(2000,approx(AvB$dose, AvB$p, 2000)$y, col="red")
text(2000,approx(AvB$dose, AvB$p, 2000)$y, col="red", labels="V20Gy (p=5.347e-05)",pos=4)


###################################################
### code chunk number 32: RadOnc.Rnw:336-337 (eval = FALSE)
###################################################
## data <- read.DICOM.RT(path="<<DICOM directory>>", verbose=TRUE)


###################################################
### code chunk number 33: RadOnc.Rnw:341-342
###################################################
data("RadOnc")


###################################################
### code chunk number 34: RadOnc.Rnw:349-351
###################################################
teeth[1:2]
c(cord, mandible)


###################################################
### code chunk number 35: RadOnc.Rnw:353-354
###################################################
temp <- teeth


###################################################
### code chunk number 36: RadOnc.Rnw:356-359
###################################################
teeth[[1]]
teeth[[1]] <- teeth[["Tooth #3"]]
teeth


###################################################
### code chunk number 37: RadOnc.Rnw:361-362
###################################################
teeth <- temp


###################################################
### code chunk number 38: RadOnc.Rnw:368-372
###################################################
names(teeth) <- c("Larry", "Curly", "Moe")
names(rev(teeth[1:3]))
length(teeth)
lapply(teeth, function(tooth) { range(tooth) })


###################################################
### code chunk number 39: RadOnc.Rnw:374-375
###################################################
teeth <- temp


###################################################
### code chunk number 40: fig9 (eval = FALSE)
###################################################
## plot(mandible)


###################################################
### code chunk number 41: fig10 (eval = FALSE)
###################################################
## plot(cord)


###################################################
### code chunk number 42: RadOnc.Rnw:413-414 (eval = FALSE)
###################################################
## compareStructures(teeth, method="grid", plot=TRUE)


###################################################
### code chunk number 43: fig11
###################################################
compareStructures(teeth, method="grid", plot=TRUE, pixels=40)


###################################################
### code chunk number 44: RadOnc.Rnw:426-427
###################################################
teeth <- teeth[c(1,3)]


###################################################
### code chunk number 45: RadOnc.Rnw:429-430
###################################################
compareStructures(teeth, method="hausdorff", hausdorff.method="mean")


###################################################
### code chunk number 46: RadOnc.Rnw:432-433
###################################################
teeth <- temp


###################################################
### code chunk number 47: RadOnc.Rnw:457-458
###################################################
news(package="RadOnc")


