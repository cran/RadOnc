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
### code chunk number 7: RadOnc.Rnw:103-104 (eval = FALSE)
###################################################
## DVHs <- read.DVH(file=c("Jane_Doe.dvh", "John_Doe.dvh"), type="aria10")


###################################################
### code chunk number 8: RadOnc.Rnw:106-108
###################################################
DVHs <- read.DVH(file=system.file(paste("extdata", c("Jane_Doe.dvh", "John_Doe.dvh"), sep="/"), package="RadOnc"), type="aria10")
DVHs


###################################################
### code chunk number 9: fig1
###################################################
plot(as(janedoe.RTdata$structures,"DVH.list"),lwd=2.5)
plot(janedoe[c(3,6:7)],new=FALSE,col="red",lwd=1.25)


###################################################
### code chunk number 10: RadOnc.Rnw:132-134
###################################################
janedoe[1:4]
c(janedoe["PTV"], johndoe[c("CTV", "DUODENUM")])


###################################################
### code chunk number 11: RadOnc.Rnw:136-137
###################################################
temp <- janedoe


###################################################
### code chunk number 12: RadOnc.Rnw:139-142
###################################################
johndoe[["CTV"]]
janedoe[[1]] <- johndoe[["CTV"]]
janedoe[1:4]


###################################################
### code chunk number 13: RadOnc.Rnw:144-145
###################################################
janedoe <- temp


###################################################
### code chunk number 14: RadOnc.Rnw:149-151
###################################################
janedoe["KIDNEY$"]
janedoe[c(2,"IGHT.*")]


###################################################
### code chunk number 15: RadOnc.Rnw:156-160
###################################################
names(janedoe)[1:4] <- c("A1", "B2", "C3", "D4")
names(rev(janedoe[1:4]))
length(johndoe)
lapply(johndoe, function(DVH) { DVH[c("DMIN", "D50%", "DMAX", "V20%")] })


###################################################
### code chunk number 16: RadOnc.Rnw:162-163
###################################################
janedoe <- temp


###################################################
### code chunk number 17: RadOnc.Rnw:167-169
###################################################
janedoe[1:2]$patients
janedoe[3:4]$ID


###################################################
### code chunk number 18: RadOnc.Rnw:176-180
###################################################
johndoe[["DUODENUM"]]["V20Gy"]
johndoe[["DUODENUM"]]["D2.5%"]
johndoe[["DUODENUM"]]["volume"] * 0.025
johndoe[["DUODENUM"]]["D2.3286cc"]


###################################################
### code chunk number 19: RadOnc.Rnw:183-185
###################################################
johndoe[["DUODENUM"]][c("V5%", "V20Gy", "D2.5%", "D2.3286cc", "Dmax")]
johndoe[1:4]$"V20Gy,Dmax"


###################################################
### code chunk number 20: RadOnc.Rnw:188-189
###################################################
johndoe[["DUODENUM"]][c("V5", "VGy", "volume", 2.5, "", "Dmax%")]


###################################################
### code chunk number 21: RadOnc.Rnw:192-193
###################################################
johndoe[["LIVER"]][c("V10Gy(%)","D25%","D25%(Gy)")]


###################################################
### code chunk number 22: RadOnc.Rnw:196-198
###################################################
johndoe[["LIVER"]][c("Dintegral","Dintegral(>0cGy)")]
johndoe[["LIVER"]][c("Dintegral(<20Gy)","Dintegral(10-20Gy)")]


###################################################
### code chunk number 23: RadOnc.Rnw:206-210
###################################################
gEUD(janedoe[1:3], 6:8)
gEUD(janedoe[1:3], 1) == unlist(janedoe[1:3]$"Dmean")
gEUD(janedoe[1:3], Inf) == unlist(janedoe[1:3]$"Dmax")
gEUD(janedoe[1:3], -Inf) == unlist(janedoe[1:3]$"Dmin")


###################################################
### code chunk number 24: RadOnc.Rnw:217-218
###################################################
LQE(c(4500, 5500, 6000), aB=3, fractions=c(20, 30))


###################################################
### code chunk number 25: fig2
###################################################
plot(janedoe[[3]], volume="relative", dose="absolute", type="cumulative")


###################################################
### code chunk number 26: fig3
###################################################
plot(janedoe[1:3], plot.type="i", col=c("red", "green", "blue"), 
legend="topright", legend.labels=names(janedoe[1:3]))


###################################################
### code chunk number 27: fig4
###################################################
plot(c(johndoe["STOMACH"],janedoe["STOMACH"]), #group 1
c(janedoe["LIVER"],johndoe["LIVER"]), #group 2
c(johndoe["DUODENUM"],janedoe["DUODENUM"]), #group 3
plot.type="g", dose="relative", col=c("blue", "red", "green"), 
lwd=2, lty="dashed", fill.lty="solid", fill.transparency=0.3)


###################################################
### code chunk number 28: fig5
###################################################
group1 <- c("CTV", "PTV")
group2 <- c("LIVER", "STOMACH", "SMALL_BOWEL")
plot(c(johndoe[group1],janedoe[group1]), 
c(janedoe[group2],johndoe[group2]),
plot.type="t", main="Target v. OAR t-Test", alpha=0.001, 
col=c("red", "blue"), lty="dashed", fill.lty="solid")


###################################################
### code chunk number 29: fig6
###################################################
plot(janedoe[2:9], plot.type="b", volume="abs", dose="rel")


###################################################
### code chunk number 30: fig7
###################################################
plot(janedoe)
plot(median(janedoe), new=FALSE, col="red", lwd=2)
plot(mean(janedoe), new=FALSE, col="blue", lwd=2, lty="dashed")


###################################################
### code chunk number 31: RadOnc.Rnw:309-312
###################################################
L.kidney <- janedoe[["LEFT_KIDNEY"]]
R.kidney <- janedoe[["RIGHT_KIDNEY"]]
total.kidney <- sum(janedoe[c("LEFT_KIDNEY", "RIGHT_KIDNEY")])


###################################################
### code chunk number 32: RadOnc.Rnw:314-320 (eval = FALSE)
###################################################
## L.kidney <- janedoe[["LEFT_KIDNEY"]]
## R.kidney <- janedoe[["RIGHT_KIDNEY"]]
## total.kidney <- sum(janedoe[c("LEFT_KIDNEY", "RIGHT_KIDNEY")])
## plot(total.kidney, type="diff", volume="abs")
## plot(L.kidney, new=FALSE, type="diff", volume="abs", col="red")
## plot(R.kidney, new=FALSE, type="diff", volume="abs", col="blue")


###################################################
### code chunk number 33: fig8
###################################################
plot(total.kidney, type="diff", volume="abs")
plot(L.kidney, new=FALSE, type="diff", volume="abs", col="red")
plot(R.kidney, new=FALSE, type="diff", volume="abs", col="blue")


###################################################
### code chunk number 34: RadOnc.Rnw:336-339
###################################################
groupA <- janedoe[c("LIVER","LEFT_KIDNEY","RIGHT_KIDNEY","CORD")]
groupB <- janedoe[c("CTV", "PTV")]
t.test(unlist(groupA$"V20Gy"), unlist(groupB$"V20Gy"))


###################################################
### code chunk number 35: RadOnc.Rnw:344-346 (eval = FALSE)
###################################################
## AvB <- t.test(groupA, groupB)
## plot(AvB$dose, AvB$p, type="l", log="y", xlab="Dose (cGy)", ylab="p-value")


###################################################
### code chunk number 36: fig9
###################################################
AvB <- t.test(groupA, groupB)
plot(AvB$dose, AvB$p, type="l", log="y", xlab="Dose (cGy)", ylab="p-value")
abline(v=2000,col="gray", lty="dashed")
points(2000,approx(AvB$dose, AvB$p, 2000)$y, col="red")
text(2000,approx(AvB$dose, AvB$p, 2000)$y, col="red", labels="V20Gy (p=5.347e-05)",pos=4)


###################################################
### code chunk number 37: RadOnc.Rnw:371-372 (eval = FALSE)
###################################################
## data <- read.DICOM.RT(path="<<DICOM directory>>", verbose=TRUE)


###################################################
### code chunk number 38: RadOnc.Rnw:376-377
###################################################
data("RadOnc")


###################################################
### code chunk number 39: RadOnc.Rnw:384-386
###################################################
teeth[1:2]
c(cord, mandible)


###################################################
### code chunk number 40: RadOnc.Rnw:388-389
###################################################
temp <- teeth


###################################################
### code chunk number 41: RadOnc.Rnw:391-394
###################################################
teeth[[1]]
teeth[[1]] <- teeth[["Tooth #3"]]
teeth


###################################################
### code chunk number 42: RadOnc.Rnw:396-397
###################################################
teeth <- temp


###################################################
### code chunk number 43: RadOnc.Rnw:401-402
###################################################
teeth["Tooth.*"]


###################################################
### code chunk number 44: RadOnc.Rnw:407-411
###################################################
names(teeth) <- c("Larry", "Curly", "Moe")
names(rev(teeth[1:3]))
length(teeth)
lapply(teeth, function(tooth) { range(tooth) })


###################################################
### code chunk number 45: RadOnc.Rnw:413-414
###################################################
teeth <- temp


###################################################
### code chunk number 46: fig10 (eval = FALSE)
###################################################
## plot(mandible)


###################################################
### code chunk number 47: fig11 (eval = FALSE)
###################################################
## plot(cord)


###################################################
### code chunk number 48: RadOnc.Rnw:452-453 (eval = FALSE)
###################################################
## compareStructures(teeth, method="axial", plot=TRUE)


###################################################
### code chunk number 49: fig12
###################################################
compareStructures(teeth, method="axial", plot=TRUE, pixels=40)


###################################################
### code chunk number 50: RadOnc.Rnw:465-466
###################################################
teeth <- teeth[c(1,3)]


###################################################
### code chunk number 51: RadOnc.Rnw:468-469
###################################################
compareStructures(teeth, method="hausdorff", hausdorff.method="mean")


###################################################
### code chunk number 52: RadOnc.Rnw:471-472
###################################################
teeth <- temp


###################################################
### code chunk number 53: RadOnc.Rnw:499-500
###################################################
news(package="RadOnc")


