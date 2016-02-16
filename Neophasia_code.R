# Neophasia project
# Code contributions by Chris Hamm (KU)

library("geomorph")
library("MASS")
library("gdata")


set.seed(2342351)
setwd("~/Desktop/Projects/Neophasia/")

# save(list = ls(), file = "Neophasia_data.R")
# load("Neophasia_data.RData")
(SesInf <- sessionInfo())

##### Remeasurement data to test for error or bias
Neo.meas <- read.csv("Data/nwo_43_re_measure.csv", sep = ",", header = TRUE)
head(Neo.meas)
dim(Neo.meas)
str(Neo.meas)

x.vars <- t(Neo.meas[, c(seq(2, 40, 2))])
str(x.vars)
y.vars <- t(Neo.meas[, c(seq(3, 41, 2))])

plot(x.vars, y.vars, pch = 19, col = rgb(0, 0, 0, 0.2), las = 1, ylab = "Y axis", xlab = "X axis")

xy.vars <- matrix(data = NA, nrow = 20, ncol = 24)
xy.vars[, c(seq(1, 23, 2))] <- x.vars
xy.vars[, c(seq(2, 24, 2))] <- y.vars
# xy.vars

rownames(xy.vars) <- paste("Sample", seq(1:20), sep = "")
colnames(xy.vars) <- paste("Coord", seq(1:24), sep = "")
head(xy.vars)
str(xy.vars)

meas.array <- arrayspecs(A = xy.vars, p = 12, k = 2)
test.gpa <- gpagen(meas.array, ProcD = TRUE, ShowPlot = TRUE)


dimnames(test.gpa$coords)[[3]] <- paste(rep("sample", 20), rep(c(1:2), each = 2), sep = "")
indivs <- dimnames(test.gpa$coords)[[3]]

err <- procD.lm(test.gpa$coords~ factor(indivs), iter = 1e3) # no evidence of between samples 
err


#####
##### Main Neophasia data set - Import and check data 
#####

Neophasia.raw <- read.delim("Data/Neophasia_raw_data.txt", sep = "\t", header = TRUE)
str(Neophasia.raw)
Neophasia.raw$Id <- as.factor(Neophasia.raw$Id)

Neop.taxa <- Neophasia.raw[, 1]
Neop.meta <- Neophasia.raw[, 1:2]
Neophasia.raw <- Neophasia.raw[, 3:26]
head(Neophasia.raw)
summary(Neop.meta)

# Prepare data for GPA using geomorph
Neop.array <- arrayspecs(A = Neophasia.raw, p = 12, k = 2)
Neop.array[, , 1]
dim(Neop.array) # 12 lms, 222 samples
dimnames(Neop.array)[[3]] <- Neop.taxa
Neop.array[1, 1, 1]

Neop.gpa <- gpagen(Neop.array, ProcD = TRUE, ShowPlot = TRUE)
plotOutliers(Neop.gpa$coords) # potential outliers to double check 
# wo_43, ge_486, la_15, me_164

Neop.2d <- two.d.array(Neop.gpa$coords)
head(Neop.2d)


#####
##### Conduct the LDA
#####

Neop.pc <- prcomp(Neop.2d)

summary(Neop.pc)$importance # PC1 explains ~25% of variance, PC2 ~15.8%

Neop.pc.shape <- cbind(Neop.meta, Neop.pc$x[, 1:20]) # remove the last four dimensions because they are empty
dim(Neop.pc.shape)

Neop.lda <- lda(as.matrix(Neop.pc.shape[, 3:22]), Neop.pc.shape$Population, method = "mle")

# plot(Neop.lda)
Neop.lda.scores <- as.matrix(Neop.pc.shape[, 3:22]) %*% as.matrix(Neop.lda$scaling)

Neop.pc.lda <- cbind(Neop.pc.shape, Neop.lda.scores)
unique(Neop.pc.shape$Population) # 8 populations


colors2 <- matrix(Neop.meta$Population, dimnames = list(Neop.meta$Population))
colors2[Neop.meta$Population == "dp"] <- "#e31a1c"
colors2[Neop.meta$Population == "ge"] <- "#a6cee3"
colors2[Neop.meta$Population == "gl"] <- "#1f78b4"
colors2[Neop.meta$Population == "wo"] <- "tomato"
colors2[Neop.meta$Population == "me"] <- "#b2df8a"
colors2[Neop.meta$Population == "ml"] <- "#33a02c"
colors2[Neop.meta$Population == "la"] <- "#ff7f00"
colors2[Neop.meta$Population == "or"] <- "#6a3d9a"

colors <- c("#e31a1c", "#1f78b4", "#a6cee3", "tomato", "#b2df8a", "#33a02c", "#ff7f00", "#6a3d9a")


# pdf(file = "Images/LDA_plot.pdf", bg = "white")
plot(Neop.pc.lda$LD1, Neop.pc.lda$LD2, col = colors2, pch = 19, las = 1, ylim = c(-4.5, 4), xlim = c(-4.5, 4.5), xlab = expression(paste("LD"[1], " 69%")), ylab = expression(paste("LD"[2], " 20%")), cex = 1.3)
legend("bottomleft", legend = c("Donner Pass", "Goat - late", "Goat - early", "Woodfords", "Mendocino - early", "Mendocino - late", "Lang", "Oregon"), col = colors, pch = 19, bty = "n", pt.cex = 1.3)
# dev.off()


# Highlight the points for Mendocino
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "me"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "me"], col = "springgreen2")
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "ml"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "ml"], col = "springgreen4")

# Plot only Mendocino early and late
plot(Neop.pc.lda$LD1[Neop.pc.lda$Population == "me"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "me"], col = "springgreen2", pch = 19, las = 1, cex = 1.5, ylim = c(-5, 5.0), xlim = c(-5, 5), ylab = expression(paste("LD"[2])), xlab = expression(paste("LD"[1])))
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "ml"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "ml"], col = "springgreen4", pch = 19, cex = 1.5)
legend("topleft", legend = c("Mendocino early", "Mendocino late"), pch = 19, col = c("springgreen2", "springgreen4"), bty = "n", pt.cex = 1.5)

# Plot only points from Goat
plot(Neop.pc.lda$LD1[Neop.pc.lda$Population == "ge"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "ge"], col = "plum2", pch = 19, las = 1, cex = 1.5, ylim = c(-3, 4.0), xlim = c(-5, 5), xlab = expression(paste("LD"[1])), ylab = expression(paste("LD"[2])))
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "gl"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "gl"], col = "purple4", pch = 19, cex = 1.5)
legend("topleft", legend = c("Goat early", "Goat late"), col = c("plum2", "purple4"), pch = 19, pt.cex = 1.5, bty = "n")


# plot the PCA
summary(Neop.pc.shape)$importance
plot(Neop.pc.shape[, 3], Neop.pc.shape[,4], pch = 19, col = colors2, ylim = c(-0.07, 0.07), xlim = c(-0.07, 0.07))
legend("bottomleft", legend = c("Donner Pass", "Goat - late", "Goat - early", "Woodfords", "Mendocino - early", "Mendocino - late", "Lang", "Oregon"), col = colors, pch = 19, bty = "n", pt.cex = 1.3) # The differentiation is much clearer in the LDA


###
### Discriminant function
###

# Asks if, given the mean shape for a population, how often is a member of that population correctly assigned to that population. Assignment is based on that individual's shape being closest to the mean shape for its home population. 
Neop.df <- predict(Neop.lda)

Neop.pc.df <- cbind(Neop.pc.lda, Neop.df$posterior)
names(Neop.pc.df)
summary(Neop.pc.df)

dp <- Neop.pc.df[Neop.pc.df$Population == "dp", ]
ge <- Neop.pc.df[Neop.pc.df$Population == "ge", ]
gl <- Neop.pc.df[Neop.pc.df$Population == "gl", ]
la <- Neop.pc.df[Neop.pc.df$Population == "la", ]
me <- Neop.pc.df[Neop.pc.df$Population == "me", ]
ml <- Neop.pc.df[Neop.pc.df$Population == "ml", ]
or <- Neop.pc.df[Neop.pc.df$Population == "or", ]
wo <- Neop.pc.df[Neop.pc.df$Population == "wo", ]

par(mfrow = c(4, 2))
hist(dp$dp, las = 1, col = "skyblue1", main = "Donner Pass")
hist(ge$ge, las = 1, col = "plum2", main = "Goat Mtn. - early")
hist(gl$gl, las = 1, col = "purple4", main = "Goat Mtn. - late")
hist(la$la, las = 1 , col = "midnightblue", main = "Lang")
hist(me$me, las = 1, col = "springgreen2", main = "Mendocino - early")
hist(ml$ml, las = 1, col = "springgreen4", main = "Mendocino - late")
hist(or$or, las = 1, col = "tomato", main = "Oregon")
hist(wo$wo, las = 1, col = "dodgerblue3", main = "Woodfords")
par(mfrow = c(1, 1))


#####
##### Classification problem - linear models 
#####

# Procrustes Distance ANOVA

# I think we should just proceed with the unbalanced samples and report results with a mention that it could be an issue. The randomization procedure used by RRPP is solid.

# Using the full data set to ask if shape varies among populations
pd0 <- lm(Neop.2d ~ 1)
pd0.1 <- lm(Neop.2d ~ Neop.meta$Population)
advanced.procD.lm(pd0, pd0.1)

apD0 <- advanced.procD.lm(pd0, pd0.1, groups = ~ Neop.meta$Population, iter = 1e4)






#####
##### Other covariate data
#####

# read in the covariate data
# Read and merge the covariate data (wings were digitized on left side)
Neophasia.covs <- read.csv("Data/Neophasia_wings.csv", header = TRUE)
head(Neophasia.covs)
summary(Neophasia.covs)
# WAR = wing area left
# WAR = wing area right
# MTL = melanization total left
# MTR = melanization total right
summary(Neophasia.covs)
which(Neophasia.covs$WAL >= 800) # We won't use WAL because of an error

Neophasia.raw2 <- read.delim("Data/Neophasia_raw_data.txt", sep = "\t", header = TRUE)
str(Neophasia.raw2)
intersect(Neophasia.raw2$Id, Neophasia.covs$Id)
Neophasia.merged <- merge(x = Neophasia.raw2, y = Neophasia.covs, by = "Id")
summary(Neophasia.merged)

Neop.taxa2 <- Neophasia.merged[, 1]
Neop.meta2 <- Neophasia.merged[, c(1:2, 27:30)]
head(Neop.meta2) # create a matched covariate data set
summary(Neop.meta2)

Neophasia.raw2 <- Neophasia.merged[, 3:26] # just the raw coordinates of the matched data
head(Neophasia.raw2)


# Prepare data for GPA using geomorph
Neop.array2 <- arrayspecs(A = Neophasia.raw2, p = 12, k = 2)
Neop.array2[, , 1]
dim(Neop.array2) # 12 lms, 180 samples
dimnames(Neop.array2)[[3]] <- Neop.taxa2
Neop.array2[1, 1, 1]

Neop.gpa2 <- gpagen(Neop.array2, ProcD = TRUE, ShowPlot = TRUE)
Neop.2d2 <- two.d.array(Neop.gpa2$coords)
head(Neop.2d2)

# sanity check
plot(Neop.meta2$WAR, log(Neop.gpa2$Csize), pch = 19) # Wing area and log centroid size very highly correlated. 



#####
##### Using melanization levels
#####

colors3 <- matrix(Neophasia.merged$Population, dimnames = list(Neophasia.merged$Population))
colors3[Neophasia.merged$Population == "dp"] <- "skyblue1"
colors3[Neophasia.merged$Population == "ge"] <- "plum2"
colors3[Neophasia.merged$Population == "gl"] <- "purple4"
colors3[Neophasia.merged$Population == "wo"] <- "dodgerblue3"
colors3[Neophasia.merged$Population == "me"] <- "springgreen2"
colors3[Neophasia.merged$Population == "ml"] <- "springgreen4"
colors3[Neophasia.merged$Population == "la"] <- "midnightblue"
colors3[Neophasia.merged$Population == "or"] <- "tomato"

# pdf(file = "Images/Mel-box.pdf", bg = "white")
boxplot(Neophasia.merged$MTR ~ sort(Neophasia.merged$Population, decreasing = FALSE), col = unique(colors3), outline = FALSE, ylab = "Melanized area", xlab = "Population", varwidth = FALSE, las = 1, staplewex = 0.95, ylim = c(35, 125))
text(x = 1, y = 113, "D")
text(x = 2, y = 79, "A")
text(x = 3, y = 120, "D")
text(x = 4, y = 104, "CD")
text(x = 5, y = 85, "A")
text(x = 6, y = 117, "BC")
text(x = 7, y = 89, "AB")
text(x = 8, y = 87, "A")
# dev.off()



# Does population explain melanization? This generates the same associations as Figure 7, but without using the model residuals. To me it is more intuitive and has the added benefit of making a population level distance matrix for melanization levels. 
pd1 <- lm(Neop.meta2$MTL ~ 1)
pd1.1 <- lm(Neop.meta2$MTL ~ Neop.meta2$Population)

apD1 <- advanced.procD.lm(pd1, pd1.1, groups = ~ Neop.meta2$Population, iter = 1e4)




# association between pairwise shape distance and genetic distance (Gst)
Gst <- read.csv("Gst-dist.csv", header = TRUE, row.names = 1)
Gst

rownames(Gst)
colnames(Gst)
rownames(apD0$Means.dist)
colnames(apD0$Means.dist)

Sdist.diag <- lowerTriangle(apD0$Means.dist)
Gst.diag <- lowerTriangle(Gst)
Mel.diag <- lowerTriangle(apD1$Means.dist)


plot(Gst.diag, Sdist.diag, pch = 19, las = 1, ylab = "Shape distance", xlab = "Genetic distance", xlim = c(0.025, 0.08), ylim = c(0.015, 0.045))
GSlm <- lm(Sdist.diag ~ Gst.diag)
summary(GSlm)
abline(GSlm)

plot(Gst.diag, Mel.diag, pch = 19, las = 1, ylab = "Melanization distance", xlab = "Genetic distance", xlim = c(0.025, 0.08), ylim = c(0, 40))
GMellm <- lm(Mel.diag ~ Gst.diag)
summary(GMellm)
abline(GMellm)

plot(Sdist.diag, Mel.diag, pch = 19, las = 1, ylab = "Melanization distance", xlab = "Shape distance", xlim = c(0.015, 0.045), ylim = c(0, 40))
SMellm <- lm(Mel.diag ~ Sdist.diag)
summary(SMellm)
abline(SMellm)

# What about making the point that the shape differences between the sympatric populations are on the scale of those we observe for the really geographically distant ones? Is it secondary contact or in situ change? 




#####
##### Plot just the sympatric pairs with warp grids
#####

p <- dim(Neop.gpa$coords)[1]
k <- dim(Neop.gpa$coords)[2]
group <- Neop.meta[, 2]
Y <- array(NA, dim = c(p, k, length(levels(group))))
dimnames(Y)[[3]] <- levels(group)

for(i in 1:length(levels(group))){
	grp <- Neop.2d[which(group == levels(group)[i]), ]
	foo <- arrayspecs(grp, p, k)
	Y[, , i] <- mshape(foo)
}
Y[,, 8] # wo
# ge is 2, gl 3, me 5, ml 6.

reference <- mshape(Neop.gpa$coords)
reference <- reference %*% matrix(c(0, 1, -1, 0), ncol = 2, byrow = TRUE) # rotate the object 180 degrees

Y.me2 <- Y[, , 5] %*% matrix(c(0, 1, -1, 0), ncol = 2, byrow = TRUE)
Y.ml2 <- Y[, , 6] %*% matrix(c(0, 1, -1, 0), ncol = 2, byrow = TRUE)
Y.ge2 <- Y[, , 2] %*% matrix(c(0, 1, -1, 0), ncol = 2, byrow = TRUE)
Y.gl2 <- Y[, , 3] %*% matrix(c(0, 1, -1, 0), ncol = 2, byrow = TRUE)


plotRefToTarget(reference, Y.me2, method = "TPS", mag = 1)


# pdf(file = "Images/Pairs-plot2.pdf", bg = "white")
plot(Neop.pc.lda$LD1[Neop.pc.lda$Population == "me"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "me"], col = "#b2df8a", pch = 19, las = 1, cex = 1.5, ylim = c(-6, 6.0), xlim = c(-6, 6), ylab = expression(paste("LD"[2])), xlab = expression(paste("LD"[1])))
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "ml"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "ml"], col = "#33a02c", pch = 19, cex = 1.5)
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "ge"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "ge"], col = "#a6cee3", pch = 19, las = 1, cex = 1.5)
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "gl"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "gl"], col = "#1f78b4", pch = 19, cex = 1.5)
legend("top", legend = c("Mendocino - early", "Mendocino - late", "Goat - early", "Goat - late"), pch = 19, col = c("#b2df8a", "#33a02c", "#a6cee3", "#1f78b4"), bty = "n", pt.cex = 1.5)



# lower right, me
par(fig = c(0.6, 1, 0.1, 0.4), new = TRUE)
plot.new()
par(mar = c(1, 1, 1, 1))
me.pars <- gridPar(tar.pt.bg = "#b2df8a", tar.pt.size = 1.5)
plotRefToTarget(reference, Y.me2, method = "TPS", mag = 5, gridPars = me.pars)

#lower left, ml
par(fig = c(0.05, 0.45, 0.1, 0.4), new = TRUE)
plot.new()
par(mar = c(rep(1, 4)))
ml.pars <- gridPar(tar.pt.bg = "#33a02c", tar.pt.size = 1.5)
plotRefToTarget(reference, Y.ml2, method = "TPS", mag = 5, gridPars = ml.pars)

# upper right, ge
par(fig = c(0.6, 1, 0.63, 0.93), new = TRUE)
plot.new()
par(mar = c(rep(1, 4)))
ge.pars <- gridPar(tar.pt.bg = "#a6cee3", tar.pt.size = 1.5)
plotRefToTarget(reference, Y.ge2, method = "TPS", mag = 5, gridPars = ge.pars)

#upper left, gl
par(fig = c(0.05, 0.45, 0.63, 0.93), new = TRUE)
plot.new()
par(mar = c(rep(1, 4)))
gl.pars <- gridPar(tar.pt.bg = "#1f78b4", tar.pt.size = 1.5)
plotRefToTarget(reference, Y.gl2, method = "TPS", mag = 5, gridPars = gl.pars)

# dev.off()

