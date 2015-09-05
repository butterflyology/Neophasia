# Neophasia project
# Code contributions by Chris Hamm (KU)

library("geomorph")
library("MASS")
library("vegan")
(SesInf <- sessionInfo())


set.seed(2342351)
setwd("~/Desktop/Projects/Neophasia/")

# save(list = ls(), file = "Neophasia_data.R")
# load("Neophasia_data.R")

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
colors2[Neop.meta$Population == "dp"] <- "goldenrod"
colors2[Neop.meta$Population == "ge"] <- "dark blue"
colors2[Neop.meta$Population == "gl"] <- "dodgerblue"
colors2[Neop.meta$Population == "wo"] <- "dark red"
colors2[Neop.meta$Population == "me"] <- "dark green"
colors2[Neop.meta$Population == "ml"] <- "dark grey"
colors2[Neop.meta$Population == "la"] <- "purple"
colors2[Neop.meta$Population == "or"] <- "red"

colors <- c("goldenrod", "dodgerblue", "dark blue", "dark red", "dark green", "dark grey", "purple", "red")


# pdf(file = "LDA_plot.pdf", bg = "white")
plot(Neop.pc.lda$LD1, Neop.pc.lda$LD2, col = colors2, pch = 19, las = 1, ylim = c(-4.5, 4), xlim = c(-4.5, 4.5), xlab = expression(paste("LD"[1], " 69%")), ylab = expression(paste("LD"[2], " 20%")), cex = 1.3)
legend("bottomleft", legend = c("Donner Pass", "Goat - late", "Goat - early", "Woodfords", "Mendocino - early", "Mendocino - late", "Lang", "Oregon"), col = colors, pch = 19, bty = "n", pt.cex = 1.3)
# dev.off()


# Highlight the points for Mendocino
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "me"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "me"], col = "red")
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "ml"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "ml"], col = "yellow")

# Plot only Mendocino early and late
plot(Neop.pc.lda$LD1[Neop.pc.lda$Population == "me"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "me"], col = "dark green", pch = 19, las = 1, cex = 1.5, ylim = c(-5, 5.0), xlim = c(-5, 5), ylab = expression(paste("LD"[2])), xlab = expression(paste("LD"[1])))
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "ml"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "ml"], col = "dark grey", pch = 19, cex = 1.5)
legend("topleft", legend = c("Mendocino early", "Mendocino late"), pch = 19, col = c("dark green", "dark grey"), bty = "n", pt.cex = 1.5)

# Plot only points from Goat
plot(Neop.pc.lda$LD1[Neop.pc.lda$Population == "ge"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "ge"], col = "dark blue", pch = 19, las = 1, cex = 1.5, ylim = c(-3, 4.0), xlim = c(-5, 5), xlab = expression(paste("LD"[1])), ylab = expression(paste("LD"[2])))
points(Neop.pc.lda$LD1[Neop.pc.lda$Population == "gl"], Neop.pc.lda$LD2[Neop.pc.lda$Population == "gl"], col = "dodgerblue", pch = 19, cex = 1.5)
legend("topleft", legend = c("Goat early", "Goat late"), col = c("dark blue", "dodgerblue"), pch = 19, pt.cex = 1.5, bty = "n")


# plot the PCA
summary(Neop.pc.shape)$importance
plot(Neop.pc.shape[, 3], Neop.pc.shape[,4], pch = 19, col = colors2, ylim = c(-0.07, 0.07), xlim = c(-0.07, 0.07))
legend("bottomleft", legend = c("Donner Pass", "Goat - late", "Goat - early", "Woodfords", "Mendocino - early", "Mendocino - late", "Lang", "Oregon"), col = colors, pch = 19, bty = "n", pt.cex = 1.3)


###
### Discriminant function
###

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
hist(dp$dp, las = 1, col = "goldenrod", main = "Donner Pass")
hist(ge$ge, las = 1, col = "dodgerblue", main = "Goat Mtn. - early")
hist(gl$gl, las = 1, col = "dark red", main = "Goat Mtn. - late")
hist(la$la, las = 1 , col = "dark green", main = "Lang")
hist(me$me, las = 1, col = "dark grey", main = "Mendocino - early")
hist(ml$ml, las = 1, col = "dark blue", main = "Mendocino - late")
hist(or$or, las = 1, col = "purple", main = "Oregon")
hist(wo$wo, las = 1, col = "red", main = "Woodfords")



#####
##### Classification problem - linear models 
#####
# Write a linear model asking if wing shape is explained by Population - using the reduced dimensionality of the PCA because sample size is not ginormous. nb, the first 13 PCs account for 95% of the variance.

lm1 <- lm(as.matrix(Neop.pc$x[, 1:5]) ~ Neop.meta$Population)
summary(lm1) # let's just consider the first five PCs worth
car::Anova(lm1, test = "Wilks", type = "III") # yes, the populations are different morphologically (but type III SS can be wonky because we may want the sequential test of the Type-I flavor)
# pairwise tests might not be cool because the sample sizes are different


# exploring how other "standard" tests handle the data
adon1 <- adonis(Neop.2d ~ Neop.meta$Population, method = "euclidean", permutations = 1e3)
adon1

lm2 <- lm(Neop.pc$x[, 1:5] ~ Neop.meta$Population)
summary(manova(lm2), test = "Wilks")

# using the geomorph function, we will probably want to report this one. For single factor designs, the two RRPP approaches are the same.
pD1 <- procD.lm(Neop.2d ~ Neop.meta$Population, iter = 1e3, RRPP = FALSE) # models Population differences on shape



pD2 <- procD.lm(Neop.2d ~ Neop.meta$Population + log(Neop.gpa$Csize), iter = 1e3, RRPP = TRUE) # models population differences and log(CS)

pD3 <- procD.lm(Neop.2d ~ Neop.meta$Population * log(Neop.gpa$Csize), iter = 1e3, RRPP = FALSE) # models size variation in populations

pD3r <- procD.lm(Neop.2d ~ Neop.meta$Population * log(Neop.gpa$Csize), iter = 1e3, RRPP = TRUE) # models size variation in populations with residual randomization permutation procedure







# read in the covariate data
# Read and merge the covariate data (wings were digitized on left side)
Neophasia.covs <- read.csv("Data/Neophasia_wings.csv", header = TRUE)
head(Neophasia.covs)
summary(Neophasia.covs)
# WAR = wing area left
# WAR = wing area right
# MTL = melanization total left
# MTR = melanization total right

str(Neophasia.covs)
head(Neophasia.covs)

summary(Neophasia.covs)
which(Neophasia.covs$WAL >= 800) # We won't use WAL because of an error


# Create a 2nd data set for the merged taxa (few individuals were landmarked but no covariates were taken from them)
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


# asking the question: is wing shape explained by population, melanization, and the interaction of the two? 


prlm0 <- procD.lm(Neop.2d2 ~ log(Neop.gpa2$Csize) + (Neop.meta2$Population * Neop.meta2$MTL) + log(Neop.gpa2$Csize) * (Neop.meta2$Population * Neop.meta2$MTL), RRPP = TRUE, iter = 1e4)
# a non significant log(CS) * (Pop * Mel) would indicate a common population by melanization level has common size-shape allometry


prlm1 <- procD.lm(Neop.2d2 ~ log(Neop.gpa2$Csize) + Neop.meta2$Population * Neop.meta2$MTL, RRPP = TRUE, iter = 1e4) 



prlm2 <- procD.lm(Neop.2d2 ~ Neop.meta2$Population * Neop.meta2$MTL, RRPP = TRUE, iter = 1e4)


pdlm1 <- procD.lm(Neop.2d2 ~ Neop.meta2$Population + Neop.meta2$MTL + Neop.meta2$Population * Neop.meta2$MTL, iter = 1e3, RRPP = TRUE)
pdlm1 # Shape varies by population, and the interaction of population and melanization level significantly interact as well. 


aDlm1 <- advanced.procD.lm(Neop.2d2 ~ Neop.meta2$Population + Neop.meta2$MTL,  ~ Neop.meta2$Population + Neop.meta2$MTL + Neop.meta2$Population * Neop.meta2$MTL, groups = ~ Neop.meta2$Population, iter = 1e3)

aDlm1a <- advanced.procD.lm(Neop.2d2 ~ Neop.meta2$Population + Neop.meta2$MTL,  ~ Neop.meta2$Population + Neop.meta2$MTL + Neop.meta2$Population * Neop.meta2$MTL, groups = ~ Neop.meta2$Population, iter = 1e3, slope = ~ Neop.meta2$MTL)


aDlm3 <- advanced.procD.lm(Neop.2d2 ~ log(Neop.gpa2$Csize) + Neop.meta2$Population + Neop.meta2$Population * log(Neop.gpa$Csize), ~ )







pdlm2 <- procD.lm(Neop.2d2 ~ Neop.meta2$MTL, iter = 1e3)
pdlm2 # melanization is strongly associated with shape

pdlm3 <- procD.lm(Neop.2d2 ~ Neop.meta2$Population, iter = 1e3)
pdlm3

# post-foc test 
adlm3 <- advanced.procD.lm(Neop.2d2 ~ Neop.meta2$Population, iter = 10)


pdlm4 <- procD.lm(Neop.2d2 ~ Neop.meta2$Population * Neop.meta2$MTL, iter = 1e3, RRPP = TRUE)
pdlm4






colors3 <- matrix(Neophasia.merged$Population, dimnames = list(Neophasia.merged$Population))
colors3[Neophasia.merged$Population == "dp"] <- "goldenrod"
colors3[Neophasia.merged$Population == "ge"] <- "dark blue"
colors3[Neophasia.merged$Population == "gl"] <- "dodgerblue"
colors3[Neophasia.merged$Population == "wo"] <- "dark red"
colors3[Neophasia.merged$Population == "me"] <- "dark green"
colors3[Neophasia.merged$Population == "ml"] <- "dark grey"
colors3[Neophasia.merged$Population == "la"] <- "purple"
colors3[Neophasia.merged$Population == "or"] <- "red"


summary(lma3 <- lm(Neophasia.merged$MTR ~ Neophasia.merged$WAR)) # 0.42 R^2 # melanization is highly correlated between wings, no surprise
plot(Neophasia.merged$WAR, Neophasia.merged$MTR, pch = 19, col = colors3, cex = 1.25, ylab = "Wing Melanization", xlab = "Wing Area") # super high correlation between melanization area and wing total wing area
abline(lma3, lwd = 2)
legend("topleft", legend = c("Donner Pass", "Goat - late", "Goat - early", "Woodfords", "Mendocino - early", "Mendocino - late", "Lang", "Oregon"), col = colors, pch = 19, bty = "n", pt.cex = 1.3)

boxplot(Neophasia.merged$MTR ~ sort(Neophasia.merged$Population, decreasing = FALSE), col = "grey", outline = FALSE, ylab = "Melanized area", xlab = "Population")

# consider plotting by latitude or altitude (sort, decreasing... etc.)
