# Neophasia project

library("geomorph")
library("MASS")

setwd("~/Desktop/Projects/Neophasia/")

# save(list = ls(), file = "Neophasia_data.R")
# load("Neophasia_data.R")

# Import and chech data 
Neophasia.raw <- read.delim("Neophasia_raw_data.txt", sep = "\t", header = TRUE)
Neop.taxa <- Neophasia.raw[, 1]
Neop.meta <- Neophasia.raw[, 1:2]
Neophasia.raw <- Neophasia.raw[, 3:26]
head(Neophasia.raw)


# Prepare data for GPA using geomorph
Neop.array <- arrayspecs(A = Neophasia.raw, p = 12, k = 2)
Neop.array[, , 1]
dim(Neop.array) # 12 lms, 222 samples
dimnames(Neop.array)[[3]] <- Neop.taxa
Neop.array[1, 1, 1]

Neop.2d <- two.d.array(Neop.array)
head(Neop.2d)

Neop.gpa <- gpagen(Neop.array, ProcD = TRUE, ShowPlot = TRUE)
plotOutliers(Neop.gpa$coords) # potential outliers to double check 
# wo_43, ge_486, la_15, me_164

# Now the LDA
Neop.pc <- prcomp(Neop.2d)
summary(Neop.pc)$importance # PC1 explains ~69% of variance, PC2 ~20%

Neop.pc.shape <- cbind(Neop.meta, Neop.pc$x[, 1:20]) # remove the last four dimensions because they are empty
dim(Neop.pc.shape)

Neop.lda <- lda(as.matrix(Neop.pc.shape[, 3:22]), Neop.pc.shape$Population, method = "mle")

# plot(Neop.lda)
Neop.lda.scores <- as.matrix(Neop.pc.shape[, 3:22]) %*% as.matrix(Neop.lda$scaling)

Neop.pc.lda <- cbind(Neop.pc.shape, Neop.lda.scores)
unique(Neop.pc.shape$Population) # 8 populations

colors <- c("goldenrod", "dodgerblue", "dark red", "dark green", "dark grey", "dark blue", "purple", "red")

plot(Neop.pc.lda$LD1, Neop.pc.lda$LD2, col = colors, pch = 19, las = 1, ylim = c(-5, 5), xlim = c(-5, 5), xlab = expression(paste("LD"[1], " 69%")), ylab = expression(paste("LD"[2], " 20%")), cex = 1.3)
legend("topleft", legend = c("Donner Pass", "Goat - late", "Woodfords", "Mendocino - early", "Mendocino - late", "Goat - early", "Lang", "Oregon"), col = colors, pch = 19, bty = "n", pt.cex = 1.5)

# still need to ask about early / late populations

# df
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
hist(dp$dp, las = 1, col = "goldenrod")
hist(ge$ge, las = 1, col = "dodgerblue")
hist(gl$gl, las = 1, col = "dark red")
hist(la$la, las = 1 , col = "dark green")
hist(me$me, las = 1, col = "grey")
hist(ml$ml, las = 1, col = "dark blue")
hist(or$or, las = 1, col = "purple")
hist(wo$wo, las = 1, col = "red")