library("vegan")

# pdf("Images/LDA-ellipse2.pdf", bg = "white")
plot(Neop.pc.lda$LD1, Neop.pc.lda$LD2, col = colors2, pch = 19, las = 1, ylim = c(-4.5, 4), xlim = c(-4.5, 4.5), xlab = expression(paste("LD"[1], " 69%")), ylab = expression(paste("LD"[2], " 20%")), cex = 1.3)
legend("bottomleft", legend = c("Donner Pass", "Goat - late", "Goat - early", "Woodfords", "Mendocino - early", "Mendocino - late", "Lang", "Oregon"), col = colors, pch = 19, bty = "n", pt.cex = 1.3)

# ordiellipse(cbind(Neop.pc.lda$LD1, Neop.pc.lda$LD2), groups = Neop.meta$Population, col = colors2[, 1], draw = "lines", kind = "sd", show.groups = "dp")


ordiellipse(cbind(Neop.pc.lda$LD1, Neop.pc.lda$LD2), groups = Neop.meta$Population, col = "#e31a1c", draw = "lines", kind = "sd", show.groups = "dp")



ordiellipse(cbind(Neop.pc.lda$LD1, Neop.pc.lda$LD2), groups = Neop.meta$Population, col = "#1f78b4", draw = "lines", kind = "sd", show.groups = "gl")
ordiellipse(cbind(Neop.pc.lda$LD1, Neop.pc.lda$LD2), groups = Neop.meta$Population, col = "#a6cee3", draw = "lines", kind = "sd", show.groups = "ge")
ordiellipse(cbind(Neop.pc.lda$LD1, Neop.pc.lda$LD2), groups = Neop.meta$Population, col = "#ff7f00", draw = "lines", kind = "sd", show.groups = "la")
ordiellipse(cbind(Neop.pc.lda$LD1, Neop.pc.lda$LD2), groups = Neop.meta$Population, col = "#1f78b4", draw = "lines", kind = "sd", show.groups = "gl")
ordiellipse(cbind(Neop.pc.lda$LD1, Neop.pc.lda$LD2), groups = Neop.meta$Population, col = "#b2df8a", draw = "lines", kind = "sd", show.groups = "me")
ordiellipse(cbind(Neop.pc.lda$LD1, Neop.pc.lda$LD2), groups = Neop.meta$Population, col = "#33a02c", draw = "lines", kind = "sd", show.groups = "ml")
ordiellipse(cbind(Neop.pc.lda$LD1, Neop.pc.lda$LD2), groups = Neop.meta$Population, col = "tomato", draw = "lines", kind = "sd", show.groups = "wo")
ordiellipse(cbind(Neop.pc.lda$LD1, Neop.pc.lda$LD2), groups = Neop.meta$Population, col = "#6a3d9a", draw = "lines", kind = "sd", show.groups = "or")
 dev.off()