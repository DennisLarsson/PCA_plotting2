#version 2.2
library(adegenet)
# install using: install.packages("adegenet")

setwd("/path/to/workDirectory/")
infile <- "spicatum.stru"
outputfile_name <- "spicatum"

# give the name of the popmap if it is in the work directory
# or the full path if it is somewhere else
popmap <- "/path/to/popmap"
number_indv <- 110
number_loci <- 2544

# where the legend will be placed in the PC1-2 plot:
#  topright, topleft, bottomright, bottomleft
ld_pos_pc12 <- "topright"
# where the legend will be placed in the PC1-3 plot:
#  topright, topleft, bottomright, bottomleft
ld_pos_pc13 <- "topright"

pop <- read.delim(popmap, header = FALSE, as.is = TRUE)

# make sure to edit n.ind to the number of individuals your
# dataset has and n.loc to the number of loci you have.
object1 <- read.structure(infile,
                          n.ind = number_indv,
                          n.loc = number_loci,
                          onerowperind = FALSE,
                          col.lab = 1,
                          col.pop = 2,
                          row.marknames = 1,
                          NA.char = "-9",
                          ask = FALSE)

data_matrix <- tab(object1, freq = TRUE, NA.method = "mean")

pca1 <- dudi.pca(data_matrix, scale = FALSE, scannf = FALSE, nf = 3)

power_value_1 <- round((pca1$eig[1] / sum(pca1$eig)) * 100, digits = 2)
power_value_2 <- round((pca1$eig[2] / sum(pca1$eig)) * 100, digits = 2)
power_value_3 <- round((pca1$eig[3] / sum(pca1$eig)) * 100, digits = 2)
val <- c(power_value_1, power_value_2, power_value_3)

# introduce function which goes through the popmap and indexes each uniq
# occurence of popname, then add this number in a new list
uniq_pop <- unique(pop[, 2])
pop_shape <- data.frame(shape = 1:length(uniq_pop), row.names = uniq_pop)
shape <- 1
i <- 1
while (i <= length(uniq_pop)) {
  pop_shape[i, 1] <- shape
  i <- i + 1
  shape <- shape + 1
  if (shape > 18) {
    shape <- 1
  }
}

pop_features <- data.frame(color = 1:length(pop[,2]),
                           shape = 1:length(pop[,2]),
                           row.names = pop[[1]])
pop_chp <- 1
indv <- 1
for (i in 1:length(pop[,2])) {
  pop_match <- match(pop[i, 2], uniq_pop)
  pop_features[i, 1] <- pop_match
  pop_features[i, 2] <- pop_shape[pop_match, 1]
}

plot.PCA <- function(x, y, pos, plot.cex, plot.pt.cex, plot.ncol) {
  plot(pca1$li[1, x], pca1$li[1, y],
       col = pop_features[1, 1], pch = pop_features[1, 2],
       xlab = paste("PC", x, " - ", val[x], "%", sep = ""),
       ylab = paste("PC", y, " - ", val[y], "%", sep = ""),
       xlim = c(min(pca1$li[, x]), max(pca1$li[, x])),
       ylim = c(min(pca1$li[, y]), max(pca1$li[, y])),
       main = paste("PCA plot of PC", x, " vs PC", y, sep = ""))

  abline(0, 0, lty = 5)
  abline(0, 180, lty = 5)

  for (sample in 2:number_indv) {
    points(pca1$li[sample, x], pca1$li[sample, y],
           col = pop_features[sample, 1], pch = pop_features[sample, 2])
  }

  legend(x = pos,
         legend = uniq_pop,
         col = unique(pop_features[[1]]),
         pch = unique(pop_features[[2]]),
         pt.cex = plot.pt.cex,
         cex = plot.cex,
         ncol = plot.ncol,
         y.intersp = 1)
}

pdf(file = paste("PCA_", outputfile_name, ".pdf", sep = ""),
    height = 8,
    width = 8,
    title = outputfile_name)
barplot(pca1$eig[1:10], xlab = "component", ylab = "eigen value")

plot.PCA(1, 2, ld_pos_pc12, plot.cex = 0.7, plot.pt.cex = 1, plot.ncol = 2)
plot.PCA(1, 3, ld_pos_pc13, plot.cex = 0.7, plot.pt.cex = 1, plot.ncol = 2)

indvnames <- as.vector(labels(pca1$li)[[1]])

plot.default(x = pca1$li[, 1], y = pca1$li[, 2],
             xlab = "PC1", ylab = "PC2",
             xlim = c(min(pca1$li[, 1]), max(pca1$li[, 1])),
             ylim = c(min(pca1$li[, 2]), max(pca1$li[, 2])),
             main = "PCA plot of PC1 vs PC2 text labels")

text(pca1$li[, 1], pca1$li[, 2], labels = indvnames, cex = 0.4, pos = 4)

plot.default(x = pca1$li[, 1], y = pca1$li[, 3],
             xlab = "PC1", ylab = "PC3",
             xlim = c(min(pca1$li[, 1]), max(pca1$li[, 1])),
             ylim = c(min(pca1$li[, 3]), max(pca1$li[, 3])),
             main = "PCA plot of PC1 vs PC3 text labels")

text(pca1$li[, 1], pca1$li[, 3], labels = indvnames, cex = 0.4, pos = 4)
dev.off()
