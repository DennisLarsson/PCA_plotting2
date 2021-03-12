library(adegenet)
#if you dont have adegenet installed, install using: install.packages("adegenet")

setwd("/path/to/workDirectory/")
infile<-"spicatum.stru"
outputfile.name <-"spicatum"
popmap="/path/to/popmap"      #give the name of the popmap if it is in the work directory or the full path if it is somewhere else
number.indv=110
number.loci=2544
ld.pos.PC12 = "topright" #where the legend will be placed in the PC1-2 plot: topright, topleft, bottomright, bottomleft
ld.pos.PC13 = "topright" #where the legend will be placed in the PC1-3 plot: topright, topleft, bottomright, bottomleft

pop <-read.delim(popmap, header = FALSE, as.is=T)
#make sure to edit n.ind to the number of individuals your dataset has and n.loc to the number of loci you have.
object1 <- read.structure(infile, n.ind = number.indv, n.loc = number.loci, onerowperind = FALSE, col.lab = 1, col.pop = 2, row.marknames = 1, 
                          NA.char = "-9", ask = FALSE)

X <- tab(object1, freq = TRUE, NA.method = "mean")

pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)

powerVal1 = round((pca1$eig[1]/sum(pca1$eig))*100,digits = 2)
powerVal2 = round((pca1$eig[2]/sum(pca1$eig))*100,digits = 2)
powerVal3 = round((pca1$eig[3]/sum(pca1$eig))*100,digits = 2)
val<-c(powerVal1,powerVal2,powerVal3)


# introduce function which goes through the popmap and indexes each uniq occurence of popname, then add this number in a new list
uniq.pop <- unique(pop[,2])
pop.shape <- data.frame(shape=1:length(uniq.pop),row.names = uniq.pop)
shape=1
i=1
while (i <= length(uniq.pop)) {
  pop.shape[i,1]=shape
  i=i+1
  shape=shape+1
  if (shape > 18) {
    shape=1
  }
}

pop.features <- data.frame(color=1:length(pop[,2]),shape=1:length(pop[,2]),row.names = pop[[1]])
pop.chp=1
indv=1
for (i in 1:length(pop[,2])) {
  pop.match=match(pop[i,2],uniq.pop)
  pop.features[i,1]=pop.match
  pop.features[i,2]=pop.shape[pop.match,1]
}

plot.PCA <- function(x,y,pos,plot.cex,plot.pt.cex,plot.ncol) {
  plot(pca1$li[1,x], pca1$li[1,y], col=pop.features[1,1], pch=pop.features[1,2], xlab=paste("PC",x," - ",val[x],"%",sep=""), ylab=paste("PC",y," - ",val[y],"%",sep=""), 
       xlim=c(min(pca1$li[,x]),max(pca1$li[,x])), ylim=c(min(pca1$li[,y]),max(pca1$li[,y])), main = paste("PCA plot of PC",x," vs PC",y,sep=""))
  abline(0,0, lty=5)
  abline(0,180, lty=5)
  for (sample in 2:number.indv) {
    points(pca1$li[sample,x], pca1$li[sample,y], col=pop.features[sample,1], pch=pop.features[sample,2])
  }
  legend(x=pos, legend=uniq.pop, col = unique(pop.features[[1]]), pch = unique(pop.features[[2]]), pt.cex=plot.pt.cex, cex=plot.cex, ncol=plot.ncol,y.intersp=1)
}

pdf(file=paste("PCA_",outputfile.name,".pdf",sep=""), height = 8, width = 8, title = outputfile.name)
barplot(pca1$eig[1:10], xlab="component", ylab="eigen value")

plot.PCA(1,2,ld.pos.PC12,plot.cex=0.7,plot.pt.cex=1,plot.ncol=2)
plot.PCA(1,3,ld.pos.PC13,plot.cex=0.7,plot.pt.cex=1,plot.ncol=2)
dev.off()

indvnames <- as.vector(labels(pca1$li)[[1]])

png(filename = paste("PCA_PC12_label_",outputfile.name,".png",sep=""),width = 1200, height = 1200)
plot.default(x=pca1$li[,1], y=pca1$li[,2], xlab="PC1", ylab="PC2",xlim=c(min(pca1$li,1),max(pca1$li,1)), ylim=c(min(pca1$li,2),max(pca1$li,2)), main = "PCA plot of PC1 vs PC2 text labels")
text(pca1$li[,1], pca1$li[,2], labels=indvnames, cex= 0.9, pos = 4)
dev.off()

png(filename = paste("PCA_PC13_label_",outputfile.name,".png",sep=""),width = 1200, height = 1200)
plot.default(x=pca1$li[,1], y=pca1$li[,3], xlab="PC1", ylab="PC3",xlim=c(min(pca1$li,1),max(pca1$li,1)), ylim=c(min(pca1$li,3),max(pca1$li,3)), main = "PCA plot of PC1 vs PC3 text labels")
text(pca1$li[,1], pca1$li[,3], labels=indvnames, cex= 0.9, pos = 4)
dev.off()
