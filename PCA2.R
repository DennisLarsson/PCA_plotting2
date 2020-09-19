library(adegenet)
#Version 2.0
#if you dont have adegenet installed, install using: install.packages("adegenet")

setwd("/path/to/workDirectory/")
infile<-"spicatum.stru"
outputfile.name <-"spicatum"
popmap="/path/to/popmap"      #give the name of the popmap if it is in the work directory or the full path if it is somewhere else
number.indv=110
number.loci=2544

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

plot.PCA <- function(x,y,pos,plot.cex,plot.pt.cex,plot.ncol) {
  pca.names<-names(object1@tab[,1])
  pop.feat=1
  pop.nr=1
  legend.names=c(pop[1,2])
  legend.feat=c(1)
  
  plot(pca1$li[1,x], pca1$li[1,y], col=object1$pop[1], pch=pop.feat, xlab=paste("PC",x," - ",val[x],"%",sep=""), ylab=paste("PC",y," - ",val[y],"%",sep=""), 
       xlim=c(min(pca1$li[,x]),max(pca1$li[,x])), ylim=c(min(pca1$li[,y]),max(pca1$li[,y])), main = paste("PCA plot of PC",x," vs PC",y,sep=""))
  
  abline(0,0, lty=5)
  abline(0,180, lty=5)
  
  for (indv in 2:number.indv) {
    points(pca1$li[indv,x], pca1$li[indv,y], col=object1$pop[indv], pch=pop.feat)
    if (indv+1>number.indv) {
      break
    }
    else if (as.character(pop[indv,2]) != as.character(pop[indv+1,2])) {
      pop.feat=pop.feat+1
      pop.nr=pop.nr+1
      legend.names[pop.nr]=pop[indv+1,2]
      legend.feat[pop.nr]=pop.feat
      if (pop.feat>18){
        pop.feat=1
        legend.feat[pop.nr]=pop.feat
      }
    }
  }
  legend(x=pos, legend=legend.names, col = unique(object1$pop), pch = legend.feat, pt.cex=plot.pt.cex, cex=plot.cex, ncol=plot.ncol,y.intersp=1)
}

pdf(file=paste("PCA_",outputfile.name,".pdf",sep=""), height = 8, width = 8, title = outputfile.name)
barplot(pca1$eig[1:10], xlab="component", ylab="eigen value")

plot.PCA(1,2,"bottomleft",plot.cex=0.7,plot.pt.cex=1,plot.ncol=2)
plot.PCA(1,3,"bottomleft",plot.cex=0.7,plot.pt.cex=1,plot.ncol=2)
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
