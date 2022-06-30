#choix du bon répertoire d'analyse
setwd("C:/Users/omics/Desktop/DU/DU_Poster/DU_Poster")
#lecture du fichier de couverture pour le séquençage Target
coverage_target <- read.table("coverage-barcode01.bedgraph", sep = "\t")
#fichier couverture annoté avec les régions du transgène
coverage_target <- read.table("barcode.bedgraph", sep = "\t")
#ajout de nom de colonnes
colnames(coverage_target)<-c("vector","start","end","cov")
#même commande pour le fichier de couverture WGS
coverage_WGS <- read.table("coverage-R-19-153-5.bedgraph", sep = "\t")
colnames(coverage_WGS)<-c("vector","start","end","cov")
#Importation des annotations qui sont dans un fichier csv
annot<-read.csv("CD4_annot.csv", sep=";",header=TRUE)
#Controle de la couverture avec un graphique simple
plot(coverage_WGS$start,coverage_WGS$cov)
plot(coverage_target$start,coverage_target$cov)
#appel de la librairie
library(circlize)
#Nettoyage de la zone graphique
circos.clear()
annot_size<-annot$end-annot$start
#modification de la zone de graphique (recommandé par un message d'erreur)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
#Initialization du circos avec les annotations
circos.initialize(factors=annot$gene,xlim=matrix(c(rep(0, 9), annot_size), ncol=2))
#Affichage de la construction
circos.track(ylim = c(0, 1),panel.fun = function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr,cex=0.5)
})
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h="top",labels.cex=0.4)
})
circos.genomicTrack(data=coverage_target, panel.fun=function(region, value, ...) {
  circos.genomicLines(region, value, type="l", col="grey50", lwd=0.6)
  circos.segments(x0=0, x1=max(annot_size), y0=100, y1=100, lwd=0.6, lty="11", col="grey90")
  circos.segments(x0=0, x1=max(annot_size), y0=300, y1=300, lwd=0.6, lty="11", col="grey90")
  #circos.segments(x0=0, x1=max(ref$V2), y0=500, y1=500, lwd=0.6, lty="11", col="grey90")
}, track.height=0.08, bg.border=F)
circos.yaxis(at=c(100, 300), labels.cex=0.25, lwd=0, tick.length=0, col="#FFFFFF")