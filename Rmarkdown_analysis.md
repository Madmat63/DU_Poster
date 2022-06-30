---
title: "DU_Poster"
author: "Mathieu Chicard"
date: "28/06/2022"
output:  
  html_document:
    toc: yes
    toc_float: yes
    keep_md: yes
    fig_caption: yes
  pdf_document:
    toc: yes
editor_options:
  chunk_output_type: inline
---



Ce script nécessite l'installion des packages suivants: 
<ul>
<li>circlize</li>
<li>RColorBrewer</li>
</ul>

Utilisation des sites :
<ul>
<li><a href="https://jokergoo.github.io/circlize_book/book/">Circular Visualization in R</a></li>
<li><a href="https://www.royfrancis.com/beautiful-circos-plots-in-r/">Beautiful circos plots in R</a></li>
</ul>

## Importation des données

Dans un premier temps les données de WGS et Target provenant de séquençage Nanopore vont être importés. Un fichier d'annotation des séquences est également importés. La couverture est dans un premier temps vérifier dans un simple <I>plot</I>


```r
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
```

![](Rmarkdown_analysis_files/figure-html/Importation des données-1.png)<!-- -->

```r
plot(coverage_target$start,coverage_target$cov)
```

![](Rmarkdown_analysis_files/figure-html/Importation des données-2.png)<!-- -->

## Création d'un circos

Dans un premier temps, il faut créer le circos avec les annotations des régions. Le découpage des régions doit être initialisé par <i>circos.initialize()</i> puis afficher avec <i>circos.track()</i> qui permet également de jouer avec les paramêtres graphiques



```r
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
```

![](Rmarkdown_analysis_files/figure-html/Creation d un circos-1.png)<!-- -->

## Ajout de coordonnées
Pour mieux se répérer et situé des "événements", des coordonnées peuvent être ajouter une nouvelle fois par la fonction <i>circos.track()<i> 


```r
#Création à nouveau du circos
library(circlize)
circos.clear()
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors=annot$gene,xlim=matrix(c(rep(0, 9), annot_size), ncol=2))
circos.track(ylim = c(0, 1),panel.fun = function(x, y) {
chr=CELL_META$sector.index
xlim=CELL_META$xlim
ylim=CELL_META$ylim
circos.text(mean(xlim), mean(ylim), chr,cex=0.5)
})

#Ajout de coordonnées
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
circos.axis(h="top",labels.cex=0.4)
})
```

![](Rmarkdown_analysis_files/figure-html/Customisation du circos-1.png)<!-- -->

## Ajout de la couverture
La couverture provenant des sorties <b>bedtools</b> peut maintenant être ajouter avec la fonction <i>circos.genomicTrack</i> et sous forme de ligne en ajoutant <i>circos.genomicLines</i>. Il faut également associé les valeurs de couvertures à l'axe Y avec <i>circos.yaxis()</i>




```r
#Création à nouveau du circos
library(circlize)
circos.clear()
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors=annot$gene,xlim=matrix(c(rep(0, 9), annot_size), ncol=2))
circos.track(ylim = c(0, 1),panel.fun = function(x, y) {
chr=CELL_META$sector.index
xlim=CELL_META$xlim
ylim=CELL_META$ylim
circos.text(mean(xlim), mean(ylim), chr,cex=0.5)
})
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
circos.axis(h="top",labels.cex=0.4)
})
#ajout de la couverture  dans les paramêtres du circos
circos.genomicTrack(data=coverage_target, panel.fun=function(region, value, ...) {
circos.genomicLines(region, value)
})
#Choix de l'axe d'affichage
circos.yaxis()
```

![](Rmarkdown_analysis_files/figure-html/Ajout de la couverture-1.png)<!-- -->

## Ravalement de façade
Pour terminer, les textes apparaisant par défaut avec le tracé de couverture vont être supprimé ainsi que les bordures des zones de traçage pour plus de lisibilité


```r
#Création à nouveau du circos
library(circlize)
circos.clear()
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors=annot$gene,xlim=matrix(c(rep(0, 9), annot_size), ncol=2))
circos.track(ylim = c(0, 1),panel.fun = function(x, y) {
chr=CELL_META$sector.index
xlim=CELL_META$xlim
ylim=CELL_META$ylim
circos.text(mean(xlim), mean(ylim), chr,cex=0.5)
})
circos.genomicTrack(data=coverage_target, panel.fun=function(region, value, ...) {
circos.genomicLines(region, value, type="l", col="grey50", lwd=0.6)
circos.segments(x0=0, x1=max(annot_size), y0=100, y1=100, lwd=0.6, lty="11", col="grey90")
circos.segments(x0=0, x1=max(annot_size), y0=300, y1=300, lwd=0.6, lty="11", col="grey90")
#circos.segments(x0=0, x1=max(ref$V2), y0=500, y1=500, lwd=0.6, lty="11", col="grey90")
}, track.height=0.08, bg.border=F)
circos.yaxis(at=c(100, 300), labels.cex=0.25, lwd=0, tick.length=0, col="#FFFFFF")
```

![](Rmarkdown_analysis_files/figure-html/Embelissement-1.png)<!-- -->

## à venir
Le circos est fonctionnelle et pourra permettre de comparer plusieurs échantillons entre eux. Ils manquent également un peu de couleurs.Pour patienter, voici quelques couleurs


```r
library("RColorBrewer")
```

```
## Warning: le package 'RColorBrewer' a été compilé avec la version R 4.1.3
```

```r
display.brewer.all()
```

![](Rmarkdown_analysis_files/figure-html/Color-1.png)<!-- -->

