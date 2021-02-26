#############################################################################
# giptools                                                                  #
#                                                                           #
# Authors: Giovanni Bussotti                                                #
# Copyright (c) 2021  Institut Pasteur                                      #
#                                                                           #
#                                                                           #
# This program is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as                   #
# published by the Free Software Foundation, either version 3 of the        #
# License, or (at your option) any later version.                           #
#                                                                           #
# This program is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.    #
#                                                                           #
#############################################################################

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='Example: Rscript karyoplotCovPerGe.R --covPerGe lsdOut/P135_C7BDUACXX.covPerGe.gz --covPerBin lsdOut/P135_C7BDUACXX.covPerBin.gz --chrSize /pasteur/projets/policy01/BioIT/Giovanni/datasets/assemblies/leishmanias/peterMyler/Leishmania_donovani_sudanese.chrSize --CHRS 1 2 3 4 5 6 7 8 9 0 11 12 13 14 15 --REPS /pasteur/projets/policy01/BioIT/Giovanni/datasets/assemblies/leishmanias/peterMyler/repeatMasker/Leishmania_donovani_sudanese/Leishmania_donovani_sudanese.fa.out.gff --amplificationThreshold 2 --depletionThreshold 0.5')
parser$add_argument("--covPerGe"  , required=TRUE , help="covPerGe file name")
parser$add_argument("--covPerBin" , required=TRUE , help="covPerGe file name")
parser$add_argument("--chrSize"   , required=TRUE , help="chr size file")
parser$add_argument("--REPS"   , required=TRUE , help="repeat masker out.gff file")
parser$add_argument("--CHRS" , nargs="+" ,  type="character" , help="space separated chromosomes identifiers to consider [default %(default)s]" , default= "NA")
parser$add_argument("--significant" ,  help="significant genes file returned by sigPeaks_mixture.R [default %(default)s]" , default="NA")
parser$add_argument("--amplificationThreshold"  , help="genes with normalized coverage > than amplificationThreshold are labeled as amplified (e.g. 2). Ignored if significant file is provided, or if depletionThreshold is NA [default %(default)s]" , default= "NA" )
parser$add_argument("--depletionThreshold" , help="genes with normalized coverage < than depletionThreshold are labeled as depleted (e.g. 0.5). Ignored if significant file is provided, or if amplificationThreshold is NA [default %(default)s]" , default= "NA" )
parser$add_argument("--outDir" ,  help="out dir [default %(default)s]" , default="covPerGeAnalysisOut")
parser$add_argument("--repeatRange" , type="integer" ,  help="repeats within repeatRange from amplified or depleted genes will be showed [default %(default)s]" , default=1000)
parser$add_argument("--minMAPQ" , type="integer"   , help="genes with MAPQ < this value are labeled [default %(default)s]", default=50 )
parser$add_argument("--outName" ,  help="out name [default %(default)s]" , default="outPlotCovPerBin")
parser$add_argument("--geneFunction" , help="file listing geneID<TAB>function [default %(default)s]" , default= "NA" )
parser$add_argument("--debug"  , action="store_true" , help="dump session and quit [default %(default)s]" , default=FALSE)
args <- parser$parse_args()

#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }

if(! is.na(amplificationThreshold)){amplificationThreshold <- as.numeric(amplificationThreshold)}
if(! is.na(depletionThreshold)){depletionThreshold <- as.numeric(depletionThreshold)}

library(karyoploteR)

if(debug){library(session);save.session("session_DEBUG_karyoplotCovPerGe");quit();}

######
#READ#
######
#covPerBin
covPerBinDf <- read.table(covPerBin,header=T,stringsAsFactors=F,sep="\t")
covPerBinDf$chromosome <- as.character(covPerBinDf$chromosome)
covPerBinDf$mid <- covPerBinDf$start + ( (covPerBinDf$end - covPerBinDf$start) /2)
#covPerGe
covPerGeDf <- read.table(covPerGe,stringsAsFactors=F,header=T)
covPerGeDf$chr   <- gsub(x=covPerGeDf$locus , pattern=":.+" ,replacement="")
covPerGeDf$start <- as.numeric(gsub(x=covPerGeDf$locus , pattern=".+:(.+)-.+" ,replacement="\\1"))
covPerGeDf$end   <- as.numeric(gsub(x=covPerGeDf$locus , pattern=".+:.+-(.+)" ,replacement="\\1"))
covPerGeDf$mid <- covPerGeDf$start + ( (covPerGeDf$end - covPerGeDf$start) /2)
#chrSize
chrSizeDf <- read.table(chrSize,stringsAsFactors=F,header=F,col.names=c("chr","size"))
#reps
repsDf <- read.table(REPS,stringsAsFactors=F)[,c(1,4,5,7,10)]
names(repsDf) <- c("chr","start","end","strand","labels")

################
#FILTER by CHRS#
################
if(! is.na(CHRS[1])){
 covPerBinDf <- covPerBinDf[covPerBinDf$chromosome %in% CHRS,]
 covPerGeDf  <- covPerGeDf[covPerGeDf$chr %in% CHRS,]
 chrSizeDf   <- chrSizeDf[chrSizeDf$chr %in% CHRS,]
 repsDf      <- repsDf[repsDf$chr %in% CHRS,]
}

#####################
#LABEL AMP/DEP GENES#
#####################
covPerGeDf$status <- "background"
covPerGeDf$color  <- "black"
amplifiedGe <- c()
depletedGe <- c()
if(! is.na(significant)){
  significantDf <- read.table(significant,header=T,stringsAsFactors=F) 
  significantDf$direction <- "depleted"
  significantDf[ significantDf$probOfTheValueOrLower > significantDf$probOfTheValueOrHigher , "direction" ] <- "amplified"
  amplifiedGe <- row.names(significantDf[ significantDf$direction == "amplified" , ])
  depletedGe  <- row.names(significantDf[ significantDf$direction == "depleted"  , ])
} else if (! is.na(amplificationThreshold) && ! is.na(depletionThreshold)) {
  amplifiedGe <- covPerGeDf[ covPerGeDf$normalizedMeanCoverage > amplificationThreshold , "gene_id" ] 
  depletedGe <- covPerGeDf[ covPerGeDf$normalizedMeanCoverage < depletionThreshold , "gene_id" ]
}
covPerGeDf[ match(amplifiedGe , covPerGeDf$gene_id) , "status"] <- "amplified"
covPerGeDf[ match(amplifiedGe , covPerGeDf$gene_id) , "color"]  <- "orange"
covPerGeDf[ match(depletedGe , covPerGeDf$gene_id) , "status"]  <- "depleted"
covPerGeDf[ match(depletedGe , covPerGeDf$gene_id) , "color"]   <- "blue"

######################
#LABEL LOW MAPQ GENES#
######################
lowmapqGe <- covPerGeDf[covPerGeDf$MAPQ < minMAPQ,"gene_id"]
covPerGeDf[ match( lowmapqGe , covPerGeDf$gene_id) , "status"] <- "lowMAPQ"
covPerGeDf[ match( lowmapqGe , covPerGeDf$gene_id) , "color"] <- "gray"

#########
#GRANGES#
#########
genomeGr <- toGRanges(data.frame(chr=chrSizeDf$chr, start=1, end=chrSizeDf$size))
repsGr <- toGRanges(repsDf)
strand(repsGr) <- repsGr$strand
mcols(repsGr)$strand <- NULL


######
#PLOT#
######
plotKP <- function(chr,ymax,repeatRange){
 #ymax should be equal to the value you have in the top tick
 system(paste("mkdir -p", outDir))
 png(paste0(outDir,"/",chr,".png") , width = 3000 , height = 1500 , type='cairo')#, width = 2500, height = 1000 , res=700

 #ideogram
 pp <- getDefaultPlotParams(plot.type=3)
 pp$data2height <- 20
 pp$data1height <- 300
 #pp$leftmargin  <- 0.1
 kp <- plotKaryotype(genome = genomeGr , plot.type=3, main="" , chromosomes=c(chr) , plot.params = pp)
 #kpDataBackground(kp, data.panel = 1)

 #select data
 tmpGeDf  <- covPerGeDf[covPerGeDf$chr == chr , ]
 tmpBinDf <- covPerBinDf[covPerBinDf$chromosome == chr , ]
 #saturate
 tmpGeDf[tmpGeDf$normalizedMeanCoverage > ymax,"normalizedMeanCoverage"]=ymax
 tmpBinDf[tmpBinDf$normalizedMeanCoverage > ymax ,"normalizedMeanCoverage"]=ymax
 #subsets
 geToMarkDf     <- tmpGeDf[ tmpGeDf$status %in% c("amplified","depleted") , ]
 backgroundGeDf <- tmpGeDf[ tmpGeDf$status == "background" , ]
 lowmapqGeDf <- tmpGeDf[ tmpGeDf$status == "lowMAPQ" , ]

 #genomic ranges
 geToMarkGr <- GRanges( seqnames=geToMarkDf$chr , ranges = IRanges(geToMarkDf$start, end = geToMarkDf$end, names = geToMarkDf$gene_id ))
 tmpBinGr <- GRanges( seqnames=tmpBinDf$chromosome , ranges = IRanges(tmpBinDf$start, end = tmpBinDf$end), mid=tmpBinDf$mid , normalizedMeanCoverage=tmpBinDf$normalizedMeanCoverage)
 tmpGeGr <- GRanges( seqnames=tmpGeDf$chr , ranges = IRanges(tmpGeDf$start, end = tmpGeDf$end, names = tmpGeDf$gene_id ), mid=tmpGeDf$mid, normalizedMeanCoverage=tmpGeDf$normalizedMeanCoverage)

 #remove bin overlapping genes
 ov <- findOverlaps(query=tmpBinGr, subject=tmpGeGr)
 tmpBinGr <- tmpBinGr[-queryHits(ov),]

 #add to binGr the genes
 tmpBinGr <- c(tmpBinGr,tmpGeGr)
 tmpBinGr <- tmpBinGr[order(tmpBinGr$mid),]

 #repeats to show must be in repeatRange  
 geToMarkGr <- geToMarkGr + repeatRange
 repsToShow <- repsGr[subjectHits(findOverlaps(geToMarkGr,repsGr)),]
 geToMarkGr <- geToMarkGr - repeatRange
 #and just the closest preceding and following the significant gene
 repsToShowIndex <- c(precede(geToMarkGr,repsToShow) , follow(geToMarkGr,repsToShow))
 repsToShowIndex <- repsToShowIndex[! is.na(repsToShowIndex)]
 if (length(repsToShowIndex) > 0) {
  repsToShow <- repsToShow[repsToShowIndex,]
 }

 #plot
 points.top <- 0.8
 kpAxis(kp, data.panel=1 , cex=2 , lwd=2 , numticks = 4, tick.pos = c(0, 0.25, 0.5, 0.75 , 1), labels = c(0, ymax/4, ymax/2, ymax/2 + ymax/4, ymax) , r1=points.top)
 kpAddLabels(kp, cex=2 ,labels="normalized coverage", srt=90, pos=3 , label.margin = 0.04 , data.panel = 1 , r1=points.top)
 #kpLines(kp, chr=chr, x=tmpBinDf$mid, y=tmpBinDf$normalizedMeanCoverage,col="gray",ymax=ymax , r1=points.top)
 kpArea(kp, chr=chr, x=tmpBinGr$mid, y=tmpBinGr$normalizedMeanCoverage,col="gray",ymax=ymax , r1=points.top )
 kpPoints(kp, chr=chr, x=backgroundGeDf$mid, y=backgroundGeDf$normalizedMeanCoverage,pch=19,col=backgroundGeDf$color,lwd=3,ymax=ymax , r1=points.top)
 kpPoints(kp, chr=chr, x=geToMarkDf$mid, y=geToMarkDf$normalizedMeanCoverage,pch=19,col=geToMarkDf$color, lwd=8 , ymax=ymax , r1=points.top)
 kpPoints(kp, chr=chr, x=lowmapqGeDf$mid, y=lowmapqGeDf$normalizedMeanCoverage,pch=19,col=lowmapqGeDf$color, lwd=3 , ymax=ymax , r1=points.top)
 if (length(geToMarkDf[,1]) > 0) { 
  kpPlotMarkers(kp, chr=geToMarkDf$chr, x=geToMarkDf$mid, labels=geToMarkDf$gene_id , r0=points.top , r1=0.9 ,
  	marker.parts = c(0.01,0.01, 0.01), text.orientation = "vertical" , max.iter=1000 , label.dist = 0.005, ignore.chromosome.ends=TRUE , cex=2)
  kpSegments(kp, chr=chr, x0=geToMarkDf$mid, x1=geToMarkDf$mid, y0=geToMarkDf$normalizedMeanCoverage, y1=ymax, ymax=ymax,  r1=points.top)
 }
 if (length(repsToShow) > 0) { 
  kpPlotRegions(kp, repsToShow, col="black", avoid.overlapping=F, r0=0 , r1=0.1 , data.panel=2 )
  kpPlotMarkers(kp, repsToShow , r0=0.1  , r1=0.3 , data.panel=2 , 
  	            max.iter=1000 , label.dist = 0.004 , ignore.chromosome.ends=TRUE, cex=1.5)
 }
 dev.off()
}
for (chr in unique(covPerGeDf$chr)){
  print(chr)
  plotKP(chr,4,repeatRange)
}

#######
#WRITE#
#######
covPerGeDf[,c("chr" , "start" , "end" , "mid")] <- NULL
if(! is.na(geneFunction)){
 geFunDf <- read.delim(geneFunction,header=F,stringsAsFactors=F ,sep="\t",col.names=c("geId","annotation"))
 covPerGeDf$annotation <- geFunDf[ match(covPerGeDf$gene_id , geFunDf$geId) , "annotation" ]
}
write.table(x=covPerGeDf,file=paste0(outDir,"/covPerGePlotData.tsv"),col.names=T,row.names=F, sep="\t")


