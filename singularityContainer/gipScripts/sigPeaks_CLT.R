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
parser <- ArgumentParser()
parser$add_argument("--input"    , type="character" , help="bin file, with coverage values normalized by median chr coverage (and GC content) [default %(default)s]")
parser$add_argument("--inFormat" , type="character" , help="format of the input bin file. Supperted are chr start end score (chrStartEndScore) or covPerBin [default %(default)s]" , default="covPerBin")
parser$add_argument("--pThresh"  , type="double"    , help="p-value threshold [default %(default)s]" , default=1 )
parser$add_argument("--padjust"  , type="character" , help="p-value multiple testing p.adjust correction method. You could use benjamini hochberg (BH), but since there is a dependency between bins (adjacent bins have more chances to be both amplified) then you should use Benjamini & Yekutieli (BY) [default %(default)s]", default="BY")
parser$add_argument("--coverageThresholds" , nargs="+" , type="character" , help="OPTION: Provide two numbers. Bins (and segments) will undergo an additional filter. Just bins (and segments) > num1 or < num2 will be selected [default %(default)s]" , default= "NA" )
parser$add_argument("--minMAPQ"      , type="integer"   , help="OPTION: if --inFormat is covPerBin filter out bins MAPQ < --minMAPQ (recommended 50) [default %(default)s]", default=0 )
parser$add_argument("--minLen"      , type="integer"   , help="OPTION: selected segments need to be at least --minLen nucleotides long [default %(default)s]", default=0 )
parser$add_argument("--outName"  , type="character" , help="out name [default %(default)s]" , default="sigPeaksOut")
parser$add_argument("--debug"  , action="store_true" , help="dump session and quit [default %(default)s]" , default=FALSE)
args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }

if (! is.na(coverageThresholds[1])){
	coverageThresholds <- as.double(coverageThresholds)
}
library(GenomicRanges)

if(debug){library(session);save.session("session_sigPeaksCLT");quit()}

############
#read input#
############
if(inFormat == "chrStartEndScore"){
	df        <- read.table(input,header=F, stringsAsFactors=F)
	names(df) <- c("chr","start","end","score")
} else if (inFormat == "covPerBin"){
	df <- read.table(input,header=T, stringsAsFactors=F)
	df <- df[df$MAPQ >= minMAPQ ,] #filter
	df <- df[,c("chromosome" , "start" , "end" , "normalizedMeanCoverage")]
	names(df) <- c("chr","start","end","score")
    df$chr <- as.character(df$chr)
} else {
	stop("input inFormat not recognized")
	quit(save = "no", status = 1, runLast = FALSE)
}
n = df[1,"end"] - df[1,"start"] + 1

########
#mu, se#
########
mu <- mean(df$score)
SD <- sd(df$score)
######
#pval#
######
p.lowerTail  <- pnorm(df$score , mean=mu , sd=SD , lower.tail=TRUE)
p.higherTail <- 1 - p.lowerTail
leftTailObs  <- df$score < mu
df$pvalue    <- p.higherTail
df$pvalue[leftTailObs] = p.lowerTail[leftTailObs]
df$direction <- "amplification"
df$direction[leftTailObs] = "depletion"
##########
#pval adj#
##########
df$pvalueAdj <- p.adjust(df$pvalue,method=padjust)
########
#filter#
########
significant <- df[df$pvalueAdj < pThresh , ]
############################################################
#collapse adjacent significant bins with the same direction#
#then average their score#
############################################################
collapse <- function (inDf){
 gr = GRanges(inDf$chr, IRanges(c(inDf$start),c(inDf$end)),strand="*")
 mcols(gr)$score     <- inDf$score
 mcols(gr)$direction <- inDf$direction
 grRed <- reduce(gr)
 mcols(gr)$segment <- subjectHits(findOverlaps(gr, grRed))
 d <- as.data.frame(gr, stringsAsFactors=F)
 d <- subset( d, select = -c(strand, width) )
 names(d)<- gsub(x=names(d),pattern="seqnames",replacement="chr")
 s <- split(d,d$segment)
 redDs <- lapply(s,function(x){
    tmp <- data.frame(chr=as.character(x$chr[1]),start=min(x$start), end=max(x$end) , direction=x$direction[1] ,stringsAsFactors=F)
    tmp$score <- round(mean(x$score),digits=2) 
    return(tmp)
 })
 return (do.call(rbind,redDs))
}
ampBins <- significant[significant$direction == "amplification" , ]
depBins <- significant[significant$direction == "depletion" , ]
clpsAmp <- NULL
clpsDep <- NULL
if(length(ampBins[,1]) > 0){clpsAmp <- collapse(ampBins)}
if(length(depBins[,1]) > 0){clpsDep <- collapse(depBins)}
redD <- rbind(clpsAmp , clpsDep)
redD <- redD[with(redD, order(chr, start)), ]

################################
#user defined segment filtering#
################################
if (! is.na(coverageThresholds[1])){
	redD <- redD[redD$score > coverageThresholds[1] | redD$score < coverageThresholds[2] , ]
        significant <- significant[significant$score > coverageThresholds[1] | significant$score < coverageThresholds[2] , ]
}
redD <- redD[(redD$end - redD$start)+1 > minLen , ]

#add segment id
if(length(redD[,1]) > 0){
  redD$segmentID <- paste0("segment_" , 1:length(redD[,1]))
} else {
  redD$segmentID <- character(0)
}

write.table(x=redD,quote=F,sep="\t",col.names=T,row.names=F,file=paste0(outName,".segments.tsv"))
write.table(x=significant,quote=F,sep="\t",col.names=T,row.names=F,file=paste0(outName,".bins.tsv"))
system(paste0("gzip ",outName,".segments.tsv"))
system(paste0("gzip ",outName,".bins.tsv"))

#write stats
write.table("##DISTRIBUTIONS FEATURES",quote=F,col.names=F,sep="\t",append=F,row.names=F,file=paste0(outName,".stats"))
write.table(paste("bin size", n),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(paste("distribution mean", mu)      ,quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(paste("distribution standard deviation" , SD),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table("##FILTERING AND FALSE DISCOVERY RATE",quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(paste("significant bins with minMAPQ",minMAPQ, ", adjusted p-val threshold", pThresh ,"and correction strategy", padjust," that were retained after possible coverageThresholds filtering are", length(significant[,1])),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(paste("collapsed adjacent bins (i.e. segments) that are retained after possible coverageThresholds and a min segment length of",minLen,"are", length(redD[,1])),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
pdf(paste0(outName,".pdf"))
hist(df$pvalue,breaks=100,col="orange",main="p-value distribution")
hist(df$pvalueAdj,breaks=100,col="orange",main="p-value distribution\nmultiple test correction")
dev.off()
