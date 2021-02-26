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
parser$add_argument("--covPerBinFiles"  , required=TRUE , help="file listing covPerBin file names")
parser$add_argument("--chrSizeFile"     , required=TRUE , help="chr size file")
parser$add_argument("--highLowPeakThresh" , nargs="+" , type="character" , help="OPTION: Provide two numbers. Bins > num1 or < num2 will be labeled. [default %(default)s]" , default= "NA" )
parser$add_argument("--minMAPQ" , type="integer"   , help="OPTION: if --inFormat is covPerBin filter out bins MAPQ < --minMAPQ (recommended 50) [default %(default)s]", default=0 )
parser$add_argument("--outName" ,  help="out name [default %(default)s]" , default="outPlotCovPerBin")
parser$add_argument("--debug"  , action="store_true" , help="dump session and quit [default %(default)s]" , default=FALSE)
parser$add_argument("--chrsToConsider" , nargs="+" , type="character" , help="list of chromosomes to be considered in the analysis [default %(default)s]" , default= "NA" )
args <- parser$parse_args()

#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }
if (! is.na(highLowPeakThresh[1])){highLowPeakThresh <- as.double(highLowPeakThresh)}

library(data.table)
library(fields)
if(debug){library(session);save.session("session_DEBUG");quit()}


#chr size cumulative sum
chrSize <- read.table(chrSizeFile,header=F,stringsAsFactors=F,col.names=c("chr","size"))
chrSize$chr <- as.character(chrSize$chr)
chrSize <- chrSize[chrSize$chr %in% chrsToConsider,]
chrSizeCumSum <- cumsum(chrSize$size)
v <- c(0,chrSizeCumSum)
chrSizeCumSum <- v[-length(v)]
names(chrSizeCumSum) <- chrSize$chr

#chr mid points
chrMids <- c()
for (i in 1:length(chrSizeCumSum)-1) {chrMids = c(chrMids , (chrSizeCumSum[i] + chrSizeCumSum[i + 1]) / 2)}
chrMids =c(chrMids ,  ( chrSizeCumSum[length(chrSizeCumSum)] + sum(chrSize$size) ) / 2)

#read covPerBin
covPerBinFilesV <- read.table(covPerBinFiles,stringsAsFactors=F)$V1
covPerBinList <- list()
for (covPerBin in covPerBinFilesV){
  print(covPerBin)
  s <-  gsub(x=basename(covPerBin),pattern=".covPerBin.gz$",replacement="")
  covPerBinList[[s]] <- fread(cmd=paste("gunzip -c", covPerBin ), select=c("chromosome","start","normalizedMeanCoverage","MAPQ"))
  positionCorrection <- chrSizeCumSum[match(covPerBinList[[s]]$chromosome, names(chrSizeCumSum))]
  covPerBinList[[s]]$startFixed <- covPerBinList[[s]]$start + positionCorrection
  covPerBinList[[s]]$sample <- s
}
covPerBinDf <- do.call(rbind,covPerBinList)
df <- covPerBinDf[with(covPerBinDf, order(startFixed)), ]

#filter
df <- df[ df$normalizedMeanCoverage <= 0.01 | (df$normalizedMeanCoverage > 0.01 & df$MAPQ > minMAPQ) ,]
df <- as.data.frame(df)

#add pseudocount just to deletions
df$normalizedMeanCoverage[df$normalizedMeanCoverage <= 0.1 ] = 0.1

#materializedPoints
materializedPoints <- df[df$normalizedMeanCoverage <= 0.5 | df$normalizedMeanCoverage >= 1.5,]
materializedPoints<- unique(materializedPoints)
maxPointsNum <- min(length(materializedPoints[,1]) , 50000)
set.seed(321)
materializedPoints <- materializedPoints[sample(nrow(materializedPoints), maxPointsNum), ]

#legend for smoothScatter
fudgeit <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  fields::image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)# , horizontal=T)
}

#plot smooth
png(paste0(outName,".smooth.png"),type='cairo' , width = 2000, height = 800,)
par(font.axis = 2 , mar = c(5,4,4,5) + .1)
smoothScatter(x=df$startFixed , y=log(df$normalizedMeanCoverage) , bandwidth=c(0.1,0.05) , nbin=400 , nrpoints=0 , xlab="chromosomes" , ylab="genomic bin coverage",xaxt = 'n' , postPlotHook = fudgeit , xaxs="i" , yaxs="i" , font=2)
points(materializedPoints$startFixed, log(materializedPoints$normalizedMeanCoverage), pch=19, col="black", cex=0.2)
segments(x0=0 , y0=log(1.5) , x1=sum(chrSize$size) , y1=log(1.5) , lty=3, col='red', lwd=1)
segments(x0=0 , y0=0 , x1=sum(chrSize$size) , y1=0 , lty=2, col='black', lwd=1)
segments(x0=0 , y0=log(0.5) , x1=sum(chrSize$size) , y1=log(0.5) , lty=3, col='red', lwd=1)
abline(v=chrSizeCumSum[-1], lty=1, col='gray', lwd=2)
axis(1, at=chrMids, labels=names(chrMids) , font=2)
dev.off()





