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
parser$add_argument("--testFile" , nargs="+", help="test gcovbin file. The file must be gzip or bgzip compressed  [default %(default)s]" )
parser$add_argument("--referenceFile" , nargs="+", help="reference gcovbin file. The file must be gzip or bgzip compressed [default %(default)s]" )
parser$add_argument("--testName" , help="test name to use internally [default %(default)s]")
parser$add_argument("--referenceName" , help="reference name to use internally [default %(default)s]")
parser$add_argument("--outName" , help="out name [default %(default)s]" , default="_outTestGcovbin_")
parser$add_argument("--ylim" , help="graphical parameter. max ratio value before saturation. Rations above ylim will be replaced with ylim in the plot (while they will be unaltered in the output dataframe) [default %(default)s]" , default="3")
parser$add_argument("--chrs" , nargs="+", help="List of chromosome names to plot. This also defines the plotting order[default %(default)s]" )
parser$add_argument("--minMAPQ"  , help="points with MAPQ < minMAPQ won't be shown. It also influence the output stats table (recommended 50) [default %(default)s]" , default="-1")
parser$add_argument("--minWindowNormMean", help="points with both meanCoverageReference and meanCoverageTest < minWindowNormMean won't be shown. At least one of the two must be > minWindowNormMean. It also influence the output stats table (recommended 0.1) [default %(default)s]" , default="-1")
parser$add_argument("--minWindowNormMedian", help="points with both medianCoverageReference and medianCoverageTest < minWindowNormMedian won't be shown. At least one of the two must be > minWindowNormMedian. It also influence the output stats table (recommended 0.1) [default %(default)s]" , default="-1")
parser$add_argument("--highRatio" , help="graphical parameter. points with ratio above highRatio will be colored differently. It also influence the output stats table [default %(default)s]" , default="1.5")
parser$add_argument("--lowRatio"  , help="graphical parameter. points with ratio below lowRatio  will be colored differently. It also influence the output stats table [default %(default)s]" , default="0.67")
parser$add_argument("--segmentation"  , help="graphical parameter. add unbalhaar segmentation to the plot [default %(default)s]" , default="no")
parser$add_argument("--divideByOne" , action="store_true" , help="set this flag to divide by one insted by the reference. It is a trick to visualize just the test sample coverage [default %(default)s]" , default=FALSE)

args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }

#correct type
ylim                <- as.numeric(ylim)
minMAPQ             <- as.numeric(minMAPQ)
chrs                <- as.character(chrs)
highRatio           <- as.numeric(highRatio)
lowRatio            <- as.numeric(lowRatio)
minWindowNormMean   <- as.numeric(minWindowNormMean)
minWindowNormMedian <- as.numeric(minWindowNormMedian)

#function
printStatsOnTable <- function(feature,df){
	#print stats of highRatio and lowRatio genomic windows and segments 
	totGenomicWindow <- length(df$coverageRatio)
	totGenomicWindowHR <- length(which(df$coverageRatio > highRatio))
	totGenomicWindowLR <- length(which(df$coverageRatio < lowRatio))
	totSegment <- length(unique(df$segmentID))
	totSegmentHR <- length(df[df$segment > highRatio,"segmentID"])
	totSegmentLR <- length(df[df$segment < lowRatio,"segmentID"])
	maxGenomicWindowRatio <- max(df$coverageRatio)
	minGenomicWindowRatio <- min(df$coverageRatio)
	maxSegment			  <- max(df$segment)
	minSegment			  <- min(df$segment)
	c(feature , totGenomicWindow, totGenomicWindowHR, totGenomicWindowLR, totSegment, totSegmentHR , totSegmentLR , maxGenomicWindowRatio , minGenomicWindowRatio , maxSegment , minSegment)
}
printBed <- function(df){
	#print the bed coverage ratio above or below highRatio and lowRatio
	#genomic windows 
	allBed <- NULL
	if(any(df$coverageRatio > highRatio)){
		genomicWindowHRbed <- df[which(df$coverageRatio > highRatio),][,c("chromosome","start","end")]
		genomicWindowHRbed$id <- paste("genomicWindowHR",row.names(genomicWindowHRbed),sep="_")
		genomicWindowHRbed$covRatio <- round(df[which(df$coverageRatio > highRatio),][,c("coverageRatio")],digits=2)
		allBed <- rbind(genomicWindowHRbed , allBed)
	}
	if(any(df$coverageRatio < lowRatio)){
		genomicWindowLRbed <- df[which(df$coverageRatio < lowRatio),][,c("chromosome","start","end")]
		genomicWindowLRbed$id <- paste("genomicWindowLR",row.names(genomicWindowLRbed),sep="_")
		genomicWindowLRbed$covRatio <- round(df[which(df$coverageRatio < lowRatio),][,c("coverageRatio")],digits=2)
		allBed <- rbind(genomicWindowLRbed , allBed)
	}
	#segments
	segments <- split(df,f=df$segmentID)
	for (seg in names(segments)){   
		if(segments[[seg]]$segment > highRatio){
			chromosome = segments[[seg]]$chromosome[1]
			start = min(segments[[seg]]$start)
			end   = max(segments[[seg]]$end)
			segmentHRbed <- data.frame(chromosome=chromosome,start=start,end=end,id=paste("segmentHR",seg,sep="_"),covRatio=round(segments[[seg]]$segment,digits=2))
			allBed <- rbind(segmentHRbed , allBed)
		}
	
		if(segments[[seg]]$segment < lowRatio){
			chromosome = segments[[seg]]$chromosome[1]
			start = min(segments[[seg]]$start)
			end   = max(segments[[seg]]$end)
			segmentLRbed <- data.frame(chromosome=chromosome,start=start,end=end,id=paste("segmentLR",seg,sep="_"),covRatio=round(segments[[seg]]$segment,digits=2))
			allBed <- rbind(segmentLRbed , allBed)
		}
	}
	
	if (! is.null(allBed)){
		allBed <- as.data.frame(allBed[order(allBed[,1], allBed[,2]) , ])
		write.table(x=allBed,quote=F,sep="\t",append=F,col.names=F,row.names=F,file=paste(outName,".extremeRatio.bed",sep=""))
		system( paste("gzip ",outName,".extremeRatio.bed",sep="") )
	}
}
addSegmentId <- function(df,chrs){
	df$segment <- uh(coverageRatio)
	segmentID <-c() ;  
	for (chr in chrs){
		dfChr <- subset(df,chromosome == chr);   
		totSegments <- length(rle(dfChr$segment)$lengths) ; 
		chrSegNames <- paste("chr",chr,"seg",1:totSegments,sep="")
		segLengths  <- rle(dfChr$segment)$lengths
		segmentID   <- c(segmentID , rep(chrSegNames , segLengths))
	}
	df$segmentID <- factor(segmentID)
	return(df)
}

#library(session)
#save.session("session")
#quit()

library(data.table)
library(ggplot2)
library(unbalhaar)
options(datatable.fread.input.cmd.message=FALSE)

#built coverageRation df
test      <- fread(paste("gunzip -c", testFile),stringsAsFactors=FALSE , header=T)
reference <- fread(paste("gunzip -c", referenceFile),stringsAsFactors=FALSE , header=T)
names(test)      <- c("chromosome","start","end","meanCoverageTest","medianCoverageTest","MAPQTest")
names(reference) <- c("chromosome","start","end","meanCoverageReference","medianCoverageReference","MAPQReference")

#divideByOne specific configurations
ylab <- paste("coverage ratio\n",testName,"/",referenceName)
if(divideByOne){
	reference$meanCoverageReference   <- 1
	reference$medianCoverageReference <- 1
	referenceName <- "noReference"
	ylab <-  paste("coverage",testName)
}

test$testName           <- testName
reference$referenceName <- referenceName
test      <- as.data.frame(test , stringsAsFactors = FALSE)
reference <- as.data.frame(reference , stringsAsFactors = FALSE)
row.names(test)      <- paste(test$chromosome, test$start, test$end,sep="_")
row.names(reference) <- paste(reference$chromosome, reference$start, reference$end,sep="_")
df <- cbind(test, reference[row.names(test),c("meanCoverageReference","medianCoverageReference","referenceName","MAPQReference")])
#select
df <- df[df$chromosome %in% chrs, ]
df <- df[df$MAPQReference >= minMAPQ & df$MAPQTest >= minMAPQ,]
df <- df[df$meanCoverageTest   >= minWindowNormMean  | df$meanCoverageReference >= minWindowNormMean,]
df <- df[df$medianCoverageTest >= minWindowNormMedian | df$meanCoverageReference >= minWindowNormMedian,]
#reformat
newDf<-NULL
for (chr in chrs){rbind(newDf,df[df$chromosome == chr,]) -> newDf}; 
df <- newDf
df$chromosome <- factor(df$chromosome , levels=chrs)
#coverage and position
df$coverageRatio <- df$meanCoverageTest / df$meanCoverageReference
df$windowMidPoint <- (df$start) + (((df$end - df$start)+1) /2)
#add color field
df$color <- "black"
df$color[df$coverageRatio > highRatio]="#E69F00"
df$color[df$coverageRatio < lowRatio]="#56B4E9"
df$color <-factor(df$color,levels=c("black","#E69F00","#56B4E9"))
df$status <- "normal"
df$status[df$coverageRatio > highRatio]="amplified"
df$status[df$coverageRatio < lowRatio]="depleted"
df$status <-factor(df$status,levels=c("normal","amplified","depleted"))
#coverage ratio
coverageRatio = df$coverageRatio
coverageRatio[is.na(coverageRatio)]=1
coverageRatio[coverageRatio == Inf]=ylim
#add segmentation field
df <- addSegmentId(df,chrs)

write.table(append=F,col.names=T,file=paste(outName,".df",sep=""),quote=F,sep="\t",row.names=T,x=df)
system( paste("gzip ",outName,".df",sep="") )
printBed(df)

#plot all together
dfAll <- data.frame(coverageRatio=df$coverageRatio , window=1:length(df$coverageRatio))
dfAll$coverageRatio[dfAll$coverageRatio > ylim] = ylim
segment_data = data.frame( x = cumsum(rle(as.character(df$chromosome))[["lengths"]]) , xend = cumsum(rle(as.character(df$chromosome))[["lengths"]]) , y=rep(0, length(unique(df$chromosome))) , yend=rep(ylim, length(unique(df$chromosome))))
#pdf(paste(outName,"_all.pdf",sep=""),height=10, width=20) 
png(paste0(outName,".all.png") ,type='cairo')
p <- ggplot(dfAll, aes(window,coverageRatio)) + geom_point(colour="black", size = 2, alpha = 0.1) + coord_cartesian(ylim =c(0,ylim)) + theme(legend.position="none") + xlab("genomic bin index") + ggtitle("all chromosomes") + theme_bw() + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=segment_data) + theme(legend.position="none") + ylab(ylab)
print(p)
dev.off()

#plot by chromosome
tableStats <- NULL
pdf(paste(outName,"_byChr.pdf",sep="")) 
for (chr in chrs){
	dfChr <- subset(df,chromosome == chr)
	dfChr$coverageRatio[dfChr$coverageRatio > ylim]=ylim
	p <- ggplot(dfChr, aes(windowMidPoint,coverageRatio)) + geom_point(colour="grey50", size = 2, alpha = 0.1) + coord_cartesian(ylim =c(0,ylim)) + theme(legend.position="none") + xlab("position") + ggtitle(paste("chromosome",chr)) + theme_bw() + ylab(ylab)
	if(segmentation == "yes"){
		p = p + geom_point(size=2,alpha = 0.3,colour="black")  + geom_point(aes(x=windowMidPoint,y=segment),colour="red",size=1,alpha=1)
	}	else {
		p = p + geom_point(size=1.5,alpha = 0.3,aes(colour=status)) + scale_colour_manual(values = c("normal" = "black","amplified" = "#E69F00","depleted" = "#56B4E9"),name="coverage ratio")
	}
	print(p)

	rbind(printStatsOnTable(chr,dfChr),tableStats) -> tableStats
}
dev.off()

#plot faceting
png(paste0(outName,".faceting.png"),height=400, width=1000,type='cairo')   #PNG works fine and is lighter than PDF, but it does not work on the Singularity implementation, so I use PDF for consistency
#pdf(paste(outName,"_faceting.pdf",sep=""),height=10, width=20)
dfFaceting <- df
dfFaceting$coverageRatio[dfFaceting$coverageRatio > ylim]=ylim
p <- ggplot(dfFaceting, aes(windowMidPoint,coverageRatio)) + geom_point(colour="grey50", size = 1, alpha = 0.1) +  facet_wrap(~ chromosome,nrow=4,scales="free") + coord_cartesian(ylim =c(0,ylim)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("position") + ylab(ylab) + scale_x_continuous(labels = function(x) format(x, scientific = FALSE))
p = p + geom_point(size=0.5,alpha=0.5,aes(colour=status)) + scale_colour_manual(values = c("normal" = "black","amplified" = "#E69F00","depleted" = "#56B4E9"))
print(p)
dev.off()

#print table stats and bed
rbind(printStatsOnTable("all",df),tableStats) -> tableStats
tableStats  <- as.data.frame(tableStats)
names(tableStats) <- c("chromosome","totGenomicWindow","totGenomicWindowHR","totGenomicWindowLR","totSegment","totSegmentHR","totSegmentLR","maxGenomicWindowRatio","minGenomicWindowRatio","maxSegment","minSegment")
write.table(x=tableStats,append=F,quote=F,sep="\t",row.names=F,file=paste(outName,".extremeRatio.stats",sep=""))



###############
#OLD FUNCTIONS#
###############

#plot faceting. Perfectly working. The only snag is that the faceting pdf are heavy.
#pdf(paste(outName,"_faceting.pdf",sep=""),height=10, width=20)
#dfFaceting <- df
#dfFaceting$coverageRatio[dfFaceting$coverageRatio > ylim]=ylim
#p <- ggplot(dfFaceting, aes(windowMidPoint,coverageRatio)) + geom_point(colour="grey50", size = 1, alpha = 0.1) +  facet_wrap(~ chromosome,nrow=4) + coord_cartesian(ylim =c(0,ylim)) + theme(legend.position="none") + xlab("windows position") + ylab(paste("coverage ratio\n",testName,"/",referenceName)) + ggtitle("all chromosomes") + theme_bw()
#p = p + geom_point(size=0.5,alpha=0.5,aes(colour=color)) + scale_colour_manual(values=levels(df$color),name="coverage ratio",labels=c("background", "enriched", "depleted")  )
#print(p)
#dev.off()




