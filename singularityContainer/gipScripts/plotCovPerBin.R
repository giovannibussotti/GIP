suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--covPerBin" , required="TRUE" , help="covPerBin file. The file must be gzip or bgzip compressed  [default %(default)s]" )
parser$add_argument("--chrSizeFile"     , required=TRUE , help="chr size file")
parser$add_argument("--outName" , help="out name [default %(default)s]" , default="_outTestGcovbin_")
parser$add_argument("--ylim" , type="double" , nargs="+" , help="min and max value shown on the y-axis. Values > ylim will be replaced with ylim in the plot (while they will be unaltered in the output dataframe) (default: 0 3)  [default %(default)s]" , default="NA")
parser$add_argument("--chrs" , nargs="+", help="List of chromosome names to plot. This also defines the plotting order[default %(default)s]" )
parser$add_argument("--minMAPQ" , type="integer" , help="points with MAPQ < minMAPQ will be shown in gray [default %(default)s]" , default="0")
parser$add_argument("--binOverviewSize" , nargs="+" ,type="integer" , help="height and width of all chromosomes plots [default %(default)s]" , default=c("400","1000"))
parser$add_argument("--significant" ,  help="significant bin from sigPeaks_CLT.R [default %(default)s]" , default="NA")
parser$add_argument("--coverageColorLimits" , nargs="+" , type="double" , help="OPTION: Provide two numbers. Bins > num1 or < num2 will be colored differently (default: 1.5 0.5) [default %(default)s]" , default= "NA" )
parser$add_argument("--debug"  , action="store_true" , help="dump session and quit [default %(default)s]" , default=FALSE)
args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }

if (is.na(coverageColorLimits[1])){ 
  coverageColorLimits <- c(1.5 , 0.5)
}
if (is.na(ylim[1])){
  ylim = c(0,3)
}


library(data.table)
library(ggplot2)
options(datatable.fread.input.cmd.message=FALSE)

if(debug){library(session);save.session("session_plotCovPerBin");quit()}


#chr size cumulative sum
chrSize <- read.table(chrSizeFile,header=F,stringsAsFactors=F,col.names=c("chr","size"))
chrSize$chr <- as.character(chrSize$chr)
chrSize <- chrSize[chrSize$chr %in% chrs,]
chrSizeCumSum <- cumsum(as.numeric(chrSize$size))
v <- c(0,chrSizeCumSum)
chrSizeCumSum <- v[-length(v)]
names(chrSizeCumSum) <- chrSize$chr

#chr mid points
chrMids <- c()
for (i in 1:length(chrSizeCumSum)-1) {chrMids = c(chrMids , (chrSizeCumSum[i] + chrSizeCumSum[i + 1]) / 2)}
chrMids =c(chrMids ,  ( chrSizeCumSum[length(chrSizeCumSum)] + sum(chrSize$size) ) / 2)

#read covPerBin
df        <- fread(cmd=paste("gunzip -c", covPerBin),stringsAsFactors=FALSE , header=T)
names(df) <- c("chromosome","start","end","meanCoverage","normalizedMeanCoverage","MAPQ")
df      <- as.data.frame(df , stringsAsFactors = FALSE)
row.names(df) <- paste(df$chromosome, df$start, df$end, sep="_")
df$chromosome <- as.character(df$chromosome)

#filter
df <- df[df$chromosome %in% chrs, ]

#reformat
newDf<-NULL
for (chr in chrs){rbind(newDf,df[df$chromosome == chr,]) -> newDf};
df <- newDf
df$binMidPoint <- (df$start) + (((df$end - df$start)+1) /2)

#color extreme bins
df$color  <- "black"
df$status <- "background"
if(! is.na(significant)){
  sigDf  <- read.table(significant,header=T,stringsAsFactors=F,sep="\t")
  row.names(sigDf) <- paste(sigDf$chr , sigDf$start , sigDf$end , sep='_')
  depletedBins <- row.names(sigDf[sigDf$direction == "depletion",])
  amplifiedBins <- row.names(sigDf[sigDf$direction == "amplification",])  
  df[amplifiedBins,"status"]="amplified"
  df[depletedBins,"status"]="depleted"
  df[amplifiedBins,"color"]="orange"
  df[depletedBins,"color"]="blue" 
} else if (! is.na(coverageColorLimits[1])){
  df$color[df$normalizedMeanCoverage > coverageColorLimits[1]] = "orange"
  df$color[df$normalizedMeanCoverage < coverageColorLimits[2]] = "blue"
  df$status[df$normalizedMeanCoverage > coverageColorLimits[1]]="enriched"
  df$status[df$normalizedMeanCoverage < coverageColorLimits[2]]="depleted"
}
df[df$MAPQ < minMAPQ , "status"] = "lowMAPQ"
df[df$MAPQ < minMAPQ , "color"]  = "gray"


#plot all together
dfAll <- data.frame(chromosome=df$chromosome , start=df$start , normalizedMeanCoverage=df$normalizedMeanCoverage , color=df$color, stringsAsFactors=F)
dfAll$normalizedMeanCoverage[dfAll$normalizedMeanCoverage > ylim[2]] = ylim[2]
positionCorrection <- chrSizeCumSum[match(dfAll$chromosome, names(chrSizeCumSum))]
dfAll$startFixed <- dfAll$start + positionCorrection
png(paste0(outName,".all.png") , height=binOverviewSize[1], width=binOverviewSize[2] , type='cairo')
plot(x=dfAll$startFixed , y=dfAll$normalizedMeanCoverage , xlab="chromosomes" , ylab="genomic bin coverage", col=dfAll$color , xaxt = 'n', xaxs="i" , yaxs="i" , font=2 , pch=19 , ylim=c(ylim[1],ylim[2]))
segments(x0=0 , y0=coverageColorLimits[1]  , x1=sum(chrSize$size) , y1=1.5 , lty=3, col='red', lwd=1)
segments(x0=0 , y0=1   , x1=sum(chrSize$size) , y1=1 , lty=2, col='black', lwd=1)
segments(x0=0 , y0=coverageColorLimits[2]  , x1=sum(chrSize$size) , y1=0.5 , lty=3, col='red', lwd=1)
abline(v=chrSizeCumSum[-1], lty=1, col='gray', lwd=2)
axis(1, at=chrMids, labels=names(chrMids) , font=2)
axis(4, at=coverageColorLimits, labels=coverageColorLimits , font=2)
dev.off()

#factorize for ggplot
df$chromosome <- factor(df$chromosome , levels=chrs)

#plot by chromosome
pdf(paste0(outName,".byChr.pdf")) 
for (chr in chrs){
  dfChr <- subset(df,chromosome == chr)
  dfChr$normalizedMeanCoverage[dfChr$normalizedMeanCoverage > ylim[2]]=ylim[2]
  p <- ggplot(dfChr, aes(binMidPoint,normalizedMeanCoverage)) + geom_point(colour="grey50", size = 2, alpha = 0.1) + coord_cartesian(ylim =c(ylim[1],ylim[2])) + theme(legend.position="none") + xlab("genomic bin position") + ggtitle(paste("chromosome",chr)) + theme_bw() + ylab("genomic bin coverage")
  p <- p + geom_point(size=1.5,alpha = 0.3,aes(colour=status)) + scale_colour_manual(values = c("background" = "black","amplified" = "orange","depleted" = "blue","lowMAPQ" = "gray"),name="copy number")
  print(p)
}
dev.off()

#plot faceting
png(paste0(outName,".faceting.png"),height=binOverviewSize[1], width=binOverviewSize[2] ,type='cairo')   
dfFaceting <- df
dfFaceting$normalizedMeanCoverage[dfFaceting$normalizedMeanCoverage > ylim[2]]=ylim[2]
p <- ggplot(dfFaceting, aes(binMidPoint,normalizedMeanCoverage)) + geom_point(colour="grey50", size = 1, alpha = 0.1) + facet_wrap(~ chromosome,nrow=4,scales="free") + coord_cartesian(ylim =c(ylim[1],ylim[2])) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("genomic bin position") + ylab("genomic bin coverage") + scale_x_continuous(labels = function(x) format(x, scientific = FALSE))
p = p + geom_point(size=0.5,alpha=0.5,aes(colour=status)) + scale_colour_manual(values = c("background" = "black","amplified" = "orange","depleted" = "blue","lowMAPQ" = "gray"),name="copy number")
print(p)
dev.off()

df[,c("meanCoverage" , "binMidPoint")] <- NULL
write.table(append=F,col.names=T,file=paste0(outName,".tsv"),quote=F,sep="\t",row.names=T,x=df)
system( paste0("gzip -f ",outName,".tsv") )
