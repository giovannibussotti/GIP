#Read coverage per bin data and returns the CNVs statistically significant
#########
#CONTEXT#
#########
#You generated single nucleotide coverage distribution (SNCD) #e.g.
#1	1	55
#1	2	56
#and you binned it in bins of size n (covPerBin or chrStartEndScore). 
#The bin score is approximated to the mean coverage in the bin, normalized by median chr coverage (to account for chr ploidy) and GC corrected.
#This is the script INPUT#e.g.
#1      1       300     50
#1      301     600     60
#Single Nucleotide coverage distribution after binning (SNCDab) is derived by re-expanding the bins into single nucleotides, with each nucleotide taking the expected normalized coverage of the bin
#e.g. SNCDab
#1	1	50
#1	1	50
#..
#1	301	60
#For the central limit theorem (CLT), regardless the shape of SNCDab, the sampling distribution of the sample means (SDSM) will be gaussian. In other words, if we sample many times SNCDab, and we compute the mean of the sample each time, we get a gaussian distribution
#e.g. SDSM
#50,60,55,70,52...
#Another important property is that The mean (mu) and the standard error (se) of SNCDab are the mean (mu) and the standard deviation (sd) of SDSM with sample size equal n
#So we do not need to do the sampling to generate SDSM to know its mu and sd.
#Given a bin of size n of asjacent bases in SNCDab and its mean score (i.e. the covPerBin score), my null-hypothesis is that we can get a bin with the same or higher score just by randomly selecting bins of n nucleotides from SNCDab. 
#The competing hypothesis is that there is something special about that specific bin.
#Because SDSM is gaussian, I can compute the P-value of each covPerBin bin by measuring how many se away each bin score is from the SNCDab mu.
##############
#SCRIPT-STEPS#
##############
#1)read covPerBin or chrStartEndScore file measuring the bin size n
#2)expand the bins generating SNCDab
#3)mu and se of SNCDab are measured to get mu and sd of SDSM
#4)A P-value is given to each input bin	
#5)The p-value is corrected by multiple testing
#6)Bins that pass the p-value filter are selected
#7)Adjacent successfull bins are merged in larger areas (segments), and the coverage score merged
#8)if --coverageThresholds the segments are selected to have a score above or below the two threshold indicated
#########
#WARNING#
#########
#In covPerBin or chrStartEndScore the bin size is regular, i.e. always the same. In covPerGe the bin size is the size of each gene, and it is not constant. So you should not use this script but sigPeaks_mixture.R

#NOTE
#if you are comparing coverage in two conditions you can first subtract the two samples (use subtractBins.R) and then run sigPeaks_CLT.R on that
#The the subtracted distribution score will be centered on 0 (not 1) so you should shift eventual --coverageThresholds filters (e.g. using -0.5 and 1) 

#see also
#sigPeaks_mixture.R, subtractBins.R 

#e.g. Rscript sigPeaks_CLT.R --input ../../../../sp-ama_Ht0_5749.covPerBin.gz --inFormat covPerBin --pThresh 0.0001 --coverageThresholds 2 0.5
#################################################
#               CONFIGURATION
#################################################
suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify our desired options # by default ArgumentParser will add an help option
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
	df <- df[,c("chromosome" , "start" , "end" , "median")]
	names(df) <- c("chr","start","end","score")
        df$chr <- as.character(df$chr)
} else {
	stop("input inFormat not recognized")
	quit(save = "no", status = 1, runLast = FALSE)
}
n = df[1,"end"] - df[1,"start"] + 1
########
#expand#
########
SNCDab <- rep(df$score,each=n)
########
#mu, se#
########
mu <- mean(SNCDab)
SNCDab_se <- sd(SNCDab)/sqrt(n)
######
#pval#
######
p.lowerTail  <- pnorm(df$score , mean=mu , sd=SNCDab_se , lower.tail=TRUE)
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
redD <- rbind(collapse(ampBins) , collapse(depBins))
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
write.table(paste("SNCDab mean", mu)      ,quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(paste("SNCDab standard error" , SNCDab_se),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table("##FILTERING AND FALSE DISCOVERY RATE",quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(paste("significant bins with minMAPQ",minMAPQ, ", adjusted p-val threshold", pThresh ,"and correction strategy", padjust," that were retained after possible coverageThresholds filtering are", length(significant[,1])),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(paste("collapsed adjacent bins (i.e. segments) that are retained after possible coverageThresholds and a min segment length of",minLen,"are", length(redD[,1])),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
pdf(paste0(outName,".pdf"))
hist(df$pvalue,breaks=100,col="orange",main="p-value distribution")
hist(df$pvalueAdj,breaks=100,col="orange",main="p-value distribution\nmultiple test correction")
dev.off()
