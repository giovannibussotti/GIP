suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--zsThresh"    , type="integer"      , help="z-score threshold (if --zscoreFilter) [default %(default)s]", default=0 )
parser$add_argument("--pThresh"     , type="double"       , help="p-value threshold [default %(default)s]" , default=1 )
parser$add_argument("--zscoreFilter", action="store_true" , help="flag. filter bins by z-score instead of pvalue [default %(default)s]" , default=FALSE )
parser$add_argument("--outName"     , type="character" , help="out name [default %(default)s]" , default="sigPeaksOut")
parser$add_argument("--input"       , type="character" , help="input bin file, with coverage values normalized by median chromosome coverage [default %(default)s]")
parser$add_argument("--padjust"     , type="character" , help="if --pFilter, corrects the p-value for multiple testing. This option is recommended. Select the p.adjust method. If you have independent bins, like in covPerGe, you can use benjamini hochberg (BH), otherwise if there is a dependency between bins (like in covPerBin, adjacent bins have more chances to be both amplified) then you should use Benjamini & Yekutieli (BY)  [default %(default)s]", default="NA")
parser$add_argument("--inFormat"    , type="character" , help="format of the input bin file. Supperted are chr start end score (chrStartEndScore), covPerBin and covPerGe [default %(default)s]" , default = "covPerGe")
parser$add_argument("--minMAPQ"      , type="integer"   , help="OPTION: if --inFormat is covPerBin or covPerGe filter out bins MAPQ < --minMAPQ (recommended 50) [default %(default)s]", default=0 )
parser$add_argument("--coverageThresholds" , nargs="+" , type="character" , help="OPTION: Provide two numbers, and selected peaks will undergo an additional filter. Just peaks > num1 or < num2 will be selected [default %(default)s]" , default= "NA" )
parser$add_argument("--minLen"      , type="integer"   , help="OPTION: selected peaks need to be at least --minLen nucleotides long. Not recommended for covPerGe inputs because you risk to discard small genes [default %(default)s]", default=0 )
parser$add_argument("--debug"  , action="store_true" , help="dump session and quit [default %(default)s]" , default=FALSE)
args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }

if(debug){library(session);save.session("session_DEBUG_sigPeaks_mixture");quit();}

pFilter=TRUE
if (zscoreFilter){
  pFilter=FALSE
}

if (! is.na(coverageThresholds[1])){
        coverageThresholds <- as.double(coverageThresholds)
}

library(mixtools)
library(GenomicRanges)
library(fitdistrplus)

#read input
if(inFormat == "chrStartEndScore"){
	df        <- read.table(input,header=F, stringsAsFactors=F)
	names(df) <- c("chr","start","end","score")
} else if (inFormat == "covPerBin"){
	df <- read.table(input,header=T, stringsAsFactors=F)
	df <- df[df$MAPQ >= minMAPQ ,] #filter
	df <- df[,c("chromosome" , "start" , "end" , "median")]
	names(df) <- c("chr","start","end","score")
} else if (inFormat == "covPerGe") {
	df       <- read.table(input,header=T,stringsAsFactors=F,sep="\t")
	df       <- df[df$MAPQ >= minMAPQ ,] #filter
	df$chr   <- gsub(x=df$locus,pattern="(.+):(.+)-(.+)$",replacement="\\1")
	df$start <- as.numeric(gsub(x=df$locus,pattern="(.+):(.+)-(.+)$",replacement="\\2"))
	df$end   <- as.numeric(gsub(x=df$locus,pattern="(.+):(.+)-(.+)$",replacement="\\3"))
	row.names(df) <- df$gene_id
        df       <- df[,c("chr" , "start" , "end" , "normalizedMeanCoverage")]
	names(df) <- c("chr","start","end","score")
} else {
	stop("input inFormat not recognized")
	quit(save = "no", status = 1, runLast = FALSE)
}
#admixure gaussian model
mixmdl = normalmixEM(df$score)
pdf(paste0(outName,".pdf"))
plot(mixmdl,which=2)
lines(density(df$score), lty=2, lwd=2)
#the central distribution is defined as the one with smaller sigma
indexCentral=0;
indexOutlier=0; 
if(mixmdl$sigma[1] < mixmdl$sigma[2]){indexCentral=1;indexOutlier=2}else{indexCentral=2;indexOutlier=1}
#estimate z-score and p-value
df$zScore <- (df$score - mixmdl$mu[indexCentral])/mixmdl$sigma[indexCentral]
df$probOfTheValueOrLower <- pnorm(df$score,mixmdl$mu[indexCentral],mixmdl$sigma[indexCentral])
df$probOfTheValueOrHigher <- 1 - df$probOfTheValueOrLower
leftTail <- df$score < mixmdl$mu[indexCentral]
df$probabilityTailed <-  df$probOfTheValueOrHigher
df$probabilityTailed[leftTail] = df$probOfTheValueOrLower[leftTail]
hist(df$probabilityTailed,breaks=100,col="orange",main="p-value distribution")
if(! is.na(padjust)){
	df$probabilityTailed <- p.adjust(df$probabilityTailed,method= padjust)
	hist(df$probabilityTailed,breaks=100,col="orange",main="p-value distribution\nmultiple test correction")
}
#filter
if(zscoreFilter){
	significant <- df[df$zScore > zsThresh | df$zScore < -zsThresh ,]
} else if(pFilter) {
	significant <- df[df$probabilityTailed < pThresh,]
} else {
	stop("filter not recognized")
	quit(save = "no", status = 1, runLast = FALSE)
}
#check normality of central distribution
centralDist <- df[df$zScore > -2 & df$zScore < 2 ,]
st <- shapiro.test(sample(centralDist$score,size=min(length(centralDist$score) , 5000) ))
descdist(centralDist$score, discrete = FALSE)
plot(fitdist(centralDist$score, "norm"))
dev.off()
#collapse adjacent bins averaging the score
redD <- NULL
if (inFormat != "covPerGe"){
 gr = GRanges(significant$chr, IRanges(c(significant$start),c(significant$end)),strand="*")
 mcols(gr)$score <- significant$score
 grRed <- reduce(gr)
 mcols(gr)$peak <- subjectHits(findOverlaps(gr, grRed))
 d <- as.data.frame(gr)
 d <- subset( d, select = -c(strand, width) )
 names(d)<- gsub(x=names(d),pattern="seqnames",replacement="chr")
 s <- split(d,d$peak)
 redDs <- lapply(split(d,d$peak),function(x){
	tmp<- data.frame(chr=x$chr[1],start=min(x$start), end=max(x$end))
	mergedValuesDf <- cbind( tmp , t(as.data.frame( round (apply( subset( x, select = -c(chr, start,end,peak) ) , 2,mean),digits=2) )) )
	row.names(mergedValuesDf)=""
	return(mergedValuesDf)
 })
 redD <- do.call(rbind,redDs)
} else {
 redD = significant
}
##########################
#arbitrary peak filtering#
##########################
if (! is.na(coverageThresholds[1])){
        redD <- redD[redD$score > coverageThresholds[1] | redD$score < coverageThresholds[2] , ]
}
redD <- redD[(redD$end - redD$start)+1 > minLen , ]
write.table(x=redD,quote=F,sep="\t",col.names=T,row.names=T,file=paste0(outName,".tsv"))
#write stats
distrs<- data.frame(distributionName=c("central distribution","outlier distribution"),mu=c(mixmdl$mu[indexCentral],mixmdl$mu[indexOutlier]) ,sigma=c(mixmdl$sigma[indexCentral],mixmdl$sigma[indexOutlier]),lambda=c(mixmdl$lambda[indexCentral],mixmdl$lambda[indexOutlier]))
write.table("##FITTING DESCRIPTION",quote=F,col.names=F,sep="\t",append=F,row.names=F,file=paste0(outName,".stats"))
write.table(paste("model",mixmdl$ft),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(paste("number of iterations to converge",length(mixmdl$all.loglik)-1),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(paste("loglik",mixmdl$loglik),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table("##DISTRIBUTIONS FEATURES",quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(distrs,quote=F,col.names=T,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table("##TESTING NORMALITY OF CENTRAL DISTRIBUTION",quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
write.table(as.matrix(unlist(st)),quote=F,col.names=F,sep="\t",append=T,file=paste0(outName,".stats"))
write.table("##FILTERING AND FALSE DISCOVERY RATE",quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
if (inFormat != "covPerGe"){
  write.table(paste("significant bins with minMAPQ",minMAPQ,"are", length(significant[,1])),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
  write.table(paste("collapsed adjacent bins (i.e. CNVs) that survived possible coverageThresholds and a min peak length of",minLen,"are", length(redD[,1])),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
} else {
  write.table(paste("significant genes with a MAPQ > than ",minMAPQ,"and possible user defined coverage or length thresholds are", length(significant[,1])),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
}

if(zscoreFilter){
	write.table(paste("z-score threshold", zsThresh),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
	FD=length(df[,1]) - (pnorm(zsThresh) * length(df[,1]))
	write.table(paste("number of expected falsely significant bins, passing the z-score threshold just by chance", FD),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
}
if(pFilter){
	write.table(paste("p-value correction used", padjust),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
	FD=length(df[,1]) - (pThresh * length(df[,1]))
	if(is.na(padjust)){
		write.table(paste("since no P-value correction is adopted, this is the number of expected falsely significant bins one would expect just by chance", FD),quote=F,col.names=F,sep="\t",append=T,row.names=F,file=paste0(outName,".stats"))
	}
}
