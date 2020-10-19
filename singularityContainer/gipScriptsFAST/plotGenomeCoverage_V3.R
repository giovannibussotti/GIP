#First generate genome coverage files (.covPerNt) where the sequencing depth of each nucleotide is normalised by the median genomic coverage 
#Then given a list of gzipped .covPerNt files (--file), this script downsample the covPerNt scores (e.g. by binning), then convert the bin scores to somy scores which will be used for the boxplots and statistics 

#The somy score is defined as the downsampled score (e.g. bin score) divided by the median downsampled score of a chromosome deemed to be disomic and multiplied by 2
#If no disomic chromosome is deemed, the somy score is just the downsampled score multiplied by 2
#To assess what chromosome is disomic you can do a dry run of the script and pick a chromsome that has median cenetered on 2 and does not vary much across the various samples

#To test whether the chromosome coverage varies in different conditions in theory one could use a wilcoxon (default), ks or aov test. 
#In practice, since there are many many observations (nucleotide) the distributions of the coverages will look always statistically significant (even if from the boxplot sometimes you can see they are very very similar). 
#Since the coverage vectors are big the statistical test have an hard time (wilcoxon or ks tests return mostly pvalues of 0 and 1). 
#that is why it is needed to adopt a downsampling approach. This script explores two alternative strategies.

#One is to use the resample function (--fft yes) which transforms the coverage profile to the Fourier Space and resize the signal before performing the inverse transform. It is explained here: https://support.bioconductor.org/p/66313/
#The size of the coverage distribution after shrinking is regulated by --maxShrinkedLength. The script will automatically estimate a shrinking factor to have a final chromosome coverage vector shorter than --maxShrinkedLength.
#It works OK for short chromosomes but it takes ages for big chromosomes. Also the fourier transform can introduce negative values (that do not make sense since it is a coverage), and the reported pvalues are not stable at all when changing --maxShrinkedLength

#The other, which is recommended, simply bins the chromosomes and estimate the median coverage of each bin (--window).
#It works much better than the fourier approach.
#As an option (--selectQuantiles) out of each distribution of bins it considers just those between the 10th and the 90th quantile, so to strip out outliers bins. This is because outliers, hotpost loci can also artificially return into statistically significant comparisons, while in fact the genomic coverage distributions are very similar (except for these outliers)

#the script runs Wilcoxon, Kolmogorov-Smirnov and AOV tests on the downsampled data for each comparison and returns one table for each test.
#it also return the difference between the median downsampled somy score for each chromosome and for each comparison

#To compare coverage distribution shapes the script generates qqplots (explained here: http://www.r-bloggers.com/exploratory-data-analysis-quantile-quantile-plots-for-new-yorks-ozone-pollution-data/), comparing the quantiles of the two downsampled distributions. If the dots (quantiles) are on the diagonal, then they distributions have the same shape.

#use --pcMapqFiles to remove the bases where the percent of mapping reads with good MAPQ is < 50%

#Tested with R version 3.2.2 , ggplot2_1.0.1 , data.table_1.9.6
#WARNING: since the coverage vectors of a chromosome in two conditions have the same length, in fact the qqplots is equivalent to sort the two vectors and plot one against the other. So the coordinates x,y of the points are the coverage score in the two conditions (you have as many quantiles as nucleotides, not a fixed quantile interval, like the 25%)
#WARNING: Version1 of this script wotks perfectly, but uses much more memory. Version2 applys the downsampling immediatly when reading the files. Version3 is a minor improvement adding the --geom and --pooled options, and the ridges plot. Also --ylim now accepts two values, the min and the max coverage range.
#WARNING: Converting covPerNt scores to somy scores first, and then downsample (e.g. binning) would results in slightly differet results (medians not exactly centered on 2)

#Example:
#Rscript plotGenomeCoverage_V3.R --DIR ../../../pipeOut/leishield_01072017/WP2/lsdOut --files infMon1_ZK27p51_P5rep1.covPerNt.gz ZK25.covPerNt.gz ZK28.covPerNt.gz infMon1_ZK47p51_P5rep1.covPerNt.gz LIPA60.covPerNt.gz ZK43.covPerNt.gz LIPA68.covPerNt.gz --NAMES ZK27 ZK25 ZK28 ZK47 LIPA60 ZK43 LIPA68 --outName karyotype --ylim 6 --chrs 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 --window 2500 --disomicChr 36 --geom violin

options(bitmapType='cairo')

#################################################
#		CONFIGURATION
#################################################
suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify our desired options # by default ArgumentParser will add an help option
parser$add_argument("--files" , nargs="+", required=TRUE, help="List of genome coverage files to load. Files must be gzipped [default %(default)s]" )
parser$add_argument("--NAMES" , nargs="+", required=TRUE, help="name of the genome coverage files loaded. This determines also the plotting order [default %(default)s]" )
parser$add_argument("--DIR" , required=TRUE , help="Directory containing the genome coverage files[default %(default)s]" )
parser$add_argument("--outName" , help="out name [default %(default)s]" , default="genomeCoveragePlot")
parser$add_argument("--ylim" , nargs="+" , type="double" , help="min and max ylim used in the plot. default 0 10 [default %(default)s]" , default="NA")
parser$add_argument("--chrs" , nargs="+", required=TRUE,  help="List of chromosome names to plot. This also defines the plotting order[default %(default)s]" )
parser$add_argument("--pcMapqFiles" , nargs="+", help="list of gzipped files of the same length of --files and --NAMES generated with the pcMapqPerNt function. In these files each base has the percent of reads with good MAPQ score. If specified --pcMapqFiles remove the bases where the percent of mapping reads with good MAPQ is < 50 [default %(default)s]", default="NA" )
parser$add_argument("--makeQqplots" , action="store_true" , help="for all samples combinations, and for each chromosome it computes qQplots [default %(default)s]" , default=FALSE)
parser$add_argument("--fft"  , help="run fourier transform to reduce the size of chromosome coverage vectore [default %(default)s]" , default="no")
parser$add_argument("--maxShrinkedLength" , type="integer" , help="if some chromosomes are very big the script dynamically increases the shrinkingFactor till the shrinked chromosome length is below --maxShrinkedLength  [default %(default)s]" , default=100)
parser$add_argument("--window" , type="integer" , help="split each chr in chunks of this size, then take mean coverage out of each window to compare ditributions (DEPENDENCY: --compareDist)  [default %(default)s]" , default=2500)
parser$add_argument("--selectQuantiles" , help="strip out bins <10 and >90 quantiles before comparing coverage distributions (DEPENDENCY: --window)  [default %(default)s]" , default="no")
parser$add_argument("--disomicChr"  , help="normalize by this chromosome [default %(default)s]", default="NA" )
parser$add_argument("--customColors", help="provide a file with header \"SAMPLE HEX\" as first two columns, specifying the color for each sample [default %(default)s]", default="NA" )
parser$add_argument("--geom" , help="select boxplot or violin [default %(default)s]" , default="boxplot")
parser$add_argument("--pooled" , action="store_true" , help="flag. pool together all the samples (i.e. one box per chromosome representing the coverage values of all samples) [default %(default)s]" , default=FALSE)

args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }

if (is.na(ylim[1])){
  ylim = c(0,10)
}

#sanity check
if (!is.na(disomicChr) & ! (disomicChr %in% chrs)) {
        stop("ERROR. disomicChr needs to be included in --chrs")
        quit(save = "no", status = 1, runLast = FALSE)  
}
#functions
resample <- function(x, length.out=length(x))
{
    stopifnot(is.numeric(x))
    if (length.out == length(x))
        return(x)
    y <- fft(x)
    n1 <- (min(length.out, length(x)) + 1L) %/% 2L
    n3 <- n1 - 1L
    n2 <- length.out - n1 - n3
    y1 <- head(y, n=n1)
    y2 <- complex(n2)
    y3 <- tail(y, n=n3)
    padded_y <- c(y1, y2, y3)
    suppressWarnings(as.numeric(fft(padded_y, inverse=TRUE))) / length(x)
}
wilcoxonTest <- function (contrasts){
  allPvalsWilcox <- NULL
  for (i in 1:length(contrasts$V1)){ 
    x <- as.character(contrasts[i,])
    contrastPvalsWilcox <- NULL
    for (chr in chrs) {
      #select the score for each chr in the two samples
      a <- subset(allSamples[[x[1]]], allSamples[[x[1]]]$chromosome == chr)$score
      b <- subset(allSamples[[x[2]]], allSamples[[x[2]]]$chromosome == chr)$score
      chrPvalW <- wilcox.test(a,b)$p.value
      contrastPvalsWilcox <- c(contrastPvalsWilcox , chrPvalW)
    }
    #combine contrast wilcox pvals
    contrastPvalsWilcox <- t(as.data.frame(contrastPvalsWilcox))
    row.names(contrastPvalsWilcox) <- paste(x[1],"_VS_",x[2],sep="")
    as.data.frame(rbind(allPvalsWilcox , contrastPvalsWilcox)) -> allPvalsWilcox
  }
  #write wilcox
  names(allPvalsWilcox) <- chrs
  allPvalsWilcox$sample <- row.names(allPvalsWilcox) ; allPvalsWilcox <- allPvalsWilcox[,c("sample",chrs)]
  write.table(x=allPvalsWilcox,append=F,quote=F,sep="\t",file=paste(outName,"_pvaluesWilcox.tsv",sep=""),col.names=T,row.names=F)
}
KSTest <- function (contrasts){
  allPvalsKS <- NULL
  for (i in 1:length(contrasts$V1)){ 
    x <- as.character(contrasts[i,])
    contrastPvalsKS <- NULL
    for (chr in chrs) {
      #select the score for each chr in the two samples
      a <- subset(allSamples[[x[1]]], allSamples[[x[1]]]$chromosome == chr)$score
      b <- subset(allSamples[[x[2]]], allSamples[[x[2]]]$chromosome == chr)$score
      chrPvalK <- ks.test(a,b)$p.value
      contrastPvalsKS     <- c(contrastPvalsKS     , chrPvalK)
    }
    #combine contrast KS pvals
    contrastPvalsKS <- t(as.data.frame(contrastPvalsKS))
    row.names(contrastPvalsKS) <- paste(x[1],"_VS_",x[2],sep="")
    as.data.frame(rbind(allPvalsKS , contrastPvalsKS)) -> allPvalsKS
  }
  #write KS
  names(allPvalsKS) <- chrs
  allPvalsKS$sample <- row.names(allPvalsKS) ; allPvalsKS <- allPvalsKS[,c("sample",chrs)]
  write.table(x=allPvalsKS,append=F,quote=F,sep="\t",file=paste(outName,"_pvaluesKS.tsv",sep=""),col.names=T,row.names=F)
}
AOVTest <- function (contrasts){
  allPvalsAOV <- NULL
  for (i in 1:length(contrasts$V1)){ 
    x <- as.character(contrasts[i,])
    contrastPvalsAOV <- NULL
    for (chr in chrs) {
      #select the score for each chr in the two samples
      a <- subset(allSamples[[x[1]]], allSamples[[x[1]]]$chromosome == chr)$score
      b <- subset(allSamples[[x[2]]], allSamples[[x[2]]]$chromosome == chr)$score
      tmpDf    <- data.frame( score=c(a , b) , sample=c(rep("a",length(a))  , rep("b",length(b))) , stringsAsFactors=F )
      chrPvalA <- summary(aov(score ~ sample , tmpDf))[[1]]$Pr[1]
      contrastPvalsAOV    <- c(contrastPvalsAOV    , chrPvalA)
    }
    #combine contrast AOV pvals
    contrastPvalsAOV <- t(as.data.frame(contrastPvalsAOV))
    row.names(contrastPvalsAOV) <- paste(x[1],"_VS_",x[2],sep="")
    as.data.frame(rbind(allPvalsAOV , contrastPvalsAOV)) -> allPvalsAOV
  }
  #write AOV
  names(allPvalsAOV) <- chrs
  allPvalsAOV$sample <- row.names(allPvalsAOV) ; allPvalsAOV <- allPvalsAOV[,c("sample",chrs)]
  write.table(x=allPvalsAOV,append=F,quote=F,sep="\t",file=paste(outName,"_pvaluesAOV.tsv",sep=""),col.names=T,row.names=F)
}
deltaMedians <- function (contrasts){
  allDeltas <- NULL
  for (i in 1:length(contrasts$V1)){ 
    x <- as.character(contrasts[i,])
    contrastDeltas <- NULL
    for (chr in chrs) {
      #select the score for each chr in the two samples
      a <- subset(allSamples[[x[1]]], allSamples[[x[1]]]$chromosome == chr)$score
      b <- subset(allSamples[[x[2]]], allSamples[[x[2]]]$chromosome == chr)$score
      delta <- median(a) - median(b)
      contrastDeltas    <- c(contrastDeltas    , delta)
    }
    #combine contrast AOV pvals
    contrastDeltas <- t(as.data.frame(contrastDeltas))
    row.names(contrastDeltas) <- paste(x[1],"_VS_",x[2],sep="")
    as.data.frame(rbind(allDeltas , contrastDeltas)) -> allDeltas
  }
  #write AOV
  names(allDeltas) <- chrs
  allDeltas$sample <- row.names(allDeltas) ; allDeltas <- allDeltas[,c("sample",chrs)]
  write.table(x=allDeltas,append=F,quote=F,sep="\t",file=paste(outName,"_deltaMedians.tsv",sep=""),col.names=T,row.names=F)
}
qqplots <- function (contrasts){
  for (i in 1:length(contrasts$V1)){ 
    x <- as.character(contrasts[i,])
    for (chr in chrs) {
      #select the score for each chr in the two samples
      a <- subset(allSamples[[x[1]]], allSamples[[x[1]]]$chromosome == chr)$score
      b <- subset(allSamples[[x[2]]], allSamples[[x[2]]]$chromosome == chr)$score        
      #qqplot
      system(paste("mkdir -p qqplots_",outName,sep=""))
      pdf(paste("qqplots_",outName,"/",x[1],"_VS_",x[2],"_chr",chr,".pdf",sep=""))
      q <- qqplot(a,b,xlab=paste("quantiles",x[1]),ylab=paste("quantiles",x[2]),main=paste("qqplot chr",chr))
      abline(1,1)
      fit <- lm(q$y ~ q$x)
      abline(fit,col="red")
      dev.off()
    }
  }
}
removeLowMAPQ <- function (df) {
	f  <- pcMapqFiles[i]
	fName <- paste(DIR,"/",f,sep="");
	dfMQ <- data.frame(fread(paste("gunzip -c", fName),colClasses=list(character=1)) ,stringsAsFactors=F )
	names(dfMQ) <- c("chromosome","position","pcMapq") 
	dfMQ$tag <- paste0(dfMQ$chromosome,"_",dfMQ$position)
	df$tag <- paste0(df$chromosome,"_",df$position)
	df$pcMapq <- dfMQ[match(df$tag,dfMQ$tag),"pcMapq"]
	df <- df[df$pcMapq > 50, c("chromosome","position","score")]
	return(df)
}
compressCovPerNt <- function (df) {
    outDf <- NULL
    for (chr in chrs) {
        tmpChr <- subset(df, df$chromosome == chr)
        a <- c()
        if (fft == "yes"){
            #dynamicly choose shrinkingFactor
            shrinkingFactor = 1
            while((length(tmpChr$score)/shrinkingFactor) > maxShrinkedLength){
                shrinkingFactor = shrinkingFactor + 1
            }
            #shrink with fourier fft the coverage vector if shrinkingFactor > 1
            a <- resample(tmpChr$score,length(tmpChr$score)/shrinkingFactor)
        } else if (window > 1){
            #binning
            a <- sapply(split(tmpChr$score, ceiling(seq_along(tmpChr$score)/window)),mean)
            if (selectQuantiles == "yes"){
                a <- subset(a , ((a > quantile(a,0.1)) & (a < quantile(a,0.9))) )
            }
        }
    outDf <- rbind(outDf,data.frame(chromosome=chr,position=1:length(a),score=a))
    }
    return(outDf)
}

library("data.table")
library(ggplot2)
library(gtools)
library(ggridges)
options(datatable.fread.input.cmd.message=FALSE)

#library(session)
#save.session("session")
#quit()

#read .gcov files in a list
allSamples <- list()
for (i in 1:length(files)){
    f = files[i]
    n = NAMES[i]
    fName = paste(DIR,"/",f,sep="");
    allSamples[[n]]    <- fread(paste("gunzip -c", fName),colClasses=list(character=1))
    names(allSamples[[n]]) <- c("chromosome","position","score") 
    allSamples[[n]] <- allSamples[[n]][allSamples[[n]]$chromosome %in% chrs,] 
    if (! is.na(pcMapqFiles[1])){
      allSamples[[n]] <- removeLowMAPQ(data.frame(allSamples[[n]] , stringsAsFactors=F))
    }
    #compress
    allSamples[[n]]        <- compressCovPerNt(allSamples[[n]])
    allSamples[[n]]$sample <- n
    #convert score to somy score
    if(!is.na(disomicChr)){
        normFact  <-  median(allSamples[[n]][ allSamples[[n]]$chromosome == disomicChr ,"score"])
        allSamples[[n]]$score <- allSamples[[n]]$score / normFact
    }
    allSamples[[n]]$score <- allSamples[[n]]$score * 2
}

#reformat
AllSamplesDf <- do.call(rbind,allSamples)
AllSamplesDf$sample     <- factor(AllSamplesDf$sample , levels=NAMES )
AllSamplesDf$chromosome <- factor(AllSamplesDf$chromosome , levels=chrs)

#customColors
if(! is.na(customColors)){
  customColDf  <- read.table(customColors,header=T,stringsAsFactors=F)
  customColVec <- customColDf$HEX
  names(customColVec) <- customColDf$SAMPLE
  customColVec <- gsub(x=customColVec,pattern="^",replacement="\\#")
}

#plot
#pdf(paste(outName,".pdf",sep=""),width=16)
png(paste0(outName,".boxplot.png"), type="cairo")
p <-  ggplot(AllSamplesDf, aes(chromosome, score)) 
if (geom == "violin"){
  if (pooled){
    p <- p + geom_violin(fill="lightblue") 
  } else {
    p <- p + geom_violin(aes(fill = sample)) 
  }
} else if (geom == "boxplot") {
  if (pooled){
    p <- p + geom_boxplot(fill="lightblue"   , outlier.shape = NA) 
  } else {
    p <- p + geom_boxplot(aes(fill = sample) , outlier.shape = NA) 
  }
}
p <- p + coord_cartesian(ylim = c(ylim[1],ylim[2]) ) + ylab("somy score") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1) , legend.title = element_text(size=13,face="bold") , axis.title = element_text(size=13,face="bold") , panel.grid.major = element_line(size=1))
if(! is.na(customColors)){
  p1 <- p + scale_fill_manual(values=customColVec)
} else {
  p1 <- p + scale_fill_brewer(palette="Dark2")
}
print(p1)
dev.off()

png(paste0(outName,".ridges.png"),type="cairo")
#ridges (pooled data, including all samples)
print( ggplot(AllSamplesDf, aes(x = score, y = chromosome, fill = stat(x))) + geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +   scale_x_continuous(expand = c(0, 0), limits=c(ylim[1],ylim[2])) +  scale_y_discrete(expand = c(0, 0)) + scale_fill_viridis_c(name = "score", option = "C") + theme_ridges(font_size = 13, grid = TRUE) )
dev.off()

#compare distributions
if(length(NAMES) > 1){
  contrasts <- as.data.frame(combinations(n=length(NAMES),r=2,NAMES),stringsAsFactors=F)
  wilcoxonTest(contrasts)
  KSTest(contrasts)
  AOVTest(contrasts)
  deltaMedians(contrasts)
  if (makeQqplots){
    qqplots(contrasts)
  }
}

#table summarizing median scores per chromosome
allMedians <- list()
for (n in names(allSamples)){
    medians <- NULL
    for (chr in chrs) {
        medians <- c(medians , as.numeric(median(subset(allSamples[[n]], allSamples[[n]]$chromosome == chr)$score)))
    }
    names(medians)  <- chrs
    allMedians[[n]] <- medians
}
allMediansDf <- as.data.frame(do.call(rbind,allMedians),stringsAsFactors=F)
newnames <- names(allMediansDf)
allMediansDf$sample <- rownames(allMediansDf)
write.table(x=allMediansDf[,c("sample",newnames)],append=F,quote=F,sep="\t",file=paste0(outName,".allMedians.tsv"),col.names=T,row.names=F)


