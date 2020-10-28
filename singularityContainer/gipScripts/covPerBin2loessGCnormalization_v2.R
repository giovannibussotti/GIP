suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify our desired options # by default ArgumentParser will add an help option
parser$add_argument("--ASSEMBLY" , help="genome reference [default %(default)s]" )
parser$add_argument("--DIR" , help="dir with the covPerBin.gz file [default %(default)s]" )
parser$add_argument("--SAMPLE" , help="sample name [default %(default)s]")
parser$add_argument("--outName" , help="out name [default %(default)s]")
parser$add_argument("--nuc" , help="bedtools nuc output (of an initial bed3 format file) matching the covPerBin file [default %(default)s]")
parser$add_argument("--ylim" , type="integer" , help="[default %(default)s]" , default="10")
parser$add_argument("--sampling" , type="integer" , help="sampling size. If 0 the entire (non-NA) dataset is considered [default %(default)s]" , default="0")
parser$add_argument("--debug"  , action="store_true" , help="dump session and quit [default %(default)s]" , default=FALSE)
args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }
if(debug){library(session);save.session("session_DEBUG");quit()}

library(bisoreg)

##############################################
#extract sequence of each bin and measure %GC#
##############################################
df    <- read.table(paste0(DIR,"/",SAMPLE,".covPerBin.gz"),header=T,stringsAsFactors=F,sep="\t")
df$chromosome <- as.character(df$chromosome)
df$CorGfreq   <- read.table(nuc , header=F , stringsAsFactors=F , sep="\t")[,5]


#############################################
#loess fitting then mean coverage correction#
#############################################
set.seed(1)
cc <- df[complete.cases(df), ]
if(sampling == 0){
        sampling = nrow(cc)
}
samp1 <- cc[sample(nrow(cc), sampling ,replace=F), ]
l <- loess.wrapper(samp1$CorGfreq , samp1$meanCoverage , span.vals = seq(0.2,1,by=0.1) , folds=5)


#On the normalizedMeanCoverage 
#for each bin extract the loess y value 
df$lPredict <- predict(l, df$CorGfreq)
#subract the y of each bin by the y of the loess model (this is the normalization step, in fact is a subtraction)
df$loessNorm <- df$normalizedMeanCoverage - df$lPredict
#fit another loess on the corrected bins
samp2 <- df[sample(nrow(df), sampling ,replace=F), ]
lAfterCorrection = loess(samp2$loessNorm ~ samp2$CorGfreq)
#for each bin extract the loess (second model, the one done with corrected bins) y value
df$lPredictAfterCorrection <- predict(lAfterCorrection, df$CorGfreq)



###################################################
#plot cov VS %GC before and after loess correction#
###################################################
pdf(paste0(outName,".gcLnorm.covPerBin.pdf"))
plotSmooth <- function (cc,N,mc,Y) {
  #bandwidth has the same unit as data
  #so if the range is 1-100 bandwidth 2 will smooth over 2%-wide regions
  #the code below measures the range in each dimension and take the 2%
  bwx <- ((max(cc$CorGfreq) - min(cc$CorGfreq)) /100) *2
  bwy <- ((max(cc[[Y]]) - min(cc[[Y]])) /100) *2
  smoothScatter(cc$CorGfreq , cc[[Y]], ylim=c(-3,ylim) , xlim=c(0.2 , 0.8) , xlab="%GC",ylab="coverage",main=N, bandwidth=c(bwx,bwy))
  legend("topleft", c( paste("R=", round(cor(x=cc[[Y]] , y=cc$CorGfreq),2)) , "linear model" , "loess" ) , lty=c(0,2,1), col=c("","black","red") , bty="n" )
  f <- paste(Y , "~ CorGfreq")
  l <- lm(f , data=cc)
  abline(l,lty=2,col="black",lwd=2)
  grid()
  lines(cc$CorGfreq[order(cc$CorGfreq)], mc[order(cc$CorGfreq)],col="red")
}
par(mfrow=c(1,2))
cc <- df[complete.cases(df), ]
plotSmooth(cc,"normalizedMeanCoverage before loess correction",cc$lPredict               ,"normalizedMeanCoverage")
plotSmooth(cc,"normalizedMeanCoverage after loess correction" ,cc$lPredictAfterCorrection,"loessNorm")
dev.off()

###################################################
#generate a covPerBin file with corrected coverage#
###################################################
covPerBin <- df[,c("chromosome" , "start" ,  "end" , "meanCoverage" , "normalizedMeanCoverage" , "MAPQ")]
covPerBin$normalizedMeanCoverage <- df$loessNorm

#add back the difference between the original median and the loess subtracted median (so to center it again on 1)
covPerBin$normalizedMeanCoverage <- covPerBin$normalizedMeanCoverage + (median(df$normalizedMeanCoverage , na.rm=T) - median(covPerBin$normalizedMeanCoverage , na.rm=T))

#after correction some bins might have negative coverage. send these bins to 0
covPerBin$normalizedMeanCoverage[covPerBin$normalizedMeanCoverage < 0]=0

#bins full of N have CorGfreq NaN, resulting in corrected mean and median of NA
#anyway these bins all have median 0 and almost always mean 0 or very low. So I send these bins to 0
covPerBin$normalizedMeanCoverage[is.na(covPerBin$normalizedMeanCoverage)] = 0
#rounding-off
is.num     <- sapply(covPerBin, is.numeric)
covPerBin[is.num] <- lapply(covPerBin[is.num], round, 3)
write.table(x=covPerBin,append=F,col.names=T,row.names=F,sep="\t",quote=F,file=paste0(outName,".gcLnorm.covPerBin"))
system(paste0("gzip " ,outName,".gcLnorm.covPerBin"))
