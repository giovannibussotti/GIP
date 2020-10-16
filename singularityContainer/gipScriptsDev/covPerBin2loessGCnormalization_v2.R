#following http://stats.stackexchange.com/questions/52860/loess-and-ma-normalization-in-r
#which is the same as the withing lane (i.e. sample specific) regression normalization described in this paper: GC content normalization for RNA-seq data (PMID 22177264).

#given a covPerBin file this script plot the GC% of each bin vs the mean (and median) coverage. 
#Then if fits a loess regression on the the mean using just a sub-sample (default all points, but you can specify a number in --sampling to speed up the process, this is sampling1) using a 5 folds cross validation exploring the loess span parameter (which relates with the fraction of points used to fit the local regressions, and influence the model smoothness). 
#Then corrects the original bin mean (or median) coverage by subtracting the values on the loess model. If there is a GC bias the loess has a bump. So bins with GC values corresponnding to a loess bump will be more penalized than other bins were the loess is straight
#another loess is computed after GC correction (sampling2 and sampling3, without cross-validation) and expected to be more or less straigth
#Finally the script plots the GC% of each bin vs the mean (and median) coverage, before and after the GCcorrection, showing also the loess, a linear model and the R pearson coefficient before and after correction.

#see also covPerGe2loessGCnormalization.R 
#one limitation of this script is that to save time you fit the loess on the mean, then you use the same fit to correct also the median (not just the mean)
#also to save time the models are built on dataset samples, not the full dataset

#e.g. Rscript covPerBin2loessGCnormalization.R --ASSEMBLY /pasteur/projets/policy01/BioIT/Giovanni/datasets/assemblies/leishmanias/ENSEMBL_ProtistsRelease29/Leishmania_donovani_bpk282a1.GCA_000227135.2.29.dna_sm.toplevel.fa --DIR /pasteur/projets/policy01/BioIT/Giovanni/leish/p2p5/bwa --SAMPLE sp-ama_Ht0_5749 --ylim 10 --outName sp-ama_Ht0_5749

suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify our desired options # by default ArgumentParser will add an help option
parser$add_argument("--ASSEMBLY" , help="genome reference [default %(default)s]" )
parser$add_argument("--DIR" , help="dir with the covPerBin.gz file [default %(default)s]" )
parser$add_argument("--SAMPLE" , help="sample name [default %(default)s]")
parser$add_argument("--outName" , help="out name [default %(default)s]")
parser$add_argument("--ylim" , type="integer" , help="[default %(default)s]" , default="10")
parser$add_argument("--sampling" , type="integer" , help="sampling size. If 0 the entire (non-NA) dataset is considered [default %(default)s]" , default="0")
args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }

library(Biostrings)
library(bisoreg)


##############################################
#extract sequence of each bin and measure %GC#
##############################################
refAssembly    <- readDNAStringSet(ASSEMBLY)
df <- read.table(paste0(DIR,"/",SAMPLE,".covPerBin.gz"),header=T,stringsAsFactors=F,sep="\t")
df$chromosome <- as.character(df$chromosome)
seqs <- NULL
CorGfreq <- NULL
for (i in 1:length(df[,1])){
	chr=df[i,"chromosome"]
	start=df[i,"start"]
	end=df[i,"end"]
	seq <- subseq(refAssembly[[chr]],start, end)
	af <- alphabetFrequency(seq)
	freq <- (af[["C"]] + af[["G"]]) / (sum(af) - af[["N"]])
	seq <- toString(seq)
	rbind(seqs,seq) -> seqs  
	rbind(CorGfreq,freq) -> CorGfreq
}
#df$seqs    <- seqs
df$CorGfreq <- CorGfreq
df$CorGfreq <- as.vector(df$CorGfreq)

#############################################
#loess fitting then mean coverage correction#
#############################################
set.seed(1)
cc <- df[complete.cases(df), ]
if(sampling == 0){
	sampling = nrow(cc)
}
samp1 <- cc[sample(nrow(cc), sampling ,replace=F), ]
l <- loess.wrapper(samp1$CorGfreq , samp1$mean , span.vals = seq(0.2,1,by=0.1) , folds=5)

#library(session)
#save.session("session")
#quit()

#On the mean coverage
#for each bin extract the loess y value 
df$lmeanPredict <- predict(l, df$CorGfreq)
#subract the y of each bin by the y of the loess model (this is the normalization step, in fact is a subtraction)
df$meanLoessNorm <- df$mean - df$lmeanPredict
#fit another loess on the corrected bins
samp2 <- df[sample(nrow(df), sampling ,replace=F), ]
lmeanAfterCorrection = loess(samp2$meanLoessNorm ~ samp2$CorGfreq)
#for each bin extract the loess (second model, the one done with corrected bins) y value
df$lmeanPredictAfterCorrection <- predict(lmeanAfterCorrection, df$CorGfreq)
###############################################
#loess fitting then median coverage correction#
###############################################
df$lmedianPredict <- predict(l , df$CorGfreq)
df$medianLoessNorm <- df$median - df$lmedianPredict
samp3 <- df[sample(nrow(df), sampling ,replace=F), ]
lmedianAfterCorrection = loess(samp3$medianLoessNorm ~ samp3$CorGfreq)
df$lmedianPredictAfterCorrection <- predict(lmedianAfterCorrection, df$CorGfreq)


###################################################
#plot cov VS %GC before and after loess correction#
###################################################
pdf(paste0(outName,".gcLnorm.covPerBin.pdf"))
plotSmooth <- function (cc,N,mc,Y) {
        smoothScatter(cc$CorGfreq , cc[[Y]], ylim=c(-3,ylim) , xlim=c(0.2 , 0.8) , xlab="%GC",ylab="coverage",main=N)
        legend("topleft", c( paste("R=", round(cor(x=cc[[Y]] , y=cc$CorGfreq),2)) , "linear model" , "loess" ) , lty=c(0,2,1), col=c("","black","red") , bty="n" )
	f <- paste(Y , "~ CorGfreq")
        l = lm(f , data=cc)
        abline(l,lty=2,col="black",lwd=2)
        grid()
        lines(cc$CorGfreq[order(cc$CorGfreq)], mc[order(cc$CorGfreq)],col="red")
}
par(mfrow=c(2,2))
cc <- df[complete.cases(df), ]
plotSmooth(cc,"mean before loess correction"  ,cc$lmeanPredict                 ,"mean")
plotSmooth(cc,"mean after loess correction"   ,cc$lmeanPredictAfterCorrection  ,"meanLoessNorm")
plotSmooth(cc,"median before loess correction",cc$lmedianPredict               ,"median")
plotSmooth(cc,"median after loess correction" ,cc$lmedianPredictAfterCorrection,"medianLoessNorm")
dev.off()

###################################################
#generate a covPerBin file with corrected coverage#
###################################################
covPerBin <- df[,c("chromosome" , "start" ,  "end" , "mean" , "median" , "MAPQ")]
covPerBin$mean <- df$meanLoessNorm
covPerBin$median <- df$medianLoessNorm

#add back the difference between the original median and the loess subtracted median (so to center it again on 1)
covPerBin$mean   <- covPerBin$mean   + (median(df$mean , na.rm=T) - median(covPerBin$mean , na.rm=T))
covPerBin$median <- covPerBin$median + (median(df$median , na.rm=T) - median(covPerBin$median , na.rm=T))

##add to all bin the coverage of the lowest bin f negative, to make sure that coverage after correction is not negative
#if (min(covPerBin$mean,na.rm=T) < 0 )  { covPerBin$mean <- covPerBin$mean + abs(min(covPerBin$mean   , na.rm=T))}
#if (min(covPerBin$median,na.rm=T) < 0 ){ covPerBin$median <- covPerBin$median + abs(min(covPerBin$median   , na.rm=T))}

#after correction some bins might have negative coverage. send these bins to 0
covPerBin$mean[covPerBin$mean < 0]=0
covPerBin$median[covPerBin$median < 0]=0

#bins full of N have CorGfreq NaN, resulting in corrected mean and median of NA
#anyway these bins all have median 0 and almost always mean 0 or very low. So I send these bins to 0
covPerBin$mean[is.na(covPerBin$mean)] = 0
covPerBin$median[is.na(covPerBin$median)] = 0
#rounding-off
is.num     <- sapply(covPerBin, is.numeric)
covPerBin[is.num] <- lapply(covPerBin[is.num], round, 3)
write.table(x=covPerBin,append=F,col.names=T,row.names=F,sep="\t",quote=F,file=paste0(outName,".gcLnorm.covPerBin"))
system(paste0("gzip " ,outName,".gcLnorm.covPerBin"))


#old plot function
#plot(df$CorGfreq , df$mean , ylim=c(-3 , ylim), main="before loess correction",ylab="mean coverage")
#lines(df$CorGfreq[order(df$CorGfreq)], df$lmeanPredict[order(df$CorGfreq)],col="red") 
#plot(df$CorGfreq , df$meanLoessNorm , ylim=c(-3 , ylim), main="after loess correction",ylab="mean coverage")
#lines(df$CorGfreq[order(df$CorGfreq)], df$lmeanPredictAfterCorrection[order(df$CorGfreq)],col="blue")
#plot(df$CorGfreq , df$median , ylim=c(-3 , ylim), main="before loess correction",ylab="median coverage")
#lines(df$CorGfreq[order(df$CorGfreq)], df$lmedianPredict[order(df$CorGfreq)],col="red") 
#plot(df$CorGfreq , df$medianLoessNorm , ylim=c(-3 , ylim), main="after loess correction",ylab="median coverage")
#lines(df$CorGfreq[order(df$CorGfreq)], df$lmedianPredictAfterCorrection[order(df$CorGfreq)],col="blue")

#implement cross validation of different models, including glm and loess
#http://www.statmethods.net/advstats/glm.html
#http://www.statmethods.net/advstats/glm.html
#require(boot)
#glm.fit1 <- glm(mean ~ poly(CorGfreq , 1 , raw=TRUE) , data=cc)
#glm.fit2 <- glm(mean ~ poly(CorGfreq , 2 , raw=TRUE) , data=cc)
#glm.fit3 <- glm(mean ~ poly(CorGfreq , 3 , raw=TRUE) , data=cc)
#sp.fit1 <- smooth.spline(cc$CorGfreq[order(cc$CorGfreq)],cc$mean[order(cc$CorGfreq)],df=10)
#sp.fit2 <- smooth.spline(cc$CorGfreq[order(cc$CorGfreq)],cc$mean[order(cc$CorGfreq)],df=20)
#loess.fit = loess(mean ~ CorGfreq , data =cc)
#smoothScatter(df$CorGfreq , df$mean , ylim=c(-3 , ylim), main="before loess correction",ylab="mean coverage")
#lines(cc$CorGfreq[order(cc$CorGfreq)] ,predict(glm.fit1, type="response")[order(cc$CorGfreq)],col="black")
#lines(cc$CorGfreq[order(cc$CorGfreq)] ,predict(glm.fit2, type="response")[order(cc$CorGfreq)],col="gray")
#lines(cc$CorGfreq[order(cc$CorGfreq)] ,predict(glm.fit3, type="response")[order(cc$CorGfreq)],col="red")
#lines(sp.fit1,col="blue")
#lines(sp.fit2,col="brown")
#lines(cc$CorGfreq[order(cc$CorGfreq)] ,predict(loess.fit,cc$CorGfreq)[order(cc$CorGfreq)],col="green")
