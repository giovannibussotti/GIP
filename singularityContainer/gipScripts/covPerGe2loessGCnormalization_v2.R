#following http://stats.stackexchange.com/questions/52860/loess-and-ma-normalization-in-r
#which is the same as the withing lane (i.e. sample specific) regression normalization described in this paper: GC content normalization for RNA-seq data (PMID 22177264).

#given a covPerGe file this script plot the GC% of each gene vs the meanCoverge (and the normalizedMeanCoverage). 
#Then if fits a loess regression
#Then is correct the original gene meanCoverge (or normalizedMeanCoverage) by subtracting the values on the loess model. If there is a GC bias the loess has a bump. So genes with GC values corresponnding to a loess bump will be more penalized than other genes were the loess is straight
#Finally the script plots again the corrected genes, fitting a new loess, that now should be more or less straight  and generates a loess corrected covPerGe (meanCoverge and normalizedMeanCoverage GC corrected)

#The script repeats the same operation for both meanCoverge and normalizedMeanCoverage, correcting both for %GC bias
#see also covPerBin2loessGCnormalization.R

#NOTE: This _v2 of the script is an improvement of the previous version because the "span" parameter of the loess function is chose using a 5 folds cross validation
#e.g. Rscript covPerGe2loessGCnormalization.R --ASSEMBLY /pasteur/projets/policy01/BioIT/Giovanni/datasets/assemblies/leishmanias/ENSEMBL_ProtistsRelease29/Leishmania_donovani_bpk282a1.GCA_000227135.2.29.dna_sm.toplevel.fa --DIR /pasteur/projets/policy01/BioIT/Giovanni/leish/p2p5/bwa --SAMPLE sp-ama_Ht0_5749 --ylim 10 --outName sp-ama_Ht0_5749

suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify our desired options # by default ArgumentParser will add an help option
parser$add_argument("--ASSEMBLY" , help="genome reference [default %(default)s]" )
parser$add_argument("--DIR" , help="dir with the covPerGe.gz file [default %(default)s]" )
parser$add_argument("--SAMPLE" , help="sample name [default %(default)s]")
parser$add_argument("--outName" , help="out name [default %(default)s]")
parser$add_argument("--nuc" , help="bedtools nuc output (of an initial bed3 format file) matching the covPerGe file [default %(default)s]")
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
df            <- read.table(paste0(DIR,"/",SAMPLE,".covPerGe.gz"),header=T,stringsAsFactors=F,sep="\t")
df$chromosome <- gsub(x=df$locus,pattern="(.+):(.+)-(.+)$",replacement="\\1")
df$start      <- as.numeric(gsub(x=df$locus,pattern="(.+):(.+)-(.+)$",replacement="\\2"))
df$end        <- as.numeric(gsub(x=df$locus,pattern="(.+):(.+)-(.+)$",replacement="\\3"))
df$tag <- paste(df$chromosome,df$start,df$end,sep="_")
nucDf   <- read.table(nuc , header=F , stringsAsFactors=F , sep="\t")[,c(1,2,3,5)]
names(nucDf) <- c("chromosome","start","end","CorGfreq")
nucDf$tag <- paste(nucDf$chromosome,nucDf$start,nucDf$end,sep="_")
df <- merge(df,nucDf[,c("tag","CorGfreq")],by="tag")
df$tag <- NULL
rm(nucDf)

###############################################################
#loess fitting then normalizedMeanCoverage coverage correction#
###############################################################
set.seed(1)
#sampling
cc <- df[complete.cases(df), ]
if ((sampling == 0) || (sampling > nrow(cc)))  {
        sampling = nrow(cc)
}
samp1 <- cc[sample(nrow(cc), sampling ,replace=F), ]
#fit a loess
lnmc = loess.wrapper(samp1$CorGfreq , samp1$normalizedMeanCoverage , span.vals = seq(0.2,1,by=0.1) , folds=5)
#for each bin extract the loess y value 
df$lnmcPredict <- predict(lnmc, df$CorGfreq)
#subract the y of each bin by the y of the loess model (this is the normalization step, in fact is a subtraction)
df$nmcLoessNorm <- df$normalizedMeanCoverage - df$lnmcPredict
#fit another loess on the corrected bins
lnmcAfterCorrection = loess.wrapper(df$CorGfreq , df$nmcLoessNorm , span.vals = seq(0.2,1,by=0.1) , folds=5)
#for each bin extract the loess (second model, the one done with corrected bins) y value
df$lnmcPredictAfterCorrection <- predict(lnmcAfterCorrection, df$CorGfreq)


###################################################
#plot cov VS %GC before and after loess correction#
###################################################
pdf(paste0(outName,".gcLnorm.covPerGe.pdf"))
plot(df$CorGfreq , df$normalizedMeanCoverage , main="before loess correction",ylab="normalizedMeanCoverage",xlab="%GC")
lines(df$CorGfreq[order(df$CorGfreq)], df$lnmcPredict[order(df$CorGfreq)],col="red") 
plot(df$CorGfreq , df$nmcLoessNorm , main="after loess correction",ylab="normalizedMeanCoverage",xlab="%GC")
lines(df$CorGfreq[order(df$CorGfreq)], df$lnmcPredictAfterCorrection[order(df$CorGfreq)],col="blue")
dev.off()

##################################################
#generate a covPerGe file with corrected coverage#
##################################################
covPerGe <- df[,c("gene_id" ,  "locus" , "meanCoverage" , "normalizedMeanCoverage" , "MAPQ")]
covPerGe$normalizedMeanCoverage <- df$nmcLoessNorm

#add back the difference between the original median and the loess subtracted median (so to center it again on 1)
covPerGe$normalizedMeanCoverage <- covPerGe$normalizedMeanCoverage + (median(df$normalizedMeanCoverage , na.rm=T) - median(covPerGe$normalizedMeanCoverage , na.rm=T))

#after correction some genes might have negative coverage. send these bins to 0
covPerGe$normalizedMeanCoverage[covPerGe$normalizedMeanCoverage < 0]=0

#genes full of N have CorGfreq NaN, resulting in corrected mean and median of NA
#anyway these genes all have median 0 and almost always mean 0 or very low. So I send these bins to 0
covPerGe$normalizedMeanCoverage[is.na(covPerGe$normalizedMeanCoverage)] = 0
#rounding-off
is.num           <- sapply(covPerGe, is.numeric)
covPerGe[is.num] <- lapply(covPerGe[is.num], round, 3)
write.table(x=covPerGe,append=F,col.names=T,row.names=F,sep="\t",quote=F,file=paste0(outName,".gcLnorm.covPerGe"))
system(paste0("gzip " ,outName,".gcLnorm.covPerGe"))

