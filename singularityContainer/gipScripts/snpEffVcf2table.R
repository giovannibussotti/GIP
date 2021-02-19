suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--vcfFile"     , help="vcf file to load. [default %(default)s]" )
parser$add_argument("--outName"     , help="tsv out table. [default %(default)s]" )
parser$add_argument("--keepMultipleAlts" , action="store_true" , help="keep positions with multiple alternative variants. WARNING: This option is deprecated. I is not sure if snpEff reports the effect of the firt or the second (or third) alternative. So it is safer to discard multi alts position, as for default [default %(default)s]" , default=FALSE)
parser$add_argument("--debug"  , action="store_true" , help="dump session and quit [default %(default)s]" , default=FALSE)
args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"  ){args[[n]] <- NA  } }
for (n in names(args)){assign(n,args[[n]]) }
if(debug){library(session);save.session("session_DEBUG_snpEffVcf2Table");quit()}

library(VariantAnnotation)

vcf <- readVcf(vcfFile, "placeHolderFakeName")
if(! keepMultipleAlts){
    howmanyalleles <- elementNROWS(alt(vcf))
    vcf <- vcf[howmanyalleles == 1,]
}
howmanyalleles <- elementNROWS(alt(vcf))

ROforEachAllele      <- rep(geno(vcf)[["RO"]],howmanyalleles)
RRforEachAllele      <- unlist(geno(vcf)[["AO"]]) / (unlist(geno(vcf)[["AO"]]) + ROforEachAllele)
AOforEachAllele      <- unlist(geno(vcf)[["AO"]])
CHRforEachAllele     <- rep(as.character(seqnames(vcf)),howmanyalleles)
startForEachAllele   <- rep(start(vcf),howmanyalleles)
depthForEachAllele   <- rep(info(vcf)$DP,howmanyalleles)
QaForEachAllele      <- unlist(info(vcf)$QA)
QrForEachAllele      <- rep(info(vcf)$QR ,howmanyalleles)
MQMforEachAllele     <- unlist(info(vcf)$MQM)
MQMRforEachAllele    <- rep(info(vcf)$MQMR,howmanyalleles)
ref_altForEachAllele <- paste( as.vector(rep(unlist(ref(vcf)),howmanyalleles))   , as.vector(unlist(alt(vcf))) ,sep="_")
annList <- as.list(info(vcf)$ANN)
anns    <- sapply(annList,function(x){paste(x,collapse=",")})
rm(annList)
df <- data.frame(
      chr=CHRforEachAllele,
      position=startForEachAllele ,
      freq=RRforEachAllele ,
      AO=AOforEachAllele,
      alt=rep(howmanyalleles,howmanyalleles),
      totDepth=depthForEachAllele ,
      QA=QaForEachAllele ,
      QR=QrForEachAllele ,
      MQM=MQMforEachAllele ,
      MQMR=MQMRforEachAllele,
      ref_alt=ref_altForEachAllele,
      EFF=rep(anns,howmanyalleles),
      stringsAsFactors=FALSE
    )
write.table(x=df,row.names=F,col.names=T,quote=F,file=outName,sep="\t")
system(paste0("gzip ",outName))



