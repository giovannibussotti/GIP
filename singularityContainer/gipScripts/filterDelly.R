suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--vcfFile" , nargs="+", help="vcf file  [default %(default)s]" )
parser$add_argument("--minMAPQ" , type="integer" , help="min Median mapping quality of paired-ends supporting the SV [default %(default)s]" , default="0")
parser$add_argument("--outDir" , help="output dir [default %(default)s]", default="_outFilterDelly")
parser$add_argument("--chrSizeFile" , help="chr size file [default %(default)s]", default="NA")
parser$add_argument("--chrEndFilter" , type="integer", help="number of bases from the chr ends to consider when filtering SV overlapping the telomers. Depends on --chrSizeFile  [default %(default)s]", default="0")
parser$add_argument("--topRcInv" , help="Select top --topRcInv inversions based on RC score [default %(default)s]", default="NA")
parser$add_argument("--topRcIns" , help="Select top --topRcIns insertions based on RC score [default %(default)s]", default="NA")
parser$add_argument("--topRcDel" , help="Select top --topRcDel deletions based on RC score [default %(default)s]", default="NA")
parser$add_argument("--topRcDup" , help="Select top --topRcDup duplications based on RC score [default %(default)s]", default="NA")
parser$add_argument("--topRcBnd" , help="Select top --topRcBnd break ends based on RC score [default %(default)s]", default="NA")
parser$add_argument("--topHqCountInv" , help="Select top --topHqCountInv inversions based on DV+RV score [default %(default)s]", default="NA")
parser$add_argument("--topHqCountIns" , help="Select top --topHqCountIns insertions based on DV+RV score [default %(default)s]", default="NA")
parser$add_argument("--topHqCountDel" , help="Select top --topHqCountDel deletions based on DV+RV score [default %(default)s]", default="NA")
parser$add_argument("--topHqCountDup" , help="Select top --topHqCountDup duplications based on DV+RV score [default %(default)s]", default="NA")
parser$add_argument("--topHqCountBnd" , help="Select top --topHqCountBnd break ends based on DV+RV score [default %(default)s]", default="NA")
parser$add_argument("--topHqPercentInv" , help="Select top --topHqPercentInv inversions based on DV+RV/DV+RV+DR+RR score [default %(default)s]", default="NA")
parser$add_argument("--topHqPercentIns" , help="Select top --topHqPercentIns insertions based on DV+RV/DV+RV+DR+RR score [default %(default)s]", default="NA")
parser$add_argument("--topHqPercentDel" , help="Select top --topHqPercentDel deletions based on DV+RV/DV+RV+DR+RR score [default %(default)s]", default="NA")
parser$add_argument("--topHqPercentDup" , help="Select top --topHqPercentDup duplications based on DV+RV/DV+RV+DR+RR score [default %(default)s]", default="NA")
parser$add_argument("--topHqPercentBnd" , help="Select top --topHqPercentBnd break ends based on DV+RV/DV+RV+DR+RR score [default %(default)s]", default="NA")
parser$add_argument("--chrsToKeep" , nargs="+", help="if specified, filter out all the SVs laying on chrs not in this list  [default %(default)s]" , default="NA" )
parser$add_argument("--rmLowQual"  , action="store_true" , help="Remove LowQual delly predictions" , default=FALSE)
parser$add_argument("--rmImprecise" , action="store_true" , help="Keep just delly predictions labeled as PRECISE", default=FALSE)
parser$add_argument("--debug"  , action="store_true" , help="Dump session and quit" , default=FALSE)
args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }

if(debug){library(session); save.session("session_filterDelly");quit()}

#read
library("VariantAnnotation")
vcf = readVcf(vcfFile, 'genome')

#fix 1-counting base shift
start(vcf) <- start(vcf) -1
end(vcf) <- end(vcf) -1
start(vcf)[start(vcf) < 1]=1

#remove lowQual
if(rmLowQual){
  vcf =vcf[rowRanges(vcf)$FILTER == "PASS",]
}
#remove SV whose the median MAPQ is too low (e.g. the reads align to repeats)
vcf =vcf[info(vcf)$MAPQ >= minMAPQ,]
#remove imprecise SV
if(rmImprecise){
  vcf <- vcf[info(vcf)$PRECISE, ]
}
#set vcf end position to the END field
end(vcf) = info(vcf)$END

#remove SVs on unwanted chrs
if(!is.na(chrsToKeep[1])){
	vcf <- vcf[seqnames(vcf) %in% chrsToKeep , ]
}
#separate BND
bnd <- vcf[info(vcf)$SVTYPE == "BND",]
vcf  <- vcf[info(vcf)$SVTYPE != "BND",]

#fix end and chr filter for BND
if(!is.na(chrsToKeep[1])){
bnd <- bnd[info(bnd)$CHR2 %in% chrsToKeep , ]
}
bndC <- data.frame(chr=info(bnd)$CHR2 , 
	    start=info(bnd)$POS2 , end=info(bnd)$POS2 +1 , stringsAsFactors=F)
bndCGr <- with(bndC, GRanges(chr, IRanges(start, end)))
names(bndCGr) <- names(bnd)

#remove SV overlapping the chr ends
filterChrEnds <- function (VCF) {
	ov  <- findOverlaps(VCF,chrExtremitiesGr)
	if (length(queryHits(ov)) > 0 ){
		VCF <- VCF[-queryHits(ov)]
	}	
	return(VCF)
}
rmSVoutsideChrRange <- function (VCF) {
	ov  <- findOverlaps(VCF,chrSizesGr)
	if (length(queryHits(ov)) > 0 ){
		VCF <- VCF[queryHits(ov)]
	}	
	return(VCF)
}
if(!is.na(chrSizeFile)){
   chrSizes       <- read.table(chrSizeFile,header=F,stringsAsFactors=F)
   chrSizes$start <- 1
   names(chrSizes)<- c("chr","end","start")
   chrBegins <- data.frame(chr=chrSizes$chr, start=chrSizes$start , end=chrEndFilter)
   chrEnds   <- data.frame(chr=chrSizes$chr, start=chrSizes$end - chrEndFilter , end=chrSizes$end)
   chrExtremities <- rbind(chrBegins,chrEnds)
   chrExtremitiesGr <- with(chrExtremities, GRanges(chr, IRanges(start, end)))
   vcf <- filterChrEnds(vcf)
   bnd <- filterChrEnds(bnd)
   bndCGr <- filterChrEnds(bndCGr)
   #rm SV outside chr range
   chrSizesGr <- with(chrSizes, GRanges(chr, IRanges(start, end)))
   vcf <- rmSVoutsideChrRange(vcf)
   bnd <- rmSVoutsideChrRange(bnd)
   bndCGr <- rmSVoutsideChrRange(bndCGr)
}
bndInter <- intersect(names(bndCGr) , names(bnd))
bndCGr <- bndCGr[bndInter,]
bnd <- bnd[bndInter,]



#select top SV AND print BED for circos
system(paste("mkdir -p",outDir))
topRc        <- list(INV= topRcInv        , INS= topRcIns        , DEL= topRcDel        , DUP= topRcDup        , BND= topRcBnd )
topHqCount   <- list(INV= topHqCountInv   , INS= topHqCountIns   , DEL= topHqCountDel   , DUP= topHqCountDup   , BND= topHqCountBnd )
topHqPercent <- list(INV= topHqPercentInv , INS= topHqPercentIns , DEL= topHqPercentDel , DUP= topHqPercentDup , BND= topHqPercentBnd )
SVsel <- NULL
genoFilter <- function(VCF,top,metrics){
  #select to SV based on a metrics
  if(is.na(top) || length(VCF) == 0 ){return(VCF)}
  top <- min(length(VCF),as.integer(top))
  if (metrics == "RC"){
      tmpDf <- as.data.frame(geno(VCF)[["RC"]] , drop=FALSE) 
   } else if (metrics == "HQcount") {
      tmpDf <- as.data.frame(geno(VCF)[["DV"]] + geno(VCF)[["RV"]] , drop=FALSE) 
   } else if (metrics == "HQpercent") {
      tmpDf <- as.data.frame( (geno(VCF)[["DV"]] + geno(VCF)[["RV"]] + 1) *100 / (geno(VCF)[["DV"]] + geno(VCF)[["RV"]] + geno(VCF)[["RR"]] + geno(VCF)[["DR"]] + 1) , drop=FALSE) 
   }
  best <- row.names(tmpDf[order(-tmpDf[,1]),drop=F,])[ 1:top ]
  return(VCF[best,])
}
for(SVTYPE in c("INV" ,"DUP", "DEL", "INS")){
   tmp <- vcf[info(vcf)$SVTYPE == SVTYPE]
   tmp <- genoFilter(VCF=tmp,top=topRc[[SVTYPE]],metrics="RC")
   tmp <- genoFilter(VCF=tmp,top=topHqCount[[SVTYPE]],metrics="HQcount")
   tmp <- genoFilter(VCF=tmp,top=topHqPercent[[SVTYPE]],metrics="HQpercent")
   SVsel <- c(SVsel, names(tmp))
   if(length(tmp) > 0 ) {
     bed <- data.frame(chr=paste0("chr",seqnames(tmp)) , start=start(tmp) , end=end(tmp) , stringsAsFactors=F )# , name=names(tmp) )
     write.table(x=bed,quote=F,append=F,col.names=F,row.names=F,sep="\t",file=paste0(outDir,"/",SVTYPE,".bed"))
   } else {
   	system(paste0("touch " , outDir,"/",SVTYPE,".bed"))
   }
}
bnd <- genoFilter(VCF=bnd,top=topRc[["BND"]],metrics="RC")
bnd <- genoFilter(VCF=bnd,top=topHqCount[["BND"]],metrics="HQcount")
bnd <- genoFilter(VCF=bnd,top=topHqPercent[["BND"]],metrics="HQpercent")
SVsel <- c(SVsel, names(bnd))
if(length(bnd) > 0 ) {
    bed <- data.frame(chr=paste0("chr",seqnames(bnd)) , start=start(bnd) , end=end(bnd) , 
	              chr2=paste0("chr",info(bnd)$CHR2) , start2=info(bnd)$POS2 , end2=info(bnd)$POS2 +1 , stringsAsFactors=F)
    write.table(x=bed,quote=F,append=F,col.names=F,row.names=F,sep="\t",file=paste0(outDir,"/BND.bed"))
} else {
    system(paste0("touch " , outDir,"/",SVTYPE,".bed"))
}
#print out VCF
writeVcf(rbind(vcf,bnd)[SVsel,],filename=paste0(outDir,"/output.vcf"))


