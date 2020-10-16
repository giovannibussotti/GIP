#Given a list of bams, this script generate a table summarizing mapping statistics using picard CollectAlignmentSummaryMetrics
#Tested with CollectAlignmentSummaryMetrics Version: 1.94(1484) and R version 3.2.2
#example: Rscript stats.R --bams sample1.bam sample2.bam --dir testBam/ --assembly assembly.toplevel.fa

#################################################
#		CONFIGURATION
#################################################
suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify our desired options # by default ArgumentParser will add an help option
parser$add_argument("--bams" , nargs="+", help="List of bam files to load [default %(default)s]" )
parser$add_argument("--dir"  , help="Directory containing the bam files[default %(default)s]" )
parser$add_argument("--assembly"  , help="reference assembly multifasta file[default %(default)s]" )
parser$add_argument("--type"  , help="type of statistics to extract. PAIR FIRST_OF_PAIR or SECOND_OF_PAIR[default %(default)s]" , default="PAIR")
parser$add_argument("--outName" , help="merged out name. Give NA if you do not want to have a merged output [default %(default)s]" , default="out.stats")
parser$add_argument("--fields" , nargs="+", help="fields to extract [default %(default)s]",default=c("TOTAL_READS","PF_READS","PF_READS_ALIGNED","PCT_PF_READS_ALIGNED","PF_ALIGNED_BASES","PF_MISMATCH_RATE","PF_INDEL_RATE","MEAN_READ_LENGTH","READS_ALIGNED_IN_PAIRS","PCT_READS_ALIGNED_IN_PAIRS","STRAND_BALANCE","PCT_CHIMERAS","PCT_ADAPTER") )
parser$add_argument("--CollectAlignmentSummaryMetrics" , help="name CollectAlignmentSummaryMetrics [default %(default)s]" , default="CollectAlignmentSummaryMetrics")
parser$add_argument("--tmpDir"  , help="TMP dir to run picard CollectAlignmentSummaryMetrics [default %(default)s]", default="NA" )

args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"  ){args[[n]] <- NA  } }
for (n in names(args)){assign(n,args[[n]]) }

#####
#RUN#
#####
#CollectAlignmentSummaryMetrics
allStats <- list()
for (bam in bams){
    if (is.na(tmpDir)){
	tmpDir = paste0("_tmpDirCollectAlignmentSummaryMetrics_",bam)	
    }
    system(paste("mkdir -p",tmpDir))
    out=gsub(x=bam,pattern=".bam$",replacement=".stats")
    print(paste0(CollectAlignmentSummaryMetrics," INPUT=",dir,"/",bam," OUTPUT=",dir,"/",out," REFERENCE_SEQUENCE=",assembly," TMP_DIR=",tmpDir))
    system(paste0(CollectAlignmentSummaryMetrics," INPUT=",dir,"/",bam," OUTPUT=",dir,"/",out," REFERENCE_SEQUENCE=",assembly," TMP_DIR=",tmpDir))
    system(paste("rm -rf ",tmpDir))
    allStats[[out]] <- read.table(paste0(dir,"/",out),skip=6,header=T,fill=T,stringsAsFactors=F,row.names=1)[,fields]
}

#merge
if (! is.na (outName) ){
	df <- do.call(rbind,lapply(allStats,function(x){x[type,]}))
	row.names(df)<-gsub(row.names(df),pattern=".stats",replacement="")
	write.table(df,file=outName,sep="\t",append="F",quote=F,col.names=T,row.names=T)
}

