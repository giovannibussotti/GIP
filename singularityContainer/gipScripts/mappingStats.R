#############################################################################
# giptools                                                                  #
#                                                                           #
# Authors: Giovanni Bussotti                                                #
# Copyright (c) 2021  Institut Pasteur                                      #
#                                                                           #
#                                                                           #
# This program is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as                   #
# published by the Free Software Foundation, either version 3 of the        #
# License, or (at your option) any later version.                           #
#                                                                           #
# This program is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.    #
#                                                                           #
#############################################################################

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
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

