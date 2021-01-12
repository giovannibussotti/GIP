suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--ALLGEANN" , required=TRUE , help="tsv file listing geneID<Tab>function [default %(default)s]" )
parser$add_argument("--SNPEFFDF" , required=TRUE , help="freebayes/snpEff df file containing the EFF field [default %(default)s]" )
parser$add_argument("--ID" , help="sample ID [default %(default)s]" , default="sampleSNVs" )
parser$add_argument("--synNames" , nargs="+" , help="snpEff effects counting as synonimous substitutions. The default is synonymous_variant , stop_retained_variant , start_retained [default %(default)s]" , default="NA")
parser$add_argument("--nonSynNames" , nargs="+" , help="snpEff effects counting as non-synonimous substitutions. The default is missense_variant, start_lost, stop_gained, stop_lost, coding_sequence_variant [default %(default)s]" , default="NA")
parser$add_argument("--debug"  , action="store_true" , help="dump session and quit [default %(default)s]" , default=FALSE)
args <- parser$parse_args()
for (n in names(args)){if(args[[n]][1] == "NA"  ){args[[n]] <- NA  } }
for (n in names(args)){assign(n,args[[n]]) }

#SETUP
if(debug){library(session);save.session("session_DEBUG_dndsRatio");quit()}
library(reshape2)
library(dplyr)
if(is.na(synNames[1])){
 synNames=c("synonymous_variant" , "stop_retained_variant" , "start_retained")
}
if(is.na(nonSynNames[1])){
 nonSynNames=c("missense_variant","start_lost","stop_gained","stop_lost","coding_sequence_variant")
}

#read table and gene function
annDF <- read.delim(ALLGEANN,header=F,col.names=c("id","ann"),stringsAsFactors = F,sep="\t")
df <- read.table(SNPEFFDF,sep="\t",header=T, stringsAsFactors=F)

#########################
#list of list the effect#
#########################
# SNV - EFFECT - 16snpEff annotations
effList <- as.list(strsplit(df$EFF, ','))
names(effList) <- paste(df$chr ,  df$position ,  df$ref_alt, sep="_")
effList <- lapply(effList,function(x){strsplit(x, '\\|')})

#################
#count dN and dS#
#################
effDf <- t(as.data.frame(effList,stringsAsFactors=F))
row.names(effDf) <- c()
effDf <- effDf[effDf[,2]!="intergenic_region",]
effDf <- data.frame(gene=effDf[,5] , type=effDf[,2] , SYN="other" , stringsAsFactors=F)
effDf$SYN[effDf$type %in% synNames]="syn"
effDf$SYN[effDf$type %in% nonSynNames]="nonSyn"
#aggregate
agg <- effDf %>% count(gene, SYN)
#cast
out <- dcast(data=agg,gene~SYN,value.var="n",fill=0)
#reformat
out$other <- NULL
names(out)[names(out)=="nonSyn"]="dN"
names(out)[names(out)=="syn"]="dS"
out <- out[,c("gene","dN","dS")]
#ratio
out$ratio <- out$dN/out$dS
out$ratio <- round(out$ratio,digits=3)
out$difference <- out$dN - out$dS
out <- out[with(out, order(-difference)), ]
out$annotation <- annDF[ match(out$gene , annDF$id) , "ann" ]


#############
#write stats#
############# 
df$EFF <- gsub(x=df$EFF,pattern=".*\\|intergenic_region\\|.*",replacement="INTERGENIC")
statFile=paste0(ID,"_cleanEFF.stats")
write(file=statFile,x=paste("total SNVs",length(df[,1])))
write(file=statFile,append=T, x=paste("intergenic",table(df$EFF == "INTERGENIC")["TRUE"]))
write(file=statFile,append=T, x=paste("genic",table(df$EFF == "INTERGENIC")["FALSE"]))
write(file=statFile,append=T, x=paste("tot genes with SNVs",length(out$gene)  ))
write(file=statFile,append=T, x=paste("the mean number of SNVs per gene is ",mean(out$dN + out$dS)  ))
posRatio <- table(out$ratio > 1)["TRUE"]
if(is.na(posRatio)) {
	posRatio <- 0
}
write(file=statFile,append=T, x=paste("the number of genes with dN/dS ratio > 1 ",posRatio  ))
negRatio <- table(out$ratio < 1)["TRUE"]
if(is.na(negRatio)) {
	negRatio <- 0
}
write(file=statFile,append=T, x=paste("the number of genes with dN/dS ratio < 1 ",negRatio  ))
eqRatio <- table(out$ratio == 1)["TRUE"]
if(is.na(eqRatio)) {
	eqRatio <- 0
}
write(file=statFile,append=T, x=paste("the number of genes with dN/dS ratio == 1 ",eqRatio  ))
write(file=statFile,append=T, x=paste("tot synonimous variations",sum(out$dS)  ))
write(file=statFile,append=T, x=paste("tot non-synonimous variations",sum(out$dN)  ))


##################
#REGENERATE TABLE#
##################
#select gene Ids
effList <- sapply(effList , function(x){
	sapply(x, function(y){
		  return(y[5])	
	})
})
#extract functions
df$gene_annotation <- sapply(effList , function(x){
	paste(annDF[match(x,annDF$id),"ann"],collapse=" ; ")
})
df$gene_annotation[df$EFF == "INTERGENIC"]="NA"
df$gene_annotation[is.na(df$gene_annotation)] = "NA"
write.table(file=paste0(ID,"_cleanEFF.tsv"),col.names=T , row.names=F,quote=F,sep="\t",x=df)



#this changes the ratio field tho character. Better do this at the end just before printing the out
out$ratio [ ! is.finite(out$ratio) ] = "infinity"
write.table(file=paste0(ID,"_dNdStables.tsv"),x=out,sep="\t",quote=F,row.names=F,col.names=T)


