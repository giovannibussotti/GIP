suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--ALLGEANN" , required=TRUE , help="tsv file listing geneID<Tab>function [default %(default)s]" )
parser$add_argument("--SNPEFFDF" , required=TRUE , help="freebayes/snpEff df file containing the EFF field [default %(default)s]" )
parser$add_argument("--ID" , help="sample ID [default %(default)s]" , default="sampleSNVs" )
parser$add_argument("--synNames" , nargs="+" , help="snpEff effects counting as synonimous substitutions. The default is SYNONYMOUS_CODING and SYNONYMOUS_STOP [default %(default)s]" , default="NA")
parser$add_argument("--nonSynNames" , nargs="+" , help="snpEff effects counting as non-synonimous substitutions. The default is NON_SYNONYMOUS_CODING NON_SYNONYMOUS_START START_LOST STOP_GAINED STOP_LOST [default %(default)s]" , default="NA")
parser$add_argument("--debug"  , action="store_true" , help="dump session and quit [default %(default)s]" , default=FALSE)
args <- parser$parse_args()
for (n in names(args)){if(args[[n]][1] == "NA"  ){args[[n]] <- NA  } }
for (n in names(args)){assign(n,args[[n]]) }

if(debug){library(session);save.session("session_DEBUG");quit()}

#all gene annotations
annDF <- read.delim(ALLGEANN,header=F,col.names=c("id","ann"),stringsAsFactors = F,sep="\t")

#clean UPSTREAM and DOWNSTREAM tags
df <- read.table(SNPEFFDF,sep="\t",header=T, stringsAsFactors=F)
df$EFF[grepl(x=df$EFF,pattern="^UPSTREAM")]="INTERGENIC"
df$EFF[grepl(x=df$EFF,pattern="^DOWNSTREAM")]="INTERGENIC"
df$EFF[grepl(x=df$EFF,pattern="^INTERGENIC")]="INTERGENIC"

#build table with all the genes with SNVs
genesWithSNVs <-  unique(gsub(x=df$EFF , pattern="^.+\\([^\\|]*\\|[^\\|]*\\|[^\\|]*\\|[^\\|]*\\|[^\\|]*\\|[^\\|]*\\|[^\\|]*\\|([^\\|]*)\\|.+$" , replacement="\\1"))
out <- data.frame(gene=genesWithSNVs[ ! genesWithSNVs %in% "INTERGENIC"],dN=0,dS=0,stringsAsFactors = F)

#simplify the EFF column
df$EFFclean <- gsub(x=df$EFF,pattern="^(.+)\\(.*",replacement="\\1")
if(is.na(synNames[1])){
 synNames=c("SYNONYMOUS_CODING","SYNONYMOUS_STOP")
}
if(is.na(nonSynNames[1])){
 nonSynNames=c("NON_SYNONYMOUS_CODING","NON_SYNONYMOUS_START","START_LOST","STOP_GAINED","STOP_LOST")
}
#count dN and dS
for (gene in out$gene){
  #extract all SNVs for each gene
  tmp <- df$EFFclean[ grep(pattern=gene , x=df$EFF) ]
  out [ which(out$gene == gene) , "dN" ] = table(tmp %in% nonSynNames)["TRUE"]
  out [ which(out$gene == gene) , "dS" ] = table(tmp %in% synNames)["TRUE"]
}
#ratio
out[is.na(out)]=0
out$ratio <- out$dN/out$dS
out$ratio <- round(out$ratio,digits=3)
out$difference <- out$dN - out$dS
out <- out[with(out, order(-difference)), ]
out$annotation <- annDF[ match(out$gene , annDF$id) , "ann" ]

#wirte stats
statFile=paste0(ID,"_cleanEFF.stats")
write(file=statFile,x=paste("total SNVs",length(df[,1])))
write(file=statFile,append=T, x=paste("intergenic",table(df$EFFclean == "INTERGENIC")["TRUE"]))
write(file=statFile,append=T, x=paste("genic",table(df$EFFclean == "INTERGENIC")["FALSE"]))
write(file=statFile,append=T, x=paste("tot genes with SNVs",length(out$gene)  ))
write(file=statFile,append=T, x=paste("the mean number of SNVs per gene is ",mean(out$dN + out$dS)  ))
write(file=statFile,append=T, x=paste("the number of genes with dN/dS ratio > 1 ",table(out$ratio > 1)["TRUE"]  ))
write(file=statFile,append=T, x=paste("the number of genes with dN/dS ratio < 1 ",table(out$ratio < 1)["TRUE"]  ))
write(file=statFile,append=T, x=paste("the number of genes with dN/dS ratio == 1 ",table(out$ratio == 1)["TRUE"]  ))
write(file=statFile,append=T, x=paste("tot synonimous variations",table(df$EFFclean %in% synNames)["TRUE"]  ))
write(file=statFile,append=T, x=paste("tot non-synonimous variations",table(df$EFFclean %in% nonSynNames)["TRUE"]  ))

df$EFFclean <- NULL
df$gene_annotation[ df$EFF == "INTERGENIC" ]="NA"
write.table(file=paste0(ID,"_cleanEFF.tsv"),col.names=T , row.names=F,quote=F,sep="\t",x=df)

out$ratio [ ! is.finite(out$ratio) ] = "infinity"

write.table(file=paste0(ID,"_dNdStables.tsv"),x=out,sep="\t",quote=F,row.names=F,col.names=T)


