#Read the chrCoverageMedians files and plot a heatmap with the coverages per chromosome
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--dir" , help="directory storing the chrCoverageMedians files  [default %(default)s]" )
parser$add_argument("--outName" , help="outName  [default %(default)s]" , default="covPerChrSummary" )
args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }


library(gplots)
allList<-list() 

allFiles1 <- list.files(path=dir , pattern="chrCoverageMedians")
allNames1 <- gsub(x=allFiles1 , pattern="chrCoverageMedians_",replacement="")
for (n in allNames1){allList[[n]] <- read.table(paste0(dir,"/chrCoverageMedians_",n),header=T,stringsAsFactors=F)[,2]   }

#allFiles2 <- list.files(path="..",pattern="chrCoverageMedians")
#allNames2 <- gsub(x=allFiles2 , pattern="chrCoverageMedians_",replacement="")
#for (n in allNames2){allList[[n]] <- read.table(paste0("../chrCoverageMedians_",n),header=T,stringsAsFactors=F)[,2]   }

df <- as.data.frame(do.call(cbind,allList))
#df <- cbind(df[,c(allNames1)] , df[,c(allNames2)])
my_palette <- colorRampPalette(c("black", "red", "gold"))(n = 299)



pdf(paste0(outName,".pdf"),height=6,width=13)
#just new samples (+1 pseudocount)
heatmap.2(as.matrix(df[,c(allNames1)]+1),Rowv=NA,trace="none",density.info="none",margins =c(14,9),col=my_palette  ,  cexCol=1.1)

##old and new comparison
#heatmap.2(as.matrix(df+1),Rowv=NA,trace="none",density.info="none",margins =c(14,9),col=my_palette , ColSideColors = c(rep("gray", length(allNames1)),rep("blue", length(allNames2))) ,  cexCol=1.1)

#log 10 scale
#heatmap.2(log10(as.matrix(df+1)),Rowv=NA,trace="none",density.info="none",margins =c(14,9),col=my_palette , ColSideColors = c(rep("gray", length(allNames1)),rep("blue", length(allNames2))) , cexCol=1.1 )
heatmap.2(log10(as.matrix(df+1)),Rowv=NA,trace="none",density.info="none",margins =c(14,9),col=my_palette , cexCol=1.1 )
dev.off()

