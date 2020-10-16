#Compare genes coverage in two conditions: sample2 over sample1
#Genes found in just one sample (unique identifier) or with low avg MAPQ will be filtered
#Alternatively, if --divideByOne, then sample 2 coverage will be just divided by 1 i.e. you can visualize gene coverage of sample 2 alone

##WARNING on lowMAPQ: A gene with 0 coverage has also 0 MAPQ, so it will be filtered by minMAPQ. 
##WARNING on low MAPQ: A gene truly deleted with very few spurious reads mapping with lowMAPQ. This gene can be filtered depending on minMAPQ 

#######
#INPUT#
#######
#The default is to read cufflinks genes.fpkm_tracking files, but it can read any file as long as:
#1) The first column contains genes identifiers (regardless the header of the field)
#2) You specify the locus (--locusFieldName) and fpkm (--scoreFieldName) fields names has in the header of the input files
#3) the locus must be in the format chromosome:start-end

#########
#OPTIONS#
#########
#You can specify --repeats or --gaps to label the genes that overlap with repetitive elements or gaps. For the moment no filtering is done based on this overlap, just labelling
#you can specify a minMAPQ score to filter out genes with bad mapping reads
#you can dump in a data frame the filtered genes

########
#OUTPUT#
########
#a  dataframe recapitulating the scores per gene and 4 plots:
#PLOT1: Ratio Sample2 over Sample1. The plot shows ratio score (all genes, all chromosomes) coloring orange sample 2 enrichments, and in blue sample 2 depletions (see --highRatio and --lowRatio variables). Genes are sorted by genome position in the chromosomes, but to have sorted chromosomes you need to specify the chromosome order via the --chr variable
#PLOT2: Plot the log10 score scatterplot (all genes, all chromosomes). Genes are colored by --highRatio and --lowRatio variables as in PLOT1. Genes with log10 score > plot24_max are not shown.
#PLOT3: One panel per chromosome, score scatterplot. Genes are colored by by ratio as in PLOT1 and 2. If --scaleFree no, the  genes with extreme scores (above >plot3_max) in either condition are not shown, but you can still see them in PLOT 4 since it is log transformed 
#PLOT4: One panel per chromosome, log10 score scatterplot. Genes are colored by by ratio as in PLOT1 and 2. If --scaleFree no, genes with log10 score > plot24_max are not shown. 

#Example1: Rscript compareGeneCoverage.R --outName test.pdf --samples donovani_MWMLLN4141_P2/genes.fpkm_tracking  donovani_MWMLLN4142_P5/genes.fpkm_tracking --NAMES donovani_MWMLLN4141_P2  donovani_MWMLLN4142_P5 
#Example2: Rscript compareGeneCoverage.R --outName test.pdf --samples donovani_MWMLLN4141_P2.covPerGe  donovani_MWMLLN4142_P5.covPerGe --NAMES donovani_MWMLLN4141_P2  donovani_MWMLLN4142_P5 --chrs $CHRS --locusFieldName locus --scoreFieldName normalizedMeanCoverage --plot3_max 2 --plot3_max 5 --plot_highRatio 1.5 --plot_lowRatio 0.7 --scoreLabel coverage
#Example3: Rscript compareGeneCoverage.R  --NAMES fake donP5 --samples ${INDIR}/donovani_MWMLLN4141_P2.covPerGe.gz ${INDIR}/donovani_MWMLLN4142_P5.covPerGe.gz --outName test2 --chrs $CHRS --repeats $REPS --gaps $GAPS $OPT --divideByOne
#################################################
#		CONFIGURATION
#################################################
suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify our desired options # by default ArgumentParser will add an help option
parser$add_argument("--samples" , nargs="+", help="Two cufflinks gene FPKM outputs named genes.fpkm_tracking or any other gene coverage file (read comments)  [default %(default)s]" )
parser$add_argument("--NAMES"   , nargs="+", help="Two names of the genes.fpkm_tracking samples in the same order [default %(default)s]" )
parser$add_argument("--outName" , help="out name [default %(default)s]" , default="compareGeneCoverage")
parser$add_argument("--locusFieldName" , help="name of the locus field in the input file. Locus must be in the format chromosome:start-end [default %(default)s]" , default="locus")
parser$add_argument("--scoreFieldName" , help="name of the fpkm (or normalizedMeanCoverage) field in the input file. [default %(default)s]" , default="FPKM")
parser$add_argument("--pseudo" , type="double" , help="add a pseudo count (e.g. 0.1 FPKM or 0.1x coverage) to compute the log10 and the ratio fields of the data frame. This is to avoid genes with coverage 0 generate infinite values [default %(default)s]" , default=0.1)
parser$add_argument("--chrs"   , nargs="+", help="FILTER: List of chromosome names to be considered. This also defines the plotting order[default %(default)s]", default="NA" )
parser$add_argument("--minMAPQ" , type="integer" , help="FILTER: consider just genes with at least this MAPQ (use if there is a MAPQ score field in covPerGe.gz) (recommended 50) [default %(default)s]", default=0 )
parser$add_argument("--repeats" ,  help="repeat elements in bed3 format. genes overlapping these element will be labelled [default %(default)s]", default="NA" )
parser$add_argument("--gaps"    ,  help="gaps elements in bed3 format. genes overlapping these element will be labelled  [default %(default)s]", default="NA" )
parser$add_argument("--dumpFilteredGenes" , action="store_true" , help="set this flag to generate an extra output .df table with the filtered genes (i.e. low MAPQ, no coverage is some sample..) [default %(default)s]" , default=FALSE)
parser$add_argument("--plot_highRatio" , type="double" ,  help="graphical parameter all plots: coverage ratio above which sample 1 genes are enriched wrt sample 2 genes (colored in orange). Ignored if --significantGenes is provided [default %(default)s]" , default=2)
parser$add_argument("--plot_lowRatio"  , type="double" , help="graphical parameter all plots: coverage ratio below which sample 1 genes are depleted wrt sample 2 genes (colored in blue). Ignored if --significantGenes is provided [default %(default)s]" , default=0.5)
parser$add_argument("--plot1_ylim"    , type="double" , help="graphical parameter plot1: max ratio value before saturation [default %(default)s]" , default=5)
parser$add_argument("--saturateRatio" , help="<<OBSOLETE: plot3 and 4 now are discrete and not continuous, so there is no need to saturate>>  saturate score ratio as in plot1 for all the plots. The .df output is not affected. For plotting,  it sets score ratio > --plot1_ylim to --plot1_ylim. While in plot1 saturating the ratio is mandatory to visualize all the dots, in plots 3 and 4 the log10 of the ratio is just used to color the dots. By default points falling outside plot3_min/plot3_max and plot24_min/plot24_max are not shown and do not contribute to the color range. However points falling in the range can still have extreme ratio values and sometimes dominate the colors in the color range, so it is nice to set --saturateRatio to yes   [default %(default)s]" , default="no")
parser$add_argument("--plot3_min" , type="double" , help="graphical parameter plot3: min shown value. DEPENDENCY:--scaleFree no [default %(default)s]" , default=0)
parser$add_argument("--plot3_max" , type="double" , help="graphical parameter plot3: max shown value. DEPENDENCY:--scaleFree no [default %(default)s]" , default=100)
parser$add_argument("--plot24_min" , type="double" , help="graphical parameter plot2&4: min shown value. DEPENDENCY:--scaleFree no [default %(default)s]" , default=-1)
parser$add_argument("--plot24_max" , type="double" , help="graphical parameter plot2&4: max shown value. DEPENDENCY:--scaleFree no [default %(default)s]" , default=3)
parser$add_argument("--scoreLabel" , help="label appearing on the x and y axis of the plot refering to the coverage or the normalized gene coverage, depending on the input of this script [default %(default)s]" , default="FPKM")
parser$add_argument("--scaleFree" , help="graphical parameter plot3&4: should faceting plot be scale free? [yes|no] [default %(default)s]" , default="yes")
parser$add_argument("--divideByOne" , action="store_true" , help="set this flag to divide by one insted by the reference. It is a trick to visualize just the test sample coverage [default %(default)s]" , default=FALSE)
parser$add_argument("--significantGenes" , help="Output of sigPeaks_mixture.R. TSV file listing significant genes. If provided this will be used for the coloring and reported in the output table [default %(default)s]" , default="NA")
parser$add_argument("--geneFunctions" , help="TSV file listing geneID<Tab>function. If provided, the output will included the function field [default %(default)s]" , default="NA")
args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }

library(ggplot2)
library(IRanges)
library(GenomicRanges)

#library(session)
#save.session("session")
#quit()

#READ
allSamples <- list() 
if(divideByOne){ NAMES[1] = "__fake__" } 
for (i in  1:length(samples)){
	sample = samples[i]
	name   = NAMES[i]
	allSamples[[name]] <- as.data.frame(read.delim(sample ,header=T,stringsAsFactors=F,row.names=1))
}
allGeneIds <- unique((unlist(lapply(allSamples,row.names))))

#FILTER
#low MAPQ genes
allSamplesFiltered <- lapply(allSamples,function(x){ x[ x$MAPQ >= minMAPQ,]  })
#consider just the ids available on the both samples
ids=intersect( row.names(allSamplesFiltered[[1]]) , row.names(allSamplesFiltered[[2]])  )

#extract chromosome and start positions info
chromosome <- factor(gsub(x=allSamplesFiltered[[1]][ids,]$locus,pattern="^([^\\:]+)\\:.+",replacement="\\1",perl=T))
geneStart  <- as.numeric(gsub(x=allSamplesFiltered[[1]][ids,]$locus,pattern="^([^\\:]+)\\:(\\d+)-(\\d+)$",replacement="\\2",perl=T))
geneEnd    <- as.numeric(gsub(x=allSamplesFiltered[[1]][ids,]$locus,pattern="^([^\\:]+)\\:(\\d+)-(\\d+)$",replacement="\\3",perl=T))

#MAKE DF
if(divideByOne){
	df <- data.frame(chromosome=chromosome , geneStart=geneStart , geneEnd=geneEnd , scoreSample1=allSamplesFiltered[[2]][ids,scoreFieldName] , scoreSample2=allSamplesFiltered[[2]][ids,scoreFieldName] , ratio=allSamplesFiltered[[2]][ids,scoreFieldName] ,color="black" , log10ScoreSample1=log10(allSamplesFiltered[[2]][ids,scoreFieldName]+pseudo) , log10ScoreSample2=log10(allSamplesFiltered[[2]][ids,scoreFieldName]+pseudo) , log10ratio=log10(allSamplesFiltered[[2]][ids,scoreFieldName]) , stringsAsFactors=FALSE)
	plot1ylab=paste(scoreLabel,NAMES[2])
	plot2xlab=paste("log10 coverage",NAMES[2],"genes")
	plot3xlab=paste(scoreLabel,NAMES[2])
	plot4xlab=paste("log10",scoreLabel,NAMES[2])
} else {
	df <- data.frame(chromosome=chromosome , geneStart=geneStart , geneEnd=geneEnd , scoreSample1=allSamplesFiltered[[1]][ids,scoreFieldName] , scoreSample2=allSamplesFiltered[[2]][ids,scoreFieldName] , ratio=(allSamplesFiltered[[2]][ids,scoreFieldName]+pseudo) /  (allSamplesFiltered[[1]][ids,scoreFieldName]+pseudo) ,color="black" , log10ScoreSample1=log10(allSamplesFiltered[[1]][ids,scoreFieldName]+pseudo) , log10ScoreSample2=log10(allSamplesFiltered[[2]][ids,scoreFieldName]+pseudo) , log10ratio=log10((allSamplesFiltered[[2]][ids,scoreFieldName]+pseudo) / (allSamplesFiltered[[1]][ids,scoreFieldName]+pseudo)) , stringsAsFactors=FALSE)
	plot1ylab=paste(scoreLabel,"ratio",NAMES[2],"/",NAMES[1])
	plot2xlab=paste("log10 coverage",NAMES[1],"genes")
	plot3xlab=paste(scoreLabel,NAMES[1])
	plot4xlab=paste("log10",scoreLabel,NAMES[1])
}
row.names(df)=ids
score1 = paste0(scoreFieldName,"_",NAMES[1])
score2 = paste0(scoreFieldName,"_",NAMES[2])
log10Score1 = paste0("log10_",score1)
log10Score2 = paste0("log10_",score2)
names(df) <- c("chromosome" ,  "geneStart" , "geneEnd", score1 ,  score2 , "ratio" , "color" , log10Score1 , log10Score2 , paste("log10ratio_",scoreFieldName,sep="") )
#select chrs and sort by position
if (! is.na (chrs[1])){
	df <- df[df$chromosome %in% chrs,]
	df$chromosome <- factor(df$chromosome, levels=chrs)
}
df <- df[with(df, order(chromosome , geneStart)), ]
#color extreme genes
if(! is.na(significantGenes)){
  sigGeDf  <- read.table(significantGenes,header=T,stringsAsFactors=F,sep="\t")
  enrichedGeNames <- row.names(sigGeDf[sigGeDf$probOfTheValueOrLower > sigGeDf$probOfTheValueOrHigher,])
  depletedGeNames <- row.names(sigGeDf[sigGeDf$probOfTheValueOrLower < sigGeDf$probOfTheValueOrHigher,])
  df$status <- "background"
  df[enrichedGeNames,"status"]="enriched"
  df[depletedGeNames,"status"]="depleted"
  df[enrichedGeNames,"color"]="orange"
  df[depletedGeNames,"color"]="blue"  
} else {
  df$color[df$ratio > plot_highRatio] = "orange"
  df$color[df$ratio < plot_lowRatio]  = "blue"
  df$status <- "background"
  df$status[df$color == "orange"]="enriched"
  df$status[df$color == "blue"]="depleted"
}
df$status <-factor(df$status,levels=c("background","enriched","depleted"))
#label overlap with gaps or repeats
addOverlapAnnotations <- function(bedFile , df ){
  bed        <- read.table(bedFile,header=F,stringsAsFactors=F)
  names(bed) <- c("chr","start","end")
  bed        <- bed[bed$chr %in% df$chromosome,]
  bedGr      <- with(bed, GRanges(chr, IRanges(start, end)))
  genesGr    <- with(df, GRanges(chromosome, IRanges(geneStart, geneEnd)))
  ovAnnotations <- rep("no",length(df[,1]))
  ov  <- findOverlaps(genesGr,bedGr)
  ovAnnotations[queryHits(ov)]="yes"
  return(ovAnnotations)	
}
if (! is.na(repeats)){ df$repeatsOverlap <- addOverlapAnnotations(repeats , df) }
if (! is.na(gaps))   { df$gapsOverlap    <- addOverlapAnnotations(gaps    , df) }

png(paste0(outName,".cloud.png"),type='cairo')
#PLOT1: cloud plot (all genes)
dfSaturatedRatio = df
dfSaturatedRatio$ratio[dfSaturatedRatio$ratio >= plot1_ylim]=plot1_ylim
plot(dfSaturatedRatio$ratio , pch=19, ylim=c(0,plot1_ylim) , col=df$color , ylab=plot1ylab , xlab="gene index")
legend(fill=c("orange","blue"),x="topleft",legend=c(paste(NAMES[2],"amplified") , paste(NAMES[2],"depleted")) , cex=0.7 )
abline(h=1,col="grey",lty=2)
if (identical(saturateRatio,"yes")){
	df <- dfSaturatedRatio
}
dev.off()
pdf(paste(outName,".pdf",sep=""))
#PLOT2: scatter plot (all genes, log transformed) (genes with log10 FPKM > 3 are not shown)
plot(log10(df[[score1]]) , log10(df[[score2]]),col=df$color,pch=19 , xlab=plot2xlab , ylab=paste("log10",scoreLabel,NAMES[2],"genes") , xlim=c(plot24_min,plot24_max) , ylim=c(plot24_min,plot24_max) )
abline(0,1 , lty=2,col="gray") 
grid()
legend(fill=c("orange","blue"),x="topleft",legend=c(paste(NAMES[2],"enriched") , paste(NAMES[2],"depleted")) , cex=0.7 )
#PLOT3: one scatterplot per chromosme (extreme FPKMs over plot3_max are excluded, and do not contribute to the legend colorbar)
names(df) <- gsub(names(df) , pattern="(\\-|\\+)",replacement="_") #patch: for aes_string compatibility replace special symbols "+" or "-" in the df names with "_"
cols=c("background" = "black","enriched" = "#E69F00","depleted" = "#56B4E9")
scatterplotNormal <- ggplot(df, aes_string(x=names(df)[4], y=names(df)[5])) + geom_abline(intercept=0,slope=1,size=0.3,colour="red", linetype="solid") + geom_point(size=1,aes(colour=status)) + scale_color_manual(values=cols) + xlab(plot3xlab) + ylab(paste(scoreLabel,NAMES[2])) + theme_bw() + theme(legend.title=element_blank() , legend.text = element_text(size = 8), legend.position="bottom" , axis.text.x = element_text(angle = 90, hjust = 1) , panel.grid.major = element_line(size=0.6,colour="gray") ) ##+ scale_colour_gradient(low="orange" , name = paste("log10",scoreLabel,"ratio")) # geom_point(size=1,aes(colour=log10(ratio))) +  scale_colour_gradientn(colours = terrain.colors(10)) #  + coord_cartesian(ylim=c(plot3_min,plot3_max)
if(scaleFree == "yes"){scatterplotNormal <- scatterplotNormal + facet_wrap(~ chromosome , scales="free") } else {scatterplotNormal <- scatterplotNormal + facet_wrap(~ chromosome) + ylim(plot3_min,plot3_max) + xlim(plot3_min,plot3_max)}
print(scatterplotNormal)
#PLOT4: one scatterplot per chromosme (log transformed) (genes with log10 FPKM > plot4_max are not shown, and do not contribute to the legend colorbar)
scatterplotLogged <- ggplot(df, aes_string(x=names(df)[8], y=names(df)[9]))  + geom_abline(intercept=0,slope=1,size=0.3,colour="red", linetype="solid") + geom_point(size=1,aes(colour=status)) + scale_color_manual(values=cols) + xlab(plot4xlab) + ylab(paste("log10",scoreLabel,NAMES[2])) + theme_bw()  +  theme(legend.title=element_blank() , legend.text = element_text(size = 8), legend.position="bottom" , axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(size=0.6,colour="gray") ) #+ scale_colour_gradient(low="orange" , name = paste("log10",scoreLabel,"ratio")) #geom_point(size=1,aes(colour=log10(ratio))) + coord_cartesian(ylim=c(plot24_min,plot24_max) 
if(scaleFree == "yes"){ scatterplotLogged <- scatterplotLogged + facet_wrap(~ chromosome , scales="free") } else {scatterplotLogged <- scatterplotLogged + facet_wrap(~ chromosome) + ylim(plot24_min,plot24_max) + xlim(plot24_min,plot24_max)}
print(scatterplotLogged)
dev.off()

#OUTPUT TABLE
#rounding-off
is.num     <- sapply(df, is.numeric)
df[is.num] <- lapply(df[is.num], round, 3)
#drop meaningless columns if --divideByOne
if(divideByOne){ df <- subset(df, select = -c(normalizedMeanCoverage___fake__,ratio,log10_normalizedMeanCoverage___fake__,log10ratio_normalizedMeanCoverage) )}
if (dumpFilteredGenes){
  #loci of filtered genes
  ids = allGeneIds[! allGeneIds %in% row.names(df)]
  if (length(ids) == 0){
    system( paste0("touch ",outName,".filtered.df ; gzip ",outName,".filtered.df" ) )
  } else {
    loci <- NULL ; for (id in ids){for (n in names(allSamples)){  locus = allSamples[[n]][id,"locus"]; if(! is.na(locus)){ rbind(loci,c(id,locus))->loci   }      }  }
    loci <- data.frame(unique(loci))
    names(loci) <- c("id","locus")
    dfF <- loci
    row.names(dfF) <- loci$id
    dfF$chromosome <- factor(gsub(x=loci$locus,pattern="^([^\\:]+)\\:.+",replacement="\\1",perl=T))
    dfF$geneStart  <- as.numeric(gsub(x=loci$locus,pattern="^([^\\:]+)\\:(\\d+)-(\\d+)$",replacement="\\2",perl=T))
    dfF$geneEnd    <- as.numeric(gsub(x=loci$locus,pattern="^([^\\:]+)\\:(\\d+)-(\\d+)$",replacement="\\3",perl=T))
    dfF[[paste0("normalizedMeanCoverage_",NAMES[1])]]<- allSamples[[NAMES[1]]][as.character(dfF$id),"normalizedMeanCoverage"]
    dfF[[paste0("normalizedMeanCoverage_",NAMES[2])]]<- allSamples[[NAMES[2]]][as.character(dfF$id),"normalizedMeanCoverage"]
    dfF[[paste0("MAPQ_",NAMES[1])]]<- allSamples[[NAMES[1]]][as.character(dfF$id),"MAPQ"]
    dfF[[paste0("MAPQ_",NAMES[2])]]<- allSamples[[NAMES[2]]][as.character(dfF$id),"MAPQ"]
    dfF <- dfF[,c(-1,-2)]
    if(divideByOne){
     dfF[ ,c('normalizedMeanCoverage___fake__', 'MAPQ___fake__')] <- list(NULL)
    }
    write.table(append=F,col.names=T,file=paste(outName,".filtered.df",sep=""),quote=F,sep="\t",row.names=T,x=dfF)
    system( paste("gzip ",outName,".filtered.df",sep="") )
  }
}

if(! is.na(geneFunctions) ){
  geFunDf <- read.delim(geneFunctions , header=F , stringsAsFactors=F, col.names=c("geId","geneFunction",sep="\t"))
  geFunDf$geneFunction[is.na(geFunDf$geneFunction)]="NA"
  df$geneFunction <- geFunDf[ match(row.names(df),geFunDf$geId) , "geneFunction"]
}

write.table(append=F,col.names=T,file=paste0(outName,".df"),quote=F,sep="\t",row.names=T,x=df)
system( paste("gzip ",outName,".df",sep="") )


#############################################
#OLD PLOT FUNCTIONS WITH CONTINUOUS COLORING#
#############################################
##PLOT3: one scatterplot per chromosme (extreme FPKMs over plot3_max are excluded, and do not contribute to the legend colorbar)
#scatterplotNormal <- ggplot(df, aes_string(x=score1, y=score2)) + geom_abline(intercept=0,slope=1,size=0.3,colour="gray", linetype="dotted") + geom_point(size=1,aes(colour=log10(ratio))) + ylim(plot3_min,plot3_max) + xlim(plot3_min,plot3_max) + facet_wrap(~ chromosome) + xlab(paste(scoreLabel,NAMES[1])) + ylab(paste(scoreLabel,NAMES[2])) + theme_bw()  +  theme(legend.title = element_text(size=10) , axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradient(low="orange" , name = paste("log10",scoreLabel,"ratio")) 
##PLOT4: one scatterplot per chromosme (log transformed) (genes with log10 FPKM > plot4_max are not shown, and do not contribute to the legend colorbar)
#scatteplotLogged <- ggplot(df, aes_string(x=log10Score1, y=log10Score2))  + geom_abline(intercept=0,slope=1,size=0.3,colour="gray", linetype="dotted") + geom_point(size=1,aes(colour=log10(ratio))) + ylim(plot24_min,plot24_max) + xlim(plot24_min,plot24_max) + facet_wrap(~ chromosome) + xlab(paste("log10",scoreLabel,NAMES[1])) + ylab(paste("log10",scoreLabel,NAMES[2])) + theme_bw()  +  theme(legend.title = element_text(size=10)) + scale_colour_gradient(low="orange" , name = paste("log10",scoreLabel,"ratio")) 

