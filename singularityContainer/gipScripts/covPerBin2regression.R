#Given one or more covPerBin files (genomic coverage binned) this script plot the fitted chromosome coverage profiles
#the script:
#1) read the covPerBin and applies filters (filterChrsNames, MAPQ, chromosome ends ..)
#2) scale (and center) all bin mean coverage (so you can plot all the profiles together)
#3) resize the vector of bin mean scores to a common value (resizeValue). This step is necessary to have all the coverage profiles of the same length (so you can plot together chromosomes of different sizes), and avoid coverage profiles too long (shring long chromosomes with many thousands bins). This step is done with the stretch function found here: https://support.bioconductor.org/p/66313/ 
#4) Fit a non-linear regression. For the moment just the smooth.splines is implemented. For the future I might try to implement loess() too. 

#######
#INPUT#
#######
#DIR: path to the covPerBin files
#SAMPLES: the covPerBin files
#NAMES:internal sample names
#filterChrsNames: list of chromosome names to be considered
#degreesOfFreedom: degrees of freedom to fit the non linear function (e.g. smooth splines, or loess)
#outName: output name
#cv: use leave one out cross validation to define smoothness
#maxCov: max scaled coverage. Values above this are sent to this value

########
#OUTPUT#
########
#outName.pdf: fitted chromosome coverage profiles

#Example: Rscript covPerBin2regression.R --DIR ./p2p5/bwa/genomecov/allPositions/covPerBin/ --SAMPLES infantum_02A_P2.covPerBin.gz infantum_03A_P5.covPerBin.gz --NAMES Linf_02A_P2 Linf_02A_P5 --filterChrsNames 1 2 3 4 5 6 7 8 9 10 --degreesOfFreedom 10 --outName test --lineColor black --resizeValue 1000 --filterChrEnds 1
#################################################
#		CONFIGURATION
#################################################
suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify our desired options # by default ArgumentParser will add an help option
parser$add_argument("--DIR"              , type="character" , help="path to the covPerBin folder [default %(default)s]")
parser$add_argument("--SAMPLES"          , type="character" , nargs="+", help="covPerBin files [default %(default)s]" )
parser$add_argument("--NAMES"            , type="character" , nargs="+", help="internal covPerBin file names [default %(default)s]" )
parser$add_argument("--cv"               , action="store_true" , help="set this flag to use leave one out Cross Validation defining the optimal degrees of freedom [default %(default)s]" , default=FALSE)
parser$add_argument("--degreesOfFreedom" , type="integer"   , help="manually specify spline degrees of fredom. Used if cross validation is false (i.e. --cv is unset) [default %(default)s]", default=10)
parser$add_argument("--outName"          , type="character" , help="out name [default %(default)s]" , default="covPerBinFitted")
parser$add_argument("--lineColor"        , type="character" , help="unique colour for all chromosomes of all samples [default %(default)s]" , default="NA")
parser$add_argument("--resizeValue"      , type="integer"   , help="stretch/compress coverage profiles to this value [default %(default)s]", default=1000)
parser$add_argument("--maxCov"                              , help="max scaled coverage score [default %(default)s]", default="NA")
parser$add_argument("--alpha"            , type="integer"   , help="GRAPHIC: color intensity [1 to 255] [default %(default)s]", default=255)
parser$add_argument("--filterChrsNames"  , type="character" , nargs="+", help="FILTER: List of chromosome names to be considered [default %(default)s]")
parser$add_argument("--filterChrEnds"    , type="integer"   , help="FILTER: bins to remove from chr extremities before fitting. Useful if telomeres have very high coverage. Removing them improves the overall fitting  (recommended 1) [default %(default)s]", default=0)
parser$add_argument("--filterMAPQ"       , type="integer"   , help="FILTER:in covPerBin each bin has a average MAPQ score. Use this filter to discard lowMAPQ bins which therefore won't participate to the fitting [default %(default)s]", default=0)
args <- parser$parse_args()
#patch NA
for (n in names(args)){if(args[[n]][1] == "NA"){args[[n]] <- NA}  }
for (n in names(args)){assign(n,args[[n]]) }


library(IRanges)
library(RColorBrewer)
#stretch <- function(r, length.out=length(x))
#{
#    stopifnot(length.out > 0)
#    r <- as(r, "Rle")
#    length.in <- length(r)
#    runLength(r) <- length.out * runLength(r)
#    Rle(mean(Views(r, breakInChunks(length(r), length.in))))
#}
stretch <- function(x, length.out) {
 length.in <- length(x)
 Rle(colMeans(matrix(rep(x, each=length.out), length.in)))
}

addTransparency = function(colors,alpha){
	newcolors = NULL 
	for(col in colors){
		z = col2rgb(col) ; 
		newcolors = c(newcolors, rgb(z[1], z[2], z[3], alpha, maxColorValue=255)) 
	}
	return(newcolors)
}

#library(session)
#save.session("session")
#quit()
######
#READ#
######
allSamplesList      <- list()
for (i in 1:length(SAMPLES)){
	s <- SAMPLES[i] 
	n <- NAMES[i]
	allSamplesList[[n]] <- read.table(paste0(DIR,"/",s),header=T,stringsAsFactors=F)
	#filters
	allSamplesList[[n]] <- allSamplesList[[n]][allSamplesList[[n]]$chromosome %in% filterChrsNames,]
	allSamplesList[[n]]$sample <- n
	if (filterChrEnds > 0){
		allSamplesList[[n]] <- allSamplesList[[n]][-filterChrEnds:-1,]
		r                   <- nrow(allSamplesList[[n]])
		allSamplesList[[n]] <- allSamplesList[[n]][-(r-(filterChrEnds-1)):-r,]
	}
	allSamplesList[[n]] <- allSamplesList[[n]][ allSamplesList[[n]]$MAPQ >= filterMAPQ, ]
}

#######
#SCALE#
#######
allSamples      <- do.call(rbind,allSamplesList)
allSamples$mean <- scale(allSamples$mean)
if(! is.na(maxCov)){ maxCov=as.numeric(maxCov); allSamples$mean [allSamples$mean > maxCov] = maxCov  }
allSamplesList  <- split(allSamples,f=allSamples$sample)

#############
#REGRESSIONS#
#############
allRegressions <- list()
for (n in NAMES){
	for (chr in filterChrsNames){
		tmp       <- allSamplesList[[n]][allSamplesList[[n]]$chromosome == chr,]
		if (cv){
			allRegressions[[n]][[chr]][["regression"]] <- smooth.spline(tmp$start,tmp$mean,cv=TRUE)
			allRegressions[[n]][[chr]][["regressionStretched"]] <- smooth.spline(1:resizeValue,stretch(tmp$mean,resizeValue),cv=TRUE)
		} else {
			allRegressions[[n]][[chr]][["regression"]] <- smooth.spline(tmp$start,tmp$mean,df=degreesOfFreedom)
                        allRegressions[[n]][[chr]][["regressionStretched"]] <- smooth.spline(1:resizeValue,stretch(tmp$mean,resizeValue),df=degreesOfFreedom)
		}		
		allRegressions[[n]][[chr]][["derivative"]] <- predict(allRegressions[[n]][[chr]][["regression"]],deriv=1)
		allRegressions[[n]][[chr]][["derivativeStretched"]] <- predict(allRegressions[[n]][[chr]][["regressionStretched"]],deriv=1)
	}
}

######
#PLOT#
######
#max and min
maxs=c();
mins=c(); 
for (n in NAMES) {
	for (chr in filterChrsNames) { 
		maxs=c(maxs,max(allRegressions[[n]][[chr]][["regressionStretched"]]$y)) 
		mins=c(mins,min(allRegressions[[n]][[chr]][["regressionStretched"]]$y))  
	}  
}
ma = max(maxs)
mi = min(mins)

#color
if (is.na(lineColor)){
	myColors <- brewer.pal(8,"Dark2")
} else {
	myColors <- rep(lineColor, length(NAMES))
}
myColors <- addTransparency(myColors,alpha)

pdf(paste0(outName,".pdf"))
plot(1, type="n", ,xlim=c(1,resizeValue),ylim=c(mi,ma),ylab="sequencing coverage")
for (i in 1:length(NAMES)) { 
		n   <- NAMES[i] 
		col <- myColors[i]  
		for (chr in filterChrsNames){ 
			lines(allRegressions[[n]][[chr]][["regressionStretched"]], col=col) 
		}
}
dev.off()

