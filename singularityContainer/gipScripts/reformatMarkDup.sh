#!/bin/bash
S=$1

perl -e 'open(F,"<'$S'.MarkDup.log"); while(<F>){if($_=>/\#\# METRICS CLASS/){$s=1;next;}if($_=>/^\s*$/){$s=0;} if($s==1){print}}' > ${S}.MarkDup.table
perl -e 'open(F,"<'$S'.MarkDup.log"); while(<F>){if($_=>/\#\# HISTOGRAM/){$s=1;next;}if($_=>/^\s*$/){$s=0;} if($s==1){print}}'     > ${S}.MarkDup.histData


R --slave << EOF 
df=read.table(list.files(pattern=".MarkDup.histData"),header=T)
names(df)<-c("sequencingFold","returnOfInvestment")
sampleName=gsub(x=list.files(pattern=".MarkDup.histData"),pattern=".MarkDup.histData",replacement="")
png(paste0(sampleName,".MarkDup.hist.png"),type='cairo')
plot(df[["sequencingFold"]] , df[["returnOfInvestment"]] , xlab="sequencing fold" , ylab="return of investment" , pch=19,col='black',lwd=3)
grid()
points(df[["sequencingFold"]] , df[["returnOfInvestment"]] , pch=19,col='red',lwd=1)
dev.off()
EOF

