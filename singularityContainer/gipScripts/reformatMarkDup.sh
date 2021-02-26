#!/bin/bash
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

