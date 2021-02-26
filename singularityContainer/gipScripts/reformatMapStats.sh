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

perl -e 'open(F,"<'$S'.alignmentMetrics"); while(<F>){if($_=>/\#\# METRICS CLASS/){$s=1;next;}if($_=>/^\s*$/){$s=0;} if($s==1){print}}' > ${S}.alignmentMetrics.table
perl -e 'open(F,"<'$S'.insertSize.metrics"); while(<F>){if($_=>/\#\# METRICS CLASS/){$s=1;next;}if($_=>/^\s*$/){$s=0;} if($s==1){print}}' > ${S}.insertSize.table
perl -e 'open(F,"<'$S'.insertSize.metrics"); while(<F>){if($_=>/\#\# HISTOGRAM/){$s=1;next;}if($_=>/^\s*$/){$s=0;} if($s==1){print}}' > ${S}.insertSize.histData


R --slave << EOF 
library(ggplot2)
library(reshape2)
df=read.table(list.files(pattern=".insertSize.histData"),header=T)
names(df) <- gsub(x=names(df),pattern="All_Reads.fr_count",replacement="FR")
names(df) <- gsub(x=names(df),pattern="All_Reads.rf_count",replacement="RF")
names(df) <- gsub(x=names(df),pattern="All_Reads.tandem_count",replacement="TANDEM")
mdf <- melt(df,id.vars="insert_size")
sampleName=gsub(x=list.files(pattern=".insertSize.histData"),pattern=".insertSize.histData",replacement="")
options(scipen=10000)
png(paste0(sampleName,".insertSize.hist.png"),type='cairo',width = 2000, height = 500 , res=300); 
p <- ggplot(mdf,aes(x=insert_size , y=value , fill=variable)) + geom_area(alpha=1)  +  scale_fill_manual(" ", values=c(FR="#E69F00", RF="#999999", TANDEM="#000000")) + labs(x="insert size",y="read count") + theme_minimal() +  scale_x_log10()
print(p);dev.off()
EOF
