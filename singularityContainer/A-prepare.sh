#download a local version of most of the packages needed to build the singularity container
#storing a local varsion ofthe package is good because in the future they may be not available for download anymore. This would prevent the singularity build command to work
#I am unpacking the packages without version number, so that if we change the tool version I just need to change this script, but not the singularity definition file
#The main output is the packages.tar.gz file that is loaded by the singularity definition file

SAMTOOLSv=1.8
BWAv=0.7.17
htslibv=1.8
bedtoolsv=2.29.2
picardv=2.18.9
python3v=3.7.0
Rv=3.6.0
dellyv=0.8.7
RepeatMaskerv=4.1.0
RMBlastv=2.10.0
freebayesv=1.3.2
cdhitv=4.8.1
circosv=0.69-9
mummerv=4.0.0rc1
iqtreev=2.1.2
trfv=4.09.1

#create the output dir
mkdir -p files/
cd files

#copy miniCRAN
cp -r ../Rpkgs/miniCRAN/ .
cp ../Rpkgs/Rpkgs.tsv ../Rpkgs/installRpkgs.R .

#SAMTOOLS
wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSv}/samtools-${SAMTOOLSv}.tar.bz2
tar -xjf samtools-${SAMTOOLSv}.tar.bz2
mv samtools-${SAMTOOLSv} samtools

#BWA
wget https://github.com/lh3/bwa/releases/download/v${BWAv}/bwa-${BWAv}.tar.bz2
tar jxf bwa-${BWAv}.tar.bz2
mv bwa-${BWAv} bwa

#htslib (tabix)
wget https://github.com/samtools/htslib/releases/download/${htslibv}/htslib-${htslibv}.tar.bz2
tar -xjf htslib-${htslibv}.tar.bz2
mv htslib-${htslibv} htslib

#bedtools
#precompiled binaries give problem on old centos systems
#wget https://github.com/arq5x/bedtools2/releases/download/v${bedtoolsv}/bedtools.static.binary
#mv bedtools.static.binary bedtools
#chmod a+x bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v${bedtoolsv}/bedtools-${bedtoolsv}.tar.gz
tar -xzf bedtools-${bedtoolsv}.tar.gz 

#picard
wget https://github.com/broadinstitute/picard/releases/download/${picardv}/picard.jar

#python 3
wget https://www.python.org/ftp/python/${python3v}/Python-${python3v}.tgz
tar -xvf Python-${python3v}.tgz
mv Python-${python3v} Python3

#R
wget http://cran.rstudio.com/src/base/R-3/R-${Rv}.tar.gz
tar -xvf R-${Rv}.tar.gz
mv R-${Rv} R

#freebayes the old version 1.0.1 installed in the cluster works. The static precompiled binaries work in the ubuntu machine but not in the centos 6 cluster (kernel too old). To make it work just compile it yoursef from source in the container as usual
git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes
git checkout v$freebayesv && git submodule update --recursive
cd ..

#snpEff
#version SnpEff 4.3t (build 2017-11-24 10:18)
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip 
rm -rf snpEff_latest_core.zip

#vcftools (Latest commit d0c95c5 on Sep 2, 2018)
git clone https://github.com/vcftools/vcftools.git

#delly
wget https://github.com/dellytools/delly/releases/download/v${dellyv}/delly_v${dellyv}_linux_x86_64bit
mv delly_v${dellyv}_linux_x86_64bit delly
chmod a+x delly

#UCSC
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod a+x bedGraphToBigWig

#repeatmasker
wget http://www.repeatmasker.org/RepeatMasker-${RepeatMaskerv}.tar.gz
tar -xvf RepeatMasker-${RepeatMaskerv}.tar.gz

#TRF (downloaded manually from http://tandem.bu.edu/trf/trf409.linux64.download.html)
#no need to download he legacy file. checking with ldd --version you can verify that xenial has GLIBC 2.23
wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v${trfv}/trf409.linux64
chmod a+x trf409.linux64

#RMBlast
wget http://www.repeatmasker.org/rmblast-${RMBlastv}+-x64-linux.tar.gz
tar -xvf rmblast-${RMBlastv}+-x64-linux.tar.gz
mv rmblast-${RMBlastv} rmblast

#MUMmer 
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-${mummerv}.tar.gz
tar -xzf mummer-${mummerv}.tar.gz
mv mummer-${mummerv} mummer

#cd-hit
wget https://github.com/weizhongli/cdhit/archive/V${cdhitv}.tar.gz
tar -xzf V${cdhitv}.tar.gz
mv cdhit-${cdhitv}/ cdhit

#circos
wget http://circos.ca/distribution/circos-${circosv}.tgz
tar xvfz circos-${circosv}.tgz
mv circos-$circosv circos

#Red
wget http://toolsmith.ens.utulsa.edu/red/data/DataSet2Unix64.tar.gz
tar -xzf DataSet2Unix64.tar.gz 
mv redUnix64/Red .

#iqtree
wget https://github.com/iqtree/iqtree2/releases/download/v${iqtreev}/iqtree-${iqtreev}-Linux.tar.gz
tar -xzf iqtree-${iqtreev}-Linux.tar.gz
mv iqtree-${iqtreev}-Linux/bin/iqtree2 .

#CLEAN
rm -rf bedtools-${bedtoolsv}.tar.gz iqtree-${iqtreev}-Linux iqtree-${iqtreev}-Linux.tar.gz mummer-${mummerv}.tar.gz Python-${python3v}.tgz _tmp htslib-${htslibv}.tar.bz2 samtools-${SAMTOOLSv}.tar.bz2 bwa-${BWAv}.tar.bz2 R-${Rv}.tar.gz RepeatMasker-${RepeatMaskerv}.tar.gz v${LGERTv}.tar.gz rmblast-${RMBlastv}+-x64-linux.tar.gz clinEff V${cdhitv}.tar.gz circos-${circosv}.tgz redUnix64/ DataSet2Unix64.tar.gz

#compress
mkdir -p packages
for X in `ls | grep -v .sh$ | grep -v README | grep -v packages `; do mv $X packages/ ; done
tar -czf packages.tar.gz packages
mv packages/* .
rmdir packages

cd ..



