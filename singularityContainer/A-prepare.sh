#download a local version of most of the packages needed to build the singularity container
#storing a local varsion ofthe package is good because in the future they may be not available for download anymore. This would prevent the singularity build command to work
#I am unpacking the packages without version number, so that if we change the tool version I just need to change this script, but not the singularity definition file
#here I download L-GERT (jobarray version) from github, then I clean up the bashFunctions file, I replace the MarkDuplicates and GATK commands, I export the functions to a subshell with typeset (very important, otherwise the functions are not available), I remove unneded scripts 
#The main output is the packages.tar.gz file that is loaded by the singularity definition file

#IMPORTANT!!! rpkginstall leishmaniaAndAncestralSharedRepeats.fa GenomeAnalysisTK.jar A-prepareGenome.sh must be downloaded/provided manually and available from a local folder named: manuallyDownlodedPkgs/ 
LGERTv=1.0
SAMTOOLSv=1.8
BWAv=0.7.17
htslibv=1.8
bedtoolsv=2.25.0
picardv=2.18.9
python3v=3.7.0
Rv=3.3.2
dellyv=0.6.7
RepeatMaskerv=4.1.0
RMBlastv=2.10.0
freebayesv=1.3.1

#create the output dir
mkdir -p files/
cd files

#copy manually downloded packages
cp -r ../manuallyDownlodedPkgs/* .

#L-GERT_jobArrays: download and convert 
wget https://github.com/giovannibussotti/L-GERT_jobArrays/archive/v${LGERTv}.tar.gz
tar -xzf v${LGERTv}.tar.gz
mv L-GERT-${LGERTv} L-GERT
#strip off unwanted file header
perl -e 'open(F,"L-GERT/LSD/bashFunctions.sh");$spy=0; while(<F>){if(($_=~/^function bwaMapSample/)or($spy==1)){print $_; $spy=1;}  } close F; ' > L-GERT/LSD/_tmp
mv L-GERT/LSD/_tmp L-GERT/LSD/bashFunctions.sh
#replace MarkDuplicates command
cat L-GERT/LSD/bashFunctions.sh | sed -e 's/MarkDuplicates /java -jar \/bin\/picard.jar MarkDuplicates /' > L-GERT/LSD/_tmp
mv L-GERT/LSD/_tmp L-GERT/LSD/bashFunctions.sh
#remove the GATK key parameters
cat L-GERT/LSD/bashFunctions.sh | sed -e 's/-et NO_ET -K $GATKKEY//' > L-GERT/LSD/_tmp
mv L-GERT/LSD/_tmp L-GERT/LSD/bashFunctions.sh
# mandatory to export a function to a subshell
# at the end of sourcing environment a subshell is launch by singularity
# exec /bin/bash --norc "$@"
grep ^function  L-GERT/LSD/bashFunctions.sh | sed -e 's/function/typeset -fx/' | sed -e 's/{//' > _tmp ; cat _tmp >> L-GERT/LSD/bashFunctions.sh
#clean and sort files
mv L-GERT/LSD/utility/covPerClstr.sh     L-GERT/
mv L-GERT/LSD/utility/covPerChrSummary.R L-GERT/
mv L-GERT/LSD/utility/binCoverage2cnvs.R L-GERT/
mv L-GERT/otherScripts/*                 L-GERT/
mv L-GERT/LSD/*                          L-GERT/
rm -rf L-GERT/testDeNovo.tsv L-GERT/test.tsv L-GERT/README.md L-GERT/LSDschema.png L-GERT/CHANGES.txt L-GERT/cleanLSD.sh L-GERT/test/ L-GERT/runSpades.sh L-GERT/S-recycler_noModule.sh L-GERT/S-recycler.sh L-GERT/Q-GAT.sh L-GERT/P-redundans.sh L-GERT/N-blastnContigsToRef.sh L-GERT/M-cosmidSeqMappingStats.sh L-GERT/A-collapseMultirunSamples.R L-GERT/I-spades.sh L-GERT/utility L-GERT/otherScripts/ L-GERT/sequenceOfPipelineTasks.txt L-GERT/runLSD L-GERT/runLSD L-GERT/configLSD L-GERT/description.docx L-GERT/LSD L-GERT/.gitignore
mv L-GERT/giovanniLibrary.pl L-GERT/customPerlLib.pl
#replace file paths
cat L-GERT/gtf2bed12.sh | sed -e "s/\/pasteur\/entites\/HubBioIT\/gio\/apps\/my_scripts\/UCSC\///"  > _tmp ; mv _tmp L-GERT/gtf2bed12.sh
cat L-GERT/exonGTF_2_fasta.sh | sed -e "s/\/pasteur\/entites\/HubBioIT\/gio\/apps\/my_scripts\/GTF_manipulator_package\///"  > _tmp ; mv _tmp L-GERT/exonGTF_2_fasta.sh 
cat L-GERT/exonGTF_2_transcriptGTF.pl | sed -e "s/\/pasteur\/entites\/HubBioIT\/gio\/apps\/my_scripts\/Scripts_Matthias\/Perl\/Lib\//\/bin\//"  > _tmp ; mv _tmp L-GERT/exonGTF_2_transcriptGTF.pl
cat L-GERT/exonGTF_2_geneGTF.pl       | sed -e "s/\/pasteur\/entites\/HubBioIT\/gio\/apps\/my_scripts\/Scripts_Matthias\/Perl\/Lib\//\/bin\//"  > _tmp ; mv _tmp L-GERT/exonGTF_2_geneGTF.pl
cat L-GERT/exonGTF_2_transcriptGTF.pl | sed -e "s/giovanniLibrary.pl/customPerlLib.pl/"  > _tmp ; mv _tmp L-GERT/exonGTF_2_transcriptGTF.pl
cat L-GERT/exonGTF_2_geneGTF.pl       | sed -e "s/giovanniLibrary.pl/customPerlLib.pl/"  > _tmp ; mv _tmp L-GERT/exonGTF_2_geneGTF.pl
chmod a+x L-GERT/*
#replace variable names
for X in `ls L-GERT/`; do cat L-GERT/$X | sed -e 's/\$SAMTOOLS/samtools/' |  sed -e "s/\'samtools\'/samtools/" > _tmp ; mv _tmp L-GERT/$X ; done
for X in `ls L-GERT/`; do cat L-GERT/$X | sed -e 's/\$BEDTOOLS/bedtools/' |  sed -e "s/\'bedtools\'/bedtools/" > _tmp ; mv _tmp L-GERT/$X ; done
for X in `ls L-GERT/`; do cat L-GERT/$X | sed -e 's/\$BWA/bwa/' > _tmp ; mv _tmp L-GERT/$X ; done
for X in `ls L-GERT/`; do cat L-GERT/$X | sed -e 's/\$GATK/\/bin\/GenomeAnalysisTK.jar/' > _tmp ; mv _tmp L-GERT/$X ; done
for X in `ls L-GERT/`; do cat L-GERT/$X | sed -e 's/\$TABIX/tabix/' > _tmp ; mv _tmp L-GERT/$X ; done
for X in `ls L-GERT/`; do cat L-GERT/$X | sed -e 's/\$BGZIP/bgzip/' > _tmp ; mv _tmp L-GERT/$X ; done
cat L-GERT/H-freebayes.sh | sed -e "s/\$FREEBAYES /freebayes /"  > _tmp ; mv _tmp L-GERT/H-freebayes.sh
#add the script to pre-process the input which is not part of the github L-GERT_jobArrays version
mv A-prepareGenome.sh L-GERT/

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
wget https://github.com/arq5x/bedtools2/releases/download/v${bedtoolsv}/bedtools-${bedtoolsv}.tar.gz
tar -zxvf bedtools-${bedtoolsv}.tar.gz

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

#GATK
#The newest GATK 4 does not use indelrealigner, so we need to stick to version 3
#But installing from source does not work
#wget https://github.com/broadgsa/gatk-protected/archive/3.5.tar.gz ; tar -xzf 3.5.tar.gz ; cd gatk-3.5
#Cloning the git repository and then using ant for the installation installation does not work
#git clone https://github.com/Frogee/gatk-protected.git ; cd gatk-protected ; ant
#copying the GenomeAnalysisTK.jar compiled binary does not work (you cannot wget directly the URL)
#So I copied GenomeAnalysisTK.jar to the shiny VM of the HUB with scp: http://www.hypexr.org/linux_scp_help.php
#To copy the file to the remote shiny VM
#scp GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2 gbussotti@shiny01.hosting.pasteur.fr:/home/gbussotti/tools/
#Copy the file from the remoty shiny VM to the local container (works but you have to add the password)
#scp gbussotti@shiny01.hosting.pasteur.fr:/home/gbussotti/tools/GenomeAnalysisTK.jar .
#or you can get it from the pasteur FTP (but you have to re-upload the file every 40 days)  
#wget http://dl.pasteur.fr/fop/KFC7iGSz/GenomeAnalysisTK.jar

#freebayes
wget https://github.com/ekg/freebayes/releases/download/v${freebayesv}/freebayes-v${freebayesv}
mv freebayes-v${freebayesv} freebayes
chmod +x freebayes

#snpEff
#version SnpEff 4.3t (build 2017-11-24 10:18)
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip 
rm -rf snpEff_latest_core.zip

#spades (version 3.13.0) ###NOT USED ANYMORE
#wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz
#tar -xzf SPAdes-3.13.0-Linux.tar.gz
#rm -rf SPAdes-3.13.0-Linux.tar.gz
#mv SPAdes-3.13.0-Linux SPAdes

#trimmomatics ###NOT USED ANYMORE
#wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
#unzip Trimmomatic-0.38.zip 
#mv Trimmomatic-0.38 Trimmomatic
#rm -rf Trimmomatic-0.38.zip 

#vcftools (Latest commit d0c95c5 on Sep 2, 2018)
git clone https://github.com/vcftools/vcftools.git

#delly
wget https://github.com/dellytools/delly/releases/download/v${dellyv}/delly_v${dellyv}_linux_x86_64bit
mv delly_v${dellyv}_linux_x86_64bit delly
chmod a+x delly

#UCSC
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod a+x bedGraphToBigWig genePredToBed gtfToGenePred

#repeatmasker
wget http://www.repeatmasker.org/RepeatMasker-${RepeatMaskerv}.tar.gz
tar -xvf RepeatMasker-${RepeatMaskerv}.tar.gz

#TRF (downloaded manually from http://tandem.bu.edu/trf/trf409.linux64.download.html)
#no need to download he legacy file. checking with ldd --version I can see that xenial has GLIBC 2.23

#RMBlast
wget http://www.repeatmasker.org/rmblast-${RMBlastv}+-x64-linux.tar.gz
tar -xvf rmblast-${RMBlastv}+-x64-linux.tar.gz
mv rmblast-${RMBlastv} rmblast

#look4TRs ###NOT USED ANYMORE
#git clone https://github.com/TulsaBioinformaticsToolsmith/Look4TRs.git
#mkdir -p Look4TRs/bin
#cd Look4TRs/bin
#cmake ..
#make

#Red  ###NOT USED ANYMORE

#MUMmer (version MUMmer3.23.tar.gz, manual download from https://sourceforge.net/projects/mummer/files/latest/download)

#leishmaniaAndAncestralSharedRepeats.fa is manually copied from repBase

#CLEAN
rm -rf Python-${python3v}.tgz _tmp htslib-${htslibv}.tar.bz2 samtools-${SAMTOOLSv}.tar.bz2 bwa-${BWAv}.tar.bz2 L-GERT/LSD/test v${LGERTver}.tar.gz bedtools-${bedtoolsv}.tar.gz R-${Rv}.tar.gz RepeatMasker-${RepeatMaskerv}.tar.gz v${LGERTv}.tar.gz rmblast-${RMBlastv}+-x64-linux.tar.gz clinEff

#compress
mkdir -p packages
for X in `ls | grep -v .sh$ | grep -v README | grep -v packages `; do mv $X packages/ ; done
tar -czf packages.tar.gz packages
mv packages/* .
rmdir packages

cd ..



