URLS=(
https://github.com/samtools/samtools/archive/1.8.tar.gz 
https://github.com/ekg/freebayes/archive/v1.3.2.tar.gz
https://github.com/lh3/bwa/archive/v0.7.17.tar.gz
https://github.com/samtools/htslib/archive/1.8.tar.gz
https://github.com/arq5x/bedtools2/archive/v2.29.2.tar.gz
https://www.python.org/ftp/python/3.7.0/Python-3.7.0.tgz
https://github.com/deeptools/deepTools/archive/2.4.2.tar.gz       
http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.0.tar.gz
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ncbi-blast-2.10.0+-src.tar.gz
https://github.com/Benson-Genomics-Lab/TRF/archive/v4.09.1.tar.gz
https://github.com/mummer4/mummer/archive/v4.0.0rc1.tar.gz          
https://github.com/pcingola/SnpEff/archive/v4.3t.tar.gz
https://github.com/weizhongli/cdhit/archive/V4.8.1.tar.gz            
https://cran.r-project.org/src/base/R-3/R-3.6.0.tar.gz
http://circos.ca/distribution/circos-0.69-9.tgz
https://github.com/broadinstitute/picard/archive/2.18.9.tar.gz
https://github.com/dellytools/delly/archive/v0.8.7.tar.gz             
http://toolsmith.ens.utulsa.edu/red/data/DataSet1Src.tar.gz
https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
)


mkdir -p source_code
cd source_code

for URL in "${URLS[@]}"; do
	PKG=`basename $URL`
	wget $URL
	tar -xzf $PKG
        rm $PKG 
done


for PKG in `ls`;do 
	tar -czf ${PKG}.tar.gz $PKG
        rm -rf $PKG
done
cd ..

mkdir -p licenses
cd licenses
cp ../../docs/source/software/LICENSE_* .
cp ../../docs/source/software/COPYING_* .
cp ../../docs/source/software/COPYRIGHT_mummer .
cd ..

tar -czf licenses.tar.gz licenses
rm -rf licenses

