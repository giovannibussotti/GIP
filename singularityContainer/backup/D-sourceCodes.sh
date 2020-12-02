O=termsAndConditions
mkdir -p $O/licenses
cd $O/licenses
cp ../../../docs/source/software/LICENSE_* .
cp ../../../docs/source/software/COPYING_* .
cp ../../../docs/source/software/COPYRIGHT_mummer .
cd ../..

mkdir -p $O/source_code/trf
cd $O/source_code/trf
wget https://github.com/Benson-Genomics-Lab/TRF/archive/v4.09.1.tar.gz ; cd ..
mkdir -p cdhit ; cd cdhit ; wget https://github.com/weizhongli/cdhit/archive/V4.8.1.tar.gz ; cd ..
mkdir -p bwa   ; cd bwa   ; wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz         ; cd ..
mkdir -p R     ; cd R     ; wget https://cran.r-project.org/src/base/R-3/R-3.6.0.tar.gz    ; cd ..
mkdir -p circos ; cd circos ; wget http://circos.ca/distribution/circos-0.69-9.tgz         ; cd ..
mkdir -p snpEff ; cd snpEff ; wget https://github.com/pcingola/SnpEff/archive/v4.3t.tar.gz ; cd ..
mkdir -p mummer ; cd mummer ; wget https://github.com/mummer4/mummer/archive/v4.0.0rc1.tar.gz ; cd ..
cd ../..




