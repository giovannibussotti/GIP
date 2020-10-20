#to develop contrasts first you need to run here the nextflow pipeline. E.G. nextflow gip.nf --genome ../inputData/dataset/Linf_test.fa --annotation ../inputData/dataset/Linf_test.ge.gtf -c gip.config --index index.tsv -resume 

#then in gipOut you need to convert the links to files
cd /home/gbussott/Desktop/GIP_githubRepo/nextflowImplementation/gipOut
find . -type l -exec bash -c "echo 'Replacing {} ...';  cp -LR '{}' '{}'.dereferenced;  rm '{}';  mv '{}'.dereferenced '{}'" \;

#then copy gipOut in the inputData folder, that will be mounted by singularity
cd ..
sudo rm -rf ../inputData/gipOut/
cp -r gipOut/ ../inputData/gipOut/

#Then create some fake extra samples running cmdsDogX.sh
cd ../inputData
bash cmdsDogX.sh

#then run singularity
Dir=/home/gbussott/Desktop/GIP_githubRepo/
sudo singularity shell --writable --no-home -B $Dir/inputData:/mnt $Dir/singularityContainer/gipContainerDev.simg/

#once in the container you can cd to /mnt and find the gipOut folder

#run R

#manually provide the required input (e.g. samples=c("Dog3","Dog2"))

#to ease the typing you build the script on mac with Sublime. To develop/test, Copy paste code chunks from sublime into the container R session

#to backup
#-select all the sublime script and in singularity:
#   cat > /bin/binCNV
#   cp /bin/binCNV .
#-outside singularity: mv binCNV to the files/L-GERT folder
#-on TARS: cat > /scripts/binCNV  and then backup on tars

#Once the script is working you can add the parsing option in /bin/sampleComparison.sh and test its execution
#e.g. singularity exec ../singularityContainer/gipContainerDev.simg  /bin/sampleComparison.sh binCNV --gipOut gipOut/ --samples Dog1-P2  Dog3


