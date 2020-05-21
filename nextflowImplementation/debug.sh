#for debugging often you need to
 #1) go in the nexflow work directorey that raised the error
 #2) run this script with bash ../../../debug.sh to replace the symbolic links with the real files
find . -type l -exec bash -c "echo 'Replacing {} ...';  cp -LR '{}' '{}'.dereferenced;  rm '{}';  mv '{}'.dereferenced '{}'" \;


 #3) Then access the container with sudo singularity shell -B $PWD:/mnt --writable ../../../../singularityContainer/lgertContainer.simg/
 #4) run again the command that failed: cd /mnt ; bash .command.sh 


