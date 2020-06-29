#for debugging often you need to
 #1) go in the nexflow work directorey that raised the error
 #2) run this script with bash ../../../debug.sh to replace the symbolic links with the real files
find . -type l -exec bash -c "echo 'Replacing {} ...';  cp -LR '{}' '{}'.dereferenced;  rm '{}';  mv '{}'.dereferenced '{}'" \;


 #3) Then access the container with 
 sudo singularity shell -B $PWD:/mnt --writable ../../../../singularityContainer/gipContainerDev.simg/
 
 #4) run again the command that failed: cd /mnt ; bash .command.sh 

 #NOTE: in nextflow, inside the process script part, you can access variables that are either 1) specified in the nexflow script but outside the process (e.g. the ones you define as command line params options) or 2) specified in the proces input part (e.g. file (annotation) ).
#in 1) you can just do $annotation, and if this is a file this is the path to the file.
#in 2) it copies a symbolic link of that file in the local work folder.
# To debug a process that failed it is best to specify all the needed input files as in 2). By doing that you have a local symbolic link that thanks to this debugging script is replaced with a copy or the real file. Then accessing singularity I can reproduce the error with all the files I need.
#On the contrary, in 1) I would have the problem that the file $annotation is not part of the local work folder, therefore not mounted in singularity
