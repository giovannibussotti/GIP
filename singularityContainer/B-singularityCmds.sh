#build container
sudo singularity build --disable-cache --sandbox lgertContainer.simg singularityBuild.def

#access container
sudo singularity shell --writable lgertContainer.simg/

#access container mounting just the data folder (safe)
sudo singularity shell --writable --no-home -B ../inputData:/mnt lgertContainer.simg/

#access container mounting my entire home directory (DANGEROUS, as from the container I can write/edit anything that is mounted)
#sudo singularity shell --writable --no-home -B /home/gbussott:/mnt lgertContainer.simg/


