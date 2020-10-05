#build container (for production)
sudo singularity build --disable-cache giptools singularityBuild.def >>giptools.build.log 2>&1

#build container (for development)
sudo singularity build --disable-cache --sandbox gipContainerDev.simg singularityBuild.def

#access container
sudo singularity shell --writable gipContainerDev.simg/

#access container mounting just the data folder (safe)
sudo singularity shell --writable --no-home -B ../inputData:/mnt gipContainerDev.simg/

#access container mounting my entire home directory (DANGEROUS, as from the container I can write/edit anything that is mounted)
#sudo singularity shell --writable --no-home -B /home/gbussott:/mnt gipContainerDev.simg/


