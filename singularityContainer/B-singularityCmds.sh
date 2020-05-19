#build container
sudo singularity build --disable-cache --sandbox lgertContainer.simg singularityBuild.def

#access container
sudo singularity shell --writable lgertContainer.simg/


