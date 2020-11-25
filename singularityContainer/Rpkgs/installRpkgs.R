pth <- "/opt/miniCRAN"

allPkgs <- read.table("/opt/Rpkgs.tsv",header=T,stringsAsFactors=F)
alreadyInstalled <- row.names(installed.packages())
manualInstallPkgs <- c("bisoreg","XML",alreadyInstalled)
allPkgsAuto <- allPkgs[ ! allPkgs$Package %in% manualInstallPkgs, ]

#XML is needed by other pkgs like karyoplote, so you need to install first
#install.packages("/opt/miniCRAN/src/contrib/XML_3.99-0.3.tar.gz", repos = NULL, type="source") #now it sould be fixed by addOldPackage()

install.packages(allPkgsAuto$Package, repos = paste0("file:///", pth), type = "source")

#first you need to install all the pkgs dependencies before installing bisoreg
#since you add it manually its dependencies are not defined, so you need to install it separatelly at the end
install.packages("/opt/miniCRAN/src/contrib/bisoreg_1.5.tar.gz", repos = NULL, type="source")

#install extra packages required for ggtree (added in a second moment)
#maintain this installation order to honor dependencies
install.packages("/opt/miniCRAN/src/contrib/patchwork_1.1.0.tar.gz", repos = NULL, type="source")
install.packages("/opt/miniCRAN/src/contrib/aplot_0.0.6.tar.gz", repos = NULL, type="source")
install.packages("/opt/miniCRAN/src/contrib/rvcheck_0.1.8.tar.gz", repos = NULL, type="source")
install.packages("/opt/miniCRAN/src/contrib/cpp11_0.2.4.tar.gz", repos = NULL, type="source")
install.packages("/opt/miniCRAN/src/contrib/tidyr_1.1.2.tar.gz", repos = NULL, type="source")
install.packages("/opt/miniCRAN/src/contrib/tidytree_0.3.3.tar.gz", repos = NULL, type="source")
install.packages("/opt/miniCRAN/src/contrib/treeio_1.14.3.tar.gz", repos = NULL, type="source")
install.packages("/opt/miniCRAN/src/contrib/ggtree_2.4.1.tar.gz", repos = NULL, type="source")


