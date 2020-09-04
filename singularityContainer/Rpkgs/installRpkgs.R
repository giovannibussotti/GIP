pth <- "/opt/miniCRAN"

allPkgs <- read.table("/opt/Rpkgs.tsv",header=T)
alreadyInstalled <- row.names(installed.packages())
manualInstallPkgs <- c("bisoreg","XML",alreadyInstalled)
allPkgsAuto <- allPkgs[ ! allPkgs$Package %in% manualInstallPkgs, ]
allPkgsAuto$Package <- factor(allPkgsAuto$Package, levels=allPkgsAuto$Package)

#XML is needed by other pkgs like karyoplote, so you need to install first
install.packages("/opt/miniCRAN/src/contrib/XML_3.99-0.3.tar.gz", repos = NULL, type="source")

install.packages(allPkgsAuto$Package, repos = paste0("file:///", pth), type = "source")

#first you need to install all the pkgs dependencies before installing bisoreg
install.packages("/opt/miniCRAN/src/contrib/bisoreg_1.5.tar.gz", repos = NULL, type="source")
