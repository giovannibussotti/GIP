library("miniCRAN")
library("BiocManager")
pth <- "/mnt/miniCRAN"
allPkgs <- read.table("/mnt/Rpkgs.tsv",header=T)
system(paste("mkdir -p",pth))

manualInstallPkgs <- c("bisoreg","XML")
allPkgsAuto <- allPkgs[ ! allPkgs$Package %in% manualInstallPkgs, ]

options(repos = BiocManager::repositories())
makeRepo(allPkgsAuto$Package , path = pth, type = c("source") , Rversion="3.6.0")
addOldPackage("XML", path = pth, vers = "3.99-0.3", type = "source",repos="https://pbil.univ-lyon1.fr/CRAN" , Rversion="3.6.0")
addOldPackage("foreign", path = pth, vers = "0.8-71", type = "source",repos="https://pbil.univ-lyon1.fr/CRAN" , Rversion="3.6.0")

#addOldPackage does not work for bisoreg, so I do it manually
system(paste0("wget -P ", pth , "/src/contrib/ https://cran.r-project.org/src/contrib/Archive/bisoreg/bisoreg_1.5.tar.gz"))
updateRepoIndex(pth, type = "source", Rversion = "3.6.0")


#while bisoreg package is now present in miniCRAN, it does not have listed the dependences
#thus to make sure it installs properly first you need to install all the other packages (which will inclide its dependencies), and then at the end manually install bisoreg


