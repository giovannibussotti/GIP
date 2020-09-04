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

#addOldPackage does not work, so I do it manually
system(paste0("wget -P ", pth , "/src/contrib/ https://cran.r-project.org/src/contrib/Archive/bisoreg/bisoreg_1.5.tar.gz"))
#system(paste0("wget -P ", pth , "/src/contrib/ https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz"))
updateRepoIndex(pth, type = "source", Rversion = "3.6.0")

#addOldPackage("bisoreg", path = pth, vers = "1.5", type = "source",repos="https://pbil.univ-lyon1.fr/CRAN" , Rversion="3.6.0")



