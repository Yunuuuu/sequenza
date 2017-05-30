
#setRepositories(graphics = FALSE, ind = 1:6)
#chooseCRANmirror(graphics = FALSE, ind = 2)
#lib <- "~/R/x86_64-unknown-linux-gnu-library/3.3"
#dir.create(path = lib, recursive = TRUE)
#.libPaths(c(lib, .libPaths()))
#install.packages(c("covr", "devtools"))
#library(devtools)
#install_bitbucket("sequenza_tools/sequenza@cleanup")

library(covr)
codecov(token="0ce2625b-f537-4638-818b-513085b2016f", branch="cleanup")
