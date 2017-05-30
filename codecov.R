
setRepositories(graphics = FALSE, ind = 1:6)
chooseCRANmirror(graphics = FALSE, ind = 2)
install.pachages(c("covr", "devtools"))
library(covr)
libary(devtools)
install_bitbucket("sequenza_tools/sequenza#cleanup")
codecov(token="0ce2625b-f537-4638-818b-513085b2016f", branch="cleanup")
