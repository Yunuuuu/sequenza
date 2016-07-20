# Sequenza: Copy Number Estimation from Tumor Genome Sequencing Data

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sequenza)](https://cran.r-project.org/package=sequenza)
[![CRAN_Downloads_Badge](http://cranlogs.r-pkg.org/badges/sequenza)](https://cran.r-project.org/package=sequenza)


Sequenza is a tool to analyze genomic sequencing data from paired normal-tumor samples, including cellularity and ploidy estimation; mutation and copy number (allele-specific and total copy number) detection, quantification and visualization.

## Installation

Install Bioconductor dependencies:

```R
source("https://bioconductor.org/biocLite.R")
biocLite("copynumber")
```

Get the released version from CRAN:

```R
install.packages("sequenza")
```

Or the development version from github:

```R
# install.packages("devtools")
library(devtools)
install_bitbucket("ffavero/sequenza")
```