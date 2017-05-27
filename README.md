
# Sequenza: Copy Number Estimation from Tumor Genome Sequencing Data
![Sequenza_logo](https://bytebucket.org/sequenza_tools/icons/raw/da034ccc8c1ab5f5f8e020402267bd3f2dd5d361/png/sequenza/sequenza.png)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sequenza)](https://cran.r-project.org/package=sequenza)
[![CRAN_Downloads_Badge](http://cranlogs.r-pkg.org/badges/sequenza)](https://cran.r-project.org/package=sequenza)
[![CRAN_licence](https://img.shields.io/cran/l/sequenza.svg)](https://www.gnu.org/licenses/gpl-3.0.txt)


Sequenza is a tool to analyze genomic sequencing data from paired normal-tumor samples, including cellularity and ploidy estimation; mutation and copy number (allele-specific and total copy number) detection, quantification and visualization.

## Installation

Get the released version from CRAN:

```R
install.packages("sequenza")
```

Or the development version from github:

```R
# install.packages("devtools")
library(devtools)
install_bitbucket("ffavero/sequenza@cleanup")
```

## Testing

Unit tests are written using the package `testthat`.
You need to have the package installed in order to run the tests.

Clone the latest version:

```bash
git clone https://bitbucket.org/ffavero/sequenza.git
# Change to sequenza root directory
cd sequenza
# Run R
R
```

Load devtools and run the tests suite:
```R
library(devtools)
# Run testthat suite:
devtools::test()
```
