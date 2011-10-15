\name{stepwiseCM-package}
\alias{stepwiseCM-package}
\alias{stepwiseCM}
\docType{package}
\title{
Stepwise classification of cancer samples using both clinical and molecular Data.
}

\description{
Given two types of data, it first evaluate the classification performances of two data types independently by the user define classification algorithm, 
then explore the two different data spaces using the proximity matrix from random forest algorithm. Based on the exact locations of the test samples in the clinical data space (presume clinical information of the test samples are given) and the "pseudo" locations in the molecular data space (presume molecular profiles of the test sample are not available), reclassification scores (RS) for each test sample is calculated. Large value of the RS means, sample benefits more if classify it with molecular data.
}

\details{
\tabular{ll}{
Package: \tab stepwiseCM\cr
Type: \tab Package\cr
Version: \tab 0.99.1\cr
Date: \tab 2011-05-16\cr
License: \tab GPL (http://www.gnu.org/copyleft/gpl.html)\cr
LazyLoad: \tab yes\cr
}
}

\author{
Askar Obulkasim

Maintainer: Askar Obulkasim <askar.wubulikasimu@vumc.nl>
}

\examples{
data(CNS)
train.cli <- t(CNS$cli[1:40, ])
test.cli <- t(CNS$cli[41:60, ])
train.gen <- CNS$mrna[, 1:40]
test.gen <- CNS$mrna[, 41:60]
train.label <- CNS$class[1:40]
test.label <- CNS$class[41:60]
result <- Curve.generator(train.cli, train.gen, train.label, test.cli, test.gen, test.label, 
               type = c("GLM_L1", "GLM_L2"), RStype = "both", Parallel = FALSE, 
               CVtype = "k-fold", outerkfold = 2, innerkfold = 2, N = 2, plot.it =FALSE)
names(result)

}
