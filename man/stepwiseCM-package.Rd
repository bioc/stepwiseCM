\name{stepwiseCM-package}
\alias{stepwiseCM-package}
\alias{stepwiseCM}
\docType{package}
\title{
Stepwise classification of cancer samples using high-dimensional data sets.
}

\description{
Given two types of data, this package is designed to evaluate the classification performances of two data types independently by user define classification algorithm (s), 
then explore the sample distributions in the two different data spaces. Based on the exact locations of the test samples in the data space  for which measreuments on all samples are available (presume that this data type is easy to obtain or relatively cheap) and the "pseudo" locations in the data space for which only partial measurements are available (presume that this data type is difficult to obtain or relatively expensive compared to former), the reclassification scores (RS) for each test sample is calculated without actually measuring the latter for large portion of samples. RS expresses our belief that a test sample is likely to be correctly classified if its covariates for the latter data types are measured. A large RS denotes, sample benefits more if classify it with latter data type and vice versa.
}

\details{
\tabular{ll}{
Package: \tab stepwiseCM\cr
Type: \tab Package\cr
Version: \tab 1.7.1\cr
Date: \tab 2013-05-20\cr
License: \tab GPL (>2)\cr
LazyLoad: \tab yes\cr
}
}

\author{
Askar Obulkasim

Maintainer: Askar Obulkasim <askar703@gmail.com>
}

\examples{
data(CNS)
train.cli <- t(CNS$cli[1:40,])
test.cli <- t(CNS$cli[41:60,])
train.gen <- CNS$mrna[,1:40]
test.gen <- CNS$mrna[,41:60]
train.label <- CNS$class[1:40]
test.label <- CNS$class[41:60]
pred.cli <- Classifier(train = train.cli, train.label = train.label, test = test.cli,
            type = "GLM_L1", CVtype = "k-fold", outerkfold = 2, innerkfold = 2)
pred.gen <- Classifier(train = train.gen, train.label = train.label, test = test.gen,
            type = "GLM_L2", CVtype = "k-fold", outerkfold = 2, innerkfold = 2)
prox1 <- Proximity(train.cli, train.label, test.cli, N = 2)$prox.test
prox2 <- Proximity(train.gen, train.label, NULL, N = 2)$prox.train
RS <- RS.generator(pred.cli$P.train, pred.gen$P.train, train.label, prox1, 
             prox2, type = "rank")
res <- Curve.generator(RS, pred.cli$P.test, pred.gen$P.test, test.label)
names(res)

}
