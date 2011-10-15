\name{Proximity}
\alias{Proximity}
\title{
A function to calculate proximity matrix.
}

\description{
A function to compute the proximity matrices using clinical and genomic data separately.
}

\usage{
Proximity(train.cli, train.label, test.cli, train.gen, N = 50, 
         Parallel = FALSE, ncpus = 2)
}

\arguments{
  \item{train.cli}{
A data frame or matrix of containing predictors for the training set, where columns correspond to samples and rows to features.
}
  \item{train.label}{
A vector of the class labels (0 or 1) of the training set. NOTE: response values should be numeric not factor.
}
  \item{test.cli}{
A data frame or matrix of containing predictors for the test set, where columns correspond to samples and rows to features.
}
  \item{train.gen}{
A data frame or matrix of containing predictors of the training set from genomic data, where columns correspond to samples and rows to features.
}
  \item{N}{
Number of repetition for calculating the proximity matrix, final proximity matrix is average of these repeats. We recommend setting this number high, so that more stable proximity matrix will be produced. Default is 50.
}
  \item{Parallel}{
Should proximity calculation use the parallel processing procedure? Default is FALSE.
}
  \item{ncpus}{
Number of CPUs assign to the parallel computation. Default is 2.
}
}

\details{
Proximity matrix is calculated using the random forest algorithm. Proximity values ranges from 0 (least similar) to 1 (perfect match). 
}

\value{A list object with the following components:
  \item{Prox.cli}{A matrix, size of ncol(test.cli) by ncol(train.cli), contains the proximity values between the test set and the training set in the clinical data space .}
  \item{Prox.gen}{A square matrix, size of ncol(train.gen), contains the proximity values between training set in the genomic data space. }
}

\references{
Breiman, L. (2001), \emph{Random Forest, 45}, 5-32.
}

\author{
Askar Obulkasim
    
Maintainer: Askar Obulkasim <askar.wubulikasimu@vumc.nl>

}

\examples{
data(CNS)
train.cli <- t(CNS$cli[1:40,])
test.cli <- t(CNS$cli[41:60,])
train.gen <- CNS$mrna[,1:40]
train.label <- CNS$class[1:40]
##without parallel processing procedure
Prox <- Proximity(train.cli, train.label, test.cli, train.gen, N = 2)
##with parallel processing procedure
\dontrun{Prox <- Proximity(train.cli, train.label, test.cli, train.gen, 
                N = 50, Parallel = TRUE, ncpus = 10)}
}                

