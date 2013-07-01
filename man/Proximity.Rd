\name{Proximity}
\alias{Proximity}
\title{
A function to calculate the proximity matrix.
}

\description{
This function computes the proximity matrix by Random Forest algorithm. Proximity values ranges from 0 (least similar) to 1 (perfect match).
}

\usage{
Proximity(train, train.label, test = NULL, N = 50, 
         Parallel = FALSE, ncpus = 2)
}

\arguments{
  \item{train}{
An object of class \code{\link{ExpressionSet}} or data frame or matrix contains the predictors for the training set, where columns correspond to samples and rows to features.
}
  \item{train.label}{
A vector of actual class labels (0 or 1) of the training set. Should be numeric not factor.
}
  \item{test}{
An object of class \code{\link{ExpressionSet}} or data frame or matrix of containing predictors for the test set, where columns correspond to samples and rows to features.
}
  \item{N}{
Number of repetition for calculating the proximity matrix, final proximity matrix is average of these repeats. We recommend to set a large number, so that stable proximity matrix will be produced. Default is 50.
}
  \item{Parallel}{
Should proximity calculation use the parallel processing procedure? Default is FALSE.
}
  \item{ncpus}{
Number of acores assign to the parallel computation. Default is 2.
}
}

%\details{
%Proximity matrix is calculated using the random forest algorithm.  
%}

\value{A list object with following components:
  \item{prox.train}{A square symmetric matrix contains the proximity values of the training set .}
  \item{prox.test}{A rectangular square matrix contains the proximity values between test set (rows) and training set (columns). Only returned when test set is supplied.}
}

\references{
Breiman, L. (2001), \emph{Random Forest, 45}, 5-32.
}

\author{
Askar Obulkasim
    
Maintainer: Askar Obulkasim <askar703@gmail.com>

}

\examples{
data(CNS)
train <- t(CNS$cli[1:40,])
test <- t(CNS$cli[41:60,])
train.label <- CNS$class[1:40]
##without parallel processing procedure
Prox <- Proximity(train, train.label, test, N = 2)
##with parallel processing procedure
\dontrun{Prox <- Proximity(train, train.label, test,  
                N = 50, Parallel = TRUE, ncpus = 10)}
}                

