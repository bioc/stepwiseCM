\name{Classifier}
\alias{Classifier}
\title{
A function to perform classification task.
}

\description{
Given the set of samples for training and the type of classification algorithm, this function constructs the classifier using the training set, and predicts the class labels of the test set using the trained classifier.
}

\usage{
Classifier(train, test = NULL, train.label, type = c("TSP", "GLM", "GLM_L1", 
           "GLM_L2", "PAM", "SVM", "plsrf_x", "plsrf_x_pv", "RF"), 
           CVtype = c("loocv", "k-fold"), outerkfold = 5, innerkfold = 5)
}

\arguments{
  \item{train}{
An object of class \code{\link{ExpressionSet}} or data frame or matrix contains predictors for the training set, where columns correspond to samples and rows to features.
}
  \item{test}{
An object of class \code{\link{ExpressionSet}} or data frame or matrix contains predictors for the test set (optional), where columns correspond to samples and rows to features.
}
  \item{train.label}{
A numeric vector contains the actual class labels (0 or 1) of the training set. NOTE: class labels should be numerical not factor. 
}
  \item{type}{
Type of classification algorithms used. Currently 9 well-known algorithm are available for user the choose from. They are: top scoring pair (TSP), logistic regression (GLM), GLM with L1 (lasso) penalty, GLM with L2 (ridge) penalty, prediction analysis for microarray (PAM), support vector machine (SVM), Random Forest combined with partial least square dimension reduction (plsrf_x), Random Forest combined with partial least square dimension reduction plus pre-validation (plsrf_x_pv), Random Forest (RF). NOTE: "TSP", "PAM", "plsrf_x" and "plsrf_x_pv"  are exclusively designed for high-dimensional data.
}
  \item{CVtype}{
Cross-validation type to obtain predicted labels of the training set. Must be either k-fold cross-validation (k-fold), or leave-one-out-cross-validation (loocv).
}
  \item{outerkfold}{
Number of cross-validation used in the training phase.
}
  \item{innerkfold}{
Number of cross validation used to estimate the model parameters. E.g. penalty parameter in "GLM_L1".
}
}

\value{A list object contains following components:
  \item{P.train}{predicted class labels of the training set.}
  \item{P.test}{predicted class labels of the test set if the test set is given.}

}

\references{
Aik Choo Tan and Daniel Q. Naiman and Lei Xu and Raimond L. Winslow and Donald Geman(2005). Simple Decision Rules for Classifying Human Cancers from Gene Expression Profiles(TSP). \emph{Bioinformatics, 21}, 3896-3904.

Anne-Laure Boulesteix and Christine Porzelius and Martin Daumer(2008). Microarray-based Classification and Clinical Predictors: on Combined Classifiers and Additional Predictive Value. \emph{Bioinformatics, 24}, 1698--1706.
}

\author{
Askar Obulkasim
    
Maintainer: Askar Obulkasim <askar703@gmail.com>
}


\seealso{
\code{\link{Classifier.par}}
}

\examples{
data(CNS)
train <- CNS$mrna[, 1:40]
test <- CNS$mrna[, 41:60]
train.label <- CNS$class[1:40]
Pred <- Classifier(train = train, test = test, train.label = train.label, 
        type = "GLM_L1", CVtype = "k-fold", outerkfold = 2, innerkfold = 2)
Pred$P.train
Pred$P.test
}
