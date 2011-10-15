\name{Classifier.par}
\alias{Classifier.par}
\title{
A function to obtain the predicted class labels of the test set using the parallel processing procedure.
}

\description{
Classification has been embedded inside the parallel processing procedure to speed up the computation for large dataset. 
}
\usage{
Classifier.par(train, test = c(), train.label, type = c("TSP", "GLM", "GLM_L1", 
               "GLM_L2", "PAM", "SVM", "plsrf_x", "plsrf_x_pv", "RF"), 
               CVtype = c("loocv", "k-fold"), outerkfold = 5, innerkfold = 5, 
               featurenames = NULL, ncpus = 2)
}

\arguments{
  \item{train}{
A data frame or matrix of containing predictors for the training set, where columns correspond to samples and rows to features.
}
  \item{test}{
A data frame or matrix containing predictors for the test set (optional), where columns correspond to samples and rows to features.
}
  \item{train.label}{
A vector of the class labels (0 or 1) of the training set. NOTE: response values should be numerical not factor.
}
  \item{type}{
Type of classification algorithms. Currently 9 different types of algorithm are available. They are: top scoring pair (TSP), logistic regression (GLM), GLM with L1 (lasso) penalty, GLM with L2 (ridge) penalty, prediction analysis for microarray (PAM), support vector machine (SVM), random forest method after partial least square dimension reduction (plsrf_x), random forest method after partial least square dimension reduction plus prevalidation (plsrf_x_pv), random forest (RF). NOTE: TSP, PAM, plsrf_x and plsrf_x_pv algorithms does not work with clinical data.
}
  \item{CVtype}{
Cross validation type.
}
  \item{outerkfold}{
Number of cross validation used in the training phase.
}
  \item{innerkfold}{
Number of cross validation used to estimate the model parameters.
}
  \item{featurenames}{
Feature names in molecular data (e.g. gene or probe names). If given, function also produces name of the selected feature during the training and test phases. Feature selection only works with "TSP", "GLM_L1" and "GLM_L2" algorithms. "RF" provide feature importance.
}
  \item{ncpus}{
Number of CPUs assign to the parallel computation.
}
}

\value{
A list object \emph{Pred} which contains following components:
  \item{P.train}{A vector of the predicted class labels of the training set.}
  \item{P.test}{A vector of the predicted class labels of the test set if the test set is given.}
  \item{selfeatname_tr}{A list object, size of \emph{outerkfold}, containing the name of the selected features during the training phase if the \emph{featurenames} is given.}
  \item{selfeatname_te}{A list object containing the name of the selected features during the test phase if the \emph{test} and \emph{featurenames} are given.} 

}

\references{
Aik Choo Tan and Daniel Q. Naiman and Lei Xu and Raimond L. Winslow and Donald Geman(2005). Simple Decision Rules for Classifying Human Cancers from Gene Expression Profiles(TSP). \emph{Bioinformatics, 21}, 3896-3904.

Anne-Laure Boulesteix and Christine Porzelius and Martin Daumer(2008). Microarray-based Classification and Clinical Predictors: on Combined Classifiers and Additional Predictive Value. \emph{Bioinformatics, 24}, 1698--1706.
}

\author{
Askar Obulkasim
    
Maintainer: Askar Obulkasim <askar.wubulikasimu@vumc.nl>
}

\examples{
data(CNS)
train <- CNS$mrna[, 1:40]
test <- CNS$mrna[, 41:60]
train.label <- CNS$class[1:40]
\dontrun{Pred <- Classifier.par(train = train, test = test, train.label = train.label, type = "GLM_L1", 
                        CVtype = "k-fold", outerkfold = 5, innerkfold = 5, ncpus = 5)}
}   
