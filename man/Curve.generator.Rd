\name{Curve.generator}
\alias{Curve.generator}
\title{
A function to generate accuracy with different cut points .
}

\description{
For given clinical and molecular data, this function first calculates the predicted class labels of the training sets and the proximity matrices using clinical and molecular data separately. In the second step based on the classification performances of the two data types and the sample distribution in the two spaces, it calculates the reclassification score for the test set. Then, produce a vector accuracies by passing different percentage of samples to molecular data. This curve can be used as a reference for choosing the RS threshold for incoming test samples.
}

\usage{
Curve.generator(train.cli, train.gen, train.label, test.cli, test.gen, 
    test.label, type = c("TSP", "GLM", "GLM_L1", "GLM_L2", "PAM", "SVM", 
    "plsrf_x", "plsrf_x_pv", "RF"), RStype = c("rank", "proximity", "both"), 
    Parallel = FALSE, CVtype = c("loocv", "k-fold"), outerkfold = 5, 
    innerkfold = 5, ncpus = 2, N = 50, featurenames = NULL, plot.it = TRUE)
}

\arguments{
  \item{train.cli}{
A data frame or matrix of containing predictors of the training set from clinical data, where columns correspond to samples and rows to features.
}
  \item{train.gen}{
A data frame or matrix of containing predictors of the training set from molecular data, where columns correspond to samples and rows to features.
}
  \item{train.label}{
A vector of the class labels (0 or 1) of the training set. NOTE: response values should be numeric not factor.
}
  \item{test.cli}{
A data frame or matrix of containing predictors of the test set from clinical data, where columns correspond to samples and rows to features.
}
  \item{test.gen}{
A data frame or matrix containing predictors of the test set from genomic data, where columns correspond to samples and rows to features.
}
  \item{test.label}{
A vector of the class labels (0 or 1) of the test set (optional). NOTE: response values should be numeric not factor. 
}
  \item{type}{
Type of classification algorithms. Currently 9 different types of algorithm are available. They are: top scoring pair (TSP), logistic regression (GLM), GLM with L1 (lasso) penalty, GLM with L2 (ridge) penalty, prediction analysis for microarray (PAM), support vector machine (SVM), random forest method after partial least square dimension reduction (plsrf_x), random forest method after partial least square dimension reduction plus prevalidation (plsrf_x_pv), random forest (RF). NOTE: TSP, PAM, plsrf_x and plsrf_x_pv algorithms does not work with clinical data.
}
  \item{RStype}{
Which values are used to calculate the reclassification score (RS)? There are three options available: \emph{proximity}, \emph{rank} and \emph{both}. If set to \emph{proximity}, RS will be calculated directly from the proximity value. If set to \emph{rank}, RS calculate based on the rank of the proximity values (more robust). If set to \emph{both}, both of them will be calculated. Default is rank.
}
  \item{Parallel}{
Should class prediction and proximity calculation use the parallel processing procedure? Default is FALSE.
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
  \item{ncpus}{
Number of CPUs assign to the parallel computation.
}
  \item{N}{
Number of repetition for calculating the proximity matrix, final proximity matrix is average of these repeats. We recommend setting this number high, so that more stable proximity matrix will be produced. Default is 50.
}
  \item{featurenames}{
Feature names in molecular data (e.g. gene or probe names). If given, function also produces name of the selected feature during the training and test phases. Feature selection only works with "TSP", "GLM_L1" and "GLM_L2" algorithms. "RF" provides feature importance.
}
  \item{plot.it}{
If set to TRUE, this function produces a plot in which Y axis denotes the accuracy and X denotes the percentage of samples passed to the second stage. In order to make this plot, class labels and molecular data for the test set must be given. Default is TRUE. 
}
}

\details{
This function requires the molecular information for a group of samples (called test here) which are not used in the training phase. Based on their RS a accuracy curve will be attained by setting different thresholds. If the option \emph{type}, has length 2 (two classification algorithm types are given), then first one used for the prediction using clinical data and second one used for the prediction using molecular data. If only one algorithm type is given, the same algorithm used for both data types. Note that, TSP, PAM, plsrf_x and plsrf_x_pv algorithms does not work with clinical data.
}

\value{A list object with the following components:
  \item{Predicted.cli}{A list object includes the predicted class labels of the training set and the test set if test is given. This is classification result with clinical data.}
  \item{Predicted.gen}{A list object includes the predicted class labels of the training set and the test set if test is given. This is classification result with molecular data.}
  \item{Proximity}{A list object of proximity matrices which includes the proximity matrix of test samples in the clinical data space (ncol(test.cli) by ncol(train.cli) matrix) and the proximity matrix of the training samples in the genomic data space (ncol(train.gen) by ncol(train.gen) matrix).}
  \item{RS}{If the "type" set to "rank" ("proximity"), it gives a vector of re-classification scores calculated from the ranking (proximity) approach , otherwise it gives a matrix with two columns and size of rows equal number of test samples, calculated using the both approaches.}
  \item{Accuracy}{If \emph{test.gen} is given, accuracies corresponding to different percentage of samples are classified with molecular data are produced. If \emph{RStype} is set to rank or proximity, accuracy will be a vector. If RStype is set to both, accuracy will be a matrix with two columns and eleven rows. First column corresponding to accuracies when RS is calculated using the rank of proximity, second column corresponding to accuracies when RS is calculated using the proximity.}
  \item{Param}{A list object contains the values of parameters specified by user.}
  \item{Matrices}{A list object contains the input data matrices.}
}

\author{
Askar Obulkasim
    
Maintainer: Askar Obulkasim <askar.wubulikasimu@vumc.nl>
}



\seealso{
\code{\link{Classifier}}, \code{\link{Classifier.par}}, \code{\link{Proximity}}, \code{\link{RS.generator}} 
}

\examples{
data(CNS)
tr.cli <- t(CNS$cli[1:40, ])
te.cli <- t(CNS$cli[41:60, ])
tr.gen <- CNS$mrna[, 1:40]
te.gen <- CNS$mrna[, 41:60]
tr.label <- CNS$class[1:40]
te.label <- CNS$class[41:60]
result <- Curve.generator(train.cli=tr.cli, train.gen=tr.gen, train.label=tr.label, test.cli=te.cli, 
                         test.gen=te.gen, test.label=te.label, type = c("GLM_L1", "GLM_L2"), 
                         RStype = "rank", Parallel = FALSE, CVtype = "k-fold", outerkfold = 2, 
                         innerkfold = 2, N = 2)
names(result)
}          
