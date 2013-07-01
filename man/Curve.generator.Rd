\name{Curve.generator}
\alias{Curve.generator}
\title{
A function to generate accuracy curve by passing different portion of samples to the data set used at the second stage.
}

\description{
Accuracy curve may be used as a reference for choosing the re-classification score (RS) threshold threshold for incoming samples.
}

\usage{
Curve.generator(RS, pred1.test, pred2.test, 
    test.label, plot.it = TRUE)
}

\arguments{
  \item{RS}{
A numeric vector contains the re-classification score of the test set.
}
  \item{pred1.test}{
A numeric vector of contains the predicted class labels of the test set from the data set used at the first stage. Should be numeric not factor. 
}
  \item{pred2.test}{
A numeric vector of contains the predicted class labels of the test set from the data set used at the second stage. Should be numeric not factor. 
}
  \item{test.label}{
A vector of actual class labels (0 or 1) of the test set. Should be numeric not factor. 
}
  \item{plot.it}{
If set to ``TRUE'', this function produces a plot in which Y axis denotes the accuracy and X denotes the percentage of samples passed to the second stage. In order to make this plot, class labels and molecular data for the test set must be given. Default is ``TRUE''. 
}
}

%\details{
%This function requires the molecular information for a group of samples (called test here) which are not used in the training phase. Based on their RS, a accuracy curve %will be obtained by passing different proportion of samples to the molecular data. 
%}

\value{A data frame of two columns. The first column denotes the percentage of samples passed to the data used at the second stage, and the second denotes the corresponding accuracy.}


\author{
Askar Obulkasim
    
Maintainer: Askar Obulkasim <askar703@gmail.com>
}



\seealso{
\code{\link{Classifier}}, \code{\link{Classifier.par}}, \code{\link{Proximity}}, \code{\link{RS.generator}} 
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
            type = "GLM_L1", CVtype = "k-fold", outerkfold = 2, innerkfold = 2)
prox1 <- Proximity(train.cli, train.label, test.cli, N = 2)$prox.test
prox2 <- Proximity(train.gen, train.label, NULL, N = 2)$prox.train
RS <- RS.generator(pred.cli$P.train, pred.gen$P.train, train.label, prox1, 
             prox2, type = "rank")
res <- Curve.generator(RS, pred.cli$P.test, pred.gen$P.test, test.label)
             
}          
