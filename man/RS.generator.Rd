\name{RS.generator}
\alias{RS.generator}
\title{
A function to generate the reclassification score.
}

\description{
This function calculates the reclassification score (RS) for the test set. See details.
}

\usage{
RS.generator(pred1.train, pred2.train, train.label, prox1, prox2, 
             type = c("rank", "proximity", "both"))
}

\arguments{
  \item{pred1.train}{
A numeric vector contains the predicted class labels of the training set from the data set used at the first stage. 
}
  \item{pred2.train}{
A numeric vector contains the predicted class labels of the training set from the data set used at the second stage.
}
  \item{train.label}{
A vector of actual class labels (0 or 1) of the training set. Should be numerical not factor.
}
  \item{prox1}{
A rectangular matrix contains the proximity values between the test set (rows) and the training set (columns) obtained from the data set used at the first stage.
}
  \item{prox2}{
A square matrix contains the proximity values between training set obtained from the data set used at the second stage.
}
  \item{type}{
Which values are used to construct the reclassification score (RS)? There are three options available: \emph{proximity}, \emph{rank} and \emph{both}. If set to proximity, RS will be calculated directly from the proximity value. If set to rank, RS calculate from the rank of proximity values (more robust). If set to \emph{both}, both of them will be calculated. Default is \emph{rank}.
}
}

\details{
For each test sample, RS is calculated using the given classification results from two data sets. Algorithm project each test sample onto the first stage data space to observe its neighbourhood and tries to gain some information about the test sample's "pseudo" neighbourhood in the second stage data space with the help of indirect mapping. If algorithm finds that the location of this test samples in the first stage data space are more "safe" (more neighbours are correctly classified) and the location in the second stage data space is surrounded by wrongly classified samples, then it gives this test sample lower RS score and vice versa. After obtaining the RS, user may order them in descending order and classify the top ranked certain portion of samples with the second stage data type.
}

\value{
If the argument "type" set to "rank" ("proximity"), function returns a vector of RS calculated by rank (proximity) based approach. If set to "both" returns a matrix of two columns corresponds to the RS obtained by rank and proximity based approaches. 
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
train.label <- CNS$class[1:40]
pred.cli <- Classifier(train = train.cli, train.label = train.label, type = "GLM_L1", 
            CVtype = "k-fold", outerkfold = 2, innerkfold = 2)
pred.gen <- Classifier(train = train.gen, train.label = train.label, type = "GLM_L1", 
            CVtype = "k-fold", outerkfold = 2, innerkfold = 2)
prox1 <- Proximity(train.cli, train.label, test.cli, N = 2)$prox.test
prox2 <- Proximity(train.gen, train.label, NULL, N = 2)$prox.train
RS <- RS.generator(pred.cli$P.train, pred.gen$P.train, train.label, prox1, 
             prox2, type = "both")
}            

