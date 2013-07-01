\name{Step.pred}
\alias{Step.pred}
\title{
A function to generate RS cutoff point based the given re-classification percentage.
}

\description{
Based on the specified percentage, this function finds the RS threshold and recommend which test samples may benefit by classifying with the data set at the second stage.   
}

\usage{
Step.pred(RS, percent)
}

\arguments{
  \item{RS}{
A vector of RS.
}
  \item{percent}{
Percentage of samples allow to pass to the second stage data set.
}
}

\value{A list object with following components:
  \item{RS.cut}{RS threshold corresponding to the specified re-classification percentage.}
  \item{ind}{a vector of binary values. 1 denotes sample is recommend to classify with the data set at the second stage and vice versa.}

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
            type = "GLM_L1", CVtype = "k-fold", outerkfold = 2, innerkfold = 2)
prox1 <- Proximity(train.cli, train.label, test.cli, N = 2)$prox.test
prox2 <- Proximity(train.gen, train.label, NULL, N = 2)$prox.train
RS <- RS.generator(pred.cli$P.train, pred.gen$P.train, train.label, prox1, 
             prox2, type = "rank")
res <- Step.pred(RS, 30)

}

