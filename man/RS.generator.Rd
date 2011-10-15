\name{RS.generator}
\alias{RS.generator}
\title{
A function to generate reclassification score.
}

\description{
A function to calculates the reclassification score (RS) using both clinical and molecular data.
}

\usage{
RS.generator(pred.cli, pred.gen, train.label, prox.gen, prox.cli, 
             type = c("rank", "proximity", "both"))
}

\arguments{
  \item{pred.cli}{
A vector containing the predicted labels of training set from clinical data. 
}
  \item{pred.gen}{
A vector containing the predicted labels of training set from genomic data.
}
  \item{train.label}{
A vector of the class labels (0 or 1) of the training set. NOTE: response values should be numerical not factor.
}
  \item{prox.gen}{
A square matrix, size of ncol(train.gen), contains the proximity values between training set in the genomic data space.
}

  \item{prox.cli}{
A matrix, size of ncol(test.cli) by ncol(train.cli), contains the proximity values between the test set and the training set in the clinical data space.
}

  \item{type}{
Which values are used to construct the reclassification score (RS)? There are three options available: \emph{proximity}, \emph{rank} and \emph{both}. If set to proximity, RS will be calculated directly from the proximity value. If set to rank, RS calculate from the rank of proximity values (more robust). If set to \emph{both}, both of them will be calculated. Default is \emph{rank}.
}
}

\details{
For each test sample, RS is calculated using the given classification results from clinical and genomic data. Algorithm project each test sample onto the clinical space check its neighborhood, also tries to gain some information about this test samples "pseudo" neighbors in the genomic space by the indirect mapping. If algorithm finds that the location of this test samples in the clinical space are more "safe" (more neighbors are correctly classified) and the location in the genomic space is surrounded by wrongly classified samples, then it will give this test sample high RS score and vice versa. After obtaining the RS, user can order them in descending order and pass the top ranked certain portion (decided by user) of samples to genomic data to classify.
}

\value{
If the "type" set to "rank" ("proximity"), then RS will be a vector of RS calculated from the ranking (proximity) approach , otherwise RS will be a matrix of RS, with two columns and size of rows equal number of test samples, calculated using the both approaches. 
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
pred.cli <- Classifier(train = train.cli, train.label = train.label, type = "GLM_L1", 
            CVtype = "k-fold", outerkfold = 2, innerkfold = 2)
pred.gen <- Classifier(train = train.gen, train.label = train.label, type = "GLM_L1", 
            CVtype = "k-fold", outerkfold = 2, innerkfold = 2)
prox <- Proximity(train.cli, train.label, test.cli, train.gen, N = 2)
RS.generator(pred.cli$P.train, pred.gen$P.train, train.label, prox$Prox.gen, 
             prox$Prox.cli, type = "both")
}            

