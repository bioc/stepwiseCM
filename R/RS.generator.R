RS.generator <- function(pred1.train, pred2.train, train.label, prox1, 
                 prox2, type = c("rank" ,"proximity", "both"))
{
 type <- match.arg(type)
 if (missing(type))
     type <- "rank"
 if (length(pred1.train) != length(pred2.train))
    stop("Predicted labels should be in the same length in two data sets")
 if (ncol(prox1) != ncol(prox2) | ncol(prox2) != nrow(prox2))
    stop("Proximity matrix dimensions mismatch") 
 Y.C <- as.numeric(train.label == pred1.train)
 Y.G <- as.numeric(train.label == pred2.train)
 K <- max(min((min(sum(Y.C == 0), sum(Y.G == 0)) - 1), (min(sum(Y.C == 1), 
      sum(Y.G == 1)) - 1)), 1)
 if (type != "both") {
     RS <- rep(NA, nrow(prox1))
 } else {
     RS <- matrix(NA, nrow(prox1), 2)
     colnames(RS) <- c("RS (rank)", "RS (proximity)")
 }
 for(i in 1:nrow(prox1)) {
    Or.prox <- order(prox1[i, ], decreasing = TRUE)
    G.correct <- Or.prox[Y.C[Or.prox] == 1][1:K]                 
    G.incorrect<- Or.prox[Y.C[Or.prox] == 0][1:K]
    Rank.correct <-(which(Y.C[Or.prox] == 1)[1:K])/(1:K)
    Rank.incorrect <-(which(Y.C[Or.prox] == 0)[1:K])/(1:K)
    Prox.correct <- prox1[i, Or.prox][Rank.correct]                         
    Prox.incorrect <- prox1[i, Or.prox][Rank.incorrect]
    Rank.G.right <- c()
    Rank.G.wrong <- c()
    Prox.G.right <- c()
    Prox.G.wrong <- c()
    for(n in 1:K) {
       P <- order(prox2[G.correct[n], ], decreasing = TRUE)
       Or.P <- prox2[G.correct[n], P]
       Prox.R.NN <- Or.P[which(Y.G[P] == 1)][1:K]                                 
       Prox.W.NN <- Or.P[which(Y.G[P] == 0)][1:K] 
       Prox.G.right <- c(Prox.G.right, sum(Prox.R.NN)-sum(Prox.W.NN))
       Rank.R.NN <- (which(Y.G[P] == 1)[1:K])/(1:K)
       Rank.W.NN <- (which(Y.G[P] == 0)[1:K])/(1:K)
       Rank.G.right <- c(Rank.G.right, (sum(Rank.W.NN)-sum(Rank.R.NN)))
    } 
    for(n in 1:K) {
       P <- order(prox2[G.incorrect[n] ,], decreasing = TRUE)
       Or.P <-prox2[G.incorrect[n], P]
       Prox.R.NN <- Or.P[which(Y.G[P] == 1)][1:K]                               
       Prox.W.NN <- Or.P[which(Y.G[P] == 0)][1:K] 
       Prox.G.wrong <- c(Prox.G.wrong, sum(Prox.R.NN)-sum(Prox.W.NN))
       Rank.R.NN <- (which(Y.G[P] == 1)[1:K])/(1:K)
       Rank.W.NN <- (which(Y.G[P] == 0)[1:K])/(1:K)
       Rank.G.wrong <- c(Rank.G.wrong, (sum(Rank.W.NN)-sum(Rank.R.NN)))
    }
    Rank <- sum(Rank.correct * Rank.G.right) - 
                sum(Rank.incorrect * Rank.G.wrong)        
    Prox <- sum(Prox.incorrect * Prox.G.right) - 
                sum(Prox.correct * Prox.G.wrong)                         
    if(type == "rank") {
        RS[i] <- Rank
    }else if (type == "proximity") {                            
        RS[i] <- Prox
    } else {
        RS[i, ] <- c(Rank, Prox)
    }
    
 }
 return(RS)

}

