Curve.generator <-
     function(train.cli, train.gen, train.label, test.cli, test.gen, test.label,
     type = c("TSP", "GLM", "GLM_L1", "GLM_L2", "PAM", "SVM", "plsrf_x", 
     "plsrf_x_pv", "RF"), RStype = c("rank", "proximity", "both"), 
     Parallel = FALSE, CVtype = c("loocv", "k-fold"), outerkfold = 5, 
     innerkfold = 5, ncpus = 2, N = 50, featurenames = NULL, plot.it = TRUE)
{ 
     if (is.null(test.cli))
         stop("Clinical data for the test set should not be empty")
     if ((is.null(test.gen) & plot.it) | (is.null(test.label) & plot.it))
         stop("Genomic data and class labels for the test samples are required
              for the plot")
     if (missing(RStype))
         RStype <- "rank"
     if (length(type) == 1 & any(type %in% c("TSP", "PAM",
                  "plsrf_x", "plsrf_x_pv"))) 
          stop("TSP, PAM, plsrf_x and plsrf_x_pv algorithms only work 
          for molecular data")
     if (length(type) == 1) { 
         type <- c(type, type)
     }
     if (!Parallel) {
         Pred.cli <- Classifier(train.cli, test.cli, train.label, type[1], 
                     CVtype = CVtype, outerkfold, innerkfold, featurenames)
         Pred.gen <- Classifier(train.gen, test.gen, train.label, type[2], 
                     CVtype = CVtype, outerkfold, innerkfold, featurenames)
     } else {
         Pred.cli <- Classifier.par(train.cli, test.cli, train.label, type[1],
                  CVtype = CVtype, outerkfold, innerkfold, featurenames, ncpus)
         Pred.gen <- Classifier.par(train.gen, test.gen, train.label, type[2], 
                  CVtype = CVtype, outerkfold, innerkfold, featurenames, ncpus)
     }
     Pred.train.cli <- Pred.cli$P.train
     Pred.test.cli <- Pred.cli$P.test
     Pred.train.gen <- Pred.gen$P.train
     if (!is.null(featurenames)) {
        selfeatname_tr <- Pred.gen$selfeatname_tr
        selfeatname_te <- Pred.gen$selfeatname_te
     }   
     Prox <- Proximity(train.cli, train.label, test.cli, train.gen, N, Parallel,
             ncpus)
     Prox.cli <- Prox$Prox.cli
     Prox.gen <- Prox$Prox.gen
     RS <- RS.generator(Pred.train.cli, Pred.train.gen, train.label, Prox.gen, 
                 Prox.cli, RStype)
     if (!is.null(test.gen)) {
         Y.C <- as.numeric(Pred.test.cli == test.label)
         Y.G <- as.numeric(Pred.gen$P.test == test.label)
         if (RStype == "rank" | RStype == "proximity") {
            Acc <- rep(NA, 11)
            for(i in seq(0.1, 0.9, 0.1)) {
                Index <- order(RS, decreasing = TRUE)[1:round(length(RS)*i)]
                Acc[i*10+1] <- (sum(Y.G[Index],Y.C[-Index])/ncol(test.gen)) * 100
            }
            Acc[1] <- (sum(Y.C)/ncol(test.gen)) * 100
            Acc[11] <- (sum(Y.G)/ncol(test.gen)) * 100
            if (plot.it){
                plot(seq(0,100,10), Acc, pch=1, type="b", main = c(RStype, 
                "based approach"), xlab="Percentage of samples passed to 
                genomic data (%)", ylab = "Accuracy (%)", ylim = c(0,100), 
                col = "black")
                points(seq(0,100,10), rep(Acc[1], 11), pch = 3, col = "red")
                points(seq(0,100,10), rep(Acc[11], 11), lty = 1, pch = 16, 
                col = "blue")
                legend("topright", cex=0.75, pch=c(1,3,16), col=c("black", "red", 
                "blue"), legend=c("stepwise", "clinical", "molecular"), ncol=2)
            }
         } else {
            Acc <- matrix(NA, 11, 2)
            rownames(Acc)=paste(seq(0, 100, 10), "%", sep="")
            colnames(Acc) <- c("Rank based approach (%)", 
                             "Proximity based approach (%)")
            for (j in 1:2) {
               for (i in seq(0.1, 0.9 ,0.1)) {
                  Index <- order(RS[, j], decreasing = TRUE)[1:round(nrow(RS)*i)]
                  Acc[i*10+1, j] <- (sum(Y.G[Index],
                                    Y.C[-Index])/ncol(test.gen)) * 100
               }
            }
            Acc[1, ] <- (sum(Y.C)/ncol(test.gen)) * 100
            Acc[11, ] <- (sum(Y.G)/ncol(test.gen)) * 100
            if (plot.it){
               par(mfrow=c(1, 2))
               plot(seq(0, 100, 10), Acc[, 1], pch=1, type = 'b', main = 
               "Rank based approach", xlab="Percentage of samples 
               passed to genomic data (%)", ylab = "Accuracy (%)",
               ylim = c(0,100), col = "black")
               points(seq(0,100,10), rep(Acc[1, 1], 11), pch = 3, col = "red")
               points(seq(0,100,10), rep(Acc[11, 1], 11), lty = 1, pch = 16, 
               col = "blue")
               legend("topright", cex=0.75, pch=c(1,3,16), col=c("black", "red", 
               "blue"), legend=c("stepwise", "clinical", "molecular"), ncol=2)
               plot(seq(0, 100, 10), Acc[, 2], pch=1, type = 'b', main =
               "Proximity based approach", xlab="Percentage of samples 
               passed to genomic data (%)", ylab = "Accuracy (%)", 
               ylim = c(0, 100), col = "black")
               points(seq(0,100,10), rep(Acc[1, 2], 11), pch = 3, col = "red")
               points(seq(0,100,10), rep(Acc[11, 2], 11), lty = 1, pch = 16, 
               col = "blue")
               legend("topright", cex=0.75, pch=c(1,3,16), col=c("black", "red", 
               "blue"), legend=c("stepwise", "clinical", "molecular"), ncol=2)
        
            }
         }
         Param <- list(type = type, RStype = RStype, CVtype = CVtype, outerkfold
                   = outerkfold, innerkfold = innerkfold, N = N)
         matrix <- list(train.cli = train.cli, train.gen = train.gen, 
                   train.label = train.label, test.cli = train.cli, test.gen
                   = train.gen, test.label = train.label)
         if (!is.null(featurenames)) {
            Result <- list(Pred.cli = Pred.cli, Pred.gen = Pred.gen, Proximity
                      = Prox, RS = RS, Accuracy = Acc, selfeatname_tr = 
                      selfeatname_tr, selfeatname_te = selfeatname_te, 
                      Param = Param, Matrices = matrix)
         } else {
            Result <- list(Pred.cli = Pred.cli, Pred.gen = Pred.gen, Proximity
                      = Prox, RS = RS, Accuracy = Acc, Param = Param, 
                      Matrices = matrix)
         }  
     } else {
         if (!is.null(featurenames)) {
            Result <- list(Pred.cli = Pred.cli, Pred.gen = Pred.gen, Proximity
                   = Prox, RS = RS, selfeatname_tr = selfeatname_tr, 
                   selfeatname_te = selfeatname_te, Param = Param, 
                   Matrices = matrix)
         } else {
            Result <- list(Pred.cli = Pred.cli, Pred.gen = Pred.gen, Proximity
                          = Prox, RS = RS, Param = Param, Matrices = matrix)
         }
     }
 
   return(Result)

}

