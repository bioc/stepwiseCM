Step.pred <- function(curve, test.cli, percent)
{ 
     if (percent < 0 | percent > 100)
         stop("value of percent should be between 0 and 100")
     Prox <- curve$Proximity
     Prox.gen <- Prox$Prox.gen
     Prox.cli <- 0
     param <- curve$Param
     data <- curve$Matrices
     for(n in 1:param$N) {
        RF.cli <- randomForest(x = t(data$train.cli), 
                  y = factor(data$train.label), xtest = t(test.cli), 
                  ytest = NULL, ntree = 1000, importance = FALSE, 
                  proximity = TRUE)
        Prox.cli <- Prox.cli + RF.cli$test$proximity[, 
                  (nrow(RF.cli$test$proximity)+1):ncol(RF.cli$test$proximity)]
     }
     Prox.cli <- Prox.cli/param$N
     Pred.cli <- curve$Pred.cli
     Pred.gen <- curve$Pred.gen
     RS <- curve$RS
     rs <- RS.generator(Pred.cli$P.train, Pred.gen$P.train, data$train.label, 
           Prox.gen, Prox.cli, param$RStype)
     Pred <- Classifier(data$train.cli, test.cli, data$train.label, 
             param$type[1], param$CVtype, param$outerkfold, 
             param$innerkfold)
     if (param$RStype != "both") {
        cut <- RS[order(RS, decreasing = TRUE)[round(length(RS)*percent)/100]]
        pass <- as.numeric(rs >= cut)
     } else {
        cut <- c()
        pass <- matrix(NA, ncol(test.cli), 2) 
        colnames(pass) <- c("Rank based approach", "Proximity based approach")
        cut <- c(cut, RS[order(RS[, 1], 
               decreasing = TRUE)[round(nrow(RS)*percent)/100], 1])
        cut <- c(cut, RS[order(RS[, 2], 
               decreasing = TRUE)[round(nrow(RS)*percent)/100], 2])
        pass[, 1] <- as.numeric(rs[, 1] >= cut[1])
        pass[, 2] <- as.numeric(rs[, 2] >= cut[2])
     }
    result <- list(Pred = Pred$P.test, RS = rs, Threshold = cut, Pass = pass)

 return(result)

}

