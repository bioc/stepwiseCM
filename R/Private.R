### internal function used to calculate the proximity matrix with multicores
.Wrapper.Proximity <- function(i, train, train.label, test)
{
   if (!is.null(test)) {
      RF <- randomForest(x = t(train), y = factor(train.label), xtest = t(test),
               ytest = NULL, ntree = 1000, importance = FALSE, proximity = TRUE)
   } else {
      RF <- randomForest(x = t(train), y = factor(train.label), xtest = test,
               ytest = NULL, ntree = 1000, importance = FALSE, proximity = TRUE) 
   }
   
  return(RF)

}



### internal function used to predict class labels for the training samples with multicore

.Wrapper.Classifier <- function(i, type, A, train, train.label, innerkfold)
{
    if(type == "TSP") {
          index <- A[[i]]
          train.tsp <- tspcalc(train[, -index], train.label[-index])
          P.tr <- as.numeric(predict(train.tsp, as.matrix(train[, index]))) 

    }  
    if(type == "GLM") {
          index <- A[[i]]
          cgh.glm.cv <- cv.glmpath(x = t(train[, -index]), y = train.label[-index], family = binomial, nfold = innerkfold,
                     mode = "lambda", type = "response", plot.it = FALSE, se = FALSE)
          min_cv <- cgh.glm.cv$fraction[which.min(cgh.glm.cv$cv.error)]
          cgh.glm <- glmpath(x = t(train[, -index]), y = train.label[-index], family = binomial, trace = FALSE)
          cgh.pred <- predict.glmpath(cgh.glm, newx = as.matrix(t(train[, index])), s = min_cv,
                     type = "response", mode = "lambda.fraction")
          P.tr <- as.vector(round(cgh.pred[1:length(index), ncol(cgh.pred)]))
    }
    if(type == "GLM_L1") {
          index <- A[[i]]
          tl <- factor(train.label[-index])
          datfr <- data.frame(t(train[, -index]))
          opt <- optL1(tl, penalized = datfr, fold = innerkfold, trace = FALSE)
          if (opt$lambda != Inf) {
             pen <- penalized(tl, penalized = datfr, lambda2 = opt$lambda, trace = FALSE)
          } else {
             pen <- penalized(tl, penalized = datfr, lambda2 = 1e-30, trace = FALSE)
          }
          whichsel <- which(pen@penalized !=0)
          datapred <- data.frame(t(train[, index]))
          names(datapred) <- names(datfr)
          P.tr <- as.vector(round(predict(pen, datapred)))

    }
    if(type == "GLM_L2") {
          index <- A[[i]]
          tl <- factor(train.label[-index])
          datfr <- data.frame(t(train[, -index]))
          opt <- optL2(tl, penalized = datfr, fold = innerkfold, trace = FALSE)
          if (opt$lambda != Inf) {
             pen <- penalized(tl, penalized = datfr, lambda2 = opt$lambda, trace = FALSE)
          } else {
             pen <- penalized(tl, penalized = datfr, lambda2 = 1e-30, trace = FALSE)
          }
          whichsel <- which(pen@penalized !=0)
          datapred <- data.frame(t(train[, index]))
          names(datapred) <- names(datfr)
          P.tr <- as.vector(round(predict(pen, datapred)))

    }
    if(type == "PAM") {
          if (!is.null(rownames(train))) {
             geneID <- rownames(train)
          } else {
             geneID <- as.character(1:nrow(train))
          }
          index <- A[[i]]
          subset.list <- list(x = train[, -index], y = factor(train.label[-index]), geneid = geneID, genenames = geneID)
          subset.train <- pamr.train(subset.list)
          for(d in 1:5) {
               subset.cv <- pamr.cv(subset.train, subset.list)
               errors <- subset.cv$error
               if(d == 1) {allerrors <- errors} else {allerrors <- allerrors+errors}
          }
          allerrors <- allerrors/5
          minim <- min(allerrors)
          marge <- 0.05
          withinmarge <- which(errors <= (marge+minim))
          threshold.chosen <- (subset.train$threshold)[withinmarge[length(withinmarge)]]
          subset.test <- as.matrix(train[, index])
          subset.pred <- pamr.predict(subset.train, subset.test, threshold = threshold.chosen)
          P.tr <- as.numeric(subset.pred) - 1

    }                 
    if(type == "SVM") {
          index <- A[[i]]
          x_sub <- t(train[, -index])
          y_sub <- as.factor(train.label[-index])
          wts <- length(y_sub)/table(y_sub)
          tuned <- tune.svm(x = x_sub, y = y_sub, gamma = c(0.5/nrow(train), 1/nrow(train), 
                  1.5/nrow(train), 2/nrow(train)), cost = c(1, 2, 4, 6, 8), tune.control(sampling = "cross", cross = 5))
          subset.train <- svm(x = x_sub, y = y_sub, kernel = "radial", gamma = tuned$best.parameters$gamma, 
                  cost = tuned$best.parameters$cost, class.weights = wts)
          subset_predicted <- predict(subset.train, t(train[, index]))
          P.tr <- as.numeric(subset_predicted) - 1
    }                     
    if(type == "plsrf_x") {
          x <- t(train)
          y <- train.label
          index <- A[[i]]
          if (length(index) != 1) {
             my.prediction <- plsrf_x(Xlearn = x[-index, ],Ylearn = y[-index], Xtest = x[index, ], ncomp=0:3,
                    ordered = NULL, nbgene = NULL)
             P.tr <- my.prediction$prediction
          } else {
             my.prediction <- plsrf_x(Xlearn = x[-index, ],Ylearn = y[-index], Xtest = x[c(index, index), ], ncomp=0:3,
                    ordered = NULL, nbgene = NULL)
             P.tr <- my.prediction$prediction[1]
          }
    }        
    if(type == "plsrf_x_pv") {
          x <- as.matrix(t(train))
          y <- train.label
          index <- A[[i]]
          if (length(index) != 1) {
             my.prediction <- plsrf_x_pv(Xlearn = x[-index, ],Ylearn = y[-index], Xtest = x[index, ], ncomp=0:3,
                    ordered = NULL, nbgene = NULL)
             P.tr <- my.prediction$prediction
          } else {
             my.prediction <- plsrf_x_pv(Xlearn = x[-index, ],Ylearn = y[-index], Xtest = x[c(index, index), ], ncomp=0:3,
                    ordered = NULL, nbgene = NULL)
             P.tr <- my.prediction$prediction[1]
          }
          
     }        
     if(type == "RF") {
          index <- A[[i]]
          rf <- randomForest(x = t(train[, -index]),y = as.factor(train.label[-index]), xtest = t(train[, index]),
                 ytest = NULL, ntree = 2000, importance = FALSE, proximity = FALSE)
          P.tr <- as.numeric(as.vector(rf$test$predicted))

     }          
   return(P.tr)
           
} 
