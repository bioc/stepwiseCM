Classifier <-
function(train, test = NULL, train.label, type = c("TSP", "GLM", "GLM_L1", 
  "GLM_L2", "PAM", "SVM", "plsrf_x", "plsrf_x_pv", "RF"), CVtype = c("loocv",
  "k-fold"), outerkfold = 5, innerkfold = 5)
{
    if (class(train) == "ExpressionSet" | class(train) == "eSet") {
       train <- exprs(train)
    }
    if (!is.null(test) & (class(test) == "ExpressionSet" | class(test) == "eSet")) {
       test <- exprs(test)
    }  
    if (length(table(train.label)) != 2)
       stop("Classifier was desined only for binary class problem")
    if (missing(type))
       stop("Classification algorithm type should be specified")
    if (any(is.na(train)) & type == "PAM")
     stop("PAM does not accept missing values!")
    if (!is.null(test)) {
       if (any(is.na(test)) & type == "PAM")
          stop("PAM does not accept missing values!")
    }
    if (CVtype == "loocv") {
       outerkfold = ncol(train)
    }
    if (!any(type %in% c("TSP", "GLM", "GLM_L1", "GLM_L2", "PAM", "SVM", "plsrf_x", "plsrf_x_pv", "RF")))
       stop("Selected classifier does not exist")
    A  <- split(sample(seq(ncol(train))), rep(1:outerkfold, 
               length = ncol(train)))
    P.tr <- c()
    P.te <- c()
    I <- c()
    call <- match.call()
    type <- match.arg(type)    
    Pred <- switch(type, TSP = 
            {
              if (length(A) != 1) {
                for(i in 1:length(A)) {
                   index <- A[[i]]
                   I <- c(I,index) 
                   train.tsp <- tspcalc(train[, -index], train.label[-index])
                   P.tr <- c(P.tr, as.numeric(predict(train.tsp, 
                             as.matrix(train[, index]))))
                }
                P.tr <- P.tr[order(I)]
              }
              if (!is.null(test)) {
                  train.tsp <- tspcalc(train, train.label)
                  P.te <- as.numeric(predict(train.tsp, test))
                  Pred <- list(P.train = P.tr, P.test = P.te)
              } else {
                  Pred <- list(P.train = P.tr)
              }

            }, GLM = 
            {
               if (length(A) != 1) {
                  for(i in 1:length(A)) {
                     index <- A[[i]]
                     I <- c(I, index)
                     cgh.glm.cv <- cv.glmpath(x = t(train[, -index]), y = 
                                   train.label[-index], family = binomial, 
                                   nfold = innerkfold, mode = "lambda", 
                                   type = "response", plot.it = FALSE, 
                                   se = FALSE)
                     min_cv <- 
                        cgh.glm.cv$fraction[which.min(cgh.glm.cv$cv.error)]
                     cgh.glm <- 
                        glmpath(x = t(train[, -index]), y = train.label[-index],
                        family = binomial, trace = FALSE)
                     cgh.pred <- 
                           predict.glmpath(cgh.glm, newx = as.matrix(t(train[, 
                           index])), s = min_cv, type = "response", 
                           mode = "lambda.fraction")
                     P.tr <- c(P.tr, as.vector(round(cgh.pred[1:length(index), 
                             ncol(cgh.pred)], 0)))
                  }
                  P.tr <- P.tr[order(I)]
               }
               if (!is.null(test)) {
                   cgh.glm.cv <- 
                         cv.glmpath(x = t(train), y = as.numeric(train.label),
                         family = binomial, nfold = innerkfold, mode = "lambda",
                         type = "response", plot.it = FALSE, se = FALSE)
                   min_cv <- cgh.glm.cv$fraction[which.min(cgh.glm.cv$cv.error)]
                   cgh.glm <- glmpath(x = t(train), y = train.label, family = 
                              binomial, trace = FALSE)
                   cgh.pred <- predict.glmpath(cgh.glm, newx = as.matrix(t(test)), 
                               s = min_cv, type = "response", mode = 
                               "lambda.fraction")
                   P.te <- as.vector(round(cgh.pred[1:ncol(test), 
                           ncol(cgh.pred)], 0))
                   Pred <- list(P.train = P.tr, P.test = P.te)
               } else {
                   Pred <- list(P.train = P.tr)
               }
            }, GLM_L1 =
            {
               if (length(A) != 1) {
                  for(i in 1:length(A)) { 
                     index <- A[[i]]
                     I <- c(I,index) 
                     tl <- factor(train.label[-index])
                     datfr <- data.frame(t(train[, -index]))
                     opt <- optL1(tl, penalized = datfr, fold = innerkfold,
                            trace = FALSE)
                     if (opt$lambda != Inf) {
                         pen <- penalized(tl, penalized = datfr, lambda1 = 
                                opt$lambda, trace = FALSE)
                     } else {
                         pen <- penalized(tl, penalized = datfr, lambda1 = 
                                1e-30, trace = FALSE)
                     }
                     whichsel <- which(pen@penalized !=0)
                     datapred <- data.frame(t(train[, index]))
                     names(datapred) <- names(datfr)
                     predresp <- round(predict(pen, datapred), 0)
                     P.tr <- c(P.tr, as.vector(predresp))
                  }   
                  P.tr <- P.tr[order(I)]  
               }
               if (!is.null(test)) {
                  tl <- factor(train.label)
                  datfr <- data.frame(t(train))
                  test <- data.frame(t(test))
                  opt <- optL1(tl, penalized = datfr,fold = innerkfold, 
                         trace = FALSE)
                  if (opt$lambda != Inf) {
                     pen <- penalized(tl, penalized = datfr, lambda2 = 
                            opt$lambda, trace = FALSE)
                  } else {
                     pen <- penalized(tl, penalized = datfr, lambda2 = 1e-30, 
                            trace = FALSE)
                  }
                  predresp <- round(predict(pen, test))
                  P.te <- as.vector(predresp)
                  whichsel <- which(pen@penalized !=0)
                  coeff <- pen@penalized[whichsel]
                  Pred <- list(P.train = P.tr, P.test = P.te)
               } else {
                  Pred <- list(P.train = P.tr)
               }
            }, GLM_L2 =
            {
               if (length(A) != 1) {
                  for(i in 1:length(A)) { 
                     index <- A[[i]]
                     I <- c(I,index) 
                     tl <- factor(train.label[-index])
                     datfr <- data.frame(t(train[, -index]))
                     opt <- optL2(tl, penalized = datfr, fold = innerkfold, 
                            trace = FALSE)
                     if (opt$lambda != Inf) {
                         pen <- penalized(tl, penalized = datfr, lambda2 = 
                                opt$lambda, trace = FALSE)
                     } else {
                         pen <- penalized(tl, penalized = datfr, lambda2 = 
                                1e-30, trace = FALSE)
                     }
                     whichsel <- which(pen@penalized !=0)
                     datapred <- data.frame(t(train[, index]))
                     names(datapred) <- names(datfr)
                     predresp <- round(predict(pen, datapred))
                     P.tr <- c(P.tr, as.vector(predresp))
                  }   
                  P.tr <- P.tr[order(I)]  
               }
               if (!is.null(test)) {
                  tl <- factor(train.label)
                  datfr <- data.frame(t(train))
                  test <- data.frame(t(test))
                  opt <- optL2(tl, penalized = datfr,fold = innerkfold, 
                         trace = FALSE)
                  if (opt$lambda != Inf) {
                     pen <- penalized(tl, penalized = datfr, lambda2 = 
                            opt$lambda, trace = FALSE)
                  } else {
                     pen <- penalized(tl, penalized = datfr, lambda2 = 
                            1e-30, trace = FALSE)
                  }
                  predresp <- round(predict(pen, test))
                  P.te <- as.vector(predresp)
                  whichsel <- which(pen@penalized !=0) 
                  coeff <- pen@penalized[whichsel]
                  Pred <- list(P.train = P.tr, P.test = P.te)
               } else {
                  Pred <- list(P.train = P.tr)
               }          
            }, PAM = 
            {
               if (!is.null(rownames(train))) {
                  geneID <- rownames(train)
               } else {
                  geneID <- as.character(1:nrow(train))
               }
               train <- as.matrix(train)
               if(length(A) != 1) {
                  for(i in 1:length(A)) {
                     index <- A[[i]]
                     I <- c(I, index)
                     subset.list <- list(x = train[, -index], y = 
                                    factor(train.label[-index]), 
                                    geneid = geneID, genenames = geneID)
                     subset.train <- pamr.train(subset.list)
                     if (CVtype == "loocv") {
                        subset.cv <- pamr.cv(subset.train, subset.list, 
                                     folds = as.list(1:(ncol(train)-1)))
                        allerrors <- subset.cv$error
                     } else {
                        for(d in 1:5) {
                           subset.cv <- pamr.cv(subset.train, subset.list, 
                                        nfold = innerkfold)
                           errors <- subset.cv$error
                           if (d == 1) {
                              allerrors <- errors
                           } else {
                              allerrors <- allerrors+errors
                           }
                        }
                        allerrors <- allerrors/5
                     }
                     minim <- min(allerrors)
                     marge <- 0.05
                     withinmarge <- which(allerrors <= (marge+minim))
                     threshold.chosen <- 
                      (subset.train$threshold)[withinmarge[length(withinmarge)]]
                     subset.test <- as.matrix(train[, index])
                     subset.pred <- pamr.predict(subset.train, subset.test, 
                                    threshold = threshold.chosen)
                                    whichsel = pamr.listgenes(subset.train, 
                                    subset.list, threshold = threshold.chosen)
                     P.tr <- c(P.tr, as.numeric(subset.pred)-1)
                  }
                  P.tr <- P.tr[order(I)]
               }
               if (!is.null(test)) {
                   subset.list <- list(x = train, y = factor(train.label), 
                                  geneid = geneID, genenames = geneID)
                   subset.train <- pamr.train(subset.list)
                   if (CVtype == "loocv") {
                      subset.cv <- pamr.cv(subset.train, subset.list, 
                                   folds = as.list(1:ncol(train)))
                      allerrors <- subset.cv$error
                   } else {
                      for(d in 1:5) {
                         subset.cv <- pamr.cv(subset.train, 
                                      subset.list,nfold = innerkfold)
                         errors <- subset.cv$error
                         if (d == 1) {
                            allerrors <- errors
                         } else {
                            allerrors <- allerrors+errors
                         }
                      }
                      allerrors <- allerrors/5
                   }
                   minim <- min(allerrors)
                   marge <- 0.05
                   withinmarge <- which(allerrors <= (marge+minim))
                   threshold.chosen <- 
                   (subset.train$threshold)[withinmarge[length(withinmarge)]]
                   subset.test <- as.matrix(as.matrix(test))
                   subset.pred <- pamr.predict(subset.train, subset.test, 
                                  threshold = threshold.chosen)
                                  whichsel = pamr.listgenes(subset.train, 
                                  subset.list, threshold = threshold.chosen)
                   P.te <- c(P.te, as.numeric(subset.pred)-1)
                   Pred <- list(P.train = P.tr, P.test = P.te)
               } else {
                   Pred <- list(P.train = P.tr)
               }
            }, SVM =
            {
               if(length(A) != 1) {
                  for(i in 1:length(A)) {
                     index <- A[[i]]
                     I <- c(I,index)
                     x_sub <- t(train[, -index])
                     y_sub <- as.factor(train.label[-index])
                     wts <- length(y_sub)/table(y_sub)
                     tuned <- tune.svm(x = x_sub, y = y_sub, gamma = 
                              c(0.5/nrow(train), 1/nrow(train), 
                              1.5/nrow(train), 2/nrow(train)), cost = 
                              c(1, 2, 4, 6, 8), tune.control(sampling = 
                              "cross", cross = innerkfold))
                     subset.train <- svm(x = x_sub, y = y_sub, kernel = 
                                     "radial", gamma = tuned$best.parameters$gamma, 
                                     cost = tuned$best.parameters$cost, 
                                     class.weights = wts)
                     subset_predicted <- predict(subset.train, t(train[, index]))
                     P.tr <- c(P.tr, as.numeric(subset_predicted)-1)
                  }
                  P.tr <- P.tr[order(I)]
               }
               if (!is.null(test)) {
                  x_sub <- t(train)
                  y_sub <- as.factor(train.label)
                  wts <- length(y_sub)/table(y_sub)
                  tuned <- tune.svm(x = x_sub, y = y_sub,gamma = 
                           c(0.5/nrow(train), 1/nrow(train), 1.5/nrow(train),
                           2/nrow(train)), cost = c(1, 2, 4, 6, 8), 
                           tune.control(sampling = "cross", cross = innerkfold))
                  subset.train <- svm(x = x_sub, y = y_sub, kernel = "radial", 
                                  gamma = tuned$best.parameters$gamma,
                                  cost = tuned$best.parameters$cost, 
                                  class.weights = wts)
                  subset_predicted <- predict(subset.train, t(test))
                  P.te <- as.numeric(subset_predicted)-1
                  Pred <- list(P.train = P.tr, P.test = P.te)
               } else {
                  Pred <- list(P.train = P.tr)
               }
            }, plsrf_x =
            {
               if (length(A) != 1) {
                  x <- t(train)
                  y <- train.label
                  for(i in 1:length(A)) {
                     index <- A[[i]]
                     I <- c(I, index)
                     if (CVtype != "loocv") {
                        my.prediction <- plsrf_x(Xlearn = x[-index, ], 
                                         Ylearn = y[-index], Xtest = x[index, ],
                                         ncomp = 0:3, ordered = NULL, 
                                         nbgene = NULL)
                        P.tr <- c(P.tr, my.prediction$prediction)
                     } else {
                        my.prediction <- plsrf_x(Xlearn = x[-index, ], 
                                         Ylearn = y[-index], Xtest = x[c(index, 
                                         index), ], ncomp = 0:3, ordered = NULL,
                                         nbgene = NULL)
                        P.tr <- c(P.tr, my.prediction$prediction[1])
                     } 
                  } 
                  P.tr <- P.tr[order(I)]
               }
               if (!is.null(test)) {
                  x <- t(train)
                  y <- train.label
                  my.prediction <- plsrf_x(Xlearn = x, Ylearn = y, Xtest = 
                                   t(test), ncomp = 0:3 , ordered = NULL, 
                                   nbgene = NULL)
                  P.te <- my.prediction$prediction
                  Pred <- list(P.train = P.tr, P.test = P.te)
               } else {
                  Pred <- list(P.train = P.tr)
               }
            }, plsrf_x_pv =
            {
               if (length(A) != 1) {
                  x <- t(train)
                  y <- train.label
                  for(i in 1:length(A)) {
                     index <- A[[i]]
                     I <- c(I, index)
                     if (CVtype != "loocv") {
                        my.prediction <- plsrf_x_pv(Xlearn = x[-index, ], 
                                         Ylearn = y[-index], Xtest = x[index, ],
                                         ncomp = 0:3, ordered = NULL, 
                                         nbgene = NULL)
                        P.tr <- c(P.tr, my.prediction$prediction)
                     } else {
                        my.prediction <- plsrf_x_pv(Xlearn = x[-index, ], Ylearn = 
                                         y[-index], Xtest = x[c(index, index), ],
                                         ncomp = 0:3, ordered = NULL, 
                                         nbgene = NULL)
                        P.tr <- c(P.tr, my.prediction$prediction[1])
                     } 
                  }
                  P.tr <- P.tr[order(I)]
               }
               if (!is.null(test)) {
                  x <- t(train)
                  y <- train.label
                  my.prediction <- plsrf_x_pv(Xlearn = x, Ylearn = y, Xtest = 
                                   t(test), ncomp = 0:3, ordered = NULL, 
                                   nbgene = NULL)
                  P.te <- my.prediction$prediction
                  Pred <- list(P.train = P.tr, P.test = P.te)
               } else {
                  Pred <- list(P.train = P.tr)
               }
            }, RF =
            {
               if(length(A) != 1) {
                  for(i in 1:length(A)) {
                     index <- A[[i]]
                     I <- c(I, index)
                     rf <- randomForest(x = t(train[, -index]),y = 
                           as.factor(train.label[-index]), xtest = 
                           t(train[, index]), ytest = NULL, ntree = 2000, 
                           importance = FALSE, proximity = FALSE)
                     P.tr <- c(P.tr, as.numeric(as.vector(rf$test$predicted)))
                  }
                  P.tr <- P.tr[order(I)]
               }
               if (!is.null(test)) {
                  rf <- randomForest(x = t(train),y = as.factor(train.label),
                        xtest = t(test), ytest = NULL, ntree = 2000,
                        importance = FALSE, proximity = FALSE)
                  P.te <- as.numeric(as.vector(rf$test$predicted))
                  Pred <- list(P.train = P.tr, P.test = P.te)
               } else {
                  Pred <- list(P.train = P.tr)
               }
            })

   return(Pred)
}

