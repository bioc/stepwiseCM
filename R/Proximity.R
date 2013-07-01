Proximity <- function(train, train.label, test = NULL,    
             N = 50, Parallel = FALSE, ncpus = 2)
{
    if (class(train) == "ExpressionSet" | class(train) == "eSet") {
       train <- exprs(train)
    }
    if (!is.null(test) & (class(test) == "ExpressionSet" | class(test) == "eSet")) {
       test <- exprs(test)
    } 
    if (ncol(train) != length(train.label))
        stop("Number of training samples and class labels should be the same")
    test.prox <- 0
    train.prox <- 0
    if (!Parallel) {
        for(n in 1:N) {
           if (!is.null(test)) {
               RF <- randomForest(x = t(train), y = factor(train.label), 
                      xtest = t(test), ytest = NULL, ntree = 1000, 
                      importance = FALSE, proximity = TRUE)
               train.prox <- train.prox + RF$proximity
               test.prox <- test.prox + RF$test$proximity[, 
                   (ncol(test)+1):ncol(RF$test$proximity)] 
            } else {
               RF <- randomForest(x = t(train), y = factor(train.label), 
                      xtest = test, ytest = NULL, ntree = 1000, 
                      importance = FALSE, proximity = TRUE)
               train.prox <- train.prox + RF$proximity
            }
        }
    } else {
         sfInit(parallel = TRUE, cpus = ncpus)
         sfLibrary(randomForest)
         sfLibrary(stepwiseCM) 
         sfExportAll()   
         Prox <- vector("list", N)
         Prox <- sfLapply(1:N, .Wrapper.Proximity, train, train.label, test)
         sfStop()
         for(n in 1:N) {
            train.prox <- train.prox + Prox[[n]]$proximity
            if (!is.null(test)) {
               test.prox <- test.prox + Prox[[n]]$test$proximity[, 
                           (ncol(test)+1):ncol(Prox[[n]]$test$proximity)] 
            } 
         } 
     }
     if (!is.null(test)) {
        return(list(prox.train = train.prox / N, prox.test = test.prox / N))
     } else {
        return(list(prox.train = train.prox / N))
     }  

}

