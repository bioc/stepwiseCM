Proximity <- function(train.cli, train.label, test.cli, train.gen,   
             N = 50 , Parallel = FALSE ,ncpus = 2)
{
    if (ncol(train.cli) != length(train.label) | 
       ncol(train.gen) != length(train.label))
        stop("Number of training samples and class labels should be the same")
    test.cli.prox <- 0
    train.gen.prox <- 0
    if (!Parallel) {
        for(n in 1:N) {
            RF.cli <- randomForest(x = t(train.cli), y = factor(train.label), 
                      xtest = t(test.cli), ytest = NULL, ntree = 1000, 
                      importance = FALSE, proximity = TRUE)
            RF.gen <- randomForest(x = t(train.gen), y = factor(train.label),
                      xtest = NULL, ytest = NULL, ntree = 1000, importance =
                      FALSE, proximity = TRUE)
            test.cli.prox <- test.cli.prox + RF.cli$test$proximity[, 
                   (nrow(RF.cli$test$proximity)+1):ncol(RF.cli$test$proximity)]
            train.gen.prox <- train.gen.prox + RF.gen$proximity
        }
    } else {
         sfInit(parallel = TRUE, cpus = ncpus)
         sfLibrary(randomForest)
         sfLibrary(Stepwise) 
         sfExportAll()   
         Prox <- vector("list", N)
         Prox <- sfLapply(1:N, .Wrapper.Proximity, train.cli = train.cli,
                 test.cli = test.cli, train.gen = train.gen, train.label 
                 = train.label)
         for(t in 1:N) {
            P <- matrix(Prox[[t]], ncol = ncol(train.cli))  
            test.cli.prox <- test.cli.prox + P[1:ncol(test.cli), ]
            train.gen.prox <- train.gen.prox + P[(ncol(test.cli) + 1):nrow(P), ]
         } 
     }

     return(list(Prox.cli = test.cli.prox/N, Prox.gen = train.gen.prox/N))

}

