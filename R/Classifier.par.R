Classifier.par <-
  function(train, test = NULL, train.label, type = c("TSP", "GLM", "GLM_L1", 
  "GLM_L2", "PAM", "SVM", "plsrf_x", "plsrf_x_pv", "RF"), CVtype = c("loocv",
  "k-fold"), outerkfold = 5, 
   innerkfold = 5, ncpus = 2)
{
    if (class(train) == "ExpressionSet" | class(train) == "eSet") {
       train <- exprs(train)
    }
    if (!is.null(test) & (class(test) == "ExpressionSet" | class(test) == "eSet")) {
       test <- exprs(test)
    } 
    if(CVtype == "loocv") {
       outerkfold = ncol(train)
    }
    A  <- split(sample(seq(ncol(train))), rep(1:outerkfold, length = 
          ncol(train)))
    sfInit(parallel = TRUE, cpus = ncpus)
    sfLibrary(stepwiseCM)
    sfExportAll()
    P <- vector("list", outerkfold)
    P <- sfLapply(1:outerkfold, .Wrapper.Classifier, type, A, train, train.label, innerkfold)
    sfStop()
    P.tr <- rep(0, ncol(train))
    for(i in 1:outerkfold) {
       P.tr[A[[i]]] <- P[[i]]
    }
    if (!is.null(test)) {
       P <- Classifier(train = train, test = test, train.label = train.label, 
            type = type, CVtype = "k-fold", outerkfold = 1, innerkfold = 
            innerkfold)
       Pred <- list(P.train = P.tr, P.test = P$P.test)
    } else {
       Pred <- list(P.train = P.tr)
    } 

 return(Pred)

}

