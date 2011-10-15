Classifier.par <-
  function(train, test = c(), train.label, type = c("TSP", "GLM", "GLM_L1", 
  "GLM_L2", "PAM", "SVM", "plsrf_x", "plsrf_x_pv", "RF"), CVtype = c("loocv",
  "k-fold"), outerkfold = 5, 
   innerkfold = 5, featurenames = NULL, ncpus = 2)
{
    if(CVtype == "loocv"){
       outerkfold = ncol(train)
    }
    A  <- split(sample(seq(ncol(train))), rep(1:outerkfold, length = 
          ncol(train)))
    sfInit(parallel = TRUE, cpus = ncpus)
    sfLibrary(stepwiseCM)
    sfExportAll()
    P <- vector("list", outerkfold)
    P <- sfLapply(1:outerkfold, .Wrapper.Classifier, train = train, 
         train.label = train.label, type = type, A = A, innerkfold = innerkfold,
         featurenames = featurenames)
    P.tr <- rep(0, ncol(train))
    selfeatname_tr <- list()
    for(i in 1:outerkfold) {
       P.tr[A[[i]]] <- P[[i]]$P.train
       if (!is.null(featurenames)) {
          selfeatname_tr <- c(selfeatname_tr, list(unlist(P[[i]]$selfeatname)))
       }
    }
    if (!is.null(test)) {
       P <- Classifier(train = train, test = test, train.label = train.label, 
            type = type, CVtype = "k-fold", outerkfold = 1, innerkfold = 
            innerkfold)
       if (is.null(featurenames)) {
          Pred <- list(P.train = P.tr, P.test = P$P.test)
       } else {
          Pred <- list(P.train = P.tr, selfeatname_tr = selfeatname_tr, 
                  P.test = P$P.test, selfeatname_te = P$selfeatname_te)
       } 
    } else {
       if (is.null(featurenames)) {
          Pred <- list(P.train = P.tr)
       } else {
          list(P.train = P.tr, selfeatname_tr = selfeatname_tr)
       }
    } 

 return(Pred)

}

