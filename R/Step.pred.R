Step.pred <- function(RS, percent)
{ 
     if (percent < 0 | percent > 100)
         stop("Value of the argument ``ercent'' should be between 0 and 100")
     cut <- RS[order(RS, decreasing = TRUE)[round(length(RS) * percent)/100]]
     pass <- as.numeric(RS >= cut)
     res <- list(RS.cut = cut, ind = pass)   

 return(res)

}

