Curve.generator <- function(RS, pred1.test, pred2.test, test.label, plot.it = TRUE)
{ 
         if ((length(pred1.test) != length(pred2.test)) | (length(pred1.test) != length(test.label)))
            stop("Imput dimensions mismatch")

         Y.C <- as.numeric(pred1.test == test.label)
         Y.G <- as.numeric(pred2.test == test.label)
         acc <- rep(NA, 11)
         for(i in seq(0.1, 0.9, 0.1)) {
            index <- order(RS, decreasing = TRUE)[1:round(length(RS) * i)]
            acc[i*10+1] <- (sum(Y.G[index],Y.C[-index]) / length(test.label)) * 100
         }
         acc[1] <- (sum(Y.C)/length(RS)) * 100
         acc[11] <- (sum(Y.G)/length(RS)) * 100
         if (plot.it){
            plot(seq(0,100,10), acc, pch = 1, type = "b", main = "RS reference curve",
                xlab="Percentage of samples passed to molecular data (%)", 
                ylab = "Accuracy (%)", ylim = c(0,100), col = "black")
                points(seq(0,100,10), rep(acc[1], 11), pch = 3, col = "red")
                points(seq(0,100,10), rep(acc[11], 11), lty = 1, pch = 16, 
                col = "blue")
                legend("topright", cex=0.75, pch=c(1,3,16), col=c("black", "red", 
                "blue"), legend=c("stepwise", "clinical", "molecular"), ncol = 1)
         }
         acc <- data.frame(mol.percentage = seq(0,100,10), accurary = acc)  
 
   return(acc)

}

