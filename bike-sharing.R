#Bike sharing dataset
#https://archive.ics.uci.edu/ml/datasets/Bike+Sharing+Dataset

dat <- read.csv("hour.csv")
dframe <- dat[,c("dteday","yr","mnth","holiday","weekday","workingday","weathersit","temp","atemp","hum","windspeed","cnt")]
dframe$cnt <- residuals(lm(sqrt(cnt) ~  as.factor(holiday)+as.factor(weekday), data=dat))
str(dframe)

quantilevec <- 0.005*(1:199)

cvxquantiless <- function(gamma){
  folds <- 5
  
  mat <- matrix(0,ncol=1+length(quantilevec),nrow=folds)
  loss_all <- numeric()
  for (i in 1:folds){
   
    # select subset 
    no_days <- length(unique(dframe$dteday))
    selected_days <- unique(dframe$dteday)[ceiling((i-1)*no_days/folds+1):ceiling(i*no_days/folds)]
    samples <- rep(TRUE,length(dframe$yr))
    samples[dframe$dteday %in% selected_days] <- FALSE
    
    # create artificial data
    fit <- lm( cbind(temp,atemp,hum,windspeed,cnt) ~ dteday, data=dframe[samples,])
    fit_const <- lm( cbind(temp,atemp,hum,windspeed,cnt) ~ 1 , data=dframe[samples,])
    newdata <- fit_const$fitted.values + fit$residuals + gamma*(fit$fitted.values-fit_const$fitted.values) 
    newdata <- as.data.frame(newdata)
    
    # run anchor regression
    newdata_fit <- lm(cnt ~ . , data = newdata)
    
    # compute loss on held-out data
    loss <- (predict(newdata_fit,newdata=dframe[!samples,]) - dframe$cnt[!samples])^2
    #dframe[!samples,][names(sort(loss,decreasing = TRUE)[1:100]),]
    loss_daily <- sapply(unique(dframe$dteday[!samples]), function(x) mean(loss[dframe$dteday[!samples] == x]))
    loss_all <- c(loss_all, loss_daily)
    
    mat[i,] <- c(mean(loss_daily),quantile(loss_daily,probs=quantilevec))
  }
  return(loss_all)
}

gamma <- exp(seq(-2,3,0.05))
start_time <- Sys.time()
library(clustermq)
list_of_loss <- Q(fun = cvxquantiless,gamma=gamma,n_jobs = 200, export=list(quantilevec=quantilevec,dframe=dframe))
Sys.time() - start_time


# compute mean and quantiles across ranges of gamma
cvx <- t(sapply(list_of_loss ,function(x) c(mean(x), quantile(x, quantilevec)) ))
cvr <- do.call(rbind, list_of_loss)


# create plot
# y-axis: daily average squared residuals
# x-axis: gamma

#pdf(file="bike2b.pdf",width=7,height=6)
par(mar=c(5,6,1,1))
coll <- rgb(0.2,0.2,0.2,0.8)
matplot( gamma, (cvx[,-1]), type="l",col=coll,log="x",lty=1,xlab="",ylab="",cex.lab=1.6,lwd=0.55,axes=FALSE)
mtext("daily average squared residuals (quantiles)",2, line=3,cex=1.6)
mtext( expression(gamma),1, line=3, cex=1.6)
sel <- c(which.min(abs(quantilevec-0.5 )))+1
for (c in 1:length(sel)) lines( gamma, cvx[,sel[c]],col=rgb(0.8,0.2,0.2),lwd=1.5)
axis(1)
axis(2)
#lines(gamma,cvx[,1],col=rgb(0.8,0.2,0.2,0.8),lwd=3,lty=3)
#dev.off()
