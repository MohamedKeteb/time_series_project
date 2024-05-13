require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries)


xm <- rev(zoo(serie_lt$index))
plot(xm)


#x  <- data.frame(data = ts(serie$index))
#plot(x)

acf(xm) # on affiche les auto corrélation 

# il ya une tendance linéaire croissante on peut donc différencier une fois 

xm_diff <- xm-lag(xm,-1)
plot(xm_diff)
acf(xm_diff)
pacf(xm_diff)

#  l'auto corrélation d'ordre 1 est loin de 1 dans ACF et PACF semble stationnaire 

plot(xm_diff)

# existe il des racine unité 

adf_result <- adf.test(xm_diff)
print(adf_result)

# on rejette H0 = racine unité 


pp_result <- pp.test(xm_diff)
print(adf_result)

# déterminer le ARIMA(p, q, d = 0) d = 0 car série tationnaire





y <- (xm_diff - mean(xm_diff))/ (var(xm_diff)*1/2) #centrer la série

plot(y)
acf(y,20)
pacf(y,20)

# p = 6, q = 1


arima(y,c(6,0,1))
arima302 <- arima(y,c(6,0,1)) #enregistre les r´esultats de l’estimation
Box.test(arima302$residuals, lag=8, type="Ljung-Box", fitdf=7) 



Qtests <- function(series, k, fitdf=0) {
pvals <- apply(matrix(1:k), 1, FUN=function(l) {
pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
return(c("lag"=l,"pval"=pval))
})
return(t(pvals))
}
Qtests(arima302$residuals, 20, 7)


# vérification de tous les models 


modelchoice <- function(p,q,data= y, k=20){
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}























