require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries)


xm <- rev(zoo(serie_lt$index)) # création de l'objet Zoo
length(xm) # le nombre d'observation de la série 


dates <- paste0(serie_lt$time, "-01")
dates <- as.Date(rev(dates), format="%Y-%m-%d")



plot(dates, xm, type = "l", col="blue", main="Indice de la production industrielle", xlab="Dates", ylab="Index (Base 100)", xaxt="n")

# Sélectionner seulement quelques dates à afficher
vecteur <- round(seq(from=1, to=337, length.out=10))

dates_to_display <- dates[vecteur]
labels_to_display <- format(dates_to_display, "%Y-%m")
axis(1, at=dates_to_display, labels=FALSE)
text(x = dates_to_display, y = par("usr")[3] - 3, labels = labels_to_display, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
grid()



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





y <- (xm_diff - mean(xm_diff)) #centrer la série

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


p = 6 ; q = 1



## fonction pour estimer un arima et en v´erifier l’ajustement et la validit´e



signif <- function(estim){ #fonction de test des significations individuelles des coefficients
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

## fonction pour estimer un arima et en verifier l'ajustement et la validite
modelchoice <- function(p,q,data=y, k=20){
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,20,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}

## fonction pour estimer et verifier tous les arima(p,q) avec p<=pmax et q<=max
armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}

armamodels <- armamodelchoice(pmax = 6,qmax = 1)



selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),]
selec

##



pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) #cr´ee
#une liste des ordres p et q des mod`eles candidats
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")") #renomme les ´el´ements de la liste
models <- lapply(pqs, function(pq) arima(y,c(pq[["p"]],0,pq[["q"]]))) #cr´ee une liste des mod`eles
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m))) #calcule les AIC et
#BIC des mod`eles candidats


T <- length(xm)
































