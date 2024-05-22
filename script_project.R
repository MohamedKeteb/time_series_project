

# Charger les packages nécessaires
require(MASS)
require(ellipse)
require(zoo)
require(tseries)


#----------------------------------------------------------------------------

xm <- rev(zoo(serie_lt$index)) # création de l'objet Zoo, la série temporelle
T = length(xm) # le nombre d'observations de la série 


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



acf(xm, main = "") # on affiche les autocorrélations


xm_diff <- xm-lag(xm,-1) # Différenciation de la série

plot(xm_diff, type = "l", col="blue", xlab="time", ylab="Index (Base 100)") # affichage de la série 


# ACF et PACF de la série différenciée
acf(xm_diff)
pacf(xm_diff)



# Test de racines unitaires et test de stationnarité

adf_result <- adf.test(xm_diff) # test ADF
print(adf_result)

pp_result <- pp.test(xm_diff) # test de Perron Phillips racine unitaire
print(pp_result)


kpss_result <- kpss.test(xm_diff) # test KPSS
print(kpss_result)




# Le ARIMA(p,d,q)



y <- (xm_diff - mean(xm_diff)) # centrer la série différenciée

plot(y) # afficher la série centrée

# PACF et ACF pour identifier l'ARMA

par(mfrow = c(1, 2))
acf(y,20, main = "")
pacf(y,20, main = "")



# vérification de tous les models 


p = 6 ; q = 1



# fonction pour estimer un ARIMA et en vérifier la significativité et la décorrélation des résidus


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
selec # les modèles sélectionnés



# Calculer l'AIC et BIC pour les modèles sélectionnés


pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) #cr´ee
#une liste des ordres p et q des mod`eles candidats
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")") #renomme les ´el´ements de la liste
models <- lapply(pqs, function(pq) arima(y,c(pq[["p"]],0,pq[["q"]]))) #cr´ee une liste des mod`eles
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m))) #calcule les AIC et




arima101 <- arima(xm,c(1,1,1)) # le modèle retenu 
print(arima101)


Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
Qtests(arima101$residuals, 20, 0) # les tests portmanteau du modèle sélectionné 



# Prévision


ris <- arima101$residuals # stocker les résidus du modèle 


mu = mean(xm_diff) # moyenne de la série différenciée

coefficients <- coef(arima101) # coéfficients du modèle 
coeff <- coefficients * c(1, -1) # pour avoir ma1 et pas l'opposé

phi = 0.2628741 # ar1
psi = 0.7032179 # ma1


variance_residus <- var(ris) # variance des résidus



# X = (1 - phi) * mu + (1 + phi) * xm[T] + (- phi * xm[T-1] - psi * ris[T])
# Y = (1- phi)*mu + (1+phi) * 100.0347 - phi * xm[T]




forecast_values <- predict(arima101, n.ahead=2) # caluculer les prédictions T+1 et T+2


X <- forecast_values$pred[1]
Y <- forecast_values$pred[2]

# Matrice de variance covariance

sigma <- variance_residus * matrix(c(1, (1+phi) - psi , (1+phi) - psi, 1 + ((1+phi) - psi)**2), nrow = 2, byrow = TRUE)

A_inv <- solve(sigma) # inverser la matrice sigma


# Tracer la région de confiance
plot(ellipse(A_inv, centre = c( X , Y),level = 0.95, draw = TRUE), type = 'l', xlab = expression(X[T+1]), ylab = expression(X[T+2]), main = 'Région de confiance 95%')




points(X, Y, col = 'red', pch = 19) # tracer le centre 
text(X, Y, labels = expression(hat(X)), pos = 3) 
grid()



# Discutons l'hypothèse de résidus gaussiens 

# Standardiser les résidus 

ris = (ris - mean(ris))/(var(ris)*(1/2))



# QQ- Plot 


qqnorm(ris, main = "")
qqline(ris, col = "red")
grid()





# Effectuer le test de Shapiro-Wilk
shapiro_test <- shapiro.test(ris)

# Afficher les résultats du test
print(shapiro_test)


# Test de Kolmogorov Smirnov 

ks_test <- ks.test(ris, "pnorm", mean(ris), sd(ris))

# Afficher les résultats du test
print(ks_test)









