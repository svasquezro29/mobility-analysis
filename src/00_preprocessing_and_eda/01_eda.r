# Cargar libreria forecast, fBasic, fArma, urca,
library(forecast)
library(fBasics)    
library(car)
library(fArma)   
library(urca)   
            

# Entrada de datos desde un archivo de texto
(z=ts(scan("f:/udea/Economet2/ejemplo_1.txt")))

# ETAPA DE IDENTIFICACI?N
# gr?fica de la serie
plot.ts(z)
plot.ts(z, type="o")

# An?lisis de estabilidad de la varianza: Es necesario usar la
# transformaci?n Box-Cox?

# An?lisis gr?fico
par(mfrow=c(2,1))
plot.ts(z)
plot.ts(diff(z))

# Estimaci?n de la transformaci?n usando la librer?a car 
(tBoxCox=powerTransform(z))
summary(tBoxCox)

# transformaci?n usando ra?z cuarta
tz=z^.25

# gr?fica de la serie original y transformada con Lambda=.25
par(mfrow=c(2,1))
plot.ts(z, type="o")
plot.ts(z^.25, type="o")
 
# gr?fica de la serie diferenciada original y transformada
par(mfrow=c(2,1))
plot.ts(diff(z))
plot.ts(diff(z^.25))

# se sigue el proceso de identificaci?n usando z^.25
# determinaci?n del valor de d
# gr?fica de la serie transformada
plot.ts(z^.25)
# correlogramas muestrales de la serie transformada 
par(mfrow=c(2,1))
Acf(z^.25, lag.max=30, ci=0,ylim=c(-1,1))
pacf(z^.25, lag.max=30, ci=0, ylim=c(-1,1))

# La serie es no estacionaria: Qu? tipo de serie es TS o DS?
# prueba de ra?ces unitarias usando la librer?a urca
(maxlag=floor(12*(length(z)/100)^(1/4)))
ru_tz=ur.df(z^.25, type = c("trend"), lags=maxlag, selectlags = c("BIC"))
summary(ru_tz)
ru_tz=ur.df(z^.25, type = c("trend"), lags=maxlag, selectlags = c("AIC"))
summary(ru_tz)

  # validaci?n de la ecuaci?n ADF
    resid=ru_tz@testreg$residuals             # residuales del modelo ADF
    plot(ru_tz)
    auto.arima(resid, max.p=5, max.q=5)       # busqueda "autom?tica" 
    cheq=Arima(resid, c(0,0,0), include.constant=TRUE) 
    summary(cheq)    
    tsdiag(cheq, gof.lag=15)
  # Verificacion de normalidad en los residuales
    qqnorm(resid,  xlab = "Cuantiles Te?ricos", ylab = "Cuantiles Muestrales")
    qqline(resid)
    shapiro.test(resid)                # prueba de Shapiro-Wilks
    jarqueberaTest(resid )             # prueba Jarque-Bera en libreria fBasics

# el modelo ADF satisface los supuestos b?sicos
#Conclusi?n: La serie z^.25 contiene al menos una ra?z unitaria.

# Examen de m?s ra?ces unitarias
# gr?fica de la serie diferenciada una vez
plot.ts(diff(z^.25), type="l")
# prueba de si hay ra?z unitarias en z^.25 diferenciada una vez.
ru_dif_tz=ur.df(diff(z^.25), type = c("drift"), lags=maxlag, selectlags = c("BIC"))
summary(ru_dif_tz)  
# reespecificaci?n del modelo con lags=0
ru_dif_tz=ur.df(diff(z^.25), type = c("drift"), lags=0)
summary(ru_dif_tz)   # valide la ecuaci?n usada para la prueba
#Conclusi?n: La serie z^.25  diferenciada no contiene r?iz unitaria.
# Por tanto, la serie z^.25 pertenece a la clase de modelos ARIMA con deriva y d=1.

# determinaci?n de los valores de (p, q) del modelos ARMA para (1-B)z^.25
# correlogramas muestrales para z diferenciada una vez
par(mfrow=c(2,1))
Acf(diff(z^.25), lag.max=12, ci=0, ylim=c(-1,1))
pacf(diff(z^.25), lag.max=12, ci=0, ylim=c(-1,1))

# correlogramas muestrales para z^.25 diferenciada una vez con bandas
par(mfrow=c(2,1))
Acf(diff(z^.25), lag.max=12, ci=0)
pacf(diff(z^.25), lag.max=12)

# selecci?n "autom?tica" del modelo
auto.arima(z^.25, d=1, max.p=5, max.q=5, ic=c("aic"))
auto.arima(z^.25, max.p=5, max.q=5, ic=c("bic"))

# Conclusi?n: modelo seleccionado es: (1-phi1*B)(1-B)z^.25=constante+at, 
# es decir, z^.25 sigue un modelo ARIMA(1,1,0) con deriva.
# la tendencia observada en la gr?fica de z^.25 es una mezcla de tendencia aleatoria y detemin?stica lineal 

# ETAPA DE ESTIMACI?N
# estimaci?n ML exacta con valores iniciales dados por la estimaci?n condicional
mod1_CSS_ML=Arima(z, c(1, 1, 0), include.drift=TRUE, lambda=.25, method = c("CSS-ML"))
summary(mod1_CSS_ML) 
(res1_CSS_ML=residuals(mod1_CSS_ML))

# ETAPA DE DIAGN?STICOS
# ra?ces de los polinomios
# usando la librer?a fArma
coeficAR=c(coef(mod1_CSS_ML)[1])
# vector con el coeficiente AR
armaRoots(coeficAR)                  # obtiene las ra?ces 
# otra forma: usando librer?a forecast usando el inverso de las ra?ces del polinomio
plot(mod1_CSS_ML)

# An?lisis de los residuales
tsdiag(mod1_CSS_ML)

# chequeo de observaciones at?picas extremas (no es un an?lisis completo de outliers)
res1_est=res1_CSS_ML/(mod1_CSS_ML$sigma2^.5)  # estandarizaci?n de los residuales
plot.ts(res1_est, type="o")
abline(a=-3, b=0)
abline(a=3, b=0)

# n?mero de observaciones esperadas fuera de los limites (-3, 3) bajo normalidad
(Nobs_Esp=round(length(z)*2*pnorm(-3, mean = 0, sd = 1, lower.tail = TRUE)))

# detecci?n de las observaciones at?picas
ind=(abs(res1_est)>3.0)
sum(ind)
(grupo=cbind(res1_est, ind))

# chequeo de normalidad
# gr?fico cuantil-cuantil
qqnorm(res1_est,  xlab = "Cuantiles Te?ricos", ylab = "Cuantiles Muestrales",
       xlim=c(-4,4), ylim=c(-4,4))
qqline(res1_est)
 
# histograma, densidad kernel y gr?fico normal
plot(density(res1_est))
mu<-mean(res1_est)
sigm<-sd(res1_est)
x<-seq(-4,4,length=100)
y<-dnorm(x,mu,sigm)
hist(res1_est,prob=T,ylim=c(0,.45),xlim=c(-4,4),col="yellow")
lines(density(res1_est))
lines(x,y,lwd=2,col="blue")
# conclusi?n: No se detecta alejamiento fuerte de la normalidad

# pruebas de normalidad
shapiro.test(res1_est)               # prueba de Shapiro-Wilks
normalTest(res1_est, method=("jb"))  # En la librer?a fBasics: puede realizar otras pruebas
                                     # "ks" for the Kolmogorov-Smirnov one?sample test, 
                                     # "sw" for the Shapiro-Wilk test, 
                                     # "jb" for the Jarque-Bera Test, 
                                     # "da" for the D'Agostino Test. The default value is "ks"
# conclusi?n: No se rechaza la normalidad.

# valores ajustados del modelo para z transformada
mod1_CSS_MLn=Arima(z^.25, c(1, 1, 0), include.drift=TRUE, method = c("CSS-ML"))
summary(mod1_CSS_ML)n) 
(res1_CSS_MLn=residuals(mod1_CSS_MLn))

(ajustn=z^.25-residuals(mod1_CSS_MLn))     
# gr?fico para los valores ajustados y los valores observados
ts.plot(z^.25,ajustn)   # gr?fico de las series contra el tiempo
lines(z^.25, col="black")
lines(ajustn, col="red")

plot(as.vector(z^.25),as.vector(ajustn), type="p")  # gr?fico de dispersi?n de la serie observada 
abline(0,1)                                      # contra la serie ajustada

# detecci?n de observaciones at?picas
ind=(abs(res1_est)>3.0)
sum(ind)
(grupo=cbind(res1_est, ind))

# n?mero de observaciones esperadas fuera del intervalo (-3, 3)
(num_obs=length(z)*2*pnorm(-3.0, mean = 0, sd = 1, lower.tail = TRUE))

# Evaluaci?n de Pron?sticos: se contruye el modelo usando 245 datos y se dejan
# los ?ltimos 5 para evaluar la capacidad de pron?stico del modelo 
# estimaci?n ML exacta con valores iniciales dados por la estimaci?n condicional
num_pron=5
(mod_Evalpron=Arima(z[1:(length(z)-num_pron)], c(1, 1, 0), include.drift=TRUE, lambda=.25, method = c("CSS-ML")))
summary(mod_Evalpron) 

(z_pred<-forecast(mod_Evalpron, h=num_pron, lambda=mod1_CSS_ML$lambda, fan=TRUE))
plot(forecast(z_pred))

(z_pred<-forecast(mod_Evalpron, h=num_pron, level=c(95), lambda=mod1_CSS_ML$lambda))
plot(forecast(z_pred))

# EVALUACI?N DE LOS PRON?STICOS
(real=ts(z[(length(z)-num_pron+1):length(z)]))
cbind(real, ts(z_pred$mean))

# gr?fico de reales, pron?sticos
ts.plot(real, ts(z_pred$mean),type="o")
lines(ts(z_pred$mean), col="red")

# gr?fico de reales, pron?sticos e intervalo de predicci?n
ts.plot(real, ts(z_pred$mean), ts(z_pred$lower), ts(z_pred$upper), type="o")
lines(ts(z_pred$mean), col="red")
lines(ts(z_pred$lower), col="blue")
lines(ts(z_pred$upper), col="blue")

# C?lculo de medidas de evaluaci?n de los pron?sticos
(recm=(mean((real-ts(z_pred$mean))^2))^.5)            #  rmse
(recmp=(mean(((real-ts(z_pred$mean))/real)^2))^.5)    #  rmspe

(eam=mean(abs(real-ts(z_pred$mean))))                 #  mae
(eamp=mean(abs((real-ts(z_pred$mean))/real)))         #  mape

# descomposici?n del error cuadr?tico medio
(ecm=(mean((real-ts(z_pred$mean))^2))) # c?lculo del error cuadr?tico medio

# proporci?n de sesgo
(prop_sesgo_med=(mean(ts(z_pred$mean))-mean(real))^2/ecm) 

# c?lculo de varianza sesgadas
sigmap=((length(real)-1)/length(real)*var(ts(z_pred$mean)))^.5   # c?lculo de desv. estand. de los pron?sticos
sigmar=((length(real)-1)/length(real)*var(ts(real)))^.5          # c?lculo de desv. estand. de los datos

# proporci?n de varianza
(prop_sesgo_var=(sigmap-sigmar)^2/ecm)         

# proporci?n de covarianza
(prop_covar=2*(1-cor(ts(z_pred$mean), real))*sigmap*sigmar/ecm)
prop_sesgo_med+prop_sesgo_var+prop_covar                         # verificaci?n de c?lculos

# Verdaderos Pron?sticos
(z_pron<-forecast(mod1_CSS_ML, h=10, level=c(90), lambda=mod1_CSS_ML$lambda, fan=FALSE))
plot(forecast(z_pron))


