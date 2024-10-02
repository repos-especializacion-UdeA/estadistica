#lee el archivo .txt
datos<-read.table(file.choose(), header=T )
#CARGA BASE DE DATOS Y NOMBRES DE VARIABLES EN R
attach(datos)
#muestra datos
datos

#son necesarios 4 paquetes
# es necesario instalarlos antes de cargarlos.
library(car)
library(leaps)
install.packages("remotes")
remotes::install_github("cran/perturb")
library(perturb)
library(perturbR)
#AJUSTE DEL MODELO
modelo<-lm(calidad ~ aroma+cuerpo+sabor+fuerza)

summary(modelo)
#Anova(modelo) #PRODUCE SUMA DE CUADRADOS SS2

#se debe calcular el siguiente anova, pero para el analisis se usan los
#resulatos de miAnova
anova(modelo) #PRODUCE SUMAS DE CUADRADOS SS1
#PARA OBTENER LA ANOVA DEL MODELO DE RLM CREAMOS LA SIGUIENTE FUNCI�N

miAnova<-function(modeloreg){
SSq<-unlist(anova(modeloreg)["Sum Sq"])
k<-length(SSq)-1
SSR<-sum(SSq[1:k])
SSE<-SSq[(k+1)]
MSR<-SSR/k
df.error<-unlist(anova(modeloreg)["Df"])[k+1]
MSE<-SSE/df.error
F0<-MSR/MSE
VP<-pf(F0,k,df.error,lower.tail=F)
result<-
data.frame(SumSq=c(SSR,SSE),Df=c(k,df.error),MeanSq=c(MSR,MSE),F0=c(round(F0,digits=3),' '),
P.value=c(format(VP,scientific = TRUE,digits=3),' '),row.names =c("Modelo","Error"))
cat("Tabla ANOVA Modelo de Regresi�n","\n")
result
}

miAnova(modelo)

#TABLA DE TODAS LAS REGRESIONES POSIBLES
allregtable<-function(modeloreg,respuesta){
t1<-summary(regsubsets(model.matrix(modeloreg)[,-1],respuesta,nbest=20),all.best=TRUE)
t2<-as.vector(apply(t1$which[,-1],1,sum))
t3<-apply(t1$which[,-1],1,function(x) as.character(paste(colnames(
model.matrix(modeloreg)[,-1])[x],collapse=" ")))
results<-data.frame(NoOfVars=t2,R2=round(t1$rsq,4),adjR2=round(t1$adjr2,4),
SSE=round(t1$rss,5),Cp=round(t1$cp,4),
MSE=round(t1$rss/(nrow(model.matrix(modeloreg)[,-1])-(t2+1)),5),
Variables.in.model=t3)
results
}
allregtable(modelo,calidad)

#CALCULO DE RESIDUALES ESTUDENTIZADOS
restud<-round(rstudent(modelo),4)
#GR�FICO DE RESIDUALES VS. VALORES PREDICHOS
plot(fitted(modelo),restud,xlab="Valores Ajustados",ylab="Residuales Estudentizados",
main="Residuales Estudentizados vs. Valores Ajustados")
abline(h=0,lty=2,col=2)

#PRUEBA DE NORMALIDAD DE SHAPIRO-WILK
test<-shapiro.test(restud)
test

#GR�FICO DE CUANTILES NORMALES
qqnorm(restud,cex=1.5)
qqline(restud,lty=2,lwd=2,col=2)

require(car)
qqplot(restud, las =1)


#PRUEBA DE ALEATORIEDAD
install.packages("tseries")
library(tseries)
x<-factor(sign(restud))
runs.test(x)

#Gr�fica de aleatoriedad
plot(restud, main='Aleatoriedad', xlab='Orden de experimentaci�n', ylab = 'Residuales')
lines(restud)

#GR�FICO DE RESIDUALES VS. VALORES PREDICHOS
plot(fitted(modelo),restud,xlab="Valores Ajustados",ylab="Residuales Estudentizados",
     main="Residuales Estudentizados vs. Valores Ajustados")
abline(h=0,lty=2,col=2)

#DIAGN�STICOS DE DATOS AT�PICOS, DE BALANCEO E INFLUENCIALES
t1<-predict(modelo,se.fit=T)
t2<-round(residuals(modelo),4)
t3<-round(cooks.distance(modelo),4)
t4<-round(hatvalues(modelo),4)
t5<-round(dffits(modelo),4)
data.frame(calidad,yhat=t1$fit,se.yhat=t1$se.fit,residuals=t2,res.estud=restud,Cooks.D=t3,
hii.value=t4,Dffits=t5)

##PARA EVALUAR MULTICOLINEALIDAD##
##################################

#MATRIZ DE DISPERSI�N
plot(datos)
#MATRIZ DE CORRELACIONES PARA VARIABLES
cor(datos)
#para ensayar no tener en cuenta una columna que no contiene numeros
#DEPARTAMENTO=which(colnames(datos)=="DEPARTAMENTO")
#cor(datos[,-DEPARTAMENTOS])


#CREANDO FUNCI�N PARA EXTRAER COEFICIENTES ESTIMADOS, SUS IC DEL 95%, VIF'S Y COEFICIENTES ESTANDARIZADOS
miscoeficientes=function(modeloreg,datosreg){
coefi=coef(modeloreg)
datos2=as.data.frame(scale(datosreg))
coef.std=c(0,coef(lm(update(formula(modeloreg),~.+0),datos2)))
limites=confint(modeloreg,level=0.95)
vifs=c(0,vif(modeloreg))
resul=data.frame(Estimación=coefi,Limites=limites,Vif=vifs,Coef.Std=coef.std)
cat("Coeficientes estimados, sus I.C, Vifs y Coeficientes estimados estandarizados","\n")
resul
}
#CREANDO FUNCION PARA EXTRAER RESULTADOS PARA DIAGN�STICOS DE MULTICOLINEALIDAD
misDiagnostcolin=function(modeloreg,centrar=F){
if(centrar==F){
X=model.matrix(modeloreg)
val.prop=prcomp(X,center=FALSE,scale=TRUE)$sdev^2
Ind=colldiag(modeloreg)
resul=data.frame(Val.propio=val.prop,Ind.Cond=Ind$condindx,Pi=Ind$pi)
cat("Diagnósticos Multicolinealidad - Intercepto incluido","\n",
"índices de Condición y Proporciones de Varianza","\n")
}
else{
X=model.matrix(modeloreg)[,-1]
val.prop=prcomp(X,center=TRUE,scale=TRUE)$sdev^2
Ind=colldiag(modeloreg,center=TRUE,scale=TRUE)
resul=data.frame(Val.propio=val.prop,Ind.Cond=Ind$condindx,Pi=Ind$pi)
cat("Diagnósticos Multicolinealidad - Intercepto ajustado","\n",
"índices de Condición y Proporciones de Varianza","\n")
}
resul
}

miscoeficientes(modelo,datos)
misDiagnostcolin(modelo)
misDiagnostcolin(modelo,centrar=T)

11#C�LCULO DE ESTAD�STICOS PARA PUNTOS DE PREDICCI�N [h00,y0hat y se(y0hat)]
x01<-c(1,4.5,5.4,6.0,3.6)
x02<-c(1,3,5,6,10)
xpred<-rbind(x01,x02)
colnames(xpred)<-colnames(model.matrix(modelo))
A<-model.matrix(modelo)
hvalues<-diag(xpred%*%solve(t(A)%*%A)%*%t(xpred))
prednew<-predict(modelo,newdata=data.frame(xpred[,-1]),se.fit=T, level = 0.95, interval = "confidence")
data.frame(h00.value=hvalues,y0hat=prednew$fit,se.y0hat=prednew$se.fit )

