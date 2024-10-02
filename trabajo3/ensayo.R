#install.packages("perturb")
#install.packages("leaps")
#install.packages("car")

library(leaps)
library(car)
datos <- read.csv("Datos.csv", encoding = "UTF-8")
mod = lm(Y ~ X1 + X2 + X3 + X4 + X5,data=datos)
summary(mod)
