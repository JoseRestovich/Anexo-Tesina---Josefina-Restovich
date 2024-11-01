########### CASO 2 ###########
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
set.seed(15)

# Parámetros de la distribución normal
media2<- 41
desvio2 <- 0.70
cantidad_muestras <- 5000

#Valores teóricos
p1__2<-pnorm(43,mean = 40,sd=desvio2)-pnorm(37,mean = 40,sd=desvio2)
plow__2<-pnorm(media2,mean = media2,sd=desvio2)-pnorm(37,mean =media2 ,sd=desvio2)    
pupp__2<-pnorm(43,mean = media2,sd=desvio2)-pnorm(media2,mean = media2,sd=desvio2)
p__2 <- 2 * min(plow__2, pupp__2)
p2__2<-pnorm(43,mean = media2,sd=desvio2)-pnorm(37,mean = media2,sd=desvio2)
sigmat__2<-sqrt((desvio2)^2+(media2-40)^2)
p3__2<-pnorm(43,mean = 40,sd=sigmat__2)-pnorm(37,mean = 40,sd=sigmat__2)

Cp_p1__2<-(-1/3*qnorm((1-p1__2)/2))
Cpk_p__2<-(-1/3*qnorm((1-p__2)/2))
Cpk_p2__2<-(-1/3*qnorm((1-p2__2)/2))
Cpm_p3__2<-(-1/3*qnorm((1-p3__2)/2))

valores_teoricos_2 <- data.frame(
  Variable = c("p1__2", "plow__2", "pupp__2", "p__2", "p2__2", "sigmat__2", "p3__2", "Cp_p1__2", "Cpk_p__2", "Cpk_p2__2", "Cpm_p3__2"),
  Valor = c(p1__2, plow__2, pupp__2, p__2, p2__2, sigmat__2, p3__2, Cp_p1__2, Cpk_p__2, Cpk_p2__2, Cpm_p3__2)
)


########## n=20 ##########
tamano_muestra <- 20

# Generar las muestras
muestras2_20 <- matrix(rnorm(tamano_muestra * cantidad_muestras, mean = media2, sd = desvio2), ncol = tamano_muestra)
#En cada fila hay una muestra

# Calcular las estimaciones de mu y sigma en cada muestra
desvios2_20 <- apply(muestras2_20, 1, sd)
medias2_20<-apply(muestras2_20,1,mean)
datos2_20<-cbind(Medias=medias2_20,Desvios=desvios2_20)
datos2_20<-as.data.frame(datos2_20)

#Calculo de SigmaT
datos2_20<-mutate(datos2_20,
                  sigmat=(sqrt((datos2_20$Desvios)^2+(datos2_20$Medias-40)^2)))

#Calculo de todas las probabilidades
datos2_20<-mutate(datos2_20,p1=(pnorm(43,mean = 40,sd=datos2_20$Desvios)-pnorm(37,mean = 40,sd=datos2_20$Desvios)),
                  plow=(pnorm(datos2_20$Medias,mean = datos2_20$Medias,sd=datos2_20$Desvios)-pnorm(37,mean = datos2_20$Medias,sd=datos2_20$Desvios)),
                  pupp=(pnorm(43,mean = datos2_20$Medias,sd=datos2_20$Desvios)-pnorm(datos2_20$Medias,mean = datos2_20$Medias,sd=datos2_20$Desvios)),
                  p2=(pnorm(43,mean = datos2_20$Medias,sd=datos2_20$Desvios)-pnorm(37,mean = datos2_20$Medias,sd=datos2_20$Desvios)),
                  p3=(pnorm(43,mean = 40,sd=datos2_20$sigmat)-pnorm(37,mean = 40,sd=datos2_20$sigmat)))

# Calcular la probabilidad p para cada muestra
datos2_20$p <- apply(datos2_20, 1, function(row) {
  p <- 2 * min(row['plow'], row['pupp'])
  return(p)
})


#Calculo de todos los índices
datos2_20 <- mutate(datos2_20, 
                    Cp_p1=(-1/3*qnorm((1-datos2_20$p1)/2)),
                    Cpk_p=(-1/3*qnorm((1-datos2_20$p)/2)),
                    Cpk_p2=(-1/3*qnorm((1-datos2_20$p2)/2)),
                    Cpm_p3=(-1/3*qnorm((1-datos2_20$p3)/2)))


# Crear un dataframe para almacenar los percentiles
intervalos2_20 <- data.frame(Estadistica = c("p1", "p", "p2", "p3", "Cp_p1", "Cpk_p", "Cpk_p2", "Cpm_p3"),
                              Percentil_2.5 = numeric(8),
                              Percentil_97.5 = numeric(8),
                              Medias = numeric(8),
                              Desvios = numeric(8))

# Calcular los percentiles 2.5 y 97.5 para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos2_20$Estadistica[i])
  
  # Calcular los percentiles 2.5 y 97.5
  percentil_2.5 <- quantile(datos2_20[, estadistica], probs = 0.025)
  percentil_97.5 <- quantile(datos2_20[, estadistica], probs = 0.975)
  
  intervalos2_20$Percentil_2.5[i] <- percentil_2.5
  intervalos2_20$Percentil_97.5[i] <- percentil_97.5
  intervalos2_20$Medias[i] <- mean(datos2_20[, estadistica])
  intervalos2_20$Desvios[i] <- sd(datos2_20[, estadistica])
}



sesgo_p1_r <- (mean(datos2_20$p1) - p1__2)/p1__2
sesgo_p_r <- (mean(datos2_20$p) - p__2)/p__2
sesgo_p2_r <- (mean(datos2_20$p2) - p2__2)/p2__2
sesgo_p3_r <- (mean(datos2_20$p3) - p3__2)/p3__2
sesgo_cp_p1_r <- (mean(datos2_20$Cp_p1) - Cp_p1__2)/Cp_p1__2
sesgo_cpk_p_r <- (mean(datos2_20$Cpk_p) - Cpk_p__2)/Cpk_p__2
sesgo_cpk_p2_r <- (mean(datos2_20$Cpk_p2) - Cpk_p2__2)/Cpk_p2__2
sesgo_cpm_p3_r <- (mean(datos2_20$Cpm_p3) - Cpm_p3__2)/Cpm_p3__2

sesgo_p1 <- (mean(datos2_20$p1) - p1__2)
sesgo_p <- (mean(datos2_20$p) - p__2)
sesgo_p2 <- (mean(datos2_20$p2) - p2__2)
sesgo_p3 <- (mean(datos2_20$p3) - p3__2)
sesgo_cp_p1 <- (mean(datos2_20$Cp_p1) - Cp_p1__2)
sesgo_cpk_p <- (mean(datos2_20$Cpk_p) - Cpk_p__2)
sesgo_cpk_p2 <- (mean(datos2_20$Cpk_p2) - Cpk_p2__2)
sesgo_cpm_p3 <- (mean(datos2_20$Cpm_p3) - Cpm_p3__2)

sesgos_2_20<-c(sesgo_p1,sesgo_p,sesgo_p2,sesgo_p3,sesgo_cp_p1,sesgo_cpk_p,sesgo_cpk_p2,sesgo_cpm_p3)
sesgos_2_20_r<-c(sesgo_p1_r,sesgo_p_r,sesgo_p2_r,sesgo_p3_r,sesgo_cp_p1_r,sesgo_cpk_p_r,sesgo_cpk_p2_r,sesgo_cpm_p3_r)

intervalos2_20$SesgoRelativo<-sesgos_2_20_r
intervalos2_20$Sesgo<-sesgos_2_20
intervalos2_20$ECM<-(intervalos2_20$Desvios)^2+(intervalos2_20$Sesgo)^2

# Imprimir los resultados
print(intervalos2_20)



########## n=25 ##########
tamano_muestra <- 25
set.seed(15)
# Generar las muestras
muestras2_25 <- matrix(rnorm(tamano_muestra * cantidad_muestras, mean = media2, sd = desvio2), ncol = tamano_muestra)
#En cada fila hay una muestra

# Calcular las estimaciones de mu y sigma en cada muestra
desvios2_25 <- apply(muestras2_25, 1, sd)
medias2_25<-apply(muestras2_25,1,mean)
datos2_25<-cbind(Medias=medias2_25,Desvios=desvios2_25)
datos2_25<-as.data.frame(datos2_25)

#Calculo de SigmaT
datos2_25<-mutate(datos2_25,
                  sigmat=(sqrt((datos2_25$Desvios)^2+(datos2_25$Medias-40)^2)))

#Calculo de todas las probabilidades
datos2_25<-mutate(datos2_25,p1=(pnorm(43,mean = 40,sd=datos2_25$Desvios)-pnorm(37,mean = 40,sd=datos2_25$Desvios)),
                  plow=(pnorm(datos2_25$Medias,mean = datos2_25$Medias,sd=datos2_25$Desvios)-pnorm(37,mean = datos2_25$Medias,sd=datos2_25$Desvios)),
                  pupp=(pnorm(43,mean = datos2_25$Medias,sd=datos2_25$Desvios)-pnorm(datos2_25$Medias,mean = datos2_25$Medias,sd=datos2_25$Desvios)),
                  p2=(pnorm(43,mean = datos2_25$Medias,sd=datos2_25$Desvios)-pnorm(37,mean = datos2_25$Medias,sd=datos2_25$Desvios)),
                  p3=(pnorm(43,mean = 40,sd=datos2_25$sigmat)-pnorm(37,mean = 40,sd=datos2_25$sigmat)))

# Calcular la probabilidad p para cada muestra
datos2_25$p <- apply(datos2_25, 1, function(row) {
  p <- 2 * min(row['plow'], row['pupp'])
  return(p)
})


#Calculo de todos los índices
datos2_25 <- mutate(datos2_25, 
                    Cp_p1=(-1/3*qnorm((1-datos2_25$p1)/2)),
                    Cpk_p=(-1/3*qnorm((1-datos2_25$p)/2)),
                    Cpk_p2=(-1/3*qnorm((1-datos2_25$p2)/2)),
                    Cpm_p3=(-1/3*qnorm((1-datos2_25$p3)/2)))


# Crear un dataframe para almacenar los percentiles
intervalos2_25 <- data.frame(Estadistica = c("p1", "p", "p2", "p3", "Cp_p1", "Cpk_p", "Cpk_p2", "Cpm_p3"),
                              Percentil_2.5 = numeric(8),
                              Percentil_97.5 = numeric(8),
                              Medias = numeric(8),
                              Desvios = numeric(8))

# Calcular los percentiles 2.5 y 97.5 para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos2_25$Estadistica[i])
  
  # Calcular los percentiles 2.5 y 97.5
  percentil_2.5 <- quantile(datos2_25[, estadistica], probs = 0.025)
  percentil_97.5 <- quantile(datos2_25[, estadistica], probs = 0.975)
  intervalos2_25$Percentil_2.5[i] <- percentil_2.5
  intervalos2_25$Percentil_97.5[i] <- percentil_97.5
  intervalos2_25$Medias[i] <- mean(datos2_25[, estadistica])
  intervalos2_25$Desvios[i] <- sd(datos2_25[, estadistica])
}


sesgo_p1_r <- (mean(datos2_25$p1) - p1__2)/p1__2
sesgo_p_r <- (mean(datos2_25$p) - p__2)/p__2
sesgo_p2_r <- (mean(datos2_25$p2) - p2__2)/p2__2
sesgo_p3_r <- (mean(datos2_25$p3) - p3__2)/p3__2
sesgo_cp_p1_r <- (mean(datos2_25$Cp_p1) - Cp_p1__2)/Cp_p1__2
sesgo_cpk_p_r <- (mean(datos2_25$Cpk_p) - Cpk_p__2)/Cpk_p__2
sesgo_cpk_p2_r <- (mean(datos2_25$Cpk_p2) - Cpk_p2__2)/Cpk_p2__2
sesgo_cpm_p3_r <- (mean(datos2_25$Cpm_p3) - Cpm_p3__2)/Cpm_p3__2

sesgo_p1 <- (mean(datos2_25$p1) - p1__2)
sesgo_p <- (mean(datos2_25$p) - p__2)
sesgo_p2 <- (mean(datos2_25$p2) - p2__2)
sesgo_p3 <- (mean(datos2_25$p3) - p3__2)
sesgo_cp_p1 <- (mean(datos2_25$Cp_p1) - Cp_p1__2)
sesgo_cpk_p <- (mean(datos2_25$Cpk_p) - Cpk_p__2)
sesgo_cpk_p2 <- (mean(datos2_25$Cpk_p2) - Cpk_p2__2)
sesgo_cpm_p3 <- (mean(datos2_25$Cpm_p3) - Cpm_p3__2)

sesgos_2_25<-c(sesgo_p1,sesgo_p,sesgo_p2,sesgo_p3,sesgo_cp_p1,sesgo_cpk_p,sesgo_cpk_p2,sesgo_cpm_p3)
sesgos_2_25_r<-c(sesgo_p1_r,sesgo_p_r,sesgo_p2_r,sesgo_p3_r,sesgo_cp_p1_r,sesgo_cpk_p_r,sesgo_cpk_p2_r,sesgo_cpm_p3_r)

intervalos2_25$SesgoRelativo<-sesgos_2_25
intervalos2_25$Sesgo<-sesgos_2_25
intervalos2_25$ECM<-(intervalos2_25$Desvios)^2+(intervalos2_25$Sesgo)^2

# Imprimir los resultados
print(intervalos2_25)



########## n=50 ##########
tamano_muestra <- 50
set.seed(15)
# Generar las muestras
muestras2_50 <- matrix(rnorm(tamano_muestra * cantidad_muestras, mean = media2, sd = desvio2), ncol = tamano_muestra)
#En cada fila hay una muestra

# Calcular las estimaciones de mu y sigma en cada muestra
desvios2_50 <- apply(muestras2_50, 1, sd)
medias2_50<-apply(muestras2_50,1,mean)
datos2_50<-cbind(Medias=medias2_50,Desvios=desvios2_50)
datos2_50<-as.data.frame(datos2_50)

#Calculo de SigmaT
datos2_50<-mutate(datos2_50,
                  sigmat=(sqrt((datos2_50$Desvios)^2+(datos2_50$Medias-40)^2)))

#Calculo de todas las probabilidades
datos2_50<-mutate(datos2_50,p1=(pnorm(43,mean = 40,sd=datos2_50$Desvios)-pnorm(37,mean = 40,sd=datos2_50$Desvios)),
                  plow=(pnorm(datos2_50$Medias,mean = datos2_50$Medias,sd=datos2_50$Desvios)-pnorm(37,mean = datos2_50$Medias,sd=datos2_50$Desvios)),
                  pupp=(pnorm(43,mean = datos2_50$Medias,sd=datos2_50$Desvios)-pnorm(datos2_50$Medias,mean = datos2_50$Medias,sd=datos2_50$Desvios)),
                  p2=(pnorm(43,mean = datos2_50$Medias,sd=datos2_50$Desvios)-pnorm(37,mean = datos2_50$Medias,sd=datos2_50$Desvios)),
                  p3=(pnorm(43,mean = 40,sd=datos2_50$sigmat)-pnorm(37,mean = 40,sd=datos2_50$sigmat)))

# Calcular la probabilidad p para cada muestra
datos2_50$p <- apply(datos2_50, 1, function(row) {
  p <- 2 * min(row['plow'], row['pupp'])
  return(p)
})


#Calculo de todos los índices
datos2_50 <- mutate(datos2_50, 
                    Cp_p1=(-1/3*qnorm((1-datos2_50$p1)/2)),
                    Cpk_p=(-1/3*qnorm((1-datos2_50$p)/2)),
                    Cpk_p2=(-1/3*qnorm((1-datos2_50$p2)/2)),
                    Cpm_p3=(-1/3*qnorm((1-datos2_50$p3)/2)))


# Crear un dataframe para almacenar los percentiles
intervalos2_50 <- data.frame(Estadistica = c("p1", "p", "p2", "p3", "Cp_p1", "Cpk_p", "Cpk_p2", "Cpm_p3"),
                              Percentil_2.5 = numeric(8),
                              Percentil_97.5 = numeric(8),
                              Medias = numeric(8),
                              Desvios = numeric(8))

# Calcular los percentiles 2.5 y 97.5 para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos2_50$Estadistica[i])
  
  # Calcular los percentiles 2.5 y 97.5
  percentil_2.5 <- quantile(datos2_50[, estadistica], probs = 0.025)
  percentil_97.5 <- quantile(datos2_50[, estadistica], probs = 0.975)
  intervalos2_50$Percentil_2.5[i] <- percentil_2.5
  intervalos2_50$Percentil_97.5[i] <- percentil_97.5
  intervalos2_50$Medias[i] <- mean(datos2_50[, estadistica])
  intervalos2_50$Desvios[i] <- sd(datos2_50[, estadistica])
}


sesgo_p1_r <- (mean(datos2_50$p1) - p1__2)/p1__2
sesgo_p_r <- (mean(datos2_50$p) - p__2)/p__2
sesgo_p2_r <- (mean(datos2_50$p2) - p2__2)/p2__2
sesgo_p3_r <- (mean(datos2_50$p3) - p3__2)/p3__2
sesgo_cp_p1_r <- (mean(datos2_50$Cp_p1) - Cp_p1__2)/Cp_p1__2
sesgo_cpk_p_r <- (mean(datos2_50$Cpk_p) - Cpk_p__2)/Cpk_p__2
sesgo_cpk_p2_r <- (mean(datos2_50$Cpk_p2) - Cpk_p2__2)/Cpk_p2__2
sesgo_cpm_p3_r <- (mean(datos2_50$Cpm_p3) - Cpm_p3__2)/Cpm_p3__2

sesgo_p1 <- (mean(datos2_50$p1) - p1__2)
sesgo_p <- (mean(datos2_50$p) - p__2)
sesgo_p2 <- (mean(datos2_50$p2) - p2__2)
sesgo_p3 <- (mean(datos2_50$p3) - p3__2)
sesgo_cp_p1 <- (mean(datos2_50$Cp_p1) - Cp_p1__2)
sesgo_cpk_p <- (mean(datos2_50$Cpk_p) - Cpk_p__2)
sesgo_cpk_p2 <- (mean(datos2_50$Cpk_p2) - Cpk_p2__2)
sesgo_cpm_p3 <- (mean(datos2_50$Cpm_p3) - Cpm_p3__2)

sesgos_2_50_r<-c(sesgo_p1_r,sesgo_p_r,sesgo_p2_r,sesgo_p3_r,sesgo_cp_p1_r,sesgo_cpk_p_r,sesgo_cpk_p2_r,sesgo_cpm_p3_r)
sesgos_2_50<-c(sesgo_p1,sesgo_p,sesgo_p2,sesgo_p3,sesgo_cp_p1,sesgo_cpk_p,sesgo_cpk_p2,sesgo_cpm_p3)

intervalos2_50$SesgoRelativo<-sesgos_2_50_r
intervalos2_50$Sesgo<-sesgos_2_50
intervalos2_50$ECM<-(intervalos2_50$Desvios)^2+(intervalos2_50$Sesgo)^2

# Imprimir los resultados
print(intervalos2_50)



########## n=100 ##########
tamano_muestra <- 100
set.seed(15)
# Generar las muestras
muestras2_100 <- matrix(rnorm(tamano_muestra * cantidad_muestras, mean = media2, sd = desvio2), ncol = tamano_muestra)
#En cada fila hay una muestra

# Calcular las estimaciones de mu y sigma en cada muestra
desvios2_100 <- apply(muestras2_100, 1, sd)
medias2_100<-apply(muestras2_100,1,mean)
datos2_100<-cbind(Medias=medias2_100,Desvios=desvios2_100)
datos2_100<-as.data.frame(datos2_100)

#Calculo de SigmaT
datos2_100<-mutate(datos2_100,
                   sigmat=(sqrt((datos2_100$Desvios)^2+(datos2_100$Medias-40)^2)))

#Calculo de todas las probabilidades
datos2_100<-mutate(datos2_100,p1=(pnorm(43,mean = 40,sd=datos2_100$Desvios)-pnorm(37,mean = 40,sd=datos2_100$Desvios)),
                   plow=(pnorm(datos2_100$Medias,mean = datos2_100$Medias,sd=datos2_100$Desvios)-pnorm(37,mean = datos2_100$Medias,sd=datos2_100$Desvios)),
                   pupp=(pnorm(43,mean = datos2_100$Medias,sd=datos2_100$Desvios)-pnorm(datos2_100$Medias,mean = datos2_100$Medias,sd=datos2_100$Desvios)),
                   p2=(pnorm(43,mean = datos2_100$Medias,sd=datos2_100$Desvios)-pnorm(37,mean = datos2_100$Medias,sd=datos2_100$Desvios)),
                   p3=(pnorm(43,mean = 40,sd=datos2_100$sigmat)-pnorm(37,mean = 40,sd=datos2_100$sigmat)))

# Calcular la probabilidad p para cada muestra
datos2_100$p <- apply(datos2_100, 1, function(row) {
  p <- 2 * min(row['plow'], row['pupp'])
  return(p)
})


#Calculo de todos los índices
datos2_100 <- mutate(datos2_100, 
                     Cp_p1=(-1/3*qnorm((1-datos2_100$p1)/2)),
                     Cpk_p=(-1/3*qnorm((1-datos2_100$p)/2)),
                     Cpk_p2=(-1/3*qnorm((1-datos2_100$p2)/2)),
                     Cpm_p3=(-1/3*qnorm((1-datos2_100$p3)/2)))


# Crear un dataframe para almacenar los percentiles
intervalos2_100 <- data.frame(Estadistica = c("p1", "p", "p2", "p3", "Cp_p1", "Cpk_p", "Cpk_p2", "Cpm_p3"),
                              Percentil_2.5 = numeric(8),
                              Percentil_97.5 = numeric(8),
                              Medias = numeric(8),
                              Desvios = numeric(8))

# Calcular los percentiles 2.5 y 97.5 para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos2_100$Estadistica[i])
  
  # Calcular los percentiles 2.5 y 97.5
  percentil_2.5 <- quantile(datos2_100[, estadistica], probs = 0.025)
  percentil_97.5 <- quantile(datos2_100[, estadistica], probs = 0.975)
  intervalos2_100$Percentil_2.5[i] <- percentil_2.5
  intervalos2_100$Percentil_97.5[i] <- percentil_97.5
  intervalos2_100$Medias[i] <- mean(datos2_100[, estadistica])
  intervalos2_100$Desvios[i] <- sd(datos2_100[, estadistica])
}



sesgo_p1_r <- (mean(datos2_100$p1) - p1__2)/p1__2
sesgo_p_r <- (mean(datos2_100$p) - p__2)/p__2
sesgo_p2_r <- (mean(datos2_100$p2) - p2__2)/p2__2
sesgo_p3_r <- (mean(datos2_100$p3) - p3__2)/p3__2
sesgo_cp_p1_r <- (mean(datos2_100$Cp_p1) - Cp_p1__2)/Cp_p1__2
sesgo_cpk_p_r <- (mean(datos2_100$Cpk_p) - Cpk_p__2)/Cpk_p__2
sesgo_cpk_p2_r <- (mean(datos2_100$Cpk_p2) - Cpk_p2__2)/Cpk_p2__2
sesgo_cpm_p3_r <- (mean(datos2_100$Cpm_p3) - Cpm_p3__2)/Cpm_p3__2

sesgo_p1 <- (mean(datos2_100$p1) - p1__2)
sesgo_p <- (mean(datos2_100$p) - p__2)
sesgo_p2 <- (mean(datos2_100$p2) - p2__2)
sesgo_p3 <- (mean(datos2_100$p3) - p3__2)
sesgo_cp_p1 <- (mean(datos2_100$Cp_p1) - Cp_p1__2)
sesgo_cpk_p <- (mean(datos2_100$Cpk_p) - Cpk_p__2)
sesgo_cpk_p2 <- (mean(datos2_100$Cpk_p2) - Cpk_p2__2)
sesgo_cpm_p3 <- (mean(datos2_100$Cpm_p3) - Cpm_p3__2)

sesgos_2_100_r<-c(sesgo_p1_r,sesgo_p_r,sesgo_p2_r,sesgo_p3_r,sesgo_cp_p1_r,sesgo_cpk_p_r,sesgo_cpk_p2_r,sesgo_cpm_p3_r)
sesgos_2_100<-c(sesgo_p1,sesgo_p,sesgo_p2,sesgo_p3,sesgo_cp_p1,sesgo_cpk_p,sesgo_cpk_p2,sesgo_cpm_p3)

intervalos2_100$SesgoRelativo<-sesgos_2_100_r
intervalos2_100$Sesgo<-sesgos_2_100
intervalos2_100$ECM<-(intervalos2_100$Desvios)^2+(intervalos2_100$Sesgo)^2

# Imprimir los resultados
print(intervalos2_100)




############GRÁFICOs####################
datos2_20$n<-"20"
datos2_25$n<-"25"
datos2_50$n<-"50"
datos2_100$n<-"100"


##############Cp_p1
# Definir los dataframes individuales
g20_Cp_p1 <- datos2_20[, c("Cp_p1", "n")]
g25_Cp_p1 <- datos2_25[, c("Cp_p1", "n")]
g50_Cp_p1 <- datos2_50[, c("Cp_p1", "n")]
g100_Cp_p1 <- datos2_100[, c("Cp_p1", "n")]

# Combinar los dataframes
df_Cp_p1 <- rbind(g20_Cp_p1, g25_Cp_p1, g50_Cp_p1, g100_Cp_p1)

# Ordenar los niveles de la variable n
df_Cp_p1$n <- factor(df_Cp_p1$n, levels = c("20", "25", "50", "100"))

##############Cpk_p
# Definir los dataframes individuales
g20_Cpk_p <- datos2_20[, c("Cpk_p", "n")]
g25_Cpk_p <- datos2_25[, c("Cpk_p", "n")]
g50_Cpk_p <- datos2_50[, c("Cpk_p", "n")]
g100_Cpk_p <- datos2_100[, c("Cpk_p", "n")]

# Combinar los dataframes
df_Cpk_p <- rbind(g20_Cpk_p, g25_Cpk_p, g50_Cpk_p, g100_Cpk_p)

# Ordenar los niveles de la variable n
df_Cpk_p$n <- factor(df_Cpk_p$n, levels = c("20", "25", "50", "100"))

##############Cpk_p2
# Definir los dataframes individuales
g20_Cpk_p2 <- datos2_20[, c("Cpk_p2", "n")]
g25_Cpk_p2 <- datos2_25[, c("Cpk_p2", "n")]
g50_Cpk_p2 <- datos2_50[, c("Cpk_p2", "n")]
g100_Cpk_p2 <- datos2_100[, c("Cpk_p2", "n")]

# Combinar los dataframes
df_Cpk_p2 <- rbind(g20_Cpk_p2, g25_Cpk_p2, g50_Cpk_p2, g100_Cpk_p2)

# Ordenar los niveles de la variable n
df_Cpk_p2$n <- factor(df_Cpk_p2$n, levels = c("20", "25", "50", "100"))

##############Cpm_p3
# Definir los dataframes individuales
g20_Cpm_p3 <- datos2_20[, c("Cpm_p3", "n")]
g25_Cpm_p3 <- datos2_25[, c("Cpm_p3", "n")]
g50_Cpm_p3 <- datos2_50[, c("Cpm_p3", "n")]
g100_Cpm_p3 <- datos2_100[, c("Cpm_p3", "n")]

# Combinar los dataframes
df_Cpm_p3 <- rbind(g20_Cpm_p3, g25_Cpm_p3, g50_Cpm_p3, g100_Cpm_p3)

# Ordenar los niveles de la variable n
df_Cpm_p3$n <- factor(df_Cpm_p3$n, levels = c("20", "25", "50", "100"))


# Crear los gráficos individuales con la leyenda activada
plot_Cp_p1 <- ggplot(df_Cp_p1, aes(x = Cp_p1, fill = n)) +
  geom_density(aes(color = n), alpha = 0.5, size = 1.5) +
  labs(title = expression("Gráfico de Densidad de " * C[p(p[1])]), x = "Valores", y = "Densidad") +
  scale_fill_manual(values = c("20" = "olivedrab3", "25" = "cornflowerblue", "50" = "turquoise", "100" = "coral"),
                    name = "Tamaño de muestra", labels = c("20", "25", "50", "100")) +
  scale_color_manual(values = c("20" = "olivedrab3", "25" = "cornflowerblue", "50" = "turquoise", "100" = "coral"),
                     name = "Tamaño de muestra", labels = c("20", "25", "50", "100")) +
  theme_bw() +
  geom_vline(xintercept = Cp_p1__2, linetype = "dashed", color = "black")

plot_Cpk_p <- ggplot(df_Cpk_p, aes(x = Cpk_p, fill = n)) +
  geom_density(aes(color = n), alpha = 0.5, size = 1.5) +
  labs(title = expression("Gráfico de Densidad de " * C[pk(p)]), x = "Valores", y = "Densidad") +
  scale_fill_manual(values = c("20" = "olivedrab3", "25" = "cornflowerblue", "50" = "turquoise", "100" = "coral"),
                    name = "Tamaño de muestra", labels = c("20", "25", "50", "100")) +
  scale_color_manual(values = c("20" = "olivedrab3", "25" = "cornflowerblue", "50" = "turquoise", "100" = "coral"),
                     name = "Tamaño de muestra", labels = c("20", "25", "50", "100")) +
  theme_bw() +
  geom_vline(xintercept = Cpk_p__2, linetype = "dashed", color = "black")

plot_Cpk_p2 <- ggplot(df_Cpk_p2, aes(x = Cpk_p2, fill = n)) +
  geom_density(aes(color = n), alpha = 0.5, size = 1.5) +
  labs(title = expression("Gráfico de Densidad de " * C[pk(p[2])]), x = "Valores", y = "Densidad") +
  scale_fill_manual(values = c("20" = "olivedrab3", "25" = "cornflowerblue", "50" = "turquoise", "100" = "coral"),
                    name = "Tamaño de muestra", labels = c("20", "25", "50", "100")) +
  scale_color_manual(values = c("20" = "olivedrab3", "25" = "cornflowerblue", "50" = "turquoise", "100" = "coral"),
                     name = "Tamaño de muestra", labels = c("20", "25", "50", "100")) +
  theme_bw() +
  geom_vline(xintercept = Cpk_p2__2, linetype = "dashed", color = "black")

plot_Cpm_p3 <- ggplot(df_Cpm_p3, aes(x = Cpm_p3, fill = n)) +
  geom_density(aes(color = n), alpha = 0.5, size = 1.5) +
  labs(title = expression("Gráfico de Densidad de " * C[pm(p[3])]), x = "Valores", y = "Densidad") +
  scale_fill_manual(values = c("20" = "olivedrab3", "25" = "cornflowerblue", "50" = "turquoise", "100" = "coral"),
                    name = "Tamaño de muestra", labels = c("20", "25", "50", "100")) +
  scale_color_manual(values = c("20" = "olivedrab3", "25" = "cornflowerblue", "50" = "turquoise", "100" = "coral"),
                     name = "Tamaño de muestra", labels = c("20", "25", "50", "100")) +
  theme_bw() +
  geom_vline(xintercept = Cpm_p3__2, linetype = "dashed", color = "black")


# Función para extraer la leyenda de un gráfico
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
# Extraer la leyenda del primer gráfico
legend <- g_legend(plot_Cp_p1)


# Quitar la leyenda de los gráficos
plot_Cp_p1 <- plot_Cp_p1 + theme(legend.position = "none")
plot_Cpk_p <- plot_Cpk_p + theme(legend.position = "none")
plot_Cpk_p2 <- plot_Cpk_p2 + theme(legend.position = "none")
plot_Cpm_p3 <- plot_Cpm_p3 + theme(legend.position = "none")

# Organizar los gráficos y la leyenda en un panel
grid.arrange(
  arrangeGrob(plot_Cp_p1, plot_Cpk_p, plot_Cpk_p2, plot_Cpm_p3, ncol = 2),
  legend, ncol = 2, widths = c(5, 1)
)
