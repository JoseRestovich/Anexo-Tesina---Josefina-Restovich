library(ggplot2)
library(tidyverse)
library(mvtnorm)
library(MASS)
library(dplyr)
library(readxl)

cantidad_muestras <- 5000
set.seed(15)

# Carga de datos
datos_multi <- read_excel("C:/Users/josefina.restovich/OneDrive - datalytics.com/Escritorio/Facultad/Tesina/Programas/Multivariado/Practica análisis de capacidad de procesos multivariado - Conjuntos de datos.xlsx", sheet = "Mecánico")
datos_multi <- datos_multi[, c(2, 3, 4, 5, 6, 7, 8)]

#################
### p=2, n=50 ###
#################
datos_2 <- datos_multi[, c(4, 4)]
#Parámetros
LIE_2 <- c(3.00, 3.00)
LSE_2 <- c(17.00, 17.00)
t_2 <- c(10, 10)
medias_2 <- colMeans(datos_2)
sigma_2 <- matrix(c(var(datos_2[,1]),
                    0.50*sqrt(var(datos_2[,1])*var(datos_2[,2])),
                    0.50*sqrt(var(datos_2[,1])*var(datos_2[,2])),
                    var(datos_2[,2])),2,2)
# Tamaño de muestra
tamano_muestra <- 50
#Cálculo de las probabilidades 
p1_2 <- pmvnorm(LIE_2, LSE_2, mean = t_2, sigma = sigma_2)
p2_2 <- pmvnorm(LIE_2, LSE_2, mean = medias_2, sigma = sigma_2)
sigmat_2 <- sigma_2 + t(t(medias_2 - t_2)) %*% t(medias_2 - t_2)
p3_2 <- pmvnorm(LIE_2, LSE_2, mean = t_2, sigma = sigmat_2)
#Cálculo de los índices
MCp_p1_2 <- (-1/3 * qnorm((1 - p1_2)/2))
MCpk_p2_2 <- (-1/3 * qnorm((1 - p2_2)/2))
MCpm_p3_2 <- (-1/3 * qnorm((1 - p3_2)/2))
# Cálculo de los valores reales de NMCp y NMCpm
cprima_2 <- min(((LSE_2[1] - t_2[1]) / sqrt(sigma_2[1, 1])), ((LSE_2[2] - t_2[2]) / sqrt(sigma_2[2, 2])))
NMCp_2 <- cprima_2 / sqrt(qchisq(1 - 0.0027, 2))
D_2 <- (1 +  (t(medias_2 - t_2) %*% solve(sigma_2) %*% (medias_2 - t_2)))^0.5
NMCpm_2 <- NMCp_2 / D_2  

#Cálculo de las medias y sigmas muestrales
medias_2_50 <- matrix(0, nrow = cantidad_muestras, ncol = 2)
sigmas_2_50 <- array(0, dim = c(2, 2, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_2, Sigma = sigma_2)
  medias_2_50[i, ] <- colMeans(muestra)
  sigmas_2_50[, , i] <- cov(muestra)
}

# Crear un dataframe para almacenar los resultados
datos_2_50 <- data.frame(p1_2_50 = numeric(0), p2_2_50 = numeric(0), p3_2_50 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_2_50[, , i]
  p1_2_50 <- pmvnorm(LIE_2, LSE_2, mean = t_2, sigma = sigma_actual)
  p2_2_50 <- pmvnorm(LIE_2, LSE_2, mean = medias_2_50[i, ], sigma = sigma_actual)
  sigmat_2_50 <- sigma_actual + (medias_2_50[i, ] - t_2) %*% t(medias_2_50[i, ] - t_2)
  p3_2_50 <- pmvnorm(LIE_2, LSE_2, mean = t_2, sigma = sigmat_2_50)
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_2[1] - t_2[1]) / sqrt(sigma_actual[1, 1])), ((LSE_2[2] - t_2[2]) / sqrt(sigma_actual[2, 2])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 2))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_2 - medias_2_50[i, ]) %*% solve(sigma_actual) %*% (t_2 - medias_2_50[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_2_50[i, ] <- c(p1_2_50, p2_2_50, p3_2_50, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_2_50 <- mutate(datos_2_50,
                     MCp_p1_2_50 = (-1 / 3 * qnorm((1 - p1_2_50) / 2)),
                     MCpk_p2_2_50 = (-1 / 3 * qnorm((1 - p2_2_50) / 2)),
                     MCpm_p3_2_50 = (-1 / 3 * qnorm((1 - p3_2_50) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_2_50 <- data.frame(Estadistica = c("p1_2_50", "p2_2_50", "p3_2_50", "MCp_p1_2_50", "MCpk_p2_2_50", "MCpm_p3_2_50", "NMCp", "NMCpm"),
                              Percentil_2.5 = numeric(8),
                              Percentil_97.5 = numeric(8),
                              Medias = numeric(8),
                              Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_2_50$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_2_50[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_2_50[[estadistica]], probs = 0.975)

  intervalos_2_50$Percentil_2.5[i] <- percentil_2.5
  intervalos_2_50$Percentil_97.5[i] <- percentil_97.5
  intervalos_2_50$Medias[i] <- mean(datos_2_50[[estadistica]])
  intervalos_2_50$Desvios[i] <- sd(datos_2_50[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_2_50 <- mean(datos_2_50$p1_2_50) - p1_2
sesgo_p2_2_50 <- mean(datos_2_50$p2_2_50) - p2_2
sesgo_p3_2_50 <- mean(datos_2_50$p3_2_50) - p3_2
sesgo_MCp_p1_2_50 <- mean(datos_2_50$MCp_p1_2_50) - MCp_p1_2
sesgo_MCpk_p2_2_50 <- mean(datos_2_50$MCpk_p2_2_50) - MCpk_p2_2
sesgo_MCpm_p3_2_50 <- mean(datos_2_50$MCpm_p3_2_50) - MCpm_p3_2
sesgo_NMCp <- mean(datos_2_50$NMCp) - NMCp_2
sesgo_NMCpm <- mean(datos_2_50$NMCpm) - NMCpm_2

sesgo_p1_2_50_r <- (mean(datos_2_50$p1_2_50) - p1_2)/p1_2
sesgo_p2_2_50_r <- (mean(datos_2_50$p2_2_50) - p2_2)/p2_2
sesgo_p3_2_50_r <- (mean(datos_2_50$p3_2_50) - p3_2)/p3_2
sesgo_MCp_p1_2_50_r <- (mean(datos_2_50$MCp_p1_2_50) - MCp_p1_2)/MCp_p1_2
sesgo_MCpk_p2_2_50_r <- (mean(datos_2_50$MCpk_p2_2_50) - MCpk_p2_2)/MCpk_p2_2
sesgo_MCpm_p3_2_50_r <- (mean(datos_2_50$MCpm_p3_2_50) - MCpm_p3_2)/MCpm_p3_2
sesgo_NMCp_r <- (mean(datos_2_50$NMCp) - NMCp_2)/NMCp_2
sesgo_NMCpm_r <- (mean(datos_2_50$NMCpm) - NMCpm_2)/NMCpm_2

# Concateno los sesgos en un vector
sesgos_2_50 <- c(sesgo_p1_2_50, sesgo_p2_2_50, sesgo_p3_2_50, sesgo_MCp_p1_2_50, sesgo_MCpk_p2_2_50, sesgo_MCpm_p3_2_50, sesgo_NMCp, sesgo_NMCpm)
sesgos_2_50_r <- c(sesgo_p1_2_50_r, sesgo_p2_2_50_r, sesgo_p3_2_50_r, sesgo_MCp_p1_2_50_r, sesgo_MCpk_p2_2_50_r, sesgo_MCpm_p3_2_50_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_2_50$Sesgo <- sesgos_2_50
intervalos_2_50$ECM <- (intervalos_2_50$Desvios)^2 + (intervalos_2_50$Sesgo)^2
intervalos_2_50$SesgoRelativo <- sesgos_2_50_r
# Mostrar los resultados
print(intervalos_2_50)

#################
### p=2, n=100 ###
#################

# Tamaño de muestra
tamano_muestra <- 100

#Cálculo de las medias y sigmas muestrales
medias_2_100 <- matrix(0, nrow = cantidad_muestras, ncol = 2)
sigmas_2_100 <- array(0, dim = c(2, 2, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_2, Sigma = sigma_2)
  medias_2_100[i, ] <- colMeans(muestra)
  sigmas_2_100[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_2_100 <- data.frame(p1_2_100 = numeric(0), p2_2_100 = numeric(0), p3_2_100 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_2_100[, , i]
  p1_2_100 <- pmvnorm(LIE_2, LSE_2, mean = t_2, sigma = sigma_actual)
  p2_2_100 <- pmvnorm(LIE_2, LSE_2, mean = medias_2_100[i, ], sigma = sigma_actual)
  sigmat_2_100 <- sigma_actual + (medias_2_100[i, ] - t_2) %*% t(medias_2_100[i, ] - t_2)
  p3_2_100 <- pmvnorm(LIE_2, LSE_2, mean = t_2, sigma = sigmat_2_100)
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_2[1] - t_2[1]) / sqrt(sigma_actual[1, 1])), ((LSE_2[2] - t_2[2]) / sqrt(sigma_actual[2, 2])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 2))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_2 - medias_2_100[i, ]) %*% solve(sigma_actual) %*% (t_2 - medias_2_100[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_2_100[i, ] <- c(p1_2_100, p2_2_100, p3_2_100, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_2_100 <- mutate(datos_2_100,
                      MCp_p1_2_100 = (-1 / 3 * qnorm((1 - p1_2_100) / 2)),
                      MCpk_p2_2_100 = (-1 / 3 * qnorm((1 - p2_2_100) / 2)),
                      MCpm_p3_2_100 = (-1 / 3 * qnorm((1 - p3_2_100) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_2_100 <- data.frame(Estadistica = c("p1_2_100", "p2_2_100", "p3_2_100", "MCp_p1_2_100", "MCpk_p2_2_100", "MCpm_p3_2_100", "NMCp", "NMCpm"),
                               Percentil_2.5 = numeric(8),
                               Percentil_97.5 = numeric(8),
                               Medias = numeric(8),
                               Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_2_100$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_2_100[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_2_100[[estadistica]], probs = 0.975)
  
  intervalos_2_100$Percentil_2.5[i] <- percentil_2.5
  intervalos_2_100$Percentil_97.5[i] <- percentil_97.5
  intervalos_2_100$Medias[i] <- mean(datos_2_100[[estadistica]])
  intervalos_2_100$Desvios[i] <- sd(datos_2_100[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_2_100 <- mean(datos_2_100$p1_2_100) - p1_2
sesgo_p2_2_100 <- mean(datos_2_100$p2_2_100) - p2_2
sesgo_p3_2_100 <- mean(datos_2_100$p3_2_100) - p3_2
sesgo_MCp_p1_2_100 <- mean(datos_2_100$MCp_p1_2_100) - MCp_p1_2
sesgo_MCpk_p2_2_100 <- mean(datos_2_100$MCpk_p2_2_100) - MCpk_p2_2
sesgo_MCpm_p3_2_100 <- mean(datos_2_100$MCpm_p3_2_100) - MCpm_p3_2
sesgo_NMCp <- mean(datos_2_100$NMCp) - NMCp_2
sesgo_NMCpm <- mean(datos_2_100$NMCpm) - NMCpm_2

sesgo_p1_2_100_r <- (mean(datos_2_100$p1_2_100) - p1_2)/p1_2
sesgo_p2_2_100_r <- (mean(datos_2_100$p2_2_100) - p2_2)/p2_2
sesgo_p3_2_100_r <- (mean(datos_2_100$p3_2_100) - p3_2)/p3_2
sesgo_MCp_p1_2_100_r <- (mean(datos_2_100$MCp_p1_2_100) - MCp_p1_2)/MCp_p1_2
sesgo_MCpk_p2_2_100_r <- (mean(datos_2_100$MCpk_p2_2_100) - MCpk_p2_2)/MCpk_p2_2
sesgo_MCpm_p3_2_100_r <- (mean(datos_2_100$MCpm_p3_2_100) - MCpm_p3_2)/MCpm_p3_2
sesgo_NMCp_r <- (mean(datos_2_100$NMCp) - NMCp_2)/NMCp_2
sesgo_NMCpm_r <- (mean(datos_2_100$NMCpm) - NMCpm_2)/NMCpm_2

# Concateno los sesgos en un vector
sesgos_2_100 <- c(sesgo_p1_2_100, sesgo_p2_2_100, sesgo_p3_2_100, sesgo_MCp_p1_2_100, sesgo_MCpk_p2_2_100, sesgo_MCpm_p3_2_100, sesgo_NMCp, sesgo_NMCpm)
sesgos_2_100_r <- c(sesgo_p1_2_100_r, sesgo_p2_2_100_r, sesgo_p3_2_100_r, sesgo_MCp_p1_2_100_r, sesgo_MCpk_p2_2_100_r, sesgo_MCpm_p3_2_100_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_2_100$Sesgo <- sesgos_2_100
intervalos_2_100$ECM <- (intervalos_2_100$Desvios)^2 + (intervalos_2_100$Sesgo)^2
intervalos_2_100$SesgoRelativo <- sesgos_2_100_r
# Mostrar los resultados
print(intervalos_2_100)

#################
### p=2, n=200 ###
#################

# Tamaño de muestra
tamano_muestra <- 200

#Cálculo de las medias y sigmas muestrales
medias_2_200 <- matrix(0, nrow = cantidad_muestras, ncol = 2)
sigmas_2_200 <- array(0, dim = c(2, 2, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_2, Sigma = sigma_2)
  medias_2_200[i, ] <- colMeans(muestra)
  sigmas_2_200[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_2_200 <- data.frame(p1_2_200 = numeric(0), p2_2_200 = numeric(0), p3_2_200 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_2_200[, , i]
  p1_2_200 <- pmvnorm(LIE_2, LSE_2, mean = t_2, sigma = sigma_actual)
  p2_2_200 <- pmvnorm(LIE_2, LSE_2, mean = medias_2_200[i, ], sigma = sigma_actual)
  sigmat_2_200 <- sigma_actual + (medias_2_200[i, ] - t_2) %*% t(medias_2_200[i, ] - t_2)
  p3_2_200 <- pmvnorm(LIE_2, LSE_2, mean = t_2, sigma = sigmat_2_200)
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_2[1] - t_2[1]) / sqrt(sigma_actual[1, 1])), ((LSE_2[2] - t_2[2]) / sqrt(sigma_actual[2, 2])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 2))
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_2 - medias_2_200[i, ]) %*% solve(sigma_actual) %*% (t_2 - medias_2_200[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_2_200[i, ] <- c(p1_2_200, p2_2_200, p3_2_200, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_2_200 <- mutate(datos_2_200,
                      MCp_p1_2_200 = (-1 / 3 * qnorm((1 - p1_2_200) / 2)),
                      MCpk_p2_2_200 = (-1 / 3 * qnorm((1 - p2_2_200) / 2)),
                      MCpm_p3_2_200 = (-1 / 3 * qnorm((1 - p3_2_200) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_2_200 <- data.frame(Estadistica = c("p1_2_200", "p2_2_200", "p3_2_200", "MCp_p1_2_200", "MCpk_p2_2_200", "MCpm_p3_2_200", "NMCp", "NMCpm"),
                               Percentil_2.5 = numeric(8),
                               Percentil_97.5 = numeric(8),
                               Medias = numeric(8),
                               Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_2_200$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_2_200[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_2_200[[estadistica]], probs = 0.975)
  
  intervalos_2_200$Percentil_2.5[i] <- percentil_2.5
  intervalos_2_200$Percentil_97.5[i] <- percentil_97.5
  intervalos_2_200$Medias[i] <- mean(datos_2_200[[estadistica]])
  intervalos_2_200$Desvios[i] <- sd(datos_2_200[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_2_200 <- mean(datos_2_200$p1_2_200) - p1_2
sesgo_p2_2_200 <- mean(datos_2_200$p2_2_200) - p2_2
sesgo_p3_2_200 <- mean(datos_2_200$p3_2_200) - p3_2
sesgo_MCp_p1_2_200 <- mean(datos_2_200$MCp_p1_2_200) - MCp_p1_2
sesgo_MCpk_p2_2_200 <- mean(datos_2_200$MCpk_p2_2_200) - MCpk_p2_2
sesgo_MCpm_p3_2_200 <- mean(datos_2_200$MCpm_p3_2_200) - MCpm_p3_2
sesgo_NMCp <- mean(datos_2_200$NMCp) - NMCp_2
sesgo_NMCpm <- mean(datos_2_200$NMCpm) - NMCpm_2

sesgo_p1_2_200_r <- (mean(datos_2_200$p1_2_200) - p1_2)/p1_2
sesgo_p2_2_200_r <- (mean(datos_2_200$p2_2_200) - p2_2)/p2_2
sesgo_p3_2_200_r <- (mean(datos_2_200$p3_2_200) - p3_2)/p3_2
sesgo_MCp_p1_2_200_r <- (mean(datos_2_200$MCp_p1_2_200) - MCp_p1_2)/MCp_p1_2
sesgo_MCpk_p2_2_200_r <- (mean(datos_2_200$MCpk_p2_2_200) - MCpk_p2_2)/MCpk_p2_2
sesgo_MCpm_p3_2_200_r <- (mean(datos_2_200$MCpm_p3_2_200) - MCpm_p3_2)/MCpm_p3_2
sesgo_NMCp_r <- (mean(datos_2_200$NMCp) - NMCp_2)/NMCp_2
sesgo_NMCpm_r <- (mean(datos_2_200$NMCpm) - NMCpm_2)/NMCpm_2

# Concateno los sesgos en un vector
sesgos_2_200 <- c(sesgo_p1_2_200, sesgo_p2_2_200, sesgo_p3_2_200, sesgo_MCp_p1_2_200, sesgo_MCpk_p2_2_200, sesgo_MCpm_p3_2_200, sesgo_NMCp, sesgo_NMCpm)
sesgos_2_200_r <- c(sesgo_p1_2_200_r, sesgo_p2_2_200_r, sesgo_p3_2_200_r, sesgo_MCp_p1_2_200_r, sesgo_MCpk_p2_2_200_r, sesgo_MCpm_p3_2_200_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_2_200$Sesgo <- sesgos_2_200
intervalos_2_200$ECM <- (intervalos_2_200$Desvios)^2 + (intervalos_2_200$Sesgo)^2
intervalos_2_200$SesgoRelativo <- sesgos_2_200_r
# Mostrar los resultados
print(intervalos_2_200)

#################
### p=3, n=50 ###
#################
datos_3 <- datos_multi[, c(4, 4, 4)]
#Parámetros
LIE_3 <- c(3.00, 3.00, 3.00)
LSE_3 <- c(17.00, 17.00, 17.00)
t_3 <- c(10, 10, 10)
medias_3 <- colMeans(datos_3)
sigma_3 <- matrix(c(var(datos_3[,1]),
                    0.50*sqrt(var(datos_3[,1])*var(datos_3[,2])),
                    0.50*sqrt(var(datos_3[,1])*var(datos_3[,3])),
                    0.50*sqrt(var(datos_3[,1])*var(datos_3[,2])),
                    var(datos_3[,2]),
                    0.50*sqrt(var(datos_3[,2])*var(datos_3[,3])),
                    0.50*sqrt(var(datos_3[,1])*var(datos_3[,3])),
                    0.50*sqrt(var(datos_3[,2])*var(datos_3[,3])),
                    var(datos_3[,3])), 3, 3)
# Tamaño de muestra
tamano_muestra <- 50
#Cálculo de las probabilidades 
p1_3 <- pmvnorm(LIE_3, LSE_3, mean = t_3, sigma = sigma_3)
p2_3 <- pmvnorm(LIE_3, LSE_3, mean = medias_3, sigma = sigma_3)
sigmat_3 <- sigma_3 + t(t(medias_3 - t_3)) %*% t(medias_3 - t_3)
p3_3 <- pmvnorm(LIE_3, LSE_3, mean = t_3, sigma = sigmat_3)
#Cálculo de los índices
MCp_p1_3 <- (-1/3 * qnorm((1 - p1_3)/2))
MCpk_p2_3 <- (-1/3 * qnorm((1 - p2_3)/2))
MCpm_p3_3 <- (-1/3 * qnorm((1 - p3_3)/2))

# Cálculo de los valores reales de NMCp y NMCpm
cprima_3 <- min(((LSE_3[1] - t_3[1]) / sqrt(sigma_3[1, 1])), ((LSE_3[2] - t_3[2]) / sqrt(sigma_3[2, 2])), ((LSE_3[3] - t_3[3]) / sqrt(sigma_3[3, 3])))
NMCp_3 <- cprima_3 / sqrt(qchisq(1 - 0.0027, 3))

# Cálculo de la distancia D
D_3 <- (1 + (t(medias_3 - t_3) %*% solve(sigma_3) %*% (medias_3 - t_3)))^0.5
NMCpm_3 <- NMCp_3 / D_3  # D_real es 1

#Cálculo de las medias y sigmas muestrales
medias_3_50 <- matrix(0, nrow = cantidad_muestras, ncol = 3)
sigmas_3_50 <- array(0, dim = c(3, 3, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_3, Sigma = sigma_3)
  medias_3_50[i, ] <- colMeans(muestra)
  sigmas_3_50[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_3_50 <- data.frame(p1_3_50 = numeric(0), p2_3_50 = numeric(0), p3_3_50 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_3_50[, , i]
  p1_3_50 <- pmvnorm(LIE_3, LSE_3, mean = t_3, sigma = sigma_actual)
  p2_3_50 <- pmvnorm(LIE_3, LSE_3, mean = medias_3_50[i, ], sigma = sigma_actual)
  sigmat_3_50 <- sigma_actual + (medias_3_50[i, ] - t_3) %*% t(medias_3_50[i, ] - t_3)
  p3_3_50 <- pmvnorm(LIE_3, LSE_3, mean = t_3, sigma = sigmat_3_50)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_3[1] - t_3[1]) / sqrt(sigma_actual[1, 1])), ((LSE_3[2] - t_3[2]) / sqrt(sigma_actual[2, 2])), ((LSE_3[3] - t_3[3]) / sqrt(sigma_actual[3, 3])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 3))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_3 - medias_3_50[i, ]) %*% solve(sigma_actual) %*% (t_3 - medias_3_50[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_3_50[i, ] <- c(p1_3_50, p2_3_50, p3_3_50, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_3_50 <- mutate(datos_3_50,
                     MCp_p1_3_50 = (-1 / 3 * qnorm((1 - p1_3_50) / 2)),
                     MCpk_p2_3_50 = (-1 / 3 * qnorm((1 - p2_3_50) / 2)),
                     MCpm_p3_3_50 = (-1 / 3 * qnorm((1 - p3_3_50) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_3_50 <- data.frame(Estadistica = c("p1_3_50", "p2_3_50", "p3_3_50", "MCp_p1_3_50", "MCpk_p2_3_50", "MCpm_p3_3_50", "NMCp", "NMCpm"),
                              Percentil_2.5 = numeric(8),
                              Percentil_97.5 = numeric(8),
                              Medias = numeric(8),
                              Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_3_50$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_3_50[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_3_50[[estadistica]], probs = 0.975)
  
  intervalos_3_50$Percentil_2.5[i] <- percentil_2.5
  intervalos_3_50$Percentil_97.5[i] <- percentil_97.5
  intervalos_3_50$Medias[i] <- mean(datos_3_50[[estadistica]])
  intervalos_3_50$Desvios[i] <- sd(datos_3_50[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_3_50 <- mean(datos_3_50$p1_3_50) - p1_3
sesgo_p2_3_50 <- mean(datos_3_50$p2_3_50) - p2_3
sesgo_p3_3_50 <- mean(datos_3_50$p3_3_50) - p3_3
sesgo_MCp_p1_3_50 <- mean(datos_3_50$MCp_p1_3_50) - MCp_p1_3
sesgo_MCpk_p2_3_50 <- mean(datos_3_50$MCpk_p2_3_50) - MCpk_p2_3
sesgo_MCpm_p3_3_50 <- mean(datos_3_50$MCpm_p3_3_50) - MCpm_p3_3
sesgo_NMCp <- mean(datos_3_50$NMCp) - NMCp_3
sesgo_NMCpm <- mean(datos_3_50$NMCpm) - NMCpm_3

sesgo_p1_3_50_r <- (mean(datos_3_50$p1_3_50) - p1_3)/p1_3
sesgo_p2_3_50_r <- (mean(datos_3_50$p2_3_50) - p2_3)/p2_3
sesgo_p3_3_50_r <- (mean(datos_3_50$p3_3_50) - p3_3)/p3_3
sesgo_MCp_p1_3_50_r <- (mean(datos_3_50$MCp_p1_3_50) - MCp_p1_3)/MCp_p1_3
sesgo_MCpk_p2_3_50_r <- (mean(datos_3_50$MCpk_p2_3_50) - MCpk_p2_3)/MCpk_p2_3
sesgo_MCpm_p3_3_50_r <- (mean(datos_3_50$MCpm_p3_3_50) - MCpm_p3_3)/MCpm_p3_3
sesgo_NMCp_r <- (mean(datos_3_50$NMCp) - NMCp_3)/NMCp_3
sesgo_NMCpm_r <- (mean(datos_3_50$NMCpm) - NMCpm_3)/NMCpm_3

# Concateno los sesgos en un vector
sesgos_3_50 <- c(sesgo_p1_3_50, sesgo_p2_3_50, sesgo_p3_3_50, sesgo_MCp_p1_3_50, sesgo_MCpk_p2_3_50, sesgo_MCpm_p3_3_50, sesgo_NMCp, sesgo_NMCpm)
sesgos_3_50_r <- c(sesgo_p1_3_50_r, sesgo_p2_3_50_r, sesgo_p3_3_50_r, sesgo_MCp_p1_3_50_r, sesgo_MCpk_p2_3_50_r, sesgo_MCpm_p3_3_50_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_3_50$Sesgo <- sesgos_3_50
intervalos_3_50$ECM <- (intervalos_3_50$Desvios)^2 + (intervalos_3_50$Sesgo)^2
intervalos_3_50$SesgoRelativo <- sesgos_3_50_r
# Mostrar los resultados
print(intervalos_3_50)

#################
### p=3, n=100 ###
#################

# Tamaño de muestra
tamano_muestra <- 100

#Cálculo de las medias y sigmas muestrales
medias_3_100 <- matrix(0, nrow = cantidad_muestras, ncol = 3)
sigmas_3_100 <- array(0, dim = c(3, 3, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_3, Sigma = sigma_3)
  medias_3_100[i, ] <- colMeans(muestra)
  sigmas_3_100[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_3_100 <- data.frame(p1_3_100 = numeric(0), p2_3_100 = numeric(0), p3_3_100 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_3_100[, , i]
  p1_3_100 <- pmvnorm(LIE_3, LSE_3, mean = t_3, sigma = sigma_actual)
  p2_3_100 <- pmvnorm(LIE_3, LSE_3, mean = medias_3_100[i, ], sigma = sigma_actual)
  sigmat_3_100 <- sigma_actual + (medias_3_100[i, ] - t_3) %*% t(medias_3_100[i, ] - t_3)
  p3_3_100 <- pmvnorm(LIE_3, LSE_3, mean = t_3, sigma = sigmat_3_100)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_3[1] - t_3[1]) / sqrt(sigma_actual[1, 1])), ((LSE_3[2] - t_3[2]) / sqrt(sigma_actual[2, 2])), ((LSE_3[3] - t_3[3]) / sqrt(sigma_actual[3, 3])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 3))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_3 - medias_3_100[i, ]) %*% solve(sigma_actual) %*% (t_3 - medias_3_100[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_3_100[i, ] <- c(p1_3_100, p2_3_100, p3_3_100, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_3_100 <- mutate(datos_3_100,
                      MCp_p1_3_100 = (-1 / 3 * qnorm((1 - p1_3_100) / 2)),
                      MCpk_p2_3_100 = (-1 / 3 * qnorm((1 - p2_3_100) / 2)),
                      MCpm_p3_3_100 = (-1 / 3 * qnorm((1 - p3_3_100) / 2)))


# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_3_100 <- data.frame(Estadistica = c("p1_3_100", "p2_3_100", "p3_3_100", "MCp_p1_3_100", "MCpk_p2_3_100", "MCpm_p3_3_100", "NMCp", "NMCpm"),
                               Percentil_2.5 = numeric(8),
                               Percentil_97.5 = numeric(8),
                               Medias = numeric(8),
                               Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_3_100$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_3_100[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_3_100[[estadistica]], probs = 0.975)
  
  intervalos_3_100$Percentil_2.5[i] <- percentil_2.5
  intervalos_3_100$Percentil_97.5[i] <- percentil_97.5
  intervalos_3_100$Medias[i] <- mean(datos_3_100[[estadistica]])
  intervalos_3_100$Desvios[i] <- sd(datos_3_100[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_3_100 <- mean(datos_3_100$p1_3_100) - p1_3
sesgo_p2_3_100 <- mean(datos_3_100$p2_3_100) - p2_3
sesgo_p3_3_100 <- mean(datos_3_100$p3_3_100) - p3_3
sesgo_MCp_p1_3_100 <- mean(datos_3_100$MCp_p1_3_100) - MCp_p1_3
sesgo_MCpk_p2_3_100 <- mean(datos_3_100$MCpk_p2_3_100) - MCpk_p2_3
sesgo_MCpm_p3_3_100 <- mean(datos_3_100$MCpm_p3_3_100) - MCpm_p3_3
sesgo_NMCp <- mean(datos_3_100$NMCp) - NMCp_3
sesgo_NMCpm <- mean(datos_3_100$NMCpm) - NMCpm_3

sesgo_p1_3_100_r <- (mean(datos_3_100$p1_3_100) - p1_3)/p1_3
sesgo_p2_3_100_r <- (mean(datos_3_100$p2_3_100) - p2_3)/p2_3
sesgo_p3_3_100_r <- (mean(datos_3_100$p3_3_100) - p3_3)/p3_3
sesgo_MCp_p1_3_100_r <- (mean(datos_3_100$MCp_p1_3_100) - MCp_p1_3)/MCp_p1_3
sesgo_MCpk_p2_3_100_r <- (mean(datos_3_100$MCpk_p2_3_100) - MCpk_p2_3)/MCpk_p2_3
sesgo_MCpm_p3_3_100_r <- (mean(datos_3_100$MCpm_p3_3_100) - MCpm_p3_3)/MCpm_p3_3
sesgo_NMCp_r <- (mean(datos_3_100$NMCp) - NMCp_3)/NMCp_3
sesgo_NMCpm_r <- (mean(datos_3_100$NMCpm) - NMCpm_3)/NMCpm_3

# Concateno los sesgos en un vector
sesgos_3_100 <- c(sesgo_p1_3_100, sesgo_p2_3_100, sesgo_p3_3_100, sesgo_MCp_p1_3_100, sesgo_MCpk_p2_3_100, sesgo_MCpm_p3_3_100, sesgo_NMCp, sesgo_NMCpm)
sesgos_3_100_r <- c(sesgo_p1_3_100_r, sesgo_p2_3_100_r, sesgo_p3_3_100_r, sesgo_MCp_p1_3_100_r, sesgo_MCpk_p2_3_100_r, sesgo_MCpm_p3_3_100_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_3_100$Sesgo <- sesgos_3_100
intervalos_3_100$ECM <- (intervalos_3_100$Desvios)^2 + (intervalos_3_100$Sesgo)^2
intervalos_3_100$SesgoRelativo <- sesgos_3_100_r
# Mostrar los resultados
print(intervalos_3_100)

#################
### p=3, n=200 ###
#################

tamano_muestra <- 200

#Cálculo de las medias y sigmas muestrales
medias_3_200 <- matrix(0, nrow = cantidad_muestras, ncol = 3)
sigmas_3_200 <- array(0, dim = c(3, 3, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_3, Sigma = sigma_3)
  medias_3_200[i, ] <- colMeans(muestra)
  sigmas_3_200[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_3_200 <- data.frame(p1_3_200 = numeric(0), p2_3_200 = numeric(0), p3_3_200 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_3_200[, , i]
  p1_3_200 <- pmvnorm(LIE_3, LSE_3, mean = t_3, sigma = sigma_actual)
  p2_3_200 <- pmvnorm(LIE_3, LSE_3, mean = medias_3_200[i, ], sigma = sigma_actual)
  sigmat_3_200 <- sigma_actual + (medias_3_200[i, ] - t_3) %*% t(medias_3_200[i, ] - t_3)
  p3_3_200 <- pmvnorm(LIE_3, LSE_3, mean = t_3, sigma = sigmat_3_200)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_3[1] - t_3[1]) / sqrt(sigma_actual[1, 1])), ((LSE_3[2] - t_3[2]) / sqrt(sigma_actual[2, 2])), ((LSE_3[3] - t_3[3]) / sqrt(sigma_actual[3, 3])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 3))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_3 - medias_3_200[i, ]) %*% solve(sigma_actual) %*% (t_3 - medias_3_200[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_3_200[i, ] <- c(p1_3_200, p2_3_200, p3_3_200, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_3_200 <- mutate(datos_3_200,
                      MCp_p1_3_200 = (-1 / 3 * qnorm((1 - p1_3_200) / 2)),
                      MCpk_p2_3_200 = (-1 / 3 * qnorm((1 - p2_3_200) / 2)),
                      MCpm_p3_3_200 = (-1 / 3 * qnorm((1 - p3_3_200) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_3_200 <- data.frame(Estadistica = c("p1_3_200", "p2_3_200", "p3_3_200", "MCp_p1_3_200", "MCpk_p2_3_200", "MCpm_p3_3_200", "NMCp", "NMCpm"),
                               Percentil_2.5 = numeric(8),
                               Percentil_97.5 = numeric(8),
                               Medias = numeric(8),
                               Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_3_200$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_3_200[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_3_200[[estadistica]], probs = 0.975)
  
  intervalos_3_200$Percentil_2.5[i] <- percentil_2.5
  intervalos_3_200$Percentil_97.5[i] <- percentil_97.5
  intervalos_3_200$Medias[i] <- mean(datos_3_200[[estadistica]])
  intervalos_3_200$Desvios[i] <- sd(datos_3_200[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_3_200 <- mean(datos_3_200$p1_3_200) - p1_3
sesgo_p2_3_200 <- mean(datos_3_200$p2_3_200) - p2_3
sesgo_p3_3_200 <- mean(datos_3_200$p3_3_200) - p3_3
sesgo_MCp_p1_3_200 <- mean(datos_3_200$MCp_p1_3_200) - MCp_p1_3
sesgo_MCpk_p2_3_200 <- mean(datos_3_200$MCpk_p2_3_200) - MCpk_p2_3
sesgo_MCpm_p3_3_200 <- mean(datos_3_200$MCpm_p3_3_200) - MCpm_p3_3
sesgo_NMCp <- mean(datos_3_200$NMCp) - NMCp_3
sesgo_NMCpm <- mean(datos_3_200$NMCpm) - NMCpm_3

sesgo_p1_3_200_r <- (mean(datos_3_200$p1_3_200) - p1_3)/p1_3
sesgo_p2_3_200_r <- (mean(datos_3_200$p2_3_200) - p2_3)/p2_3
sesgo_p3_3_200_r <- (mean(datos_3_200$p3_3_200) - p3_3)/p3_3
sesgo_MCp_p1_3_200_r <- (mean(datos_3_200$MCp_p1_3_200) - MCp_p1_3)/MCp_p1_3
sesgo_MCpk_p2_3_200_r <- (mean(datos_3_200$MCpk_p2_3_200) - MCpk_p2_3)/MCpk_p2_3
sesgo_MCpm_p3_3_200_r <- (mean(datos_3_200$MCpm_p3_3_200) - MCpm_p3_3)/MCpm_p3_3
sesgo_NMCp_r <- (mean(datos_3_200$NMCp) - NMCp_3)/NMCp_3
sesgo_NMCpm_r <- (mean(datos_3_200$NMCpm) - NMCpm_3)/NMCpm_3

# Concateno los sesgos en un vector
sesgos_3_200 <- c(sesgo_p1_3_200, sesgo_p2_3_200, sesgo_p3_3_200, sesgo_MCp_p1_3_200, sesgo_MCpk_p2_3_200, sesgo_MCpm_p3_3_200, sesgo_NMCp, sesgo_NMCpm)
sesgos_3_200_r <- c(sesgo_p1_3_200_r, sesgo_p2_3_200_r, sesgo_p3_3_200_r, sesgo_MCp_p1_3_200_r, sesgo_MCpk_p2_3_200_r, sesgo_MCpm_p3_3_200_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_3_200$Sesgo <- sesgos_3_200
intervalos_3_200$ECM <- (intervalos_3_200$Desvios)^2 + (intervalos_3_200$Sesgo)^2
intervalos_3_200$SesgoRelativo <- sesgos_3_200_r
# Mostrar los resultados
print(intervalos_3_200)


#################
### p=4, n= 50 ###
#################
datos_4 <- datos_multi[, c(4, 4, 4, 4)]
#Parámetros
LIE_4 <- c(3.00, 3.00, 3.00, 3.00)
LSE_4 <- c(17.00, 17.00, 17.00, 17.00)
t_4 <- c(10, 10, 10, 10)
medias_4 <- colMeans(datos_4)
sigma_4 <- matrix(c(var(datos_4[,1]),
                    0.50*sqrt(var(datos_4[,1])*var(datos_4[,2])),
                    0.50*sqrt(var(datos_4[,1])*var(datos_4[,3])),
                    0.50*sqrt(var(datos_4[,1])*var(datos_4[,4])),
                    0.50*sqrt(var(datos_4[,1])*var(datos_4[,2])),
                    var(datos_4[,2]),
                    0.50*sqrt(var(datos_4[,2])*var(datos_4[,3])),
                    0.50*sqrt(var(datos_4[,2])*var(datos_4[,4])),
                    0.50*sqrt(var(datos_4[,1])*var(datos_4[,3])),
                    0.50*sqrt(var(datos_4[,2])*var(datos_4[,3])),
                    var(datos_4[,3]),
                    0.50*sqrt(var(datos_4[,3])*var(datos_4[,4])),
                    0.50*sqrt(var(datos_4[,1])*var(datos_4[,4])),
                    0.50*sqrt(var(datos_4[,2])*var(datos_4[,4])),
                    0.50*sqrt(var(datos_4[,3])*var(datos_4[,4])),
                    var(datos_4[,4])), 4, 4)
# Tamaño de muestra
tamano_muestra <- 50
#Cálculo de las probabilidades 
p1_4 <- pmvnorm(LIE_4, LSE_4, mean = t_4, sigma = sigma_4)
p2_4 <- pmvnorm(LIE_4, LSE_4, mean = medias_4, sigma = sigma_4)
sigmat_4 <- sigma_4 + t(t(medias_4 - t_4)) %*% t(medias_4 - t_4)
p3_4 <- pmvnorm(LIE_4, LSE_4, mean = t_4, sigma = sigmat_4)
#Cálculo de los índices
MCp_p1_4 <- (-1/3 * qnorm((1 - p1_4)/2))
MCpk_p2_4 <- (-1/3 * qnorm((1 - p2_4)/2))
MCpm_p3_4 <- (-1/3 * qnorm((1 - p3_4)/2))

# Cálculo de los valores reales de NMCp y NMCpm
cprima_4 <- min(((LSE_4[1] - t_4[1]) / sqrt(sigma_4[1, 1])), ((LSE_4[2] - t_4[2]) / sqrt(sigma_4[2, 2])), ((LSE_4[3] - t_4[3]) / sqrt(sigma_4[3, 3])), ((LSE_4[4] - t_4[4]) / sqrt(sigma_4[4, 4])))
NMCp_4 <- cprima_4 / sqrt(qchisq(1 - 0.0027, 4))

# Cálculo de la distancia D
D_4 <- (1 + (t(medias_4 - t_4) %*% solve(sigma_4) %*% (medias_4 - t_4)))^0.5
NMCpm_4 <- NMCp_4 / D_4  # D_real es 1

#Cálculo de las medias y sigmas muestrales
medias_4_50 <- matrix(0, nrow = cantidad_muestras, ncol = 4)
sigmas_4_50 <- array(0, dim = c(4, 4, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_4, Sigma = sigma_4)
  medias_4_50[i, ] <- colMeans(muestra)
  sigmas_4_50[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_4_50 <- data.frame(p1_4_50 = numeric(0), p2_4_50 = numeric(0), p3_4_50 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_4_50[, , i]
  p1_4_50 <- pmvnorm(LIE_4, LSE_4, mean = t_4, sigma = sigma_actual)
  p2_4_50 <- pmvnorm(LIE_4, LSE_4, mean = medias_4_50[i, ], sigma = sigma_actual)
  sigmat_4_50 <- sigma_actual + (medias_4_50[i, ] - t_4) %*% t(medias_4_50[i, ] - t_4)
  p3_4_50 <- pmvnorm(LIE_4, LSE_4, mean = t_4, sigma = sigmat_4_50)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_4[1] - t_4[1]) / sqrt(sigma_actual[1, 1])), ((LSE_4[2] - t_4[2]) / sqrt(sigma_actual[2, 2])), ((LSE_4[3] - t_4[3]) / sqrt(sigma_actual[3, 3])), ((LSE_4[4] - t_4[4]) / sqrt(sigma_actual[4, 4])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 4))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_4 - medias_4_50[i, ]) %*% solve(sigma_actual) %*% (t_4 - medias_4_50[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_4_50[i, ] <- c(p1_4_50, p2_4_50, p3_4_50, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_4_50 <- mutate(datos_4_50,
                     MCp_p1_4_50 = (-1 / 3 * qnorm((1 - p1_4_50) / 2)),
                     MCpk_p2_4_50 = (-1 / 3 * qnorm((1 - p2_4_50) / 2)),
                     MCpm_p3_4_50 = (-1 / 3 * qnorm((1 - p3_4_50) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_4_50 <- data.frame(Estadistica = c("p1_4_50", "p2_4_50", "p3_4_50", "MCp_p1_4_50", "MCpk_p2_4_50", "MCpm_p3_4_50", "NMCp", "NMCpm"),
                              Percentil_2.5 = numeric(8),
                              Percentil_97.5 = numeric(8),
                              Medias = numeric(8),
                              Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_4_50$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_4_50[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_4_50[[estadistica]], probs = 0.975)
  
  intervalos_4_50$Percentil_2.5[i] <- percentil_2.5
  intervalos_4_50$Percentil_97.5[i] <- percentil_97.5
  intervalos_4_50$Medias[i] <- mean(datos_4_50[[estadistica]])
  intervalos_4_50$Desvios[i] <- sd(datos_4_50[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_4_50 <- mean(datos_4_50$p1_4_50) - p1_4
sesgo_p2_4_50 <- mean(datos_4_50$p2_4_50) - p2_4
sesgo_p3_4_50 <- mean(datos_4_50$p3_4_50) - p3_4
sesgo_MCp_p1_4_50 <- mean(datos_4_50$MCp_p1_4_50) - MCp_p1_4
sesgo_MCpk_p2_4_50 <- mean(datos_4_50$MCpk_p2_4_50) - MCpk_p2_4
sesgo_MCpm_p3_4_50 <- mean(datos_4_50$MCpm_p3_4_50) - MCpm_p3_4
sesgo_NMCp <- mean(datos_4_50$NMCp) - NMCp_4
sesgo_NMCpm <- mean(datos_4_50$NMCpm) - NMCpm_4

sesgo_p1_4_50_r <- (mean(datos_4_50$p1_4_50) - p1_4)/p1_4
sesgo_p2_4_50_r <- (mean(datos_4_50$p2_4_50) - p2_4)/p2_4
sesgo_p3_4_50_r <- (mean(datos_4_50$p3_4_50) - p3_4)/p3_4
sesgo_MCp_p1_4_50_r <- (mean(datos_4_50$MCp_p1_4_50) - MCp_p1_4)/MCp_p1_4
sesgo_MCpk_p2_4_50_r <- (mean(datos_4_50$MCpk_p2_4_50) - MCpk_p2_4)/MCpk_p2_4
sesgo_MCpm_p3_4_50_r <- (mean(datos_4_50$MCpm_p3_4_50) - MCpm_p3_4)/MCpm_p3_4
sesgo_NMCp_r <- (mean(datos_4_50$NMCp) - NMCp_4)/NMCp_4
sesgo_NMCpm_r <- (mean(datos_4_50$NMCpm) - NMCpm_4)/NMCpm_4

# Concateno los sesgos en un vector
sesgos_4_50 <- c(sesgo_p1_4_50, sesgo_p2_4_50, sesgo_p3_4_50, sesgo_MCp_p1_4_50, sesgo_MCpk_p2_4_50, sesgo_MCpm_p3_4_50, sesgo_NMCp, sesgo_NMCpm)
sesgos_4_50_r <- c(sesgo_p1_4_50_r, sesgo_p2_4_50_r, sesgo_p3_4_50_r, sesgo_MCp_p1_4_50_r, sesgo_MCpk_p2_4_50_r, sesgo_MCpm_p3_4_50_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_4_50$Sesgo <- sesgos_4_50
intervalos_4_50$ECM <- (intervalos_4_50$Desvios)^2 + (intervalos_4_50$Sesgo)^2
intervalos_4_50$SesgoRelativo <- sesgos_4_50_r
# Mostrar los resultados
print(intervalos_4_50)

#################
### p=4, n=100 ###
#################

# Tamaño de muestra
tamano_muestra <- 100

#Cálculo de las medias y sigmas muestrales
medias_4_100 <- matrix(0, nrow = cantidad_muestras, ncol = 4)
sigmas_4_100 <- array(0, dim = c(4, 4, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_4, Sigma = sigma_4)
  medias_4_100[i, ] <- colMeans(muestra)
  sigmas_4_100[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_4_100 <- data.frame(p1_4_100 = numeric(0), p2_4_100 = numeric(0), p3_4_100 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_4_100[, , i]
  p1_4_100 <- pmvnorm(LIE_4, LSE_4, mean = t_4, sigma = sigma_actual)
  p2_4_100 <- pmvnorm(LIE_4, LSE_4, mean = medias_4_100[i, ], sigma = sigma_actual)
  sigmat_4_100 <- sigma_actual + (medias_4_100[i, ] - t_4) %*% t(medias_4_100[i, ] - t_4)
  p3_4_100 <- pmvnorm(LIE_4, LSE_4, mean = t_4, sigma = sigmat_4_100)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_4[1] - t_4[1]) / sqrt(sigma_actual[1, 1])), ((LSE_4[2] - t_4[2]) / sqrt(sigma_actual[2, 2])), ((LSE_4[4] - t_4[4]) / sqrt(sigma_actual[4, 4])), ((LSE_4[4] - t_4[4]) / sqrt(sigma_actual[4, 4])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 4))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_4 - medias_4_100[i, ]) %*% solve(sigma_actual) %*% (t_4 - medias_4_100[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_4_100[i, ] <- c(p1_4_100, p2_4_100, p3_4_100, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_4_100 <- mutate(datos_4_100,
                      MCp_p1_4_100 = (-1 / 3 * qnorm((1 - p1_4_100) / 2)),
                      MCpk_p2_4_100 = (-1 / 3 * qnorm((1 - p2_4_100) / 2)),
                      MCpm_p3_4_100 = (-1 / 3 * qnorm((1 - p3_4_100) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_4_100 <- data.frame(Estadistica = c("p1_4_100", "p2_4_100", "p3_4_100", "MCp_p1_4_100", "MCpk_p2_4_100", "MCpm_p3_4_100", "NMCp", "NMCpm"),
                               Percentil_2.5 = numeric(8),
                               Percentil_97.5 = numeric(8),
                               Medias = numeric(8),
                               Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_4_100$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_4_100[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_4_100[[estadistica]], probs = 0.975)
  
  intervalos_4_100$Percentil_2.5[i] <- percentil_2.5
  intervalos_4_100$Percentil_97.5[i] <- percentil_97.5
  intervalos_4_100$Medias[i] <- mean(datos_4_100[[estadistica]])
  intervalos_4_100$Desvios[i] <- sd(datos_4_100[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_4_100 <- mean(datos_4_100$p1_4_100) - p1_4
sesgo_p2_4_100 <- mean(datos_4_100$p2_4_100) - p2_4
sesgo_p3_4_100 <- mean(datos_4_100$p3_4_100) - p3_4
sesgo_MCp_p1_4_100 <- mean(datos_4_100$MCp_p1_4_100) - MCp_p1_4
sesgo_MCpk_p2_4_100 <- mean(datos_4_100$MCpk_p2_4_100) - MCpk_p2_4
sesgo_MCpm_p3_4_100 <- mean(datos_4_100$MCpm_p3_4_100) - MCpm_p3_4
sesgo_NMCp <- mean(datos_4_100$NMCp) - NMCp_4
sesgo_NMCpm <- mean(datos_4_100$NMCpm) - NMCpm_4

sesgo_p1_4_100_r <- (mean(datos_4_100$p1_4_100) - p1_4)/p1_4
sesgo_p2_4_100_r <- (mean(datos_4_100$p2_4_100) - p2_4)/p2_4
sesgo_p3_4_100_r <- (mean(datos_4_100$p3_4_100) - p3_4)/p3_4
sesgo_MCp_p1_4_100_r <- (mean(datos_4_100$MCp_p1_4_100) - MCp_p1_4)/MCp_p1_4
sesgo_MCpk_p2_4_100_r <- (mean(datos_4_100$MCpk_p2_4_100) - MCpk_p2_4)/MCpk_p2_4
sesgo_MCpm_p3_4_100_r <- (mean(datos_4_100$MCpm_p3_4_100) - MCpm_p3_4)/MCpm_p3_4
sesgo_NMCp_r <- (mean(datos_4_100$NMCp) - NMCp_4)/NMCp_4
sesgo_NMCpm_r <- (mean(datos_4_100$NMCpm) - NMCpm_4)/NMCpm_4

# Concateno los sesgos en un vector
sesgos_4_100 <- c(sesgo_p1_4_100, sesgo_p2_4_100, sesgo_p3_4_100, sesgo_MCp_p1_4_100, sesgo_MCpk_p2_4_100, sesgo_MCpm_p3_4_100, sesgo_NMCp, sesgo_NMCpm)
sesgos_4_100_r <- c(sesgo_p1_4_100_r, sesgo_p2_4_100_r, sesgo_p3_4_100_r, sesgo_MCp_p1_4_100_r, sesgo_MCpk_p2_4_100_r, sesgo_MCpm_p3_4_100_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_4_100$Sesgo <- sesgos_4_100
intervalos_4_100$ECM <- (intervalos_4_100$Desvios)^2 + (intervalos_4_100$Sesgo)^2
intervalos_4_100$SesgoRelativo <- sesgos_4_100_r
# Mostrar los resultados
print(intervalos_4_100)

#################
### p=4, n=200 ###
#################

# Tamaño de muestra
tamano_muestra <- 200

#Cálculo de las medias y sigmas muestrales
medias_4_200 <- matrix(0, nrow = cantidad_muestras, ncol = 4)
sigmas_4_200 <- array(0, dim = c(4, 4, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_4, Sigma = sigma_4)
  medias_4_200[i, ] <- colMeans(muestra)
  sigmas_4_200[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_4_200 <- data.frame(p1_4_200 = numeric(0), p2_4_200 = numeric(0), p3_4_200 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_4_200[, , i]
  p1_4_200 <- pmvnorm(LIE_4, LSE_4, mean = t_4, sigma = sigma_actual)
  p2_4_200 <- pmvnorm(LIE_4, LSE_4, mean = medias_4_200[i, ], sigma = sigma_actual)
  sigmat_4_200 <- sigma_actual + (medias_4_200[i, ] - t_4) %*% t(medias_4_200[i, ] - t_4)
  p3_4_200 <- pmvnorm(LIE_4, LSE_4, mean = t_4, sigma = sigmat_4_200)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_4[1] - t_4[1]) / sqrt(sigma_actual[1, 1])), ((LSE_4[2] - t_4[2]) / sqrt(sigma_actual[2, 2])), ((LSE_4[4] - t_4[4]) / sqrt(sigma_actual[4, 4])), ((LSE_4[4] - t_4[4]) / sqrt(sigma_actual[4, 4])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 4))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_4 - medias_4_200[i, ]) %*% solve(sigma_actual) %*% (t_4 - medias_4_200[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_4_200[i, ] <- c(p1_4_200, p2_4_200, p3_4_200, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_4_200 <- mutate(datos_4_200,
                      MCp_p1_4_200 = (-1 / 3 * qnorm((1 - p1_4_200) / 2)),
                      MCpk_p2_4_200 = (-1 / 3 * qnorm((1 - p2_4_200) / 2)),
                      MCpm_p3_4_200 = (-1 / 3 * qnorm((1 - p3_4_200) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_4_200 <- data.frame(Estadistica = c("p1_4_200", "p2_4_200", "p3_4_200", "MCp_p1_4_200", "MCpk_p2_4_200", "MCpm_p3_4_200", "NMCp", "NMCpm"),
                               Percentil_2.5 = numeric(8),
                               Percentil_97.5 = numeric(8),
                               Medias = numeric(8),
                               Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_4_200$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_4_200[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_4_200[[estadistica]], probs = 0.975)
  
  intervalos_4_200$Percentil_2.5[i] <- percentil_2.5
  intervalos_4_200$Percentil_97.5[i] <- percentil_97.5
  intervalos_4_200$Medias[i] <- mean(datos_4_200[[estadistica]])
  intervalos_4_200$Desvios[i] <- sd(datos_4_200[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_4_200 <- mean(datos_4_200$p1_4_200) - p1_4
sesgo_p2_4_200 <- mean(datos_4_200$p2_4_200) - p2_4
sesgo_p3_4_200 <- mean(datos_4_200$p3_4_200) - p3_4
sesgo_MCp_p1_4_200 <- mean(datos_4_200$MCp_p1_4_200) - MCp_p1_4
sesgo_MCpk_p2_4_200 <- mean(datos_4_200$MCpk_p2_4_200) - MCpk_p2_4
sesgo_MCpm_p3_4_200 <- mean(datos_4_200$MCpm_p3_4_200) - MCpm_p3_4
sesgo_NMCp <- mean(datos_4_200$NMCp) - NMCp_4
sesgo_NMCpm <- mean(datos_4_200$NMCpm) - NMCpm_4

sesgo_p1_4_200_r <- (mean(datos_4_200$p1_4_200) - p1_4)/p1_4
sesgo_p2_4_200_r <- (mean(datos_4_200$p2_4_200) - p2_4)/p2_4
sesgo_p3_4_200_r <- (mean(datos_4_200$p3_4_200) - p3_4)/p3_4
sesgo_MCp_p1_4_200_r <- (mean(datos_4_200$MCp_p1_4_200) - MCp_p1_4)/MCp_p1_4
sesgo_MCpk_p2_4_200_r <- (mean(datos_4_200$MCpk_p2_4_200) - MCpk_p2_4)/MCpk_p2_4
sesgo_MCpm_p3_4_200_r <- (mean(datos_4_200$MCpm_p3_4_200) - MCpm_p3_4)/MCpm_p3_4
sesgo_NMCp_r <- (mean(datos_4_200$NMCp) - NMCp_4)/NMCp_4
sesgo_NMCpm_r <- (mean(datos_4_200$NMCpm) - NMCpm_4)/NMCpm_4

# Concateno los sesgos en un vector
sesgos_4_200 <- c(sesgo_p1_4_200, sesgo_p2_4_200, sesgo_p3_4_200, sesgo_MCp_p1_4_200, sesgo_MCpk_p2_4_200, sesgo_MCpm_p3_4_200, sesgo_NMCp, sesgo_NMCpm)
sesgos_4_200_r <- c(sesgo_p1_4_200_r, sesgo_p2_4_200_r, sesgo_p3_4_200_r, sesgo_MCp_p1_4_200_r, sesgo_MCpk_p2_4_200_r, sesgo_MCpm_p3_4_200_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_4_200$Sesgo <- sesgos_4_200
intervalos_4_200$ECM <- (intervalos_4_200$Desvios)^2 + (intervalos_4_200$Sesgo)^2
intervalos_4_200$SesgoRelativo <- sesgos_4_200_r
# Mostrar los resultados
print(intervalos_4_200)


#################
### p=5, n= 50 ###
#################
datos_5 <- datos_multi[, c(4, 4, 4, 4, 4)]
#Parámetros
LIE_5 <- c(3.00, 3.00, 3.00, 3.00, 3.00)
LSE_5 <- c(17.00, 17.00, 17.00, 17.00, 17.00)
t_5 <- c(10, 10, 10, 10, 10)
medias_5 <- colMeans(datos_5)
sigma_5 <- matrix(c(var(datos_5[,1]),
                    0.50*sqrt(var(datos_5[,1])*var(datos_5[,2])),
                    0.50*sqrt(var(datos_5[,1])*var(datos_5[,3])),
                    0.50*sqrt(var(datos_5[,1])*var(datos_5[,4])),
                    0.50*sqrt(var(datos_5[,1])*var(datos_5[,5])),
                    0.50*sqrt(var(datos_5[,1])*var(datos_5[,2])),
                    var(datos_5[,2]),
                    0.50*sqrt(var(datos_5[,2])*var(datos_5[,3])),
                    0.50*sqrt(var(datos_5[,2])*var(datos_5[,4])),
                    0.50*sqrt(var(datos_5[,2])*var(datos_5[,5])),
                    0.50*sqrt(var(datos_5[,1])*var(datos_5[,3])),
                    0.50*sqrt(var(datos_5[,2])*var(datos_5[,3])),
                    var(datos_5[,3]),
                    0.50*sqrt(var(datos_5[,3])*var(datos_5[,4])),
                    0.50*sqrt(var(datos_5[,3])*var(datos_5[,5])),
                    0.50*sqrt(var(datos_5[,1])*var(datos_5[,4])),
                    0.50*sqrt(var(datos_5[,2])*var(datos_5[,4])),
                    0.50*sqrt(var(datos_5[,3])*var(datos_5[,4])),
                    var(datos_5[,4]),
                    0.50*sqrt(var(datos_5[,4])*var(datos_5[,5])),
                    0.50*sqrt(var(datos_5[,1])*var(datos_5[,5])),
                    0.50*sqrt(var(datos_5[,2])*var(datos_5[,5])),
                    0.50*sqrt(var(datos_5[,3])*var(datos_5[,5])),
                    0.50*sqrt(var(datos_5[,4])*var(datos_5[,5])),
                    var(datos_5[,5])), 5, 5)
# Tamaño de muestra
tamano_muestra <- 50
#Cálculo de las probabilidades 
p1_5 <- pmvnorm(LIE_5, LSE_5, mean = t_5, sigma = sigma_5)
p2_5 <- pmvnorm(LIE_5, LSE_5, mean = medias_5, sigma = sigma_5)
sigmat_5 <- sigma_5 + t(t(medias_5 - t_5)) %*% t(medias_5 - t_5)
p3_5 <- pmvnorm(LIE_5, LSE_5, mean = t_5, sigma = sigmat_5)
#Cálculo de los índices
MCp_p1_5 <- (-1/3 * qnorm((1 - p1_5)/2))
MCpk_p2_5 <- (-1/3 * qnorm((1 - p2_5)/2))
MCpm_p3_5 <- (-1/3 * qnorm((1 - p3_5)/2))

# Cálculo de los valores reales de NMCp y NMCpm
cprima_5 <- min(((LSE_5[1] - t_5[1]) / sqrt(sigma_5[1, 1])), ((LSE_5[2] - t_5[2]) / sqrt(sigma_5[2, 2])), ((LSE_5[3] - t_5[3]) / sqrt(sigma_5[3, 3])), ((LSE_5[4] - t_5[4]) / sqrt(sigma_5[4, 4])), ((LSE_5[5] - t_5[5]) / sqrt(sigma_5[5, 5])))
NMCp_5 <- cprima_5 / sqrt(qchisq(1 - 0.0027, 5))

# Cálculo de la distancia D
D_5 <- (1 + (t(medias_5 - t_5) %*% solve(sigma_5) %*% (medias_5 - t_5)))^0.5
NMCpm_5 <- NMCp_5 / D_5  # D_real es 1

#Cálculo de las medias y sigmas muestrales
medias_5_50 <- matrix(0, nrow = cantidad_muestras, ncol = 5)
sigmas_5_50 <- array(0, dim = c(5, 5, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_5, Sigma = sigma_5)
  medias_5_50[i, ] <- colMeans(muestra)
  sigmas_5_50[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_5_50 <- data.frame(p1_5_50 = numeric(0), p2_5_50 = numeric(0), p3_5_50 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_5_50[, , i]
  p1_5_50 <- pmvnorm(LIE_5, LSE_5, mean = t_5, sigma = sigma_actual)
  p2_5_50 <- pmvnorm(LIE_5, LSE_5, mean = medias_5_50[i, ], sigma = sigma_actual)
  sigmat_5_50 <- sigma_actual + (medias_5_50[i, ] - t_5) %*% t(medias_5_50[i, ] - t_5)
  p3_5_50 <- pmvnorm(LIE_5, LSE_5, mean = t_5, sigma = sigmat_5_50)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_5[1] - t_5[1]) / sqrt(sigma_actual[1, 1])), ((LSE_5[2] - t_5[2]) / sqrt(sigma_actual[2, 2])), ((LSE_5[3] - t_5[3]) / sqrt(sigma_actual[3, 3])), ((LSE_5[4] - t_5[4]) / sqrt(sigma_actual[4, 4])), ((LSE_5[5] - t_5[5]) / sqrt(sigma_actual[5, 5])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 5))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_5 - medias_5_50[i, ]) %*% solve(sigma_actual) %*% (t_5 - medias_5_50[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_5_50[i, ] <- c(p1_5_50, p2_5_50, p3_5_50, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_5_50 <- mutate(datos_5_50,
                     MCp_p1_5_50 = (-1 / 3 * qnorm((1 - p1_5_50) / 2)),
                     MCpk_p2_5_50 = (-1 / 3 * qnorm((1 - p2_5_50) / 2)),
                     MCpm_p3_5_50 = (-1 / 3 * qnorm((1 - p3_5_50) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_5_50 <- data.frame(Estadistica = c("p1_5_50", "p2_5_50", "p3_5_50", "MCp_p1_5_50", "MCpk_p2_5_50", "MCpm_p3_5_50", "NMCp", "NMCpm"),
                              Percentil_2.5 = numeric(8),
                              Percentil_97.5 = numeric(8),
                              Medias = numeric(8),
                              Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_5_50$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_5_50[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_5_50[[estadistica]], probs = 0.975)
  
  intervalos_5_50$Percentil_2.5[i] <- percentil_2.5
  intervalos_5_50$Percentil_97.5[i] <- percentil_97.5
  intervalos_5_50$Medias[i] <- mean(datos_5_50[[estadistica]])
  intervalos_5_50$Desvios[i] <- sd(datos_5_50[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_5_50 <- mean(datos_5_50$p1_5_50) - p1_5
sesgo_p2_5_50 <- mean(datos_5_50$p2_5_50) - p2_5
sesgo_p3_5_50 <- mean(datos_5_50$p3_5_50) - p3_5
sesgo_MCp_p1_5_50 <- mean(datos_5_50$MCp_p1_5_50) - MCp_p1_5
sesgo_MCpk_p2_5_50 <- mean(datos_5_50$MCpk_p2_5_50) - MCpk_p2_5
sesgo_MCpm_p3_5_50 <- mean(datos_5_50$MCpm_p3_5_50) - MCpm_p3_5
sesgo_NMCp <- mean(datos_5_50$NMCp) - NMCp_5
sesgo_NMCpm <- mean(datos_5_50$NMCpm) - NMCpm_5

sesgo_p1_5_50_r <- (mean(datos_5_50$p1_5_50) - p1_5)/p1_5
sesgo_p2_5_50_r <- (mean(datos_5_50$p2_5_50) - p2_5)/p2_5
sesgo_p3_5_50_r <- (mean(datos_5_50$p3_5_50) - p3_5)/p3_5
sesgo_MCp_p1_5_50_r <- (mean(datos_5_50$MCp_p1_5_50) - MCp_p1_5)/MCp_p1_5
sesgo_MCpk_p2_5_50_r <- (mean(datos_5_50$MCpk_p2_5_50) - MCpk_p2_5)/MCpk_p2_5
sesgo_MCpm_p3_5_50_r <- (mean(datos_5_50$MCpm_p3_5_50) - MCpm_p3_5)/MCpm_p3_5
sesgo_NMCp_r <- (mean(datos_5_50$NMCp) - NMCp_5)/NMCp_5
sesgo_NMCpm_r <- (mean(datos_5_50$NMCpm) - NMCpm_5)/NMCpm_5

# Concateno los sesgos en un vector
sesgos_5_50 <- c(sesgo_p1_5_50, sesgo_p2_5_50, sesgo_p3_5_50, sesgo_MCp_p1_5_50, sesgo_MCpk_p2_5_50, sesgo_MCpm_p3_5_50, sesgo_NMCp, sesgo_NMCpm)
sesgos_5_50_r <- c(sesgo_p1_5_50_r, sesgo_p2_5_50_r, sesgo_p3_5_50_r, sesgo_MCp_p1_5_50_r, sesgo_MCpk_p2_5_50_r, sesgo_MCpm_p3_5_50_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_5_50$Sesgo <- sesgos_5_50
intervalos_5_50$ECM <- (intervalos_5_50$Desvios)^2 + (intervalos_5_50$Sesgo)^2
intervalos_5_50$SesgoRelativo <- sesgos_5_50_r
# Mostrar los resultados
print(intervalos_5_50)

#################
### p=5, n= 100 ###
#################

# Tamaño de muestra
tamano_muestra <- 100

#Cálculo de las medias y sigmas muestrales
medias_5_100 <- matrix(0, nrow = cantidad_muestras, ncol = 5)
sigmas_5_100 <- array(0, dim = c(5, 5, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_5, Sigma = sigma_5)
  medias_5_100[i, ] <- colMeans(muestra)
  sigmas_5_100[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_5_100 <- data.frame(p1_5_100 = numeric(0), p2_5_100 = numeric(0), p3_5_100 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_5_100[, , i]
  p1_5_100 <- pmvnorm(LIE_5, LSE_5, mean = t_5, sigma = sigma_actual)
  p2_5_100 <- pmvnorm(LIE_5, LSE_5, mean = medias_5_100[i, ], sigma = sigma_actual)
  sigmat_5_100 <- sigma_actual + (medias_5_100[i, ] - t_5) %*% t(medias_5_100[i, ] - t_5)
  p3_5_100 <- pmvnorm(LIE_5, LSE_5, mean = t_5, sigma = sigmat_5_100)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_5[1] - t_5[1]) / sqrt(sigma_actual[1, 1])), ((LSE_5[2] - t_5[2]) / sqrt(sigma_actual[2, 2])), ((LSE_5[3] - t_5[3]) / sqrt(sigma_actual[3, 3])), ((LSE_5[4] - t_5[4]) / sqrt(sigma_actual[4, 4])), ((LSE_5[5] - t_5[5]) / sqrt(sigma_actual[5, 5])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 5))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_5 - medias_5_100[i, ]) %*% solve(sigma_actual) %*% (t_5 - medias_5_100[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_5_100[i, ] <- c(p1_5_100, p2_5_100, p3_5_100, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_5_100 <- mutate(datos_5_100,
                      MCp_p1_5_100 = (-1 / 3 * qnorm((1 - p1_5_100) / 2)),
                      MCpk_p2_5_100 = (-1 / 3 * qnorm((1 - p2_5_100) / 2)),
                      MCpm_p3_5_100 = (-1 / 3 * qnorm((1 - p3_5_100) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_5_100 <- data.frame(Estadistica = c("p1_5_100", "p2_5_100", "p3_5_100", "MCp_p1_5_100", "MCpk_p2_5_100", "MCpm_p3_5_100", "NMCp", "NMCpm"),
                               Percentil_2.5 = numeric(8),
                               Percentil_97.5 = numeric(8),
                               Medias = numeric(8),
                               Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_5_100$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_5_100[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_5_100[[estadistica]], probs = 0.975)
  
  intervalos_5_100$Percentil_2.5[i] <- percentil_2.5
  intervalos_5_100$Percentil_97.5[i] <- percentil_97.5
  intervalos_5_100$Medias[i] <- mean(datos_5_100[[estadistica]])
  intervalos_5_100$Desvios[i] <- sd(datos_5_100[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_5_100 <- mean(datos_5_100$p1_5_100) - p1_5
sesgo_p2_5_100 <- mean(datos_5_100$p2_5_100) - p2_5
sesgo_p3_5_100 <- mean(datos_5_100$p3_5_100) - p3_5
sesgo_MCp_p1_5_100 <- mean(datos_5_100$MCp_p1_5_100) - MCp_p1_5
sesgo_MCpk_p2_5_100 <- mean(datos_5_100$MCpk_p2_5_100) - MCpk_p2_5
sesgo_MCpm_p3_5_100 <- mean(datos_5_100$MCpm_p3_5_100) - MCpm_p3_5
sesgo_NMCp <- mean(datos_5_100$NMCp) - NMCp_5
sesgo_NMCpm <- mean(datos_5_100$NMCpm) - NMCpm_5

sesgo_p1_5_100_r <- (mean(datos_5_100$p1_5_100) - p1_5)/p1_5
sesgo_p2_5_100_r <- (mean(datos_5_100$p2_5_100) - p2_5)/p2_5
sesgo_p3_5_100_r <- (mean(datos_5_100$p3_5_100) - p3_5)/p3_5
sesgo_MCp_p1_5_100_r <- (mean(datos_5_100$MCp_p1_5_100) - MCp_p1_5)/MCp_p1_5
sesgo_MCpk_p2_5_100_r <- (mean(datos_5_100$MCpk_p2_5_100) - MCpk_p2_5)/MCpk_p2_5
sesgo_MCpm_p3_5_100_r <- (mean(datos_5_100$MCpm_p3_5_100) - MCpm_p3_5)/MCpm_p3_5
sesgo_NMCp_r <- (mean(datos_5_100$NMCp) - NMCp_5)/NMCp_5
sesgo_NMCpm_r <- (mean(datos_5_100$NMCpm) - NMCpm_5)/NMCpm_5

# Concateno los sesgos en un vector
sesgos_5_100 <- c(sesgo_p1_5_100, sesgo_p2_5_100, sesgo_p3_5_100, sesgo_MCp_p1_5_100, sesgo_MCpk_p2_5_100, sesgo_MCpm_p3_5_100, sesgo_NMCp, sesgo_NMCpm)
sesgos_5_100_r <- c(sesgo_p1_5_100_r, sesgo_p2_5_100_r, sesgo_p3_5_100_r, sesgo_MCp_p1_5_100_r, sesgo_MCpk_p2_5_100_r, sesgo_MCpm_p3_5_100_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_5_100$Sesgo <- sesgos_5_100
intervalos_5_100$ECM <- (intervalos_5_100$Desvios)^2 + (intervalos_5_100$Sesgo)^2
intervalos_5_100$SesgoRelativo <- sesgos_5_100_r
# Mostrar los resultados
print(intervalos_5_100)

#################
### p=5, n= 200 ###
#################

# Tamaño de muestra
tamano_muestra <- 200

#Cálculo de las medias y sigmas muestrales
medias_5_200 <- matrix(0, nrow = cantidad_muestras, ncol = 5)
sigmas_5_200 <- array(0, dim = c(5, 5, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_5, Sigma = sigma_5)
  medias_5_200[i, ] <- colMeans(muestra)
  sigmas_5_200[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_5_200 <- data.frame(p1_5_200 = numeric(0), p2_5_200 = numeric(0), p3_5_200 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_5_200[, , i]
  p1_5_200 <- pmvnorm(LIE_5, LSE_5, mean = t_5, sigma = sigma_actual)
  p2_5_200 <- pmvnorm(LIE_5, LSE_5, mean = medias_5_200[i, ], sigma = sigma_actual)
  sigmat_5_200 <- sigma_actual + (medias_5_200[i, ] - t_5) %*% t(medias_5_200[i, ] - t_5)
  p3_5_200 <- pmvnorm(LIE_5, LSE_5, mean = t_5, sigma = sigmat_5_200)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_5[1] - t_5[1]) / sqrt(sigma_actual[1, 1])), ((LSE_5[2] - t_5[2]) / sqrt(sigma_actual[2, 2])), ((LSE_5[3] - t_5[3]) / sqrt(sigma_actual[3, 3])), ((LSE_5[4] - t_5[4]) / sqrt(sigma_actual[4, 4])), ((LSE_5[5] - t_5[5]) / sqrt(sigma_actual[5, 5])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 5))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_5 - medias_5_200[i, ]) %*% solve(sigma_actual) %*% (t_5 - medias_5_200[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_5_200[i, ] <- c(p1_5_200, p2_5_200, p3_5_200, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_5_200 <- mutate(datos_5_200,
                      MCp_p1_5_200 = (-1 / 3 * qnorm((1 - p1_5_200) / 2)),
                      MCpk_p2_5_200 = (-1 / 3 * qnorm((1 - p2_5_200) / 2)),
                      MCpm_p3_5_200 = (-1 / 3 * qnorm((1 - p3_5_200) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_5_200 <- data.frame(Estadistica = c("p1_5_200", "p2_5_200", "p3_5_200", "MCp_p1_5_200", "MCpk_p2_5_200", "MCpm_p3_5_200", "NMCp", "NMCpm"),
                               Percentil_2.5 = numeric(8),
                               Percentil_97.5 = numeric(8),
                               Medias = numeric(8),
                               Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_5_200$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_5_200[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_5_200[[estadistica]], probs = 0.975)
  
  intervalos_5_200$Percentil_2.5[i] <- percentil_2.5
  intervalos_5_200$Percentil_97.5[i] <- percentil_97.5
  intervalos_5_200$Medias[i] <- mean(datos_5_200[[estadistica]])
  intervalos_5_200$Desvios[i] <- sd(datos_5_200[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_5_200 <- mean(datos_5_200$p1_5_200) - p1_5
sesgo_p2_5_200 <- mean(datos_5_200$p2_5_200) - p2_5
sesgo_p3_5_200 <- mean(datos_5_200$p3_5_200) - p3_5
sesgo_MCp_p1_5_200 <- mean(datos_5_200$MCp_p1_5_200) - MCp_p1_5
sesgo_MCpk_p2_5_200 <- mean(datos_5_200$MCpk_p2_5_200) - MCpk_p2_5
sesgo_MCpm_p3_5_200 <- mean(datos_5_200$MCpm_p3_5_200) - MCpm_p3_5
sesgo_NMCp <- mean(datos_5_200$NMCp) - NMCp_5
sesgo_NMCpm <- mean(datos_5_200$NMCpm) - NMCpm_5

sesgo_p1_5_200_r <- (mean(datos_5_200$p1_5_200) - p1_5)/p1_5
sesgo_p2_5_200_r <- (mean(datos_5_200$p2_5_200) - p2_5)/p2_5
sesgo_p3_5_200_r <- (mean(datos_5_200$p3_5_200) - p3_5)/p3_5
sesgo_MCp_p1_5_200_r <- (mean(datos_5_200$MCp_p1_5_200) - MCp_p1_5)/MCp_p1_5
sesgo_MCpk_p2_5_200_r <- (mean(datos_5_200$MCpk_p2_5_200) - MCpk_p2_5)/MCpk_p2_5
sesgo_MCpm_p3_5_200_r <- (mean(datos_5_200$MCpm_p3_5_200) - MCpm_p3_5)/MCpm_p3_5
sesgo_NMCp_r <- (mean(datos_5_200$NMCp) - NMCp_5)/NMCp_5
sesgo_NMCpm_r <- (mean(datos_5_200$NMCpm) - NMCpm_5)/NMCpm_5

# Concateno los sesgos en un vector
sesgos_5_200 <- c(sesgo_p1_5_200, sesgo_p2_5_200, sesgo_p3_5_200, sesgo_MCp_p1_5_200, sesgo_MCpk_p2_5_200, sesgo_MCpm_p3_5_200, sesgo_NMCp, sesgo_NMCpm)
sesgos_5_200_r <- c(sesgo_p1_5_200_r, sesgo_p2_5_200_r, sesgo_p3_5_200_r, sesgo_MCp_p1_5_200_r, sesgo_MCpk_p2_5_200_r, sesgo_MCpm_p3_5_200_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_5_200$Sesgo <- sesgos_5_200
intervalos_5_200$ECM <- (intervalos_5_200$Desvios)^2 + (intervalos_5_200$Sesgo)^2
intervalos_5_200$SesgoRelativo <- sesgos_5_200_r
# Mostrar los resultados
print(intervalos_5_200)


#################
### p=6, n= 50 ###
#################
datos_6 <- datos_multi[, c(4, 4, 4, 4, 4, 4)]
#Parámetros
LIE_6 <- c(3.00, 3.00, 3.00, 3.00, 3.00, 3.00)
LSE_6 <- c(17.00, 17.00, 17.00, 17.00, 17.00, 17.00)
t_6 <- c(10, 10, 10, 10, 10, 10)
medias_6 <- colMeans(datos_6)
sigma_6 <- matrix(c(var(datos_6[,1]),
                    0.50*sqrt(var(datos_6[,1])*var(datos_6[,2])),
                    0.50*sqrt(var(datos_6[,1])*var(datos_6[,3])),
                    0.50*sqrt(var(datos_6[,1])*var(datos_6[,4])),
                    0.50*sqrt(var(datos_6[,1])*var(datos_6[,5])),
                    0.50*sqrt(var(datos_6[,1])*var(datos_6[,6])),
                    0.50*sqrt(var(datos_6[,1])*var(datos_6[,2])),
                    var(datos_6[,2]),
                    0.50*sqrt(var(datos_6[,2])*var(datos_6[,3])),
                    0.50*sqrt(var(datos_6[,2])*var(datos_6[,4])),
                    0.50*sqrt(var(datos_6[,2])*var(datos_6[,5])),
                    0.50*sqrt(var(datos_6[,2])*var(datos_6[,6])),
                    0.50*sqrt(var(datos_6[,1])*var(datos_6[,3])),
                    0.50*sqrt(var(datos_6[,2])*var(datos_6[,3])),
                    var(datos_6[,3]),
                    0.50*sqrt(var(datos_6[,3])*var(datos_6[,4])),
                    0.50*sqrt(var(datos_6[,3])*var(datos_6[,5])),
                    0.50*sqrt(var(datos_6[,3])*var(datos_6[,6])),
                    0.50*sqrt(var(datos_6[,1])*var(datos_6[,4])),
                    0.50*sqrt(var(datos_6[,2])*var(datos_6[,4])),
                    0.50*sqrt(var(datos_6[,3])*var(datos_6[,4])),
                    var(datos_6[,4]),
                    0.50*sqrt(var(datos_6[,4])*var(datos_6[,5])),
                    0.50*sqrt(var(datos_6[,4])*var(datos_6[,6])),
                    0.50*sqrt(var(datos_6[,1])*var(datos_6[,5])),
                    0.50*sqrt(var(datos_6[,2])*var(datos_6[,5])),
                    0.50*sqrt(var(datos_6[,3])*var(datos_6[,5])),
                    0.50*sqrt(var(datos_6[,4])*var(datos_6[,5])),
                    var(datos_6[,5]),
                    0.50*sqrt(var(datos_6[,5])*var(datos_6[,6])),
                    0.50*sqrt(var(datos_6[,1])*var(datos_6[,6])),
                    0.50*sqrt(var(datos_6[,2])*var(datos_6[,6])),
                    0.50*sqrt(var(datos_6[,3])*var(datos_6[,6])),
                    0.50*sqrt(var(datos_6[,4])*var(datos_6[,6])),
                    0.50*sqrt(var(datos_6[,5])*var(datos_6[,6])),
                    var(datos_6[,6])), 6, 6)

# Tamaño de muestra
tamano_muestra <- 50
#Cálculo de las probabilidades 
p1_6 <- pmvnorm(LIE_6, LSE_6, mean = t_6, sigma = sigma_6)
p2_6 <- pmvnorm(LIE_6, LSE_6, mean = medias_6, sigma = sigma_6)
sigmat_6 <- sigma_6 + t(t(medias_6 - t_6)) %*% t(medias_6 - t_6)
p3_6 <- pmvnorm(LIE_6, LSE_6, mean = t_6, sigma = sigmat_6)
#Cálculo de los índices
MCp_p1_6 <- (-1/3 * qnorm((1 - p1_6)/2))
MCpk_p2_6 <- (-1/3 * qnorm((1 - p2_6)/2))
MCpm_p3_6 <- (-1/3 * qnorm((1 - p3_6)/2))

# Cálculo de los valores reales de NMCp y NMCpm
cprima_6 <- min(((LSE_6[1] - t_6[1]) / sqrt(sigma_6[1, 1])), ((LSE_6[2] - t_6[2]) / sqrt(sigma_6[2, 2])), ((LSE_6[3] - t_6[3]) / sqrt(sigma_6[3, 3])), ((LSE_6[4] - t_6[4]) / sqrt(sigma_6[4, 4])), ((LSE_6[5] - t_6[5]) / sqrt(sigma_6[5, 5])), ((LSE_6[6] - t_6[6]) / sqrt(sigma_6[6, 6])))
NMCp_6 <- cprima_6 / sqrt(qchisq(1 - 0.0027, 6))

# Cálculo de la distancia D
D_6 <- (1 + tamano_muestra / (tamano_muestra - 1) * 
          (t(medias_6 - t_6) %*% solve(sigma_6) %*% (medias_6 - t_6)))^0.5
NMCpm_6 <- NMCp_6 / D_6  # D_real es 1

#Cálculo de las medias y sigmas muestrales
medias_6_50 <- matrix(0, nrow = cantidad_muestras, ncol = 6)
sigmas_6_50 <- array(0, dim = c(6, 6, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_6, Sigma = sigma_6)
  medias_6_50[i, ] <- colMeans(muestra)
  sigmas_6_50[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_6_50 <- data.frame(p1_6_50 = numeric(0), p2_6_50 = numeric(0), p3_6_50 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_6_50[, , i]
  p1_6_50 <- pmvnorm(LIE_6, LSE_6, mean = t_6, sigma = sigma_actual)
  p2_6_50 <- pmvnorm(LIE_6, LSE_6, mean = medias_6_50[i, ], sigma = sigma_actual)
  sigmat_6_50 <- sigma_actual + (medias_6_50[i, ] - t_6) %*% t(medias_6_50[i, ] - t_6)
  p3_6_50 <- pmvnorm(LIE_6, LSE_6, mean = t_6, sigma = sigmat_6_50)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_6[1] - t_6[1]) / sqrt(sigma_actual[1, 1])), ((LSE_6[2] - t_6[2]) / sqrt(sigma_actual[2, 2])), ((LSE_6[3] - t_6[3]) / sqrt(sigma_actual[3, 3])), ((LSE_6[4] - t_6[4]) / sqrt(sigma_actual[4, 4])), ((LSE_6[5] - t_6[5]) / sqrt(sigma_actual[5, 5])), ((LSE_6[6] - t_6[6]) / sqrt(sigma_actual[6, 6])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 6))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_6 - medias_6_50[i, ]) %*% solve(sigma_actual) %*% (t_6 - medias_6_50[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_6_50[i, ] <- c(p1_6_50, p2_6_50, p3_6_50, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_6_50 <- mutate(datos_6_50,
                     MCp_p1_6_50 = (-1 / 3 * qnorm((1 - p1_6_50) / 2)),
                     MCpk_p2_6_50 = (-1 / 3 * qnorm((1 - p2_6_50) / 2)),
                     MCpm_p3_6_50 = (-1 / 3 * qnorm((1 - p3_6_50) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_6_50 <- data.frame(Estadistica = c("p1_6_50", "p2_6_50", "p3_6_50", "MCp_p1_6_50", "MCpk_p2_6_50", "MCpm_p3_6_50", "NMCp", "NMCpm"),
                              Percentil_2.5 = numeric(8),
                              Percentil_97.5 = numeric(8),
                              Medias = numeric(8),
                              Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_6_50$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_6_50[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_6_50[[estadistica]], probs = 0.975)
  
  intervalos_6_50$Percentil_2.5[i] <- percentil_2.5
  intervalos_6_50$Percentil_97.5[i] <- percentil_97.5
  intervalos_6_50$Medias[i] <- mean(datos_6_50[[estadistica]])
  intervalos_6_50$Desvios[i] <- sd(datos_6_50[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_6_50 <- mean(datos_6_50$p1_6_50) - p1_6
sesgo_p2_6_50 <- mean(datos_6_50$p2_6_50) - p2_6
sesgo_p3_6_50 <- mean(datos_6_50$p3_6_50) - p3_6
sesgo_MCp_p1_6_50 <- mean(datos_6_50$MCp_p1_6_50) - MCp_p1_6
sesgo_MCpk_p2_6_50 <- mean(datos_6_50$MCpk_p2_6_50) - MCpk_p2_6
sesgo_MCpm_p3_6_50 <- mean(datos_6_50$MCpm_p3_6_50) - MCpm_p3_6
sesgo_NMCp <- mean(datos_6_50$NMCp) - NMCp_6
sesgo_NMCpm <- mean(datos_6_50$NMCpm) - NMCpm_6

sesgo_p1_6_50_r <- (mean(datos_6_50$p1_6_50) - p1_6)/p1_6
sesgo_p2_6_50_r <- (mean(datos_6_50$p2_6_50) - p2_6)/p2_6
sesgo_p3_6_50_r <- (mean(datos_6_50$p3_6_50) - p3_6)/p3_6
sesgo_MCp_p1_6_50_r <- (mean(datos_6_50$MCp_p1_6_50) - MCp_p1_6)/MCp_p1_6
sesgo_MCpk_p2_6_50_r <- (mean(datos_6_50$MCpk_p2_6_50) - MCpk_p2_6)/MCpk_p2_6
sesgo_MCpm_p3_6_50_r <- (mean(datos_6_50$MCpm_p3_6_50) - MCpm_p3_6)/MCpm_p3_6
sesgo_NMCp_r <- (mean(datos_6_50$NMCp) - NMCp_6)/NMCp_6
sesgo_NMCpm_r <- (mean(datos_6_50$NMCpm) - NMCpm_6)/NMCpm_6

# Concateno los sesgos en un vector
sesgos_6_50 <- c(sesgo_p1_6_50, sesgo_p2_6_50, sesgo_p3_6_50, sesgo_MCp_p1_6_50, sesgo_MCpk_p2_6_50, sesgo_MCpm_p3_6_50, sesgo_NMCp, sesgo_NMCpm)
sesgos_6_50_r <- c(sesgo_p1_6_50_r, sesgo_p2_6_50_r, sesgo_p3_6_50_r, sesgo_MCp_p1_6_50_r, sesgo_MCpk_p2_6_50_r, sesgo_MCpm_p3_6_50_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_6_50$Sesgo <- sesgos_6_50
intervalos_6_50$ECM <- (intervalos_6_50$Desvios)^2 + (intervalos_6_50$Sesgo)^2
intervalos_6_50$SesgoRelativo <- sesgos_6_50_r
# Mostrar los resultados
print(intervalos_6_50)

#################
### p=6, n= 100 ###
#################

# Tamaño de muestra
tamano_muestra <- 100

#Cálculo de las medias y sigmas muestrales
medias_6_100 <- matrix(0, nrow = cantidad_muestras, ncol = 6)
sigmas_6_100 <- array(0, dim = c(6, 6, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_6, Sigma = sigma_6)
  medias_6_100[i, ] <- colMeans(muestra)
  sigmas_6_100[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_6_100 <- data.frame(p1_6_100 = numeric(0), p2_6_100 = numeric(0), p3_6_100 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_6_100[, , i]
  p1_6_100 <- pmvnorm(LIE_6, LSE_6, mean = t_6, sigma = sigma_actual)
  p2_6_100 <- pmvnorm(LIE_6, LSE_6, mean = medias_6_100[i, ], sigma = sigma_actual)
  sigmat_6_100 <- sigma_actual + (medias_6_100[i, ] - t_6) %*% t(medias_6_100[i, ] - t_6)
  p3_6_100 <- pmvnorm(LIE_6, LSE_6, mean = t_6, sigma = sigmat_6_100)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_6[1] - t_6[1]) / sqrt(sigma_actual[1, 1])), ((LSE_6[2] - t_6[2]) / sqrt(sigma_actual[2, 2])), ((LSE_6[3] - t_6[3]) / sqrt(sigma_actual[3, 3])), ((LSE_6[4] - t_6[4]) / sqrt(sigma_actual[4, 4])), ((LSE_6[5] - t_6[5]) / sqrt(sigma_actual[5, 5])), ((LSE_6[6] - t_6[6]) / sqrt(sigma_actual[6, 6])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 6))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_6 - medias_6_100[i, ]) %*% solve(sigma_actual) %*% (t_6 - medias_6_100[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_6_100[i, ] <- c(p1_6_100, p2_6_100, p3_6_100, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_6_100 <- mutate(datos_6_100,
                      MCp_p1_6_100 = (-1 / 3 * qnorm((1 - p1_6_100) / 2)),
                      MCpk_p2_6_100 = (-1 / 3 * qnorm((1 - p2_6_100) / 2)),
                      MCpm_p3_6_100 = (-1 / 3 * qnorm((1 - p3_6_100) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_6_100 <- data.frame(Estadistica = c("p1_6_100", "p2_6_100", "p3_6_100", "MCp_p1_6_100", "MCpk_p2_6_100", "MCpm_p3_6_100", "NMCp", "NMCpm"),
                               Percentil_2.5 = numeric(8),
                               Percentil_97.5 = numeric(8),
                               Medias = numeric(8),
                               Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_6_100$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_6_100[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_6_100[[estadistica]], probs = 0.975)
  
  intervalos_6_100$Percentil_2.5[i] <- percentil_2.5
  intervalos_6_100$Percentil_97.5[i] <- percentil_97.5
  intervalos_6_100$Medias[i] <- mean(datos_6_100[[estadistica]])
  intervalos_6_100$Desvios[i] <- sd(datos_6_100[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_6_100 <- mean(datos_6_100$p1_6_100) - p1_6
sesgo_p2_6_100 <- mean(datos_6_100$p2_6_100) - p2_6
sesgo_p3_6_100 <- mean(datos_6_100$p3_6_100) - p3_6
sesgo_MCp_p1_6_100 <- mean(datos_6_100$MCp_p1_6_100) - MCp_p1_6
sesgo_MCpk_p2_6_100 <- mean(datos_6_100$MCpk_p2_6_100) - MCpk_p2_6
sesgo_MCpm_p3_6_100 <- mean(datos_6_100$MCpm_p3_6_100) - MCpm_p3_6
sesgo_NMCp <- mean(datos_6_100$NMCp) - NMCp_6
sesgo_NMCpm <- mean(datos_6_100$NMCpm) - NMCpm_6

sesgo_p1_6_100_r <- (mean(datos_6_100$p1_6_100) - p1_6)/p1_6
sesgo_p2_6_100_r <- (mean(datos_6_100$p2_6_100) - p2_6)/p2_6
sesgo_p3_6_100_r <- (mean(datos_6_100$p3_6_100) - p3_6)/p3_6
sesgo_MCp_p1_6_100_r <- (mean(datos_6_100$MCp_p1_6_100) - MCp_p1_6)/MCp_p1_6
sesgo_MCpk_p2_6_100_r <- (mean(datos_6_100$MCpk_p2_6_100) - MCpk_p2_6)/MCpk_p2_6
sesgo_MCpm_p3_6_100_r <- (mean(datos_6_100$MCpm_p3_6_100) - MCpm_p3_6)/MCpm_p3_6
sesgo_NMCp_r <- (mean(datos_6_100$NMCp) - NMCp_6)/NMCp_6
sesgo_NMCpm_r <- (mean(datos_6_100$NMCpm) - NMCpm_6)/NMCpm_6

# Concateno los sesgos en un vector
sesgos_6_100 <- c(sesgo_p1_6_100, sesgo_p2_6_100, sesgo_p3_6_100, sesgo_MCp_p1_6_100, sesgo_MCpk_p2_6_100, sesgo_MCpm_p3_6_100, sesgo_NMCp, sesgo_NMCpm)
sesgos_6_100_r <- c(sesgo_p1_6_100_r, sesgo_p2_6_100_r, sesgo_p3_6_100_r, sesgo_MCp_p1_6_100_r, sesgo_MCpk_p2_6_100_r, sesgo_MCpm_p3_6_100_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_6_100$Sesgo <- sesgos_6_100
intervalos_6_100$ECM <- (intervalos_6_100$Desvios)^2 + (intervalos_6_100$Sesgo)^2
intervalos_6_100$SesgoRelativo <- sesgos_6_100_r
# Mostrar los resultados
print(intervalos_6_100)

#################
### p=6, n= 200 ###
#################

# Tamaño de muestra
tamano_muestra <- 200

#Cálculo de las medias y sigmas muestrales
medias_6_200 <- matrix(0, nrow = cantidad_muestras, ncol = 6)
sigmas_6_200 <- array(0, dim = c(6, 6, cantidad_muestras))
# Generar las muestras con las medias y las cov
for (i in 1:cantidad_muestras) {
  muestra <- mvrnorm(n = tamano_muestra, mu = medias_6, Sigma = sigma_6)
  medias_6_200[i, ] <- colMeans(muestra)
  sigmas_6_200[, , i] <- cov(muestra)
}
# Crear un dataframe para almacenar los resultados
datos_6_200 <- data.frame(p1_6_200 = numeric(0), p2_6_200 = numeric(0), p3_6_200 = numeric(0), NMCp = numeric(0), NMCpm = numeric(0))
# Calcular las probabilidades y los índices NMCp y NMCpm para cada conjunto de parámetros
for (i in 1:cantidad_muestras) {
  sigma_actual <- sigmas_6_200[, , i]
  p1_6_200 <- pmvnorm(LIE_6, LSE_6, mean = t_6, sigma = sigma_actual)
  p2_6_200 <- pmvnorm(LIE_6, LSE_6, mean = medias_6_200[i, ], sigma = sigma_actual)
  sigmat_6_200 <- sigma_actual + (medias_6_200[i, ] - t_6) %*% t(medias_6_200[i, ] - t_6)
  p3_6_200 <- pmvnorm(LIE_6, LSE_6, mean = t_6, sigma = sigmat_6_200)
  
  # Calcular NMCp y NMCpm
  cprima <- min(((LSE_6[1] - t_6[1]) / sqrt(sigma_actual[1, 1])), ((LSE_6[2] - t_6[2]) / sqrt(sigma_actual[2, 2])), ((LSE_6[3] - t_6[3]) / sqrt(sigma_actual[3, 3])), ((LSE_6[4] - t_6[4]) / sqrt(sigma_actual[4, 4])), ((LSE_6[5] - t_6[5]) / sqrt(sigma_actual[5, 5])), ((LSE_6[6] - t_6[6]) / sqrt(sigma_actual[6, 6])))
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 6))  
  D <- (1 + tamano_muestra / (tamano_muestra - 1) * (t(t_6 - medias_6_200[i, ]) %*% solve(sigma_actual) %*% (t_6 - medias_6_200[i, ])))^0.5
  NMCpm <- NMCp / D
  
  datos_6_200[i, ] <- c(p1_6_200, p2_6_200, p3_6_200, NMCp, NMCpm)
}

# Calcular MCp, MCpk y MCpm para las muestras
datos_6_200 <- mutate(datos_6_200,
                      MCp_p1_6_200 = (-1 / 3 * qnorm((1 - p1_6_200) / 2)),
                      MCpk_p2_6_200 = (-1 / 3 * qnorm((1 - p2_6_200) / 2)),
                      MCpm_p3_6_200 = (-1 / 3 * qnorm((1 - p3_6_200) / 2)))

# Crear un dataframe para almacenar los percentiles y estadísticas
intervalos_6_200 <- data.frame(Estadistica = c("p1_6_200", "p2_6_200", "p3_6_200", "MCp_p1_6_200", "MCpk_p2_6_200", "MCpm_p3_6_200", "NMCp", "NMCpm"),
                               Percentil_2.5 = numeric(8),
                               Percentil_97.5 = numeric(8),
                               Medias = numeric(8),
                               Desvios = numeric(8))

# Calcular los percentiles, medias y desviaciones para cada estadística
for (i in 1:8) {
  estadistica <- as.character(intervalos_6_200$Estadistica[i])
  
  percentil_2.5 <- quantile(datos_6_200[[estadistica]], probs = 0.025)
  percentil_97.5 <- quantile(datos_6_200[[estadistica]], probs = 0.975)
  
  intervalos_6_200$Percentil_2.5[i] <- percentil_2.5
  intervalos_6_200$Percentil_97.5[i] <- percentil_97.5
  intervalos_6_200$Medias[i] <- mean(datos_6_200[[estadistica]])
  intervalos_6_200$Desvios[i] <- sd(datos_6_200[[estadistica]])
}

# Calcular sesgos y ECM
sesgo_p1_6_200 <- mean(datos_6_200$p1_6_200) - p1_6
sesgo_p2_6_200 <- mean(datos_6_200$p2_6_200) - p2_6
sesgo_p3_6_200 <- mean(datos_6_200$p3_6_200) - p3_6
sesgo_MCp_p1_6_200 <- mean(datos_6_200$MCp_p1_6_200) - MCp_p1_6
sesgo_MCpk_p2_6_200 <- mean(datos_6_200$MCpk_p2_6_200) - MCpk_p2_6
sesgo_MCpm_p3_6_200 <- mean(datos_6_200$MCpm_p3_6_200) - MCpm_p3_6

sesgo_NMCp <- mean(datos_6_200$NMCp) - NMCp_6
sesgo_NMCpm <- mean(datos_6_200$NMCpm) - NMCpm_6

sesgo_p1_6_200_r <- (mean(datos_6_200$p1_6_200) - p1_6)/p1_6
sesgo_p2_6_200_r <- (mean(datos_6_200$p2_6_200) - p2_6)/p2_6
sesgo_p3_6_200_r <- (mean(datos_6_200$p3_6_200) - p3_6)/p3_6
sesgo_MCp_p1_6_200_r <- (mean(datos_6_200$MCp_p1_6_200) - MCp_p1_6)/MCp_p1_6
sesgo_MCpk_p2_6_200_r <- (mean(datos_6_200$MCpk_p2_6_200) - MCpk_p2_6)/MCpk_p2_6
sesgo_MCpm_p3_6_200_r <- (mean(datos_6_200$MCpm_p3_6_200) - MCpm_p3_6)/MCpm_p3_6
sesgo_NMCp_r <- (mean(datos_6_200$NMCp) - NMCp_6)/NMCp_6
sesgo_NMCpm_r <- (mean(datos_6_200$NMCpm) - NMCpm_6)/NMCpm_6

# Concateno los sesgos en un vector
sesgos_6_200 <- c(sesgo_p1_6_200, sesgo_p2_6_200, sesgo_p3_6_200, sesgo_MCp_p1_6_200, sesgo_MCpk_p2_6_200, sesgo_MCpm_p3_6_200, sesgo_NMCp, sesgo_NMCpm)
sesgos_6_200_r <- c(sesgo_p1_6_200_r, sesgo_p2_6_200_r, sesgo_p3_6_200_r, sesgo_MCp_p1_6_200_r, sesgo_MCpk_p2_6_200_r, sesgo_MCpm_p3_6_200_r, sesgo_NMCp_r, sesgo_NMCpm_r)
intervalos_6_200$Sesgo <- sesgos_6_200
intervalos_6_200$ECM <- (intervalos_6_200$Desvios)^2 + (intervalos_6_200$Sesgo)^2
intervalos_6_200$SesgoRelativo <- sesgos_6_200_r
# Mostrar los resultados
print(intervalos_6_200)

############GRÁFICOs####################

colnames(datos_2_50)[which(colnames(datos_2_50) == "MCp_p1_2_50")] <- "MCp_p1"
colnames(datos_2_50)[which(colnames(datos_2_50) == "MCpk_p2_2_50")] <- "MCpk_p2"
colnames(datos_2_50)[which(colnames(datos_2_50) == "MCpm_p3_2_50")] <- "MCpm_p3"

colnames(datos_3_50)[which(colnames(datos_3_50) == "MCp_p1_3_50")] <- "MCp_p1"
colnames(datos_3_50)[which(colnames(datos_3_50) == "MCpk_p2_3_50")] <- "MCpk_p2"
colnames(datos_3_50)[which(colnames(datos_3_50) == "MCpm_p3_3_50")] <- "MCpm_p3"

colnames(datos_4_50)[which(colnames(datos_4_50) == "MCp_p1_4_50")] <- "MCp_p1"
colnames(datos_4_50)[which(colnames(datos_4_50) == "MCpk_p2_4_50")] <- "MCpk_p2"
colnames(datos_4_50)[which(colnames(datos_4_50) == "MCpm_p3_4_50")] <- "MCpm_p3"

colnames(datos_5_50)[which(colnames(datos_5_50) == "MCp_p1_5_50")] <- "MCp_p1"
colnames(datos_5_50)[which(colnames(datos_5_50) == "MCpk_p2_5_50")] <- "MCpk_p2"
colnames(datos_5_50)[which(colnames(datos_5_50) == "MCpm_p3_5_50")] <- "MCpm_p3"

colnames(datos_6_50)[which(colnames(datos_6_50) == "MCp_p1_6_50")] <- "MCp_p1"
colnames(datos_6_50)[which(colnames(datos_6_50) == "MCpk_p2_6_50")] <- "MCpk_p2"
colnames(datos_6_50)[which(colnames(datos_6_50) == "MCpm_p3_6_50")] <- "MCpm_p3"
########### n=50 ###########
datos_2_50$p<-"2"
datos_3_50$p<-"3"
datos_4_50$p<-"4"
datos_5_50$p<-"5"
datos_6_50$p<-"6"

##############MCp_p1
# Definir los dataframes individuales
g2_MCp_p1_50 <- datos_2_50[, c("MCp_p1", "p")]
g3_MCp_p1_50 <- datos_3_50[, c("MCp_p1", "p")]
g4_MCp_p1_50 <- datos_4_50[, c("MCp_p1", "p")]
g5_MCp_p1_50 <- datos_5_50[, c("MCp_p1", "p")]
g6_MCp_p1_50 <- datos_6_50[, c("MCp_p1", "p")]

# Combinar los data frames
df_MCp_p1_50 <- rbind(g2_MCp_p1_50, g3_MCp_p1_50, g4_MCp_p1_50, g5_MCp_p1_50, g6_MCp_p1_50)

# Ordenar los niveles de la variable n
df_MCp_p1_50$p <- factor(df_MCp_p1_50$p, levels = c("2", "3", "4", "5", "6"))
# Crear el gráfico de densidad con ggplot2
plot_df_MCp_p1_50<-ggplot(df_MCp_p1_50, aes(x = MCp_p1, fill = p)) +
  geom_density(aes(color = p), alpha = 0.5, size = 1.5, show.legend = TRUE) +
  labs(title = expression("Índice " * hat(MC)[p(p[1])] * ""), x = "Índice estimado", y = "Densidad") +
  scale_fill_manual(values = c("2" = "olivedrab3", "3" = "cornflowerblue", "4" = "turquoise", "5" = "coral", "6"="pink"),
                    name = "Cantidad de variables", labels = c("2", "3", "4", "5", "6")) +
  scale_color_manual(values = c("2" = "olivedrab3", "3" = "cornflowerblue", "4" = "turquoise", "5" = "coral", "6"="pink"),
                     name = "Cantidad de variables", labels = c("2", "3", "4", "5", "6")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 16)
  )

##############MCpk_p2
# Definir los dataframes individuales
g2_MCpk_p2_50 <- datos_2_50[, c("MCpk_p2", "p")]
g3_MCpk_p2_50 <- datos_3_50[, c("MCpk_p2", "p")]
g4_MCpk_p2_50 <- datos_4_50[, c("MCpk_p2", "p")]
g5_MCpk_p2_50 <- datos_5_50[, c("MCpk_p2", "p")]
g6_MCpk_p2_50 <- datos_6_50[, c("MCpk_p2", "p")]

# Combinar los data frames
df_MCpk_p2_50 <- rbind(g2_MCpk_p2_50, g3_MCpk_p2_50, g4_MCpk_p2_50, g5_MCpk_p2_50, g6_MCpk_p2_50)

# Ordenar los niveles de la variable n
df_MCpk_p2_50$p <- factor(df_MCpk_p2_50$p, levels = c("2", "3", "4", "5", "6"))

# Crear el gráfico de densidad con ggplot2
plot_df_MCpk_p2_50<-ggplot(df_MCpk_p2_50, aes(x = MCpk_p2, fill = p)) +
  geom_density(aes(color = p), alpha = 0.5, size = 1.5, show.legend = TRUE) +
  labs(title = expression("Índice " * hat(MC)[pk(p[2])] * ""), x = "Índice estimado", y = "Densidad") +
  scale_fill_manual(values = c("2" = "olivedrab3", "3" = "cornflowerblue", "4" = "turquoise", "5" = "coral", "6"="pink"),
                    name = "Cantidad de variables", labels = c("2", "3", "4", "5", "6")) +
  scale_color_manual(values = c("2" = "olivedrab3", "3" = "cornflowerblue", "4" = "turquoise", "5" = "coral", "6"="pink"),
                     name = "Cantidad de variables", labels = c("2", "3", "4", "5", "6")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 16)
  )

##############MCpm_p3
# Definir los dataframes individuales
g2_MCpm_p3_50 <- datos_2_50[, c("MCpm_p3", "p")]
g3_MCpm_p3_50 <- datos_3_50[, c("MCpm_p3", "p")]
g4_MCpm_p3_50 <- datos_4_50[, c("MCpm_p3", "p")]
g5_MCpm_p3_50 <- datos_5_50[, c("MCpm_p3", "p")]
g6_MCpm_p3_50 <- datos_6_50[, c("MCpm_p3", "p")]

# Combinar los data frames
df_MCpm_p3_50 <- rbind(g2_MCpm_p3_50, g3_MCpm_p3_50, g4_MCpm_p3_50, g5_MCpm_p3_50, g6_MCpm_p3_50)

# Ordenar los niveles de la variable n
df_MCpm_p3_50$p <- factor(df_MCpm_p3_50$p, levels = c("2", "3", "4", "5", "6"))

# Crear el gráfico de densidad con ggplot2
plot_df_MCpm_p3_50<-ggplot(df_MCpm_p3_50, aes(x = MCpm_p3, fill = p)) +
  geom_density(aes(color = p), alpha = 0.5, size = 1.5, show.legend = TRUE) +
  labs(title = expression("Índice " * hat(MC)[pm(p[3])] * ""), x = "Índice estimado", y = "Densidad") +
  scale_fill_manual(values = c("2" = "olivedrab3", "3" = "cornflowerblue", "4" = "turquoise", "5" = "coral", "6"="pink"),
                    name = "Cantidad de variables", labels = c("2", "3", "4", "5", "6")) +
  scale_color_manual(values = c("2" = "olivedrab3", "3" = "cornflowerblue", "4" = "turquoise", "5" = "coral", "6"="pink"),
                     name = "Cantidad de variables", labels = c("2", "3", "4", "5", "6")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 16)
  )

##############NMCp
# Definir los dataframes individuales
g2_NMCp_50 <- datos_2_50[, c("NMCp", "p")]
g3_NMCp_50 <- datos_3_50[, c("NMCp", "p")]
g4_NMCp_50 <- datos_4_50[, c("NMCp", "p")]
g5_NMCp_50 <- datos_5_50[, c("NMCp", "p")]
g6_NMCp_50 <- datos_6_50[, c("NMCp", "p")]

# Combinar los data frames
df_NMCp_50 <- rbind(g2_NMCp_50, g3_NMCp_50, g4_NMCp_50, g5_NMCp_50, g6_NMCp_50)

# Ordenar los niveles de la variable n
df_NMCp_50$p <- factor(df_NMCp_50$p, levels = c("2", "3", "4", "5", "6"))

# Crear el gráfico de densidad con ggplot2
plot_df_NMCp_50<-ggplot(df_NMCp_50, aes(x = NMCp, fill = p)) +
  geom_density(aes(color = p), alpha = 0.5, size = 1.5, show.legend = TRUE) +
  labs(title = expression("Índice " * hat(NMC)[p] * ""), x = "Índice estimado", y = "Densidad") +
  scale_fill_manual(values = c("2" = "olivedrab3", "3" = "cornflowerblue", "4" = "turquoise", "5" = "coral", "6"="pink"),
                    name = "Cantidad de variables", labels = c("2", "3", "4", "5", "6")) +
  scale_color_manual(values = c("2" = "olivedrab3", "3" = "cornflowerblue", "4" = "turquoise", "5" = "coral", "6"="pink"),
                     name = "Cantidad de variables", labels = c("2", "3", "4", "5", "6")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 16)
  )

##############NMCpm
# Definir los dataframes individuales
g2_NMCpm_50 <- datos_2_50[, c("NMCpm", "p")]
g3_NMCpm_50 <- datos_3_50[, c("NMCpm", "p")]
g4_NMCpm_50 <- datos_4_50[, c("NMCpm", "p")]
g5_NMCpm_50 <- datos_5_50[, c("NMCpm", "p")]
g6_NMCpm_50 <- datos_6_50[, c("NMCpm", "p")]

# Combinar los data frames
df_NMCpm_50 <- rbind(g2_NMCpm_50, g3_NMCpm_50, g4_NMCpm_50, g5_NMCpm_50, g6_NMCpm_50)

# Ordenar los niveles de la variable n
df_NMCpm_50$p <- factor(df_NMCpm_50$p, levels = c("2", "3", "4", "5", "6"))

# Crear el gráfico de densidad con ggplot2
plot_df_NMCpm_50<-ggplot(df_NMCpm_50, aes(x = NMCpm, fill = p)) +
  geom_density(aes(color = p), alpha = 0.5, size = 1.5, show.legend = TRUE) +
  labs(title = expression("Índice " * hat(NMC)[pm] * ""), x = "Índice estimado", y = "Densidad") +
  scale_fill_manual(values = c("2" = "olivedrab3", "3" = "cornflowerblue", "4" = "turquoise", "5" = "coral", "6"="pink"),
                    name = "Cantidad de variables", labels = c("2", "3", "4", "5", "6")) +
  scale_color_manual(values = c("2" = "olivedrab3", "3" = "cornflowerblue", "4" = "turquoise", "5" = "coral", "6"="pink"),
                     name = "Cantidad de variables", labels = c("2", "3", "4", "5", "6")) +
  theme_bw()  +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 16)
  )


# Función para extraer la leyenda de un gráfico
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
# Extraer la leyenda del primer gráfico
legend <- g_legend(plot_df_MCp_p1_50)


# Quitar la leyenda de los gráficos
plot_df_MCp_p1_50 <- plot_df_MCp_p1_50 + theme(legend.position = "none")
plot_df_MCpk_p2_50 <- plot_df_MCpk_p2_50 + theme(legend.position = "none")
plot_df_MCpm_p3_50 <- plot_df_MCpm_p3_50 + theme(legend.position = "none")
plot_df_NMCp_50 <- plot_df_NMCp_50 + theme(legend.position = "none")
plot_df_NMCpm_50 <- plot_df_NMCpm_50 + theme(legend.position = "none")


# Organizar los gráficos y la leyenda en un panel, leyenda en el último espacio
grid.arrange(
  arrangeGrob(
    plot_df_MCp_p1_50, plot_df_MCpk_p2_50, plot_df_MCpm_p3_50,
    plot_df_NMCp_50, plot_df_NMCpm_50, legend,  # Leyenda en el 6to espacio
    ncol = 3  # 3 columnas
  )
)

