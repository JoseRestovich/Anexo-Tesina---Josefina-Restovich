library(ggplot2)
library(tidyverse)
library(mvtnorm)
library(MASS)
library(dplyr)
library(readxl)
library(gridExtra)


cantidad_muestras <- 5000
set.seed(15)

# Carga de datos
datos_multi <- read_excel("C:/Users/josefina.restovich/OneDrive - datalytics.com/Escritorio/Facultad/Tesina/Programas/Multivariado/Practica análisis de capacidad de procesos multivariado - Conjuntos de datos.xlsx", sheet = "Mecánico")
datos_multi <- datos_multi[, c(2, 3, 4, 5, 6, 7, 8)]


#################
### p=6, n= 50 ###
#################
datos_6 <- datos_multi[, c(4, 7, 1, 6, 3, 2)]
#Parámetros
LIE_6 <- c(3.00, 118.00, 5.00, 37.00, 2.50, 33.00)
LSE_6 <- c(17.00, 122.00, 15.00, 43.00, 6.50, 37.00)
t_6 <- c(10, 120, 10, 40, 5, 35)
medias_6 <- colMeans(datos_6)
sigma_6 <- cov(datos_6)
# Tamaño de muestra
tamano_muestra <- 50
#Cálculo de las probabilidades 
p1_6 <- pmvnorm(LIE_6, LSE_6, mean = t_6, sigma = sigma_6)
p2_6 <- pmvnorm(LIE_6, LSE_6, mean = medias_6, sigma = sigma_6)
sigmat_6 <- sigma_6 + (medias_6 - t_6) %*% t(medias_6 - t_6)
p3_6 <- pmvnorm(LIE_6, LSE_6, mean = t_6, sigma = sigmat_6)
#Cálculo de los índices
MCp_p1_6 <- (-1/3 * qnorm((1 - p1_6)/2))
MCpk_p2_6 <- (-1/3 * qnorm((1 - p2_6)/2))
MCpm_p3_6 <- (-1/3 * qnorm((1 - p3_6)/2))

# Cálculo de los valores reales de NMCp y NMCpm
cprima_6 <- min(((LSE_6[1] - t_6[1]) / sqrt(sigma_6[1, 1])), ((LSE_6[2] - t_6[2]) / sqrt(sigma_6[2, 2])), ((LSE_6[3] - t_6[3]) / sqrt(sigma_6[3, 3])), ((LSE_6[4] - t_6[4]) / sqrt(sigma_6[4, 4])), ((LSE_6[5] - t_6[5]) / sqrt(sigma_6[5, 5])), ((LSE_6[6] - t_6[6]) / sqrt(sigma_6[6, 6])))
NMCp_6 <- cprima_6 / sqrt(qchisq(1 - 0.0027, 6))

# Cálculo de la distancia D
D_6 <- (1 + (t(medias_6 - t_6) %*% solve(sigma_6) %*% (medias_6 - t_6)))^0.5
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
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 6))  # Usar alfa = 0.0027
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
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 6))  # Usar alfa = 0.0027
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
  NMCp <- cprima / sqrt(qchisq(1 - 0.0027, 6))  # Usar alfa = 0.0027
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

############################
##########GRÁFICOS##########
############################

#########caso 6
datos_6_50$n<-"50"
datos_6_100$n<-"100"
datos_6_200$n<-"200"
colnames(datos_6_50)[which(colnames(datos_6_50) == "MCp_p1_6_50")] <- "MCp_p1"
colnames(datos_6_100)[which(colnames(datos_6_100) == "MCp_p1_6_100")] <- "MCp_p1"
colnames(datos_6_200)[which(colnames(datos_6_200) == "MCp_p1_6_200")] <- "MCp_p1"

colnames(datos_6_50)[which(colnames(datos_6_50) == "MCpk_p2_6_50")] <- "MCpk_p2"
colnames(datos_6_100)[which(colnames(datos_6_100) == "MCpk_p2_6_100")] <- "MCpk_p2"
colnames(datos_6_200)[which(colnames(datos_6_200) == "MCpk_p2_6_200")] <- "MCpk_p2"

colnames(datos_6_50)[which(colnames(datos_6_50) == "MCpm_p3_6_50")] <- "MCpm_p3"
colnames(datos_6_100)[which(colnames(datos_6_100) == "MCpm_p3_6_100")] <- "MCpm_p3"
colnames(datos_6_200)[which(colnames(datos_6_200) == "MCpm_p3_6_200")] <- "MCpm_p3"

##############MCp_p1_6
# Definir los dataframes individuales
g50_MCp_p1 <- datos_6_50[, c("MCp_p1", "n")]
g100_MCp_p1 <- datos_6_100[, c("MCp_p1", "n")]
g200_MCp_p1 <- datos_6_200[, c("MCp_p1", "n")]


# Combinar los dataframes
df_MCp_p1_6 <- rbind(g50_MCp_p1, g100_MCp_p1, g200_MCp_p1)

# Ordenar los niveles de la variable n
df_MCp_p1_6$n <- factor(df_MCp_p1_6$n, levels = c("50", "100", "200"))

##############MCpk_p2_6
# Definir los dataframes individuales
g50_MCpk_p2 <- datos_6_50[, c("MCpk_p2", "n")]
g100_MCpk_p2 <- datos_6_100[, c("MCpk_p2", "n")]
g200_MCpk_p2 <- datos_6_200[, c("MCpk_p2", "n")]


# Combinar los dataframes
df_MCpk_p2_6 <- rbind(g50_MCpk_p2, g100_MCpk_p2, g200_MCpk_p2)

# Ordenar los niveles de la variable n
df_MCpk_p2_6$n <- factor(df_MCpk_p2_6$n, levels = c("50", "100", "200"))

##############MCpm_p3_6
# Definir los dataframes individuales
g50_MCpm_p3 <- datos_6_50[, c("MCpm_p3", "n")]
g100_MCpm_p3 <- datos_6_100[, c("MCpm_p3", "n")]
g200_MCpm_p3 <- datos_6_200[, c("MCpm_p3", "n")]


# Combinar los dataframes
df_MCpm_p3_6 <- rbind(g50_MCpm_p3, g100_MCpm_p3, g200_MCpm_p3)

# Ordenar los niveles de la variable n
df_MCpm_p3_6$n <- factor(df_MCpm_p3_6$n, levels = c("50", "100", "200"))

##############NMCp_6
# Definir los dataframes individuales
g50_NMCp <- datos_6_50[, c("NMCp", "n")]
g100_NMCp <- datos_6_100[, c("NMCp", "n")]
g200_NMCp <- datos_6_200[, c("NMCp", "n")]


# Combinar los dataframes
df_NMCp_6 <- rbind(g50_NMCp, g100_NMCp, g200_NMCp)

# Ordenar los niveles de la variable n
df_NMCp_6$n <- factor(df_NMCp_6$n, levels = c("50", "100", "200"))

##############NMCpm_6
# Definir los dataframes individuales
g50_NMCpm <- datos_6_50[, c("NMCpm", "n")]
g100_NMCpm <- datos_6_100[, c("NMCpm", "n")]
g200_NMCpm <- datos_6_200[, c("NMCpm", "n")]


# Combinar los dataframes
df_NMCpm_6 <- rbind(g50_NMCpm, g100_NMCpm, g200_NMCpm)

# Ordenar los niveles de la variable n
df_NMCpm_6$n <- factor(df_NMCpm_6$n, levels = c("50", "100", "200"))

# Crear los gráficos individuales con la leyenda activada
plot_MCp_p1_6 <- ggplot(df_MCp_p1_6, aes(x = MCp_p1, fill = n)) +
  geom_density(aes(color = n), alpha = 0.5, size = 1.5) +
  labs(title = expression("Índice" * hat(MC)[p(p[1])] * ""), x = "Índice estimado", y = "Densidad") +
  scale_fill_manual(values = c("50" = "olivedrab3", "100" = "cornflowerblue", "200" = "turquoise"),
                    name = "Tamaño de muestra", labels = c("50", "100", "200")) +
  scale_color_manual(values = c("50" = "olivedrab3", "100" = "cornflowerblue", "200" = "turquoise"),
                     name = "Tamaño de muestra", labels = c("50", "100", "200")) +
  theme_bw() +
  geom_vline(xintercept = MCp_p1_6, linetype = "dashed", color = "black") +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 16)
  )

plot_MCpk_p2_6 <- ggplot(df_MCpk_p2_6, aes(x = MCpk_p2, fill = n)) +
  geom_density(aes(color = n), alpha = 0.5, size = 1.5) +
  labs(title = expression("Índice " * hat(MC)[pk(p[2])] * ""), x = "Índice estimado", y = "Densidad") +
  scale_fill_manual(values = c("50" = "olivedrab3", "100" = "cornflowerblue", "200" = "turquoise"),
                    name = "Tamaño de muestra", labels = c("50", "100", "200")) +
  scale_color_manual(values = c("50" = "olivedrab3", "100" = "cornflowerblue", "200" = "turquoise"),
                     name = "Tamaño de muestra", labels = c("50", "100", "200")) +
  theme_bw() +
  geom_vline(xintercept = MCpk_p2_6, linetype = "dashed", color = "black") +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 16)
  )

plot_MCpm_p3_6 <- ggplot(df_MCpm_p3_6, aes(x = MCpm_p3, fill = n)) +
  geom_density(aes(color = n), alpha = 0.5, size = 1.5) +
  labs(title = expression("Índice " * hat(MC)[pm(p[3])] * ""), x = "Índice estimado", y = "Densidad") +
  scale_fill_manual(values = c("50" = "olivedrab3", "100" = "cornflowerblue", "200" = "turquoise"),
                    name = "Tamaño de muestra", labels = c("50", "100", "200")) +
  scale_color_manual(values = c("50" = "olivedrab3", "100" = "cornflowerblue", "200" = "turquoise"),
                     name = "Tamaño de muestra", labels = c("50", "100", "200")) +
  theme_bw() +
  geom_vline(xintercept = MCpm_p3_6, linetype = "dashed", color = "black") +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 16)
  )

plot_NMCp_6 <- ggplot(df_NMCp_6, aes(x = NMCp, fill = n)) +
  geom_density(aes(color = n), alpha = 0.5, size = 1.5) +
  labs(title = expression("Índice " * hat(NMC)[p] * ""), x = "Índice estimado", y = "Densidad") +
  scale_fill_manual(values = c("50" = "olivedrab3", "100" = "cornflowerblue", "200" = "turquoise"),
                    name = "Tamaño de muestra", labels = c("50", "100", "200")) +
  scale_color_manual(values = c("50" = "olivedrab3", "100" = "cornflowerblue", "200" = "turquoise"),
                     name = "Tamaño de muestra", labels = c("50", "100", "200")) +
  theme_bw() +
  geom_vline(xintercept = NMCp_6, linetype = "dashed", color = "black") +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 16)
  )

plot_NMCpm_6 <- ggplot(df_NMCpm_6, aes(x = NMCpm, fill = n)) +
  geom_density(aes(color = n), alpha = 0.5, size = 1.5) +
  labs(title = expression("Índicd " * hat(NMC)[pm] * ""), x = "Índice estimado", y = "Densidad") +
  scale_fill_manual(values = c("50" = "olivedrab3", "100" = "cornflowerblue", "200" = "turquoise"),
                    name = "Tamaño de muestra", labels = c("50", "100", "200")) +
  scale_color_manual(values = c("50" = "olivedrab3", "100" = "cornflowerblue", "200" = "turquoise"),
                     name = "Tamaño de muestra", labels = c("50", "100", "200")) +
  theme_bw() +
  geom_vline(xintercept = NMCpm_6, linetype = "dashed", color = "black") +
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
legend <- g_legend(plot_MCp_p1_6)

# Quitar la leyenda de los gráficos
plot_MCp_p1_6 <- plot_MCp_p1_6 + theme(legend.position = "none")
plot_MCpk_p2_6 <- plot_MCpk_p2_6 + theme(legend.position = "none")
plot_MCpm_p3_6 <- plot_MCpm_p3_6 + theme(legend.position = "none")
plot_NMCp_6 <- plot_NMCp_6 + theme(legend.position = "none")
plot_NMCpm_6 <- plot_NMCpm_6 + theme(legend.position = "none")


# Organizar los gráficos y la leyenda en un panel, leyenda en el último espacio
grid.arrange(
  arrangeGrob(
    plot_MCp_p1_6, plot_MCpk_p2_6, plot_MCpm_p3_6,
    plot_NMCp_6, plot_NMCpm_6, legend,  # Leyenda en el 6to espacio
    ncol = 3  # 3 columnas
  )
)

