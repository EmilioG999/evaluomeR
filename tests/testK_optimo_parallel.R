library(evaluomeR)

seed = 100

# Casos para comprobar la rapidez de getOptimalKValue

# CASO 1

# Calculamos el tiempo de los índices
calculo1_indices <- system.time({
  stabData1 <- stabilityRange(ontMetrics, k.range=c(2,10), bs=100, numCores=12)
  qualData1 <- qualityRange(ontMetrics, k.range=c(2,10))
})

# Calculamos lo que se tarda en realizar el cálculo de k óptimo una vez obtenidos los índices
calculo1_kopt <- system.time({
  kopt1 <- getOptimalKValue(stabData1, qualData1, k.range=c(2,10))
})

cat("Tiempo de cálculo 1 de índices:", calculo1_indices["elapsed"], "seg\n")
cat("Tiempo de cálculo 1 de k óptimo:", calculo1_kopt["elapsed"], "seg\n")


# CASO 2

# Calculamos el tiempo de los índices
calculo2_indices <- system.time({
  stabData2 <- stabilityRange(nci60_k8, k.range=c(2,6), bs=100, numCores=12)
  qualData2 <- qualityRange(nci60_k8, k.range=c(2,6))
})

# Calculamos lo que se tarda en realizar el cálculo de k óptimo una vez obtenidos los índices
calculo2_kopt <- system.time({
  kopt2 <- getOptimalKValue(stabData2, qualData2, k.range=c(2,6))
})

cat("Tiempo de cálculo 2 de índices:", calculo2_indices["elapsed"], "seg\n")
cat("Tiempo de cálculo 2 de k óptimo:", calculo2_kopt["elapsed"], "seg\n")


# CASO 3

# Calculamos el tiempo de los índices
calculo3_indices <- system.time({
  stabData3 <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=12)
  qualData3 <- qualityRange(nci60_k8, k.range=c(3,10))
})

# Calculamos lo que se tarda en realizar el cálculo de k óptimo una vez obtenidos los índices
calculo3_kopt <- system.time({
  kopt3 <- getOptimalKValue(stabData3, qualData3, k.range=c(3,10))
})

cat("Tiempo de cálculo 3 de índices:", calculo3_indices["elapsed"], "seg\n")
cat("Tiempo de cálculo 3 de k óptimo:", calculo3_kopt["elapsed"], "seg\n")



# Comparativa cálculo completo k óptimo secuencial y paralelo

# CASO 1

# Calculo secuencial de k óptimo
kopt1_secuencial <- system.time({
  stabData1 <- stabilityRange(ontMetrics, k.range=c(2,5), bs=50)
  qualData1 <- qualityRange(ontMetrics, k.range=c(2,5))
  kopt1_sec <- getOptimalKValue(stabData1, qualData1, k.range=c(2,5))
})

# Calculo paralelo de k óptimo
kopt1_paralelo <- system.time({
  stabData1 <- stabilityRange(ontMetrics, k.range=c(2,5), bs=50, numCores=6)
  qualData1 <- qualityRange(ontMetrics, k.range=c(2,5))
  kopt1_par <- getOptimalKValue(stabData1, qualData1, k.range=c(2,5))
})

cat("Tiempo secuencial k_optimo caso 1:", kopt1_secuencial["elapsed"], "seg\n")
cat("Tiempo paralelo k_optimo caso 1:", kopt1_paralelo["elapsed"], "seg\n")


# CASO 2

# Calculo secuencial de k óptimo
kopt2_secuencial <- system.time({
  stabData2 <- stabilityRange(ontMetrics, k.range=c(2,10), bs=100)
  qualData2 <- qualityRange(ontMetrics, k.range=c(2,10))
  kopt2_sec <- getOptimalKValue(stabData2, qualData2, k.range=c(2,10))
})

# Calculo paralelo de k óptimo
kopt2_paralelo <- system.time({
  stabData2 <- stabilityRange(ontMetrics, k.range=c(2,10), bs=100, numCores=6)
  qualData2 <- qualityRange(ontMetrics, k.range=c(2,10))
  kopt2_par <- getOptimalKValue(stabData2, qualData2, k.range=c(2,10))
})

cat("Tiempo secuencial k_optimo caso 2:", kopt2_secuencial["elapsed"], "seg\n")
cat("Tiempo paralelo k_optimo caso 2:", kopt2_paralelo["elapsed"], "seg\n")


# CASO 3

# Calculo secuencial de k óptimo
kopt3_secuencial <- system.time({
  stabData3 <- stabilityRange(nci60_k8, k.range=c(2,10), bs=100)
  qualData3 <- qualityRange(nci60_k8, k.range=c(2,10))
  kopt3_sec <- getOptimalKValue(stabData3, qualData3, k.range=c(2,10))
})

# Calculo paralelo de k óptimo
kopt3_paralelo <- system.time({
  stabData3 <- stabilityRange(nci60_k8, k.range=c(2,10), bs=100, numCores=6)
  qualData3 <- qualityRange(nci60_k8, k.range=c(2,10))
  kopt3_par <- getOptimalKValue(stabData3, qualData3, k.range=c(2,10))
})

cat("Tiempo secuencial k_optimo caso 3:", kopt3_secuencial["elapsed"], "seg\n")
cat("Tiempo paralelo k_optimo caso 3:", kopt3_paralelo["elapsed"], "seg\n")


# CASO 4

# Calculo secuencial de k óptimo
kopt4_secuencial <- system.time({
  stabData4 <- stabilityRange(nci60_k8, k.range=c(2,15), bs=150)
  qualData4 <- qualityRange(nci60_k8, k.range=c(2,15))
  kopt4_sec <- getOptimalKValue(stabData4, qualData4, k.range=c(2,15))
})

# Calculo paralelo de k óptimo
kopt4_paralelo <- system.time({
  stabData4 <- stabilityRange(nci60_k8, k.range=c(2,15), bs=150, numCores=6)
  qualData4 <- qualityRange(nci60_k8, k.range=c(2,15))
  kopt4_par <- getOptimalKValue(stabData4, qualData4, k.range=c(2,15))
})

cat("Tiempo secuencial k_optimo caso 4:", kopt4_secuencial["elapsed"], "seg\n")
cat("Tiempo paralelo k_optimo caso 4:", kopt4_paralelo["elapsed"], "seg\n")


# Comprobamos si los resultados son iguales
# ident1 <- identical(kopt1_sec, kopt1_par)
# ident2 <- identical(kopt2_sec, kopt2_par)
# ident3 <- identical(kopt3_sec, kopt3_par)
# ident4 <- identical(kopt4_sec, kopt4_par)

# cat("Ident 1:", ident1)
# cat("Ident 2:", ident2)
# cat("Ident 3:", ident3)
# cat("Ident 4:", ident4)

