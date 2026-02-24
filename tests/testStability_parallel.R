library(evaluomeR)
library(RSKC)
library(sparcl)

# CASO 1

# Calculamos el tiempo secuencial
stability1_secuencial <- system.time({
  r_stability1_secuencial <- stabilityRange(ontMetrics, k.range=c(3,5), bs=50)
})

# Calculamos el tiempo paralelo
stability1_parallel <- system.time({
  r_stability1_parallel <- stabilityRange_parallel(ontMetrics, k.range=c(3,5), bs=50)
})

# Comprobar si los resultados son iguales
# ident1 <- identical(as.data.frame(assay(r_stability1_secuencial[[1]])), as.data.frame(assay(r_stability1_parallel[[1]])))

cat("Tiempo caso 1 secuencial:", stability1_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 1 paralelo:", stability1_parallel["elapsed"], "seg\n")


# CASO 2

# Calculamos el tiempo secuencial
stability2_secuencial <- system.time({
  r_stability2_secuencial <- stabilityRange(ontMetrics, k.range=c(2,10), bs=100)
})

# Calculamos el tiempo paralelo
stability2_parallel <- system.time({
  r_stability2_parallel <- stabilityRange_parallel(ontMetrics, k.range=c(2,10), bs=100)
})

# Comprobar si los resultados son iguales
# ident2 <- identical(as.data.frame(assay(r_stability2_secuencial[[1]])), as.data.frame(assay(r_stability2_parallel[[1]])))

cat("Tiempo caso 2 secuencial:", stability2_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 2 paralelo:", stability2_parallel["elapsed"], "seg\n")

# CASO 3

# Calculamos el tiempo secuencial
stability3_secuencial <- system.time({
  r_stability3_secuencial <- stabilityRange(golub, k.range=c(2,8), bs=100)
})

# Calculamos el tiempo paralelo
stability3_parallel <- system.time({
  r_stability3_parallel <- stabilityRange_parallel(golub, k.range=c(2,8), bs=100)
})

# Comprobar si los resultados son iguales
# ident3 <- identical(as.data.frame(assay(r_stability3_secuencial[[1]])), as.data.frame(assay(r_stability3_parallel[[1]])))

cat("Tiempo caso 3 secuencial:", stability3_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 3 paralelo:", stability3_parallel["elapsed"], "seg\n")


# CASO 4

# Calculamos el tiempo secuencial
stability4_secuencial <- system.time({
  r_stability4_secuencial <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100)
})

# Calculamos el tiempo paralelo
stability4_parallel <- system.time({
  r_stability4_parallel <- stabilityRange_parallel(nci60_k8, k.range=c(3,10), bs=100)
})

# Comprobar si los resultados son iguales
# ident4 <- identical(as.data.frame(assay(r_stability4_secuencial[[1]])), as.data.frame(assay(r_stability4_parallel[[1]])))

cat("Tiempo caso 4 secuencial:", stability4_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 4 paralelo:", stability4_parallel["elapsed"], "seg\n")


# CASO 5

# Calculamos el tiempo secuencial
stability5_secuencial <- system.time({
  r_stability5_secuencial <- stabilityRange(ontMetrics, k.range=c(3,5), bs=50, all_metrics=TRUE)
})

# Calculamos el tiempo paralelo
stability5_parallel <- system.time({
  r_stability5_parallel <- stabilityRange_parallel(ontMetrics, k.range=c(3,5), bs=50, all_metrics=TRUE)
})

# Comprobar si los resultados son iguales
# ident5 <- identical(as.data.frame(assay(r_stability5_secuencial[[1]])),as.data.frame(assay(r_stability5_parallel[[1]])))

cat("Tiempo caso 5 secuencial:", stability5_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 5 paralelo:", stability5_parallel["elapsed"], "seg\n")

