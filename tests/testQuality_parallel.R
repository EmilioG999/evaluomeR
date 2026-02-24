library(evaluomeR)
library(RSKC)
library(sparcl)
seed=100

# CASO 1

# Calculamos el tiempo secuencial
quality1_secuencial <- system.time({
  r_quality1_secuencial <- qualityRange(ontMetrics, k.range=c(3,5))
})

# Calculamos el tiempo paralelo
quality1_parallel <- system.time({
  r_quality1_parallel <- qualityRange_parallel(ontMetrics, k.range=c(3,5))
})

# Comprobar si los resultados son iguales
# ident1 <- identical(as.data.frame(assay(r_quality1_secuencial[[1]])), as.data.frame(assay(r_quality1_parallel[[1]])))

cat("Tiempo secuencial:", quality1_secuencial["elapsed"], "seg\n")
cat("Tiempo paralelo:", quality1_parallel["elapsed"], "seg\n")


# CASO 2

# Calculamos el tiempo secuencial
quality2_secuencial <- system.time({
  r_quality2_secuencial <- qualityRange(ontMetrics, k.range=c(2,15))
})

# Calculamos el tiempo paralelo
quality2_parallel <- system.time({
  r_quality2_parallel <- qualityRange_parallel(ontMetrics, k.range=c(2,15))
})

# Comprobar si los resultados son iguales
# ident2 <- identical(as.data.frame(assay(r_quality2_secuencial[[1]])), as.data.frame(assay(r_quality2_parallel[[1]])))

cat("Tiempo secuencial:", quality2_secuencial["elapsed"], "seg\n")
cat("Tiempo paralelo:", quality2_parallel["elapsed"], "seg\n")


# CASO 3

# Calculamos el tiempo secuencial
quality3_secuencial <- system.time({
  r_quality3_secuencial <- qualityRange(nci60_k8, k.range=c(2,8))
})

# Calculamos el tiempo paralelo
quality3_parallel <- system.time({
  r_quality3_parallel <- qualityRange_parallel(nci60_k8, k.range=c(2,8))
})

# Comprobar si los resultados son iguales
# ident3 <- identical(as.data.frame(assay(r_quality3_secuencial[[1]])), as.data.frame(assay(r_quality3_parallel[[1]])))

cat("Tiempo secuencial:", quality3_secuencial["elapsed"], "seg\n")
cat("Tiempo paralelo:", quality3_parallel["elapsed"], "seg\n")


# CASO 4

# Calculamos el tiempo secuencial
quality4_secuencial <- system.time({
  r_quality4_secuencial <- qualityRange(nci60_k8, k.range=c(3,15))
})

# Calculamos el tiempo paralelo
quality4_parallel <- system.time({
  r_quality4_parallel <- qualityRange_parallel(nci60_k8, k.range=c(3,15))
})

# Comprobar si los resultados son iguales
# ident4 <- identical(as.data.frame(assay(r_quality4_secuencial[[1]])), as.data.frame(assay(r_quality4_parallel[[1]])))

cat("Tiempo secuencial:", quality4_secuencial["elapsed"], "seg\n")
cat("Tiempo paralelo:", quality4_parallel["elapsed"], "seg\n")


# CASO 5

# Calculamos el tiempo secuencial
quality5_secuencial <- system.time({
  r_quality5_secuencial <- qualityRange(ontMetrics, k.range=c(3,5), all_metrics=TRUE)
})

# Calculamos el tiempo paralelo
quality5_parallel <- system.time({
  r_quality5_parallel <- qualityRange_parallel(ontMetrics, k.range=c(3,5), all_metrics=TRUE)
})

# Comprobar si los resultados son iguales
# ident5 <- identical(as.data.frame(assay(r_quality5_secuencial[[1]])), as.data.frame(assay(r_quality5_parallel[[1]])))

cat("Tiempo secuencial:", quality5_secuencial["elapsed"], "seg\n")
cat("Tiempo paralelo:", quality5_parallel["elapsed"], "seg\n")

