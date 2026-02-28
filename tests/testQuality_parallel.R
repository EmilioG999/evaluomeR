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
quality1_parallel_2cores <- system.time({
  r_quality1_parallel_2cores <- qualityRange(ontMetrics, k.range=c(3,5), numCores=2)
})

quality1_parallel_12cores <- system.time({
  r_quality1_parallel_12cores <- qualityRange(ontMetrics, k.range=c(3,5), numCores=12)
})

# Comprobar si los resultados son iguales
# ident1 <- identical(as.data.frame(assay(r_quality1_secuencial[[1]])), as.data.frame(assay(r_quality1_parallel_2cores[[1]])))

cat("Tiempo caso 1 secuencial:", quality1_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 1 paralelo 2 cores:", quality1_parallel_2cores["elapsed"], "seg\n")
cat("Tiempo caso 1 paralelo 12 cores:", quality1_parallel_12cores["elapsed"], "seg\n")


# CASO 2

# Calculamos el tiempo secuencial
quality2_secuencial <- system.time({
  r_quality2_secuencial <- qualityRange(ontMetrics, k.range=c(2,15))
})

# Calculamos el tiempo paralelo
quality2_parallel_2cores <- system.time({
  r_quality2_parallel_2cores <- qualityRange(ontMetrics, k.range=c(2,15), numCores=2)
})

quality2_parallel_12cores <- system.time({
  r_quality2_parallel_12cores <- qualityRange(ontMetrics, k.range=c(2,15), numCores=12)
})

# Comprobar si los resultados son iguales
# ident2 <- identical(as.data.frame(assay(r_quality2_secuencial[[1]])), as.data.frame(assay(r_quality2_parallel_2cores[[1]])))

cat("Tiempo caso 2 secuencial:", quality2_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 2 paralelo 2 cores:", quality2_parallel_2cores["elapsed"], "seg\n")
cat("Tiempo caso 2 paralelo 12 cores:", quality2_parallel_12cores["elapsed"], "seg\n")

# CASO 3

# Calculamos el tiempo secuencial
quality3_secuencial <- system.time({
  r_quality3_secuencial <- qualityRange(nci60_k8, k.range=c(2,8))
})

# Calculamos el tiempo paralelo
quality3_parallel_2cores <- system.time({
  r_quality3_parallel_2cores <- qualityRange(nci60_k8, k.range=c(2,8), numCores=2)
})

quality3_parallel_12cores <- system.time({
  r_quality3_parallel_12cores <- qualityRange(nci60_k8, k.range=c(2,8), numCores=12)
})

# Comprobar si los resultados son iguales
# ident3 <- identical(as.data.frame(assay(r_quality3_secuencial[[1]])), as.data.frame(assay(r_quality3_parallel_2cores[[1]])))

cat("Tiempo caso 3 secuencial:", quality3_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 3 paralelo 2 cores:", quality3_parallel_2cores["elapsed"], "seg\n")
cat("Tiempo caso 3 paralelo 12 cores:", quality3_parallel_12cores["elapsed"], "seg\n")


# CASO 4

# Calculamos el tiempo secuencial
quality4_secuencial <- system.time({
  r_quality4_secuencial <- qualityRange(nci60_k8, k.range=c(3,15))
})

# Calculamos el tiempo paralelo
quality4_parallel_2cores <- system.time({
  r_quality4_parallel_2cores <- qualityRange(nci60_k8, k.range=c(3,15), numCores=2)
})

quality4_parallel_12cores <- system.time({
  r_quality4_parallel_12cores <- qualityRange(nci60_k8, k.range=c(3,15), numCores=12)
})

# Comprobar si los resultados son iguales
# ident4 <- identical(as.data.frame(assay(r_quality4_secuencial[[1]])), as.data.frame(assay(r_quality4_parallel_2cores[[1]])))

cat("Tiempo caso 4 secuencial:", quality4_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 4 paralelo 2 cores:", quality4_parallel_2cores["elapsed"], "seg\n")
cat("Tiempo caso 4 paralelo 12 cores:", quality4_parallel_12cores["elapsed"], "seg\n")


# CASO 5

# Calculamos el tiempo secuencial
quality5_secuencial <- system.time({
  r_quality5_secuencial <- qualityRange(ontMetrics, k.range=c(3,5), all_metrics=TRUE)
})

# Calculamos el tiempo paralelo
quality5_parallel_2cores <- system.time({
  r_quality5_parallel_2cores <- qualityRange(ontMetrics, k.range=c(3,5), all_metrics=TRUE, numCores=2)
})

quality5_parallel_12cores <- system.time({
  r_quality5_parallel_12cores <- qualityRange(ontMetrics, k.range=c(3,5), all_metrics=TRUE, numCores=12)
})


# Comprobar si los resultados son iguales
# ident5 <- identical(as.data.frame(assay(r_quality5_secuencial[[1]])), as.data.frame(assay(r_quality5_parallel_2cores[[1]])))

cat("Tiempo caso 5 secuencial:", quality5_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 5 paralelo 2 cores:", quality5_parallel_2cores["elapsed"], "seg\n")
cat("Tiempo caso 5 paralelo 12 cores:", quality5_parallel_12cores["elapsed"], "seg\n")

