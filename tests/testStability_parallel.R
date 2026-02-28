library(evaluomeR)
library(RSKC)
library(sparcl)

# CASO 1

# Calculamos el tiempo secuencial
stability1_secuencial <- system.time({
  r_stability1_secuencial <- stabilityRange(ontMetrics, k.range=c(3,5), bs=50)
})

# Calculamos el tiempo paralelo
stability1_parallel_2cores <- system.time({
  r_stability1_parallel_2cores <- stabilityRange(ontMetrics, k.range=c(3,5), bs=50, numCores=2)
})

stability1_parallel_12cores <- system.time({
  r_stability1_parallel_12cores <- stabilityRange(ontMetrics, k.range=c(3,5), bs=50, numCores=12)
})

# Comprobar si los resultados son iguales
# ident1 <- identical(as.data.frame(assay(r_stability1_secuencial[[1]])), as.data.frame(assay(r_stability1_parallel_2cores[[1]])))

cat("Tiempo caso 1 secuencial:", stability1_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 1 paralelo 2 cores:", stability1_parallel_2cores["elapsed"], "seg\n")
cat("Tiempo caso 1 paralelo 12 cores:", stability1_parallel_12cores["elapsed"], "seg\n")


# CASO 2

# Calculamos el tiempo secuencial
stability2_secuencial <- system.time({
  r_stability2_secuencial <- stabilityRange(ontMetrics, k.range=c(2,10), bs=100)
})

# Calculamos el tiempo paralelo
stability2_parallel_2cores <- system.time({
  r_stability2_parallel_2cores <- stabilityRange(ontMetrics, k.range=c(2,10), bs=100, numCores=2)
})

stability2_parallel_12cores <- system.time({
  r_stability2_parallel_12cores <- stabilityRange(ontMetrics, k.range=c(2,10), bs=100, numCores=12)
})

# Comprobar si los resultados son iguales
# ident2 <- identical(as.data.frame(assay(r_stability2_secuencial[[1]])), as.data.frame(assay(r_stability2_parallel_2cores[[1]])))

cat("Tiempo caso 2 secuencial:", stability2_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 2 paralelo 2 cores:", stability2_parallel_2cores["elapsed"], "seg\n")
cat("Tiempo caso 2 paralelo 12 cores:", stability2_parallel_12cores["elapsed"], "seg\n")


# CASO 3

# Calculamos el tiempo secuencial
stability3_secuencial <- system.time({
  r_stability3_secuencial <- stabilityRange(golub, k.range=c(2,8), bs=100)
})

# Calculamos el tiempo paralelo
stability3_parallel_2cores <- system.time({
  r_stability3_parallel_2cores <- stabilityRange(golub, k.range=c(2,8), bs=100, numCores=2)
})

stability3_parallel_12cores <- system.time({
  r_stability3_parallel_12cores <- stabilityRange(golub, k.range=c(2,8), bs=100, numCores=12)
})

# Comprobar si los resultados son iguales
# ident3 <- identical(as.data.frame(assay(r_stability3_secuencial[[1]])), as.data.frame(assay(r_stability3_parallel_2cores[[1]])))

cat("Tiempo caso 3 secuencial:", stability3_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 3 paralelo 2 cores:", stability3_parallel_2cores["elapsed"], "seg\n")
cat("Tiempo caso 3 paralelo 12 cores:", stability3_parallel_12cores["elapsed"], "seg\n")


# CASO 4

# Calculamos el tiempo secuencial
stability4_secuencial <- system.time({
  r_stability4_secuencial <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100)
})

# Calculamos el tiempo paralelo
stability4_parallel_2cores <- system.time({
  r_stability4_parallel_2cores <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=2)
})

stability4_parallel_12cores <- system.time({
  r_stability4_parallel_12cores <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=12)
})

# Comprobar si los resultados son iguales
# ident4 <- identical(as.data.frame(assay(r_stability4_secuencial[[1]])), as.data.frame(assay(r_stability4_parallel_2cores[[1]])))

cat("Tiempo caso 4 secuencial:", stability4_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 4 paralelo 2 cores:", stability4_parallel_2cores["elapsed"], "seg\n")
cat("Tiempo caso 4 paralelo 12 cores:", stability4_parallel_12cores["elapsed"], "seg\n")


# CASO 5

# Calculamos el tiempo secuencial
stability5_secuencial <- system.time({
  r_stability5_secuencial <- stabilityRange(ontMetrics, k.range=c(3,5), bs=50, all_metrics=TRUE)
})

# Calculamos el tiempo paralelo
stability5_parallel_2cores <- system.time({
  r_stability5_parallel_2cores <- stabilityRange(ontMetrics, k.range=c(3,5), bs=50, all_metrics=TRUE, numCores=2)
})

stability5_parallel_12cores <- system.time({
  r_stability5_parallel_12cores <- stabilityRange(ontMetrics, k.range=c(3,5), bs=50, all_metrics=TRUE, numCores=12)
})

# Comprobar si los resultados son iguales
# ident5 <- identical(as.data.frame(assay(r_stability5_secuencial[[1]])),as.data.frame(assay(r_stability5_parallel_2cores[[1]])))

cat("Tiempo caso 5 secuencial:", stability5_secuencial["elapsed"], "seg\n")
cat("Tiempo caso 5 paralelo 2 cores:", stability5_parallel_2cores["elapsed"], "seg\n")
cat("Tiempo caso 5 paralelo 12 cores:", stability5_parallel_12cores["elapsed"], "seg\n")



# Tests para comprobar parámetro a parámetro la diferencia de rendimiento entre versiones


# NUMERO DE CORES

# Fijamos un dataset (nci60_k8), un número de iteraciones de bootstrapping (bs=100) y rango de k (3,10)

stability_1cores <- system.time({
  r_stability_1cores <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=1)
})

stability_2cores <- system.time({
  r_stability_2cores <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=2)
})

stability_4cores <- system.time({
  r_stability_4cores <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=4)
})

stability_6cores <- system.time({
  r_stability_6cores <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=6)
})

stability_8cores <- system.time({
  r_stability_8cores <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=8)
})

stability_10cores <- system.time({
  r_stability_10cores <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=10)
})

stability_12cores <- system.time({
  r_stability_12cores <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=12)
})

cat("Tiempo secuencial (1 core):", stability_1cores["elapsed"], "seg\n")
cat("Tiempo 2 cores:", stability_2cores["elapsed"], "seg\n")
cat("Tiempo 4 cores:", stability_4cores["elapsed"], "seg\n")
cat("Tiempo 6 cores:", stability_6cores["elapsed"], "seg\n")
cat("Tiempo 8 cores:", stability_8cores["elapsed"], "seg\n")
cat("Tiempo 10 cores:", stability_10cores["elapsed"], "seg\n")
cat("Tiempo 12 cores:", stability_12cores["elapsed"], "seg\n")


# NUMERO DE MÉTRICAS

# Fijamos un número de cores (12), un número de iteraciones de bootstrapping (bs=100) y rango de k (3,10)

stability_metrica1_sec <- system.time({
  r_stability_metrica1_sec <- stabilityRange(rnaMetrics, k.range=c(3,10), bs=100, numCores=1)
})

stability_metrica1_par <- system.time({
  r_stability_metrica1_par <- stabilityRange(rnaMetrics, k.range=c(3,10), bs=100, numCores=12)
})

stability_metrica2_sec <- system.time({
  r_stability_metrica2_sec <- stabilityRange(bioMetrics, k.range=c(3,10), bs=100, numCores=1)
})

stability_metrica2_par <- system.time({
  r_stability_metrica2_par <- stabilityRange(bioMetrics, k.range=c(3,10), bs=100, numCores=12)
})

stability_metrica3_sec <- system.time({
  r_stability_metrica3_sec <- stabilityRange(ontMetrics, k.range=c(3,10), bs=100, numCores=1)
})

stability_metrica3_par <- system.time({
  r_stability_metrica3_par <- stabilityRange(ontMetrics, k.range=c(3,10), bs=100, numCores=12)
})

stability_metrica4_sec <- system.time({
  r_stability_metrica4_sec <- stabilityRange(golub, k.range=c(3,10), bs=100, numCores=1)
})

stability_metrica4_par <- system.time({
  r_stability_metrica4_par <- stabilityRange(golub, k.range=c(3,10), bs=100, numCores=12)
})

stability_metrica5_sec <- system.time({
  r_stability_metrica5_sec <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=1)
})

stability_metrica5_par <- system.time({
  r_stability_metrica5_par <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=12)
})

stability_metrica6_sec <- system.time({
  r_stability_metrica6_sec <- stabilityRange(nci60_k10, k.range=c(3,10), bs=100, numCores=1)
})

stability_metrica6_par <- system.time({
  r_stability_metrica6_par <- stabilityRange(nci60_k10, k.range=c(3,10), bs=100, numCores=12)
})

cat("Tiempo rnaMetrics sec:", stability_metrica1_sec["elapsed"], "seg\n")
cat("Tiempo bioMetrics sec:", stability_metrica2_sec["elapsed"], "seg\n")
cat("Tiempo ontMetrics sec:", stability_metrica3_sec["elapsed"], "seg\n")
cat("Tiempo golub sec:", stability_metrica4_sec["elapsed"], "seg\n")
cat("Tiempo nci60_k8 sec:", stability_metrica5_sec["elapsed"], "seg\n")
cat("Tiempo nci60_k10 sec:", stability_metrica6_sec["elapsed"], "seg\n")

cat("Tiempo rnaMetrics par:", stability_metrica1_par["elapsed"], "seg\n")
cat("Tiempo bioMetrics par:", stability_metrica2_par["elapsed"], "seg\n")
cat("Tiempo ontMetrics par:", stability_metrica3_par["elapsed"], "seg\n")
cat("Tiempo golub par:", stability_metrica4_par["elapsed"], "seg\n")
cat("Tiempo nci60_k8 par:", stability_metrica5_par["elapsed"], "seg\n")
cat("Tiempo nci60_k10 par:", stability_metrica6_par["elapsed"], "seg\n")


# NÚMERO DE ITERACIONES DE BOOTSTRAPPING

# Fijamos un dataset (nci60_k8), un número de cores (12) y rango de k (3,10)

stability_10bs_sec <- system.time({
  r_stability_10bs_sec <- stabilityRange(nci60_k8, k.range=c(3,10), bs=10, numCores=1)
})

stability_10bs_par <- system.time({
  r_stability_10bs_par <- stabilityRange(nci60_k8, k.range=c(3,10), bs=10, numCores=12)
})

stability_20bs_sec <- system.time({
  r_stability_20bs_sec <- stabilityRange(nci60_k8, k.range=c(3,10), bs=20, numCores=1)
})

stability_20bs_par <- system.time({
  r_stability_20bs_par <- stabilityRange(nci60_k8, k.range=c(3,10), bs=20, numCores=12)
})

stability_50bs_sec <- system.time({
  r_stability_50bs_sec <- stabilityRange(nci60_k8, k.range=c(3,10), bs=50, numCores=1)
})

stability_50bs_par <- system.time({
  r_stability_50bs_par <- stabilityRange(nci60_k8, k.range=c(3,10), bs=50, numCores=12)
})

stability_100bs_sec <- system.time({
  r_stability_100bs_sec <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=1)
})

stability_100bs_par <- system.time({
  r_stability_100bs_par <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=12)
})

stability_150bs_sec <- system.time({
  r_stability_150bs_sec <- stabilityRange(nci60_k8, k.range=c(3,10), bs=150, numCores=1)
})

stability_150bs_par <- system.time({
  r_stability_150bs_par <- stabilityRange(nci60_k8, k.range=c(3,10), bs=150, numCores=12)
})

stability_200bs_sec <- system.time({
  r_stability_200bs_sec <- stabilityRange(nci60_k8, k.range=c(3,10), bs=200, numCores=1)
})

stability_200bs_par <- system.time({
  r_stability_200bs_par <- stabilityRange(nci60_k8, k.range=c(3,10), bs=200, numCores=12)
})

cat("Tiempo 10 iteraciones de bootstrapping sec:", stability_10bs_sec["elapsed"], "seg\n")
cat("Tiempo 20 iteraciones de bootstrapping sec:", stability_20bs_sec["elapsed"], "seg\n")
cat("Tiempo 50 iteraciones de bootstrapping sec:", stability_50bs_sec["elapsed"], "seg\n")
cat("Tiempo 100 iteraciones de bootstrapping sec:", stability_100bs_sec["elapsed"], "seg\n")
cat("Tiempo 150 iteraciones de bootstrapping sec:", stability_150bs_sec["elapsed"], "seg\n")
cat("Tiempo 200 iteraciones de bootstrapping sec:", stability_200bs_sec["elapsed"], "seg\n")

cat("Tiempo 10 iteraciones de bootstrapping par:", stability_10bs_par["elapsed"], "seg\n")
cat("Tiempo 20 iteraciones de bootstrapping par:", stability_20bs_par["elapsed"], "seg\n")
cat("Tiempo 50 iteraciones de bootstrapping par:", stability_50bs_par["elapsed"], "seg\n")
cat("Tiempo 100 iteraciones de bootstrapping par:", stability_100bs_par["elapsed"], "seg\n")
cat("Tiempo 150 iteraciones de bootstrapping par:", stability_150bs_par["elapsed"], "seg\n")
cat("Tiempo 200 iteraciones de bootstrapping par:", stability_200bs_par["elapsed"], "seg\n")


# RANGO DE K

# Fijamos un dataset (nci60_k8), un número de cores (12) y número de iteraciones de bootstrapping (bs=100)

stability_k2_3_sec <- system.time({
  r_stability_k2_3_sec <- stabilityRange(nci60_k8, k.range=c(2,3), bs=100, numCores=1)
})

stability_k2_3_par <- system.time({
  r_stability_k2_3_par <- stabilityRange(nci60_k8, k.range=c(2,3), bs=100, numCores=12)
})

stability_k2_5_sec <- system.time({
  r_stability_k2_5_sec <- stabilityRange(nci60_k8, k.range=c(2,5), bs=100, numCores=1)
})

stability_k2_5_par <- system.time({
  r_stability_k2_5_par <- stabilityRange(nci60_k8, k.range=c(2,5), bs=100, numCores=12)
})

stability_k2_8_sec <- system.time({
  r_stability_k2_8_sec <- stabilityRange(nci60_k8, k.range=c(2,8), bs=100, numCores=1)
})

stability_k2_8_par <- system.time({
  r_stability_k2_8_par <- stabilityRange(nci60_k8, k.range=c(2,8), bs=100, numCores=12)
})

stability_k2_10_sec <- system.time({
  r_stability_k2_10_sec <- stabilityRange(nci60_k8, k.range=c(2,10), bs=100, numCores=1)
})

stability_k2_10_par <- system.time({
  r_stability_k2_10_par <- stabilityRange(nci60_k8, k.range=c(2,10), bs=100, numCores=12)
})

stability_k2_12_sec <- system.time({
  r_stability_k2_12_sec <- stabilityRange(nci60_k8, k.range=c(2,12), bs=100, numCores=1)
})

stability_k2_12_par <- system.time({
  r_stability_k2_12_par <- stabilityRange(nci60_k8, k.range=c(2,12), bs=100, numCores=12)
})

stability_k2_15_sec <- system.time({
  r_stability_k2_15_sec <- stabilityRange(nci60_k8, k.range=c(2,15), bs=100, numCores=1)
})

stability_k2_15_par <- system.time({
  r_stability_k2_15_par <- stabilityRange(nci60_k8, k.range=c(2,15), bs=100, numCores=12)
})

cat("Tiempo k = (2,3) sec:", stability_k2_3_sec["elapsed"], "seg\n")
cat("Tiempo k = (2,5) sec:", stability_k2_5_sec["elapsed"], "seg\n")
cat("Tiempo k = (2,8) sec:", stability_k2_8_sec["elapsed"], "seg\n")
cat("Tiempo k = (2,10) sec:", stability_k2_10_sec["elapsed"], "seg\n")
cat("Tiempo k = (2,12) sec:", stability_k2_12_sec["elapsed"], "seg\n")
cat("Tiempo k = (2,15) sec:", stability_k2_15_sec["elapsed"], "seg\n")

cat("Tiempo k = (2,3) par:", stability_k2_3_par["elapsed"], "seg\n")
cat("Tiempo k = (2,5) par:", stability_k2_5_par["elapsed"], "seg\n")
cat("Tiempo k = (2,8) par:", stability_k2_8_par["elapsed"], "seg\n")
cat("Tiempo k = (2,10) par:", stability_k2_10_par["elapsed"], "seg\n")
cat("Tiempo k = (2,12) par:", stability_k2_12_par["elapsed"], "seg\n")
cat("Tiempo k = (2,15) par:", stability_k2_15_par["elapsed"], "seg\n")

