# Librería para paralelizar
library(parallel)

# Versión paralela de la función stability
stability_parallel <- function(data, k=5, bs=100, cbi="kmeans",
                      getImages=FALSE, all_metrics=FALSE, seed=NULL,
                      gold_standard=NULL,...) {
  
  data <- as.data.frame(assay(data))
  
  checkKValue(k)
  runStabilityIndex_parallel(data, k.min=k, k.max=k, bs=bs, cbi=cbi,
                    all_metrics=all_metrics, gold_standard=gold_standard, seed=seed, ...)
  stabilityDataFrame <- suppressWarnings(
    runStabilityIndexTableRange(data, k.min=k, k.max=k))
  if (getImages == TRUE) {
    suppressWarnings(
      runStabilityIndexK_IMG(bs, k.min = k, k.max = k))
  }
  se <- createSEList(stabilityDataFrame)
  return(se)
}

# Versión paralela de la función stabilityRange
stabilityRange_parallel <- function(data, k.range=c(2,15), bs=100, cbi="kmeans",
                                    getImages=FALSE, all_metrics=FALSE, seed=NULL,
                                    gold_standard=NULL,...) {
  k.range.length = length(k.range)
  if (k.range.length != 2) {
    stop("k.range length must be 2")
  }
  k.min = k.range[1]
  k.max = k.range[2]
  checkKValue(k.min)
  checkKValue(k.max)
  if (k.max < k.min) {
    stop("The first value of k.range cannot be greater than its second value")
  }
  data <- as.data.frame(SummarizedExperiment::assay(data))
  
  runStabilityIndex_parallel(data, k.min=k.min, k.max=k.max, bs, cbi, all_metrics=all_metrics, gold_standard=gold_standard, seed=seed, ...)
  
  stabilityDataFrame <- suppressWarnings(
    runStabilityIndexTableRange(data, k.min=k.min, k.max=k.max))
  
  if (getImages == TRUE) {
    suppressWarnings(
      runStabilityIndexK_IMG(bs, k.min=k.min, k.max=k.max))
    suppressWarnings(
      runStabilityIndexMetric_IMG(bs, k.min=k.min, k.max=k.max))
  }
  se <- createSEList(stabilityDataFrame)
  return(se)
}

# Versión paralela de la función stabilitySet
stabilitySet_parallel <- function(data, k.set=c(2,3), bs=100, cbi="kmeans",
                         getImages=FALSE, all_metrics=FALSE, seed=NULL,
                         gold_standard=NULL, ...) {
  k.set.length = length(k.set)
  if (k.set.length == 0) {
    stop("k.set list is empty")
  } else if (k.set.length == 1) {
    stop("k.set list contains only one element. For one K analysis use 'stability' method")
  }
  k.set = sort(k.set)
  for (k in k.set) {
    checkKValue(k)
  }
  
  data <- as.data.frame(SummarizedExperiment::assay(data))
  
  runStabilityIndex_parallel(data, k.set = k.set, bs=bs, cbi=cbi, all_metrics=all_metrics, gold_standard=gold_standard, seed=seed, ...)
  stabilityDataFrame <- suppressWarnings(
    runStabilityIndexTableRange(data, k.set = k.set))
  
  if (getImages == TRUE) {
    suppressWarnings(
      runStabilityIndexK_IMG(bs, k.set = k.set))
    suppressWarnings(
      runStabilityIndexMetric_IMG(bs, k.set = k.set))
  }
  se <- createSEList(stabilityDataFrame)
  return(se)
}

# Versión paralela de la función runStabilityIndex
# Función paralelizada que divide el número de métricas entre un número de clusters (añadimos numCores)
runStabilityIndex_parallel <- function(data, k.min=NULL, k.max=NULL, bs,
                                       cbi, all_metrics, seed, k.set=NULL,
                                       gold_standard=NULL, numCores=NULL,...) {
  
  if (is.null(seed)) {
    seed = pkg.env$seed
  }
  if (is.null(k.min) && is.null(k.max) && is.null(k.set)) {
    stop("runStabilityIndex: All k parameters are null!")
  }
  if (!is.null(gold_standard)) {
    if (all_metrics==FALSE) {
      stop("Gold standard parameter can be set only if the clustering of all the metrics is selected (all_metrics = TRUE)")
    }
    if (bs != 0) {
      message("Warning: 'gold_standard' parameter is set, argument 'bs' will be ignored.")
    }
  }
  
  data <- removeNAValues(data)
  dfStats(data)
  
  inversa=NULL
  m.stab.global = NULL
  m.stab.global.csv = NULL # To store new CSV output measures without altering legacy code
  todo.estable = NULL
  datos.bruto=data
  names.metr=names(datos.bruto)[-c(1)]
  
  pkg.env$names.metr = names.metr
  
  bs.values=c(bs)
  contador=0
  i.min=k.min
  i.max=k.max
  
  k.range = NULL
  k.range.length = NULL
  if (!is.null(k.set)) {
    k.range = k.set
    k.range.length = length(k.set)
  } else {
    k.range = i.min:i.max
    k.range.length = length(i.min:i.max)+1
  }
  
  if (all_metrics == TRUE) { # Processing all metrics as one
    num.metrics = 1
    pkg.env$names.metr = c("all_metrics")
  } else {
    num.metrics = length(names.metr)
    pkg.env$names.metr = names.metr
  }
  
  # Detectamos el número de cores, también se puede asignar manualmente aquí si se prefiere
  numCores <- detectCores() 
  message("Número de cores que participan en la paralelización: ", numCores)
  cl <- makeCluster(numCores)
  on.exit(stopCluster(cl), add = TRUE)
  
  # A partir del entorno actual de la función busca esas variables y las pasa a cl
  clusterExport(cl, c("datos.bruto", "names.metr", "cbi", "seed", "k.range", "k.range.length", "bs.values", "gold_standard", "all_metrics"), envir=environment())
  
  # Cargamos el paquete de evaluomeR
  clusterEvalQ(cl, {library(evaluomeR)})
  
  # Bucle a paralelizar, número de métricas
  resultados <- parLapply(cl, 1:num.metrics, function(i.metr){
    
    # Guardamos los resultados de una métrica específica
    stab_valores <- rep(NA, max(k.range))
    
    csv_data <- vector("list", max(k.range))
    
    
    if (all_metrics == TRUE) {
      message("Processing all metrics, 'merge', in dataframe (", length(names.metr),")")
    } else {
      message("Processing metric: ", names.metr[i.metr],"(", i.metr,")")
    }
    
    for (j.k in k.range) {
      message("\tCalculation of k = ", j.k,"")
      estable=NULL
      contador=contador+1
      i=i.metr+1
      estable$n.metric=i.metr
      estable$name.metric=names.metr[i.metr]
      estable$n.k=j.k
      estable$name.ontology=names(datos.bruto[1])
      
      if (all_metrics == TRUE) { # Processing all metrics as one
        data_to_cluster = datos.bruto[,-1] # Removing first column
      } else {
        data_to_cluster = datos.bruto[,i]
      }
      
      km5=NULL
      v.size = 0
      if (!all_metrics) {
        v.size = length(levels(as.factor(data_to_cluster)))
      } else {
        v.size = j.k
      }
      if (v.size>=j.k) {
        #km5$cluster=boot.cluster(data=datos.bruto[,i],
        #                         nk=j.k, B=bs, seed=seed)
        #km5$jac=km5$cluster$means
        
        
        
        clusterbootData = tryCatch({
          clusterbootWrapper(data=data_to_cluster, B=bs,
                             bootmethod="boot",
                             cbi=cbi,
                             gold_standard=gold_standard,
                             krange=j.k, seed=seed, ...)
        }, error = function(error_condition) {
          error_condition
        })
        
        if(inherits(clusterbootData, "error")) {
          message(paste0("\t", clusterbootData))
          message("\tWarning: Could not process data for k = ", j.k)
          km5$bspart=rep(NA,length(datos.bruto[,i]))
          km5$jac=rep(NA,j.k)
          km5$centr=rep(NA,j.k)
          km5$means=km5$bspart
          km5$bspart.or=km5$bspart
          km5$bspart.inv=km5$means
          km5$jac.or=km5$jac
          km5$jac.inv=km5$jac
          km5$partition=km5$bspart.inv
          km5$jac.stab=km5$jac.inv
          
          km5$csv = NULL
          km5$csv$cluster_partition = NULL
          km5$csv$cluster_mean = NULL
          km5$csv$cluster_centers = NULL
          km5$csv$cluster_size = NULL
          km5$csv$cluster_betweenss = NULL
          km5$csv$cluster_totss = NULL
          km5$csv$cluster_tot.withinss = NULL
          km5$csv$cluster_anova = NULL
          
          # CAMBIAMOS A NUESTRAS ESTRUCTURAS LOCALES, NO A LA GLOBAL
          
          stab_valores[j.k] = mean(km5$jac)
          csv_data[[j.k]] = km5
          next
        }
        
        km5$cluster = clusterbootData
        
        km5$jac=km5$cluster$bootmean
        km5$bspart=km5$cluster$partition
        
        km5$csv = NULL
        km5$csv$cluster_partition = km5$cluster$partition
        km5$csv$cluster_mean = km5$cluster$bootmean
        km5$csv$cluster_centers = km5$cluster$result$result$centers
        km5$csv$cluster_size = km5$cluster$result$result$size
        km5$csv$cluster_betweenss = km5$cluster$result$result$betweenss
        km5$csv$cluster_totss = km5$cluster$result$result$totss
        km5$csv$cluster_tot.withinss = km5$cluster$result$result$tot.withinss
        km5$csv$cluster_anova = fAnova(km5$cluster$result$result, j.k, num.metrics)
        
        km5$centr=km5$cluster$result$result$centers
        for (km5.i in 1:length(km5$centr)) {
          km5$means[which(km5$bspart==km5.i)]=km5$centr[km5.i]
        }
        
        #km5$bspart.or=ordered(km5$means,labels=seq(1,length(km5$centr)))
        
        #km5$bspart.inv=ordered(km5$means,labels=seq(length(km5$centr),1))
        
        #km5$jac.or=km5$jac[order(km5$centr)]
        #km5$jac.inv=km5$jac[order(km5$centr,decreasing=TRUE)]
        
        #if (any(inversa==names.metr[i.metr])) {
        #  km5$partition=km5$bspart.inv
        #  km5$jac.stab=km5$jac.inv
        #} else {
        #  km5$partition=km5$bspart.or
        #  km5$jac.stab=km5$jac.or
        #}
        #m.stab.global[[i.metr]][j.k] = mean(km5$jac.stab)
        
        # CAMBIAMOS A NUESTRAS ESTRUCTURAS LOCALES, NO A LA GLOBAL
        stab_valores[j.k] = mean(km5$jac)
        csv_data[[j.k]] = km5
      } else {
        message("\tWarning: Could not process data for k = ", j.k)
        km5$bspart=rep(NA,length(datos.bruto[,i]))
        km5$jac=rep(NA,j.k)
        km5$centr=rep(NA,j.k)
        km5$means=km5$bspart
        km5$bspart.or=km5$bspart
        km5$bspart.inv=km5$means
        km5$jac.or=km5$jac
        km5$jac.inv=km5$jac
        km5$partition=km5$bspart.inv
        km5$jac.stab=km5$jac.inv
        
        km5$csv = NULL
        km5$csv$cluster_partition = NULL
        km5$csv$cluster_mean = NULL
        km5$csv$cluster_centers = NULL
        km5$csv$cluster_size = NULL
        km5$csv$cluster_betweenss = NULL
        km5$csv$cluster_totss = NULL
        km5$csv$cluster_tot.withinss = NULL
        km5$csv$cluster_anova = NULL
        
        # m.stab.global[[i.metr]][j.k] = mean(km5$jac.stab)
        
        # CAMBIAMOS A NUESTRAS ESTRUCTURAS LOCALES, NO A LA GLOBAL
        stab_valores[j.k] = mean(km5$jac)
        csv_data[[j.k]] = km5
      }
    }
    return(list(stab_valores = stab_valores, csv_data = csv_data))
  })
  
  
  # COMBINAMOS
  m.stab.global <- list()
  m.stab.global.csv <- list()
  
  for (i.metr in 1:num.metrics) {
    m.stab.global[[i.metr]] <- resultados[[i.metr]]$stab_valores
    m.stab.global.csv[[i.metr]] <- resultados[[i.metr]]$csv_data
  }
  
  e.stab.global=NULL
  for (j.k in k.range) {
    #e.stab.global[[j.k]]=matrix(data=NA, nrow=length(names.metr), ncol=length(i.min:i.max))
    e.stab.global[[j.k]]=matrix(data=NA, nrow=k.range.length, ncol=1)
    for (i.metr in 1:num.metrics) {
      e.stab.global[[j.k]][i.metr]=m.stab.global[[i.metr]][j.k]
    }
  }
  
  
  pkg.env$m.stab.global.csv = m.stab.global.csv
  pkg.env$m.stab.global = m.stab.global
  pkg.env$e.stab.global = e.stab.global
  #pkg.env$names.metr = names.metr
  #return(stabilityDataFrame)
  #return(NULL)
}