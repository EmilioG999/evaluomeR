# Librería para paralelizar
library(parallel)

# Versión paralela de la función quality
quality_parallel <- function(data, k=5, cbi="kmeans", getImages=FALSE,
                    all_metrics=FALSE, seed=NULL, ...) {
  
  checkKValue(k)
  
  data <- as.data.frame(assay(data))
  
  suppressWarnings(
    runQualityIndicesSilhouette_parallel(data, k.min = k,
                                k.max = k, bs = 1, cbi, all_metrics, seed=seed, ...))
  silhouetteData =  suppressWarnings(
    runSilhouetteTableRange(data, k.min = k, k.max = k))
  
  if (getImages == TRUE) {
    suppressWarnings(
      runQualityIndicesSilhouetteK_IMG(k.min = k, k.max = k))
    suppressWarnings(
      runQualityIndicesSilhouetteMetric_IMG(k.min = k, k.max = k))
  }
  seList <- createSEList(silhouetteData)
  
}

# Versión paralela de la función qualityRange
qualityRange_parallel <- function(data, k.range=c(3,5), cbi="kmeans", getImages=FALSE,
                                  all_metrics=FALSE, seed=NULL, ...) {
  
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
  
  data <- as.data.frame(assay(data))
  
  suppressWarnings(
    runQualityIndicesSilhouette_parallel(data, k.min = k.min,
                                         k.max = k.max, bs = 1, cbi, all_metrics, seed=seed, ...))
  silhouetteData =  suppressWarnings(
    runSilhouetteTableRange(data, k.min = k.min, k.max = k.max))
  
  if (getImages == TRUE) {
    suppressWarnings(
      runQualityIndicesSilhouetteK_IMG(k.min = k.min, k.max = k.max))
    suppressWarnings(
      runQualityIndicesSilhouetteMetric_IMG(k.min = k.min, k.max = k.max))
  }
  seList <- createSEList(silhouetteData)
  return(seList)
}

# Versión paralela de la función qulitySet
qualitySet_parallel <- function(data, k.set=c(2,4), cbi="kmeans", all_metrics=FALSE,
                       getImages=FALSE, seed=NULL, ...) {
  
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
  
  data <- as.data.frame(assay(data))
  
  suppressWarnings(
    runQualityIndicesSilhouette_parallel(data, bs = 1, seed=seed, cbi=cbi, all_metrics=all_metrics,
                                k.set=k.set, ...))
  silhouetteData =  suppressWarnings(
    runSilhouetteTableRange(data, k.set=k.set))
  
  if (getImages == TRUE) {
    suppressWarnings(
      runQualityIndicesSilhouetteK_IMG(k.set=k.set))
    suppressWarnings(
      runQualityIndicesSilhouetteMetric_IMG(k.set=k.set))
  }
  seList <- createSEList(silhouetteData)
  return(seList)
}

# Versión paralela de la función runQualityIndicesSilhouette
# Función paralelizada que divide el número de métricas entre un número de clusters (añadimos numCores)
runQualityIndicesSilhouette_parallel2 <- function(data, k.min=NULL, k.max=NULL, bs,
                                                 cbi, all_metrics, seed=NULL, k.set=NULL, numCores=NULL,...) {
  if (is.null(seed)) {
    seed = pkg.env$seed
  }
  if (is.null(k.min) && is.null(k.max) && is.null(k.set)) {
    stop("runQualityIndicesSilhouette: All k parameters are null!")
  }
  
  data <- removeNAValues(data)
  dfStats(data)
  
  datos.bruto=data
  names.metr=names(datos.bruto)[-c(1)]
  pkg.env$names.metr = names.metr
  names.index=c("sil")
  pkg.env$names.index = names.index
  k.min=k.min
  k.max=k.max
  
  estable=NULL
  m.global=NULL
  e.global=NULL
  contador=0
  remuestreo=bs
  
  i.min=k.min
  i.max=k.max
  
  k.range = NULL
  k.range.length = NULL
  if (!is.null(k.set)) {
    k.range = k.set
    k.range.length = length(k.set)
    nrow = max(k.set)
  } else {
    k.range = i.min:i.max
    k.range.length = length(i.min:i.max)+1
    nrow = i.max
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
  clusterExport(cl, c("datos.bruto", "names.metr", "cbi", "seed", "k.range", "nrow", "num.metrics", "k.range.length", "all_metrics"), envir=environment())
  
  # Cargamos el paquete de evaluomeR
  clusterEvalQ(cl, {
    library(evaluomeR)
    library(cluster)
  })
  
  resultados <- parLapply(cl, 1:num.metrics, function(i.metr){
    
    # Guardamos los resultados de una métrica específica
    quality_valores <- matrix(NA, nrow=nrow, ncol=k.range.length)
    
    
    if (all_metrics == TRUE) {
      message("Processing all metrics, 'merge', in dataframe (", length(names.metr),")")
    } else {
      message("Processing metric: ", names.metr[i.metr],"(", i.metr,")")
    }
    
    
    for (j.k in k.range) {
      message("\tCalculation of k = ", j.k,"")
      e.res=NULL
      e.res.or=NULL
      #contador=contador+1
      i=i.metr+1
      j=j.k
      
      if (all_metrics == TRUE) { # Processing all metrics as one
        data_to_cluster = datos.bruto[,-1] # Removing first column
      } else {
        data_to_cluster = datos.bruto[,i]
      }
      
      #e.res$n=contador
      e.res$n.metric=i.metr
      e.res$name.metric=names.metr[i.metr]
      
      e.res$n.k=j.k
      e.res$name.ontology=datos.bruto$Description
      unique.values = length(unique(data_to_cluster))
      
      if (unique.values < j.k) {
        # SUSTITUIMOS POR NUESTRA ESTRUCTURA
        quality_valores[j.k,] <- NA
        message("\tWarning: Could not process data for k = ", j.k)
      } else {
        # bootClusterResult <- boot.cluster(data=datos.bruto[,i],
        #                                  nk=j.k, B=bs, seed=seed)
        bootClusterResult <- clusteringWrapper(data=data_to_cluster, cbi=cbi,
                                               krange=j.k, seed=seed, ...)
        # bootClusterResult <- clusterbootWrapper(data=datos.bruto[,i], B=bs,
        #                    bootmethod="boot",
        #                    cbi=cbi,
        #                    krange=j.k, seed=seed)
        
        e.res$kmk.dynamic.bs <- as.integer(bootClusterResult$partition)
        
        e.res.or$centr=bootClusterResult$result$centers
        #e.res.or$centr=by(datos.bruto[,i],e.res$kmk.dynamic.bs,mean)
        #for (e.res.or.i in 1:length(e.res.or$centr)) {
        #  e.res.or$means[which(e.res$kmk.dynamic.bs==e.res.or.i)]=e.res.or$centr[e.res.or.i]}
        #e.res$kmk.dynamic.bs.or=ordered(e.res.or$means,labels=seq(1,length(e.res.or$centr)))
        
        e.res$kmk.dynamic.bs.or = bootClusterResult$partition
        ## Using Silhouette width as index
        metric.onto=data_to_cluster
        # part.onto=as.numeric(e.res$kmk.dynamic.bs.or)
        part.onto = bootClusterResult$partition
        sil.w=silhouette(part.onto, dist(metric.onto))
        sil.c = NULL
        sil.c$n=length(sil.w[,1])
        sil.c$cluster.size = as.numeric(summary(sil.w)$clus.sizes)
        sil.c$cluster.pos = part.onto
        sil.c$cluster.labels = e.res$name.ontology
        sil.c$cluster.number = length(summary(sil.w)$cluster.size)
        sil.c$clus.avg.silwidths = summary(sil.w)$clus.avg.widths
        sil.c$avg.silwidths = summary(sil.w)$avg.width
        e.res$sil.w = sil.w
        e.res$sil.c = sil.c
        #estable[[contador]] = e.res
        
        # SUSTITUIMOS POR NUESTRA ESTRUCTURA
        quality_valores[j.k,] = mean(sil.w[,"sil_width"])
      }
    }
    return(list(quality_valores = quality_valores))
  })
  
  # COMBINAMOS
  m.global <- list()
  for (i.metr in 1:num.metrics) {
    m.global[[i.metr]] <- resultados[[i.metr]]$quality_valores
  }
  
  e.global <- list()
  for (j.k in k.range) {
    e.global[[j.k]]=matrix(data=NA, nrow=num.metrics, ncol=k.range.length)
    for (i.metr in 1:num.metrics) {
      e.global[[j.k]][i.metr,]=m.global[[i.metr]][j.k,]
    }
  }
  
  pkg.env$m.global = m.global
  pkg.env$e.global = e.global
  #pkg.env$estable = estable
}