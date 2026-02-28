#' @title Goodness of classifications.
#' @name quality
#' @aliases quality
#' @description
#' The goodness of the classifications are assessed by validating the clusters
#' generated. For this purpose, we use the Silhouette width as validity index.
#' This index computes and compares the quality of the clustering outputs found
#' by the different metrics, thus enabling to measure the goodness of the
#' classification for both instances and metrics. More precisely, this goodness measurement
#' provides an assessment of how similar an instance is to other instances from
#' the same cluster and dissimilar to all the other clusters. The average on all
#' the instances quantifies how appropriately the instances are clustered. Kaufman
#' and Rousseeuw suggested the interpretation of the global Silhouette width score
#' as the effectiveness of the clustering structure. The values are in the
#' range [0,1], having the following meaning:
#'
#' \itemize{
#' \item There is no substantial clustering structure: [-1, 0.25].
#' \item The clustering structure is weak and could be artificial: ]0.25, 0.50].
#' \item There is a reasonable clustering structure: ]0.50, 0.70].
#' \item A strong clustering structure has been found: ]0.70, 1].
#' }
#'
#' @inheritParams stability
#'
#' @return A \code{\link{SummarizedExperiment}} containing the silhouette width measurements and
#' cluster sizes for cluster \code{k}.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' result = quality(ontMetrics, k=4)
#'
#' @references
#' \insertRef{kaufman2009finding}{evaluomeR}
#'
quality <- function(data, k=5, cbi="kmeans", getImages=FALSE,
                    all_metrics=FALSE, seed=NULL, numCores=1,...) {

  checkKValue(k)

  data <- as.data.frame(assay(data))
  
  if(numCores>1){
    # Versión paralela
    suppressWarnings(
      runQualityIndicesSilhouette_parallel(data, k.min = k,
                                  k.max = k, bs = 1, cbi, all_metrics, seed=seed, numCores=numCores, ...))
    silhouetteData =  suppressWarnings(
      runSilhouetteTableRange(data, k.min = k, k.max = k))
  }else{
    suppressWarnings(
      runQualityIndicesSilhouette(data, k.min = k,
                                  k.max = k, bs = 1, cbi, all_metrics, seed=seed, ...))
    silhouetteData =  suppressWarnings(
      runSilhouetteTableRange(data, k.min = k, k.max = k))
  }
  
  if (getImages == TRUE) {
    suppressWarnings(
      runQualityIndicesSilhouetteK_IMG(k.min = k, k.max = k))
    suppressWarnings(
      runQualityIndicesSilhouetteMetric_IMG(k.min = k, k.max = k))
  }
  seList <- createSEList(silhouetteData)

}

#' @title Goodness of classifications for a range of k clusters.
#' @name qualityRange
#' @aliases qualityRange
#' @description
#' The goodness of the classifications are assessed by validating the clusters
#' generated for a range of k values. For this purpose, we use the Silhouette width as validity index.
#' This index computes and compares the quality of the clustering outputs found
#' by the different metrics, thus enabling to measure the goodness of the
#' classification for both instances and metrics. More precisely, this measurement
#' provides an assessment of how similar an instance is to other instances from
#' the same cluster and dissimilar to the rest of clusters. The average on all
#' the instances quantifies how the instances appropriately are clustered. Kaufman
#' and Rousseeuw suggested the interpretation of the global Silhouette width score
#' as the effectiveness of the clustering structure. The values are in the
#' range [0,1], having the following meaning:
#'
#' \itemize{
#' \item There is no substantial clustering structure: [-1, 0.25].
#' \item The clustering structure is weak and could be artificial: ]0.25, 0.50].
#' \item There is a reasonable clustering structure: ]0.50, 0.70].
#' \item A strong clustering structure has been found: ]0.70, 1].
#' }
#'
#' @inheritParams stability
#' @param k.range Concatenation of two positive integers.
#' The first value \code{k.range[1]} is considered as the lower bound of the range,
#' whilst the second one, \code{k.range[2]}, as the higher. Both values must be
#' contained in [2,15] range.
#'
#' @return A list of \code{\link{SummarizedExperiment}} containing the silhouette width measurements and
#' cluster sizes from \code{k.range[1]} to \code{k.range[2]}. The position on the list matches
#' with the k-value used in that dataframe. For instance, position 5
#' represents the dataframe with k = 5.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' # Without plotting
#' dataFrameList = qualityRange(ontMetrics, k.range=c(2,3), getImages = FALSE)
#'
#' @references
#' \insertRef{kaufman2009finding}{evaluomeR}
#'
qualityRange <- function(data, k.range=c(3,5), cbi="kmeans", getImages=FALSE,
                         all_metrics=FALSE, seed=NULL, numCores=1, ...) {

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

  if(numCores>1){
    # Versión paralela
    suppressWarnings(
      runQualityIndicesSilhouette_parallel(data, k.min = k.min,
                                  k.max = k.max, bs = 1, cbi, all_metrics, seed=seed, numCores=numCores, ...))
    silhouetteData =  suppressWarnings(
      runSilhouetteTableRange(data, k.min = k.min, k.max = k.max))
  }else{
    suppressWarnings(
      runQualityIndicesSilhouette(data, k.min = k.min,
                                  k.max = k.max, bs = 1, cbi, all_metrics, seed=seed, ...))
    silhouetteData =  suppressWarnings(
      runSilhouetteTableRange(data, k.min = k.min, k.max = k.max))
  }
  

  if (getImages == TRUE) {
    suppressWarnings(
      runQualityIndicesSilhouetteK_IMG(k.min = k.min, k.max = k.max))
    suppressWarnings(
      runQualityIndicesSilhouetteMetric_IMG(k.min = k.min, k.max = k.max))
  }
  seList <- createSEList(silhouetteData)
  return(seList)
}

#' @title Goodness of classifications for a set of k clusters.
#' @name qualitySet
#' @aliases qualitySet
#' @description
#' The goodness of the classifications are assessed by validating the clusters
#' generated for a range of k values. For this purpose, we use the Silhouette width as validity index.
#' This index computes and compares the quality of the clustering outputs found
#' by the different metrics, thus enabling to measure the goodness of the
#' classification for both instances and metrics. More precisely, this measurement
#' provides an assessment of how similar an instance is to other instances from
#' the same cluster and dissimilar to the rest of clusters. The average on all
#' the instances quantifies how the instances appropriately are clustered. Kaufman
#' and Rousseeuw suggested the interpretation of the global Silhouette width score
#' as the effectiveness of the clustering structure. The values are in the
#' range [0,1], having the following meaning:
#'
#' \itemize{
#' \item There is no substantial clustering structure: [-1, 0.25].
#' \item The clustering structure is weak and could be artificial: ]0.25, 0.50].
#' \item There is a reasonable clustering structure: ]0.50, 0.70].
#' \item A strong clustering structure has been found: ]0.70, 1].
#' }
#'
#' @inheritParams stability
#' @param k.set A list of integer values of \code{k}, as in c(2,4,8).
#' The values must be contained in [2,15] range.
#'
#' @return A list of \code{\link{SummarizedExperiment}} containing the silhouette width measurements and
#' cluster sizes from \code{k.set}.
#'
#' @examples
#' # Using example data from our package
#' data("rnaMetrics")
#' # Without plotting
#' dataFrameList = qualitySet(rnaMetrics, k.set=c(2,3), getImages = FALSE)
#'
#' @references
#' \insertRef{kaufman2009finding}{evaluomeR}
#'
qualitySet <- function(data, k.set=c(2,4), cbi="kmeans", all_metrics=FALSE,
                       getImages=FALSE, seed=NULL, numCores=1, ...) {

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

  if(numCores>1){
    suppressWarnings(
      runQualityIndicesSilhouette_parallel(data, bs = 1, seed=seed, cbi=cbi, all_metrics=all_metrics,
                                  k.set=k.set, numCores=numCores, ...))
    silhouetteData =  suppressWarnings(
      runSilhouetteTableRange(data, k.set=k.set))
  }else{
    suppressWarnings(
      runQualityIndicesSilhouette(data, bs = 1, seed=seed, cbi=cbi, all_metrics=all_metrics,
                                  k.set=k.set, ...))
    silhouetteData =  suppressWarnings(
      runSilhouetteTableRange(data, k.set=k.set))
  }
  
  if (getImages == TRUE) {
    suppressWarnings(
      runQualityIndicesSilhouetteK_IMG(k.set=k.set))
    suppressWarnings(
      runQualityIndicesSilhouetteMetric_IMG(k.set=k.set))
  }
  seList <- createSEList(silhouetteData)
  return(seList)
}

runQualityIndicesSilhouette <- function(data, k.min=NULL, k.max=NULL, bs,
                                        cbi, all_metrics, seed=NULL, k.set=NULL, ...) {
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
  } else {
    num.metrics = length(names.metr)
  }

  for (i.metr in 1:num.metrics) {
    if (all_metrics == TRUE) {
      message("Processing all metrics, 'merge', in dataframe (", length(names.metr),")")
      pkg.env$names.metr = c("all_metrics")
    } else {
      message("Processing metric: ", names.metr[i.metr],"(", i.metr,")")
    }

    m.global[[i.metr]]=matrix(data=NA, nrow=nrow, ncol=k.range.length)

    for (j.k in k.range) {
      message("\tCalculation of k = ", j.k,"")
      e.res=NULL
      e.res.or=NULL
      contador=contador+1
      i=i.metr+1
      j=j.k

      if (all_metrics == TRUE) { # Processing all metrics as one
        data_to_cluster = datos.bruto[,-1] # Removing first column
      } else {
        data_to_cluster = datos.bruto[,i]
      }

      e.res$n=contador
      e.res$n.metric=i.metr
      e.res$name.metric=pkg.env$names.metr[i.metr]

      e.res$n.k=j.k
      e.res$name.ontology=datos.bruto$Description
      unique.values = length(unique(data_to_cluster))

      if (unique.values < j.k) {
        estable[[contador]] = NA
        m.global[[i.metr]][j.k,] = NA
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
        estable[[contador]] = e.res

        m.global[[i.metr]][j.k,] = mean(sil.w[,"sil_width"])
      }
    }
  }
  for (j.k in k.range) {
    e.global[[j.k]]=matrix(data=NA, nrow=length(names.metr), ncol=k.range.length)
    for (i.metr in 1:length(pkg.env$names.metr)) {
      e.global[[j.k]][i.metr,]=m.global[[i.metr]][j.k,]
    }
  }

  pkg.env$m.global = m.global
  pkg.env$e.global = e.global
  pkg.env$estable = estable
}

# Versión paralela de la función runQualityIndicesSilhouette
# Función paralelizada que divide el número de métricas entre un número de clusters (añadimos numCores)
runQualityIndicesSilhouette_parallel <- function(data, k.min=NULL, k.max=NULL, bs,
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
  
  # Se puede también detectar el número de cores o se puede asignar manualmente aquí si se prefiere
  # numCores <- detectCores() 
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

# Silhouette width per k (x values = metrics)
runQualityIndicesSilhouetteK_IMG <- function(k.min=NULL, k.max=NULL, k.set=NULL) {
  if (is.null(k.min) && is.null(k.max) && is.null(k.set)) {
    stop("runQualityIndicesSilhouetteK_IMG: All k parameters are null!")
  }
  ancho=6
  alto=4
  escala=0.75
  escalax=escala
  escalal=0.85
  ajuste=0.5
  escalat=0.5
  escalap=0.4

  m.global = pkg.env$m.global
  e.global = pkg.env$e.global
  e.mat.global=e.global

  names.index = pkg.env$names.index
  i.min=1
  i.max=k.max-(k.min-1)
  leg.g = NULL
  names.metr = pkg.env$names.metr
  x=seq(1,length(names.metr))
  x.label="Metrics"
  x.name=xnames=as.character(names.metr)
  y.label="Silhouette avg. width"
  #Pattern: QualityIndices_K_2, ..., QualityIndices_K_N
  figurename="QualityIndices_K_"

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
  margins <- par(mar=c(5,5,3,3))
  on.exit(par(margins))
  stype <- c(1:k.range.length)
  metrics_length = length(names.metr)
  num_metrics_plot = 19
  num_iterations = round(metrics_length/num_metrics_plot)
  if (num_iterations > 0) {
    num_iterations = num_iterations - 1
  }
  for (iteration in 0:num_iterations) {
    i = 1
    labels = list()
    rangeStart = (iteration*num_metrics_plot)+1
    rangeEnd = rangeStart+num_metrics_plot
    if (rangeEnd > metrics_length) {
      rangeEnd = metrics_length
    }
    new_xnames = x.name[rangeStart:rangeEnd]
    if (!is.null(k.set)) {
      setAsStrList = paste(as.character(k.set),collapse=", ",sep="")
      g.main=paste(" Qual. Indices of the metrics for k in {", setAsStrList, "}",sep="")
    } else {
      g.main=paste(" Qual. Indices of the metrics for k in [", i.min, ",", i.max, "]",sep="")
    }

    for (m.g in k.range) {
      c.max=dim(e.mat.global[[m.g]])[2]
      ymarcas=round(seq(0,1,length.out=5),2)

      k.classes=m.g

      for (m in length(names.index)) {
        y=e.mat.global[[m.g]][,m]
        if (all(is.na(y))) { # Skip if all values are NA?
          next
        }

        y = y[rangeStart:rangeEnd]
        y.name=names.index[m]
        #leg.g[m] <- paste(y.name," avg. width",sep="")
        plot(y, main=g.main, axes=TRUE, col.axis="white",
             xlim=c(0.75,length(new_xnames)+0.25), xlab="", ylim=c(0,1),
             ylab="", col="black", type="o", lwd=1, lty=stype[i], pch=stype[i])
        #par(new=TRUE)
        labels = c(labels,(paste0("k=", m.g)))
        i = i + 1
        par(new=TRUE)
      }
      par(new=TRUE)
    }
    mtext(side=1, text=x.label,line=4)
    mtext(side=2, text=y.label,line=3)
    axis(1,at=1:length(new_xnames),labels=new_xnames,las=2,cex.axis=0.75)
    axis(2,las=3,cex.axis=0.85)
    legend("bottomright", legend=labels, inset=.01, lwd=1, lty=stype, col="black", cex=0.7, pch=stype)
    text(0.76, 0.75, "Strong", cex=0.6, col = "black")
    abline(h = 0.7, col="black", lwd=1, lty=1) # Strong clust. strct.: (0.70, 1
    text(0.76, 0.55, "Reasn.", cex=0.6, col = "black")
    abline(h = 0.5, col="black", lwd=1, lty=1) # Reasonable clust. strct.: (0.50, 0.70]
    text(0.76, 0.3, "Weak", cex=0.6, col = "black")
    abline(h = 0.25, col="black", lwd=1, lty=1) # Weak clust. strct.: (0.25, 0.50]
    text(0.77, 0.05, "No.strct", cex=0.6, col = "black")
    #abline(h = -1, col="black", lwd=1, lty=4) # No clust. strct.: [-1, 0.25]
    par(new=FALSE)

  }
}

# Silhouette width per metric (x values = k range)
runQualityIndicesSilhouetteMetric_IMG <- function(k.min=NULL, k.max=NULL, k.set=NULL) {
  if (is.null(k.min) && is.null(k.max) && is.null(k.set)) {
    stop("runQualityIndicesSilhouetteMetric_IMG: All k parameters are null!")
  }

  ancho=6
  alto=4
  escala=0.9
  escalax=escala
  escalal=0.85
  ajuste=0.5
  escalat=0.5
  escalap=0.4

  m.global = pkg.env$m.global
  m.mat.global=m.global
  names.index = pkg.env$names.index
  names.metr = pkg.env$names.metr

  y.label="Silhouette avg. width"
  #Pattern: QualityIndices__MetricX, ..., QualityIndices__MetricN
  figurename="QualityIndices_"

  i.min=k.min
  i.max=k.max

  k.range = NULL
  k.range.length = NULL

  if (!is.null(k.set)) {
    k.range = k.set
    k.range.length = length(k.set)
    x=k.set
  } else {
    k.range = i.min:i.max
    k.range.length = length(i.min:i.max)+1
  }

  x=c(k.range)
  x.name=as.character(k.range)
  x.label="K values"


  margins <- par(mar=c(5,5,3,3))
  on.exit(par(margins))
  for (m.g in 1:length(names.metr)) {
    cur.k.width = m.mat.global[[m.g]][,1]
    cur.k.width = cur.k.width[k.range]
    #cur.k.width = cur.k.width[!is.na(cur.k.width)]
    leg.g=NULL
    xmin=min(x)-0.25
    xmax=max(x)+0.25
    xleg=((xmax-xmin)*escalal)+3.2
    c.max=dim(m.mat.global[[m.g]])[2]
    ymin=min(cur.k.width)
    if (is.na(ymin)) {
      ymin = 0
    }
    ymax=1
    ymarcas=round(seq(ymin,ymax,length.out=5),2)
    yleg=ymin+((ymax-ymin)/2)*seq(c.max,1,-1)/(2*c.max)
    t.linea=seq(1,c.max)
    t.color=rep("black",c.max)

    if (!is.null(k.set)) {
      setAsStrList = paste(as.character(k.set),collapse=", ",sep="")
      g.main=paste(" Qual. Indices of '", names.metr[m.g], "' for k in {", setAsStrList, "}",sep="")
    } else {
      g.main=paste(" Qual. Indices of '", names.metr[m.g], "' for k in [",
                   i.min, ",", i.max,"]",sep="")
    }

    y=cur.k.width
    y.name=names.index[1]
    leg.g[1] <- paste(y.name," avg. width",sep="")
    plot(x,y, type="l", xaxt="n", yaxt="n", xlab="", ylab="", main=g.main, xlim=c(xmin,xmax), ylim=c(ymin,ymax), lty=t.linea[1], col=t.color[1])
    par(new=TRUE)
    plot(x,y, type="o", xaxt="n", yaxt="n", xlab="", ylab="", main=g.main, xlim=c(xmin,xmax), ylim=c(ymin,ymax), lty=t.linea[1], col=t.color[1])
    par(new=TRUE)

    mtext(side=1, text=x.label,line=3)
    mtext(side=2, text=y.label,line=3)
    axis(side=1, at=x, labels=x.name, las=1, cex.axis=escalax)
    axis(side=2, at=ymarcas, labels=ymarcas, cex.axis=escalal)
    par(new=FALSE)
  }
}

runSilhouetteIMG <- function(data, k) {
  names.metr = pkg.env$names.metr
  datos.bruto = data
  estable = pkg.env$estable

  ancho=7
  alto=6
  escala=1     #new 0.6
  escalax=0.7  #new escala
  escalal=0.75  #new 0.8
  ajuste=0.5
  escalat=0.5
  escalap=0.4
  par(new=FALSE,bg="white",fg="black", cex=1, mex=.6)
  onto.matrix=matrix(data=NA, nrow=length(datos.bruto[,1]), ncol=(length(names.metr)+1))
  onto.matrix[,1]=as.character(datos.bruto[,1])
  colnames(onto.matrix)=c("Datasets",paste(names.metr,sep="."))

  margenes=c(6,4,6,8)
  margins <- par(mar=margenes, cex=escala, mex=escalal)
  on.exit(par(margins))

  k.cl = k
  colores=c(2:(k.cl+1)) # 2 to k.cl+1, avoid number 1 since it's black and it's not pretty
  #Pattern: Silhouette_K_N_MetricX, ..., Silhouette_K_N_MetricN
  figurename="Silhouette_K_"

  for (i.metr in 1:length(names.metr)) {

    x.leyenda=0.99

    metric.onto=datos.bruto[,i.metr+1]
    metric.name=names(datos.bruto)[i.metr+1]

    i.datos=i.metr
    if (!is.list(estable[[i.datos]]) && is.na(estable[[i.datos]])) {
      next
    }
    for (estable.content in estable[[i.datos]]) {
      if (is.list(estable.content)) {
        # Could not calculate silhouette clustering for this metric
        # (Data used for horizontal bars graph)
        next
      }
    }
    if (estable[[i.datos]]$n.k==k.cl & estable[[i.datos]]$name.metric==metric.name) {
      part.onto=as.numeric(estable[[i.datos]]$kmk.dynamic.bs.or)
      onto.matrix[,(i.metr+1)]=part.onto

      sil.w = estable[[i.datos]]$sil.w
      sil.c = estable[[i.datos]]$sil.c

      estable[[i.datos]]$kmk.dynamic.bs.or.numeric=part.onto
      estable[[i.datos]]$sil.width=sil.w

      g.main=paste(metric.name,sep="")

      plot(sil.w, col=colores, main=g.main, border=NULL,
           mar=margenes, cex=escala, mex=escalal,
           cex.names = par("cex.axis"), do.n.k = TRUE, do.clus.stat = FALSE)
      t.leyenda=c(expression('j:  n'['j']), expression(' | ave'['i' %in% 'C'['j']]), expression('s'['i']))
      legend(x=x.leyenda,y=sil.c$n+1, legend=expression('j:  n'['j']), col="black",
             xjust=0, yjust=0, bty="n", xpd=TRUE, inset=c(-0.1,0), cex=escalax)
      legend(x=x.leyenda+0.03,y=sil.c$n+1, legend=expression(' | ave'['i' %in% 'C'['j']]), col="black",
             xjust=0, yjust=0, bty="n", xpd=TRUE, inset=c(-0.1,0), cex=escalax)
      legend(x=x.leyenda+0.1,y=sil.c$n+1, legend=expression(' s'['i']), col="black",
             xjust=0, yjust=0, bty="n", xpd=TRUE, inset=c(-0.1,0), cex=escalax)
      xleyenda=rep(x.leyenda,k.cl)
      #yleyenda=(sil.c$cluster.size==1)*0.6*(sil.c$n-cumsum(sil.c$cluster.size))+
      #  (sil.c$n-cumsum(sil.c$cluster.size))+sil.c$cluster.size*3/k.cl+2
      yleyenda=(sil.c$cluster.size==1)*0.6*(sil.c$n-cumsum(sil.c$cluster.size))+
        (sil.c$n-cumsum(sil.c$cluster.size))+sil.c$cluster.size/k.cl+2
      leyenda=paste(names(sil.c$clus.avg.silwidths),rep(": ",k.cl),
                    sil.c$cluster.size,"|",round(sil.c$clus.avg.silwidths,digits=2),sec="")
      for (i.leyenda in 1:k.cl){
        legend(list(x=xleyenda[i.leyenda],y=yleyenda[i.leyenda]), legend=leyenda[i.leyenda], col="black",
               xjust=0, yjust=1, bty="n", xpd=TRUE, inset=c(-0.1,0), cex=escalal)
        }
    }
  }
}

runSilhouetteTable <- function(data, k) {
  data = removeNAValues(data, verbose=FALSE)
  names.metr = pkg.env$names.metr
  datos.bruto = data
  estable = pkg.env$estable
  k.cl = k
  ##
  #  Building table header
  ##
  silhouetteData <- list()
  silhouetteData$header <- list("Metric")
  for (i in 1:k.cl) {
    header = paste("Cluster_", i, "_SilScore", sep="")
    silhouetteData$header <- c(silhouetteData$header, header)
  }
  silhouetteData$header <- c(silhouetteData$header, "Avg_Silhouette_Width")
  for (i in 1:k.cl) {
    header = paste("Cluster_", i, "_Size", sep="")
    silhouetteData$header <- c(silhouetteData$header, header)
  }

  silhouetteData$header = unlist(silhouetteData$header, use.names=FALSE)
  ##
  #  Building table header
  ##

  #onto.matrix=matrix(data=NA, nrow=length(datos.bruto[,1]), ncol=(length(names.metr)+1))
  #onto.matrix[,1]=as.character(datos.bruto[,1])
  #colnames(onto.matrix)=c("Datasets",paste(names.metr,sep="."))
  for (i.metr in 1:length(names.metr)) { # i.metr= n de metrica     i.metr=5

    #metric.onto=datos.bruto[,i.metr+1]
    metric.name=names.metr[i.metr]
    #x.leyenda=0.99
    #
    i.datos=i.metr
    if (!is.list(estable[[i.datos]]) && is.na(estable[[i.datos]])) {
      next
    }
    for (estable.content in estable[[i.datos]]) {
      if (is.null(estable.content) || !is.list(estable.content)) {
        # Could not calculate silhouette clustering for this metric
        # (Data used for horizontal bars graph)
        next
      }
    }
    if (estable[[i.datos]]$n.k==k.cl &
        estable[[i.datos]]$name.metric==metric.name) {
        part.onto=as.numeric(estable[[i.datos]]$kmk.dynamic.bs.or)
        #onto.matrix[,(i.metr+1)]=part.onto

        sil.w = estable[[i.datos]]$sil.w
        sil.c = estable[[i.datos]]$sil.c

        ## Building body rows
        silhouetteData$body[[i.metr]]=c(metric.name, sil.c$clus.avg.silwidths, mean(sil.w[,"sil_width"]), sil.c$cluster.size)
        silhouetteData$body[[i.metr]]=unlist(silhouetteData$body[[i.metr]], use.names=FALSE)
    }
  }  # end for i.metr
  # Remove empty positions in the list
  # (just in case a metric could not be processed)
  silhouetteData$body = silhouetteData$body[lapply(silhouetteData$body, length)>0]
  silhouetteDataFrame = t(data.frame(silhouetteData$body))
  colnames(silhouetteDataFrame) = silhouetteData$header
  rownames(silhouetteDataFrame) <- NULL
  return(silhouetteDataFrame)
}

runSilhouetteTableRange <- function(data, k.min=NULL, k.max=NULL, k.set=NULL) {
  if (is.null(k.min) && is.null(k.max) && is.null(k.set)) {
    stop("runSilhouetteTableRange: All k parameters are null!")
  }

  getHeader <- function(k) {
    silhouetteData$header <- list("Metric")
    for (i in 1:k) {
      header = paste("Cluster_", i, "_SilScore", sep="")
      silhouetteData$header <- c(silhouetteData$header, header)
    }
    silhouetteData$header <- c(silhouetteData$header, "Avg_Silhouette_Width")
    for (i in 1:k) {
      header = paste("Cluster_", i, "_Size", sep="")
      silhouetteData$header <- c(silhouetteData$header, header)
    }

    header = paste(c("Cluster_position","Cluster_labels"))
    silhouetteData$header <- c(silhouetteData$header, header)

    silhouetteData$header = unlist(silhouetteData$header, use.names=FALSE)
    return(silhouetteData$header)
  }

  names.metr = pkg.env$names.metr
  datos.bruto = data

  estable = pkg.env$estable
  k.min = k.min
  k.max = k.max

  k.range = NULL
  k.range.length = NULL
  if (!is.null(k.set)) {
    k.range = k.set
    k.range.length = length(k.set)
  } else {
    k.range = k.min:k.max
    k.range.length = length(k.min:k.max)+1
  }

  onto.matrix=matrix(data=NA, nrow=length(datos.bruto[,1]), ncol=(length(names.metr)+1))
  onto.matrix[,1]=as.character(datos.bruto[,1])
  colnames(onto.matrix)=c("Datasets",paste(names.metr,sep="."))
  offset = 0
  estableLength = length(estable)
  names.metrLength = length(names.metr)
  silhouetteData <- list()
  silhouetteDataIndex <- vector(mode="integer", length=k.range.length)
  for (k in k.range) {
    header <- getHeader(k = k)
    silhouetteData[[k]] <- data.frame(matrix(ncol = length(header), nrow = names.metrLength))
    colnames(silhouetteData[[k]]) = header
    rownames(silhouetteData[[k]]) <- NULL
    silhouetteDataIndex[k] = 1
  }
  # estable object stores names.metr * length(k.min:k.max) entries
  for (i.metr in 1:estableLength) {

    cur.metr = as.integer(abs(i.metr-(names.metrLength*offset)))
    cur.data = estable[[i.metr]]
    if (is.list(cur.data) && !is.null(cur.data)) {
      cur.k = cur.data$n.k
      cur.row <- list(cur.data$name.metric)
      cur.row <- c(cur.row, cur.data$sil.c$clus.avg.silwidths)
      #cur.row <- c(cur.row, cur.data$sil.c$avg.silwidth)
      cur.row <- c(cur.row, mean(cur.data$sil.w[,"sil_width"]))
      cur.row <- c(cur.row, cur.data$sil.c$cluster.size)

      cur.row <- c(cur.row, paste(cur.data$sil.c$cluster.pos, collapse  = ","))
      cur.row <- c(cur.row, paste(cur.data$sil.c$cluster.labels, collapse  = ","))

      cur.row <- unlist(cur.row, use.names = FALSE)

      index = silhouetteDataIndex[cur.k]
      silhouetteData[[cur.k]] = insertRow(silhouetteData[[cur.k]], cur.row,index)
      silhouetteDataIndex[cur.k] = index + 1
    }

    if (cur.metr == names.metrLength) { # Last metric
      offset = offset + 1
    }
  }

  ##
  # Data cleaning
  ##
  # Matrix inserts NA by default, remove them before returning the data
  emptyDataFrames = list()
  emptyDataFramesIndex = 1
  for (k in k.range) {
    silhouetteData[[k]] <- na.omit(silhouetteData[[k]])
    if (nrow(silhouetteData[[k]]) == 0) {
      emptyDataFrames[[emptyDataFramesIndex]] = k
      emptyDataFramesIndex = emptyDataFramesIndex + 1
    }
  }

  # Delete empty dfs (this occurs when no bootstrap is performed for a k)
  silhouetteData = Filter(NROW, silhouetteData)

  # Delete k if its df was empty
  k.range = k.range[!k.range %in% emptyDataFrames]
  silhouetteData[sapply(silhouetteData, is.null)] <- NULL
  names(silhouetteData) <- paste("k_", k.range, sep = "")
  return(silhouetteData)
}

checkKValue <- function(k) {
  if (k < 2 || k > 15) {
    error=paste("k value (",k,") is not in range [2,15]", sep="")
    stop(error)
  }
}