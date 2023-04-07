
response_marker_enrichment_analysis_v2 <- function(protein_data,
                                                marker_type=c("crispr","drug","rnai"),
                                                tumour_type=c("haem","solid")){

  # in protein_data:
  #     - second column = folds (log2)
  #     - third column = pvalues
  #     - protein_acc column must be present, these are the XXX_HUMAN accessions to match to marker signatures

  library(ggplot2)
  library(ggrepel)
  library(foreach)
  library(doParallel)

  # Get signatues to analyse
  signature_list <- ReMEA::list_of_signatures
  signatures_to_analyse <- foreach(m=marker_type, .combine="c")%do%{
    foreach (tt = tumour_type, .combine = "c")%do%{
      return(names(signature_list)[grepl(m,names(signature_list),fixed = T) &
                            grepl(tt,names(signature_list),fixed = T)])
    }
  }
  cl <- makeCluster(detectCores(logical = TRUE)-1)
  registerDoParallel(cl)
  df <- foreach(s = signatures_to_analyse,.combine = 'rbind')%dopar%{
    sl <- signature_list[[s]]
    if (nrow(sl)>1){
      return(ReMEA::signature_enrichment_v2(protein_data = protein_data ,
                                df.signatures = sl,
                                signature.type = s))
    }
  }
  stopCluster(cl)
  df$delta.zscore <- df$zscore.sensitivity-df$zscore.resistance
  df$total_counts <- df$n.resistance.markers+df$n.sensitivity.markers

  df$delta_rank_mean <- log2(df$av.rank.sensitivity)-log2(df$av.rank.resistance)
  df$delta_rank_median <- log2(df$med.rank.sensitivity)-log2(df$med.rank.resistance)
  df$delta_rank_geommean <- log2(df$geommean.rank.sensitivity)-log2(df$geommean.rank.resistance)
  df$qvalue <- p.adjust(df$pvalue.ks,method = 'fdr')
  df$max_delta_score <- apply(df[,c('delta.zscore','delta_rank_geommean','delta_rank_mean','delta_rank_median')],1,max,na.rm=T)

  if ("drug" %in% marker_type){
    df.druginfo <- ReMEA::drug_info_from_Shelleckchem_pharmacoDB_plus_Finland
    df$perturbagens <- gsub("-",'.',df$perturbagens,fixed = T)
    df$perturbagens <- gsub("(",'.',df$perturbagens,fixed = T)
    df$perturbagens <- gsub(")",'.',df$perturbagens,fixed = T)
    df <- merge.data.frame(df,drug_info_from_Shelleckchem_pharmacoDB_plus_Finland,
                           by.x='perturbagens',by.y = 'drug')
    df <- data.frame(drug_target=df$drug.target,df)
    df$perturbagens <- df$drug_target
  }
  return(list(enrichment_data=df))
}


signature_enrichment_v2 <- function(protein_data,
                                 df.signatures,
                                 signature.type){
  library(BWStest)
  protein_data <- na.omit(protein_data)
  protein_data[,2] <- scale(protein_data[,2])
  protein_data$protein_rank <- rank(protein_data[,2])
  mean.population <- mean(protein_data[,2])
  sd.population <- sd(protein_data[,2])

  nr <- nrow(df.signatures)

  av.rank.resistance <- numeric(nr)
  av.rank.sensitivity <- numeric(nr)

  geommean.rank.resistance <- numeric(nr)
  geommean.rank.sensitivity <- numeric(nr)

  med.rank.resistance <- numeric(nr)
  med.rank.sensitivity <- numeric(nr)


  zscore.resistance <- numeric(nr)
  zscore.sensitivity <- numeric(nr)
  pvalue.ks <- numeric(nr)
  perturbagen.signature <- character(nr)
  perturbagens <- character(nr)
  pvalue.bws <- numeric(nr)

  #q1 <- quantile(protein_data[,2])[2]
  #q3 <- quantile(protein_data[,2])[4]

  #df.prot.increased <- subset(protein_data,protein_data[,2]>q3 & protein_data[,3]<pvalue.cutoff)
  #df.prot.decreased <- subset(protein_data,protein_data[,2]<q1 & protein_data[,3]<pvalue.cutoff)



  n.resistance.markers <- numeric(nr)
  n.sensitivity.markers <- numeric(nr)
  resistance.markers <- character(nr)
  sensitivity.markers <- character(nr)

  r <- 1
  for (r in 1:nr) {
    n.res <- df.signatures[r,]$n.resistance.markers
    n.sen <- df.signatures[r,]$n.sensitivity.markers

    signature <- paste(df.signatures[r,]$perturbagen, signature.type,sep = "_")

    if(n.res>2 & n.res>2){
      resistance.signature <- unlist(strsplit(df.signatures[r,]$proteins.resistance.markers,';'))
      sensitivity.signature <- unlist(strsplit(df.signatures[r,]$proteins.sensitivity.markers,';'))
      df.prot.res <- protein_data[protein_data$protein_acc %in% resistance.signature,]
      df.prot.sen <- protein_data[protein_data$protein_acc %in% sensitivity.signature,]

      if (nrow(df.prot.res)>2 & nrow(df.prot.sen)>2){

        resistance.markers[r] <- paste(unique(df.prot.res$protein_acc),collapse = ";")
        sensitivity.markers[r] <- paste(unique( df.prot.sen$protein_acc),collapse = ";")
        ks.p <- ks.test(df.prot.res[,2],df.prot.sen[,2])
        pvalue.ks[r] <- ks.p$p.value
        zscore.resistance[r] <- (mean(df.prot.res[,2])-mean.population)/sd.population
        zscore.sensitivity[r] <- (mean(df.prot.sen[,2])-mean.population)/sd.population
        perturbagens[r] <- df.signatures[r,]$perturbagen
        perturbagen.signature[r] <- signature
        #aa <-
        #bb <-
        aa <-
        bb <-
        av.rank.resistance[r] <-  mean(df.prot.res$protein_rank)
        av.rank.sensitivity[r] <- mean(df.prot.sen$protein_rank)

        geommean.rank.resistance[r] <- ReMEA::gm_mean(df.prot.res$protein_rank)
        geommean.rank.sensitivity[r] <- ReMEA::gm_mean(df.prot.sen$protein_rank)

        med.rank.resistance[r] <-  median(df.prot.res$protein_rank)
        med.rank.sensitivity[r] <- median(df.prot.sen$protein_rank)


        n.resistance.markers[r] <- nrow(df.prot.res)
        n.sensitivity.markers[r] <- nrow(df.prot.sen)
      }
    }
  }


  df <- data.frame(perturbagens,
                   perturbagen.signature,
                   zscore.resistance,
                   zscore.sensitivity,
                   pvalue.ks,
                   pvalue.bws,
                   resistance.markers,
                   sensitivity.markers,
                   av.rank.resistance,
                   av.rank.sensitivity,
                   geommean.rank.resistance,
                   geommean.rank.sensitivity,
                   med.rank.resistance,
                   med.rank.sensitivity,
                   n.resistance.markers,
                   n.sensitivity.markers,
                   signature.type)

  df <- subset(df,df$perturbagen.signature !="")


  return(df)
}


response_marker_enrichment_analysis <- function(protein_data,
                                                marker_type=c("crispr","drug","rnai"),
                                                tumour_type=c("haem","solid"),
                                                pvalue.cutoff=0.01){

  # in protein_data:
  #     - second column = folds (log2)
  #     - third column = pvalues
  #     - protein_acc column must be present, these are the XXX_HUMAN accessions to match to marker signatures

  library(ggplot2)
  library(ggrepel)
  library(foreach)
  library(doParallel)


  #tumour_type <- 'haem'
  #marker_type <- "crispr"

  signature_list <- ReMEA::list_of_signatures



  signatures_to_analyse <- foreach(m=marker_type, .combine="c")%do%{

    names(signature_list)[grepl(m,names(signature_list),fixed = T) &
                            grepl(tumour_type,names(signature_list),fixed = T)]
  }

  cl <- makeCluster(detectCores(logical = TRUE)-1)
  registerDoParallel(cl)

  df <- foreach(s = signatures_to_analyse,.combine = 'rbind')%dopar%{

    ReMEA::signature_enrichment(protein_data = protein_data ,
                                df.signatures = signature_list[[s]],
                                signature.type = s)

  }
  stopCluster(cl)


  df$delta.zscore <- df$zscore.sensitivity-df$zscore.resistance
  x1 <- log2(df$e.senstivity_in_increased)-log2(df$e.senstivity_in_decreased)
  x2 <- log2(df$e.resistance_in_decreased)-log2(df$e.resistance_in_increased)

  df$sensitivity_score <- x1-x2

  y1 <- df$counts.senstivity_in_increased-df$counts.senstivity_in_decreased
  y2 <- df$counts.resistance_in_decreased-df$counts.resistance_in_increased
  df$sensitivity_score_by_counts <- y1/y2

  a <- df$counts.senstivity_in_increased+df$counts.resistance_in_decreased
  b <- df$counts.senstivity_in_decreased+df$counts.resistance_in_increased
  df$relative_counts <- a/b

  df$total_counts <- a

  df$delta_rank <- df$av.rank.resistance-df$av.rank.sensitivity

  if ("drug" %in% marker_type){
    df.druginfo <- ReMEA::drug_info_from_Shelleckchem_pharmacoDB_plus_Finland
    df$perturbagens <- gsub("-",'.',df$perturbagens,fixed = T)
    df$perturbagens <- gsub("(",'.',df$perturbagens,fixed = T)
    df$perturbagens <- gsub(")",'.',df$perturbagens,fixed = T)
    df <- merge.data.frame(df,drug_info_from_Shelleckchem_pharmacoDB_plus_Finland,
                           by.x='perturbagens',by.y = 'drug')
    df <- data.frame(drug_target=df$drug.target,df)
    df$perturbagens <- df$drug_target


  }

  df$qvalue <- p.adjust(df$pvalue.bws,method = 'fdr')


  return(list(enrichment_data=df))



}



signature_enrichment <- function(protein_data,
                                 df.signatures,
                                 signature.type,
                                 pvalue.cutoff=0.01){


  library(BWStest)


  protein_data <- na.omit(protein_data)
  protein_data[,2] <- scale(protein_data[,2])
  protein_data$protein_rank <- rank(-protein_data[,2])



  mean.population <- mean(protein_data[,2])
  sd.population <- sd(protein_data[,2])
  background_population <- protein_data$protein_acc
  value.cutoff <- 0#mean.population+(sd.population*1.5)

  #protein_data$protein <- paste0(";", protein_data$protein)

  nr <- nrow(df.signatures)

  e.resistance_in_increased <- numeric(nr)
  e.senstivity_in_increased <- numeric(nr)
  e.resistance_in_decreased <- numeric(nr)
  e.senstivity_in_decreased <- numeric(nr)

  p.resistance_in_increased <- numeric(nr)
  p.senstivity_in_increased <- numeric(nr)
  p.resistance_in_decreased <- numeric(nr)
  p.senstivity_in_decreased <- numeric(nr)

  counts.resistance_in_increased <- numeric(nr)
  counts.senstivity_in_increased <- numeric(nr)
  counts.resistance_in_decreased <- numeric(nr)
  counts.senstivity_in_decreased <- numeric(nr)

  av.rank.resistance <- numeric(nr)
  av.rank.sensitivity <- numeric(nr)

  zscore.resistance <- numeric(nr)
  zscore.sensitivity <- numeric(nr)
  pvalue.ks <- numeric(nr)
  perturbagen.signature <- character(nr)
  perturbagens <- character(nr)

  pvalue.bws <- numeric(nr)

  q1 <- quantile(protein_data[,2])[2]
  q3 <- quantile(protein_data[,2])[4]

  df.prot.increased <- subset(protein_data,protein_data[,2]>q3 & protein_data[,3]<pvalue.cutoff)
  df.prot.decreased <- subset(protein_data,protein_data[,2]<q1 & protein_data[,3]<pvalue.cutoff)




  resistance.markers <- character(nr)
  sensitivity.markers <- character(nr)

r <- 1
  for (r in 1:nr) {
    n.res <- df.signatures[r,]$n.resistance.markers
    n.sen <- df.signatures[r,]$n.sensitivity.markers

    signature <- paste(df.signatures[r,]$perturbagen, signature.type,sep = "_")

    if(n.res>2 & n.res>2){

    resistance.signature <- unlist(strsplit(df.signatures[r,]$proteins.resistance.markers,';'))
    sensitivity.signature <- unlist(strsplit(df.signatures[r,]$proteins.sensitivity.markers,';'))

    df.prot.res <- protein_data[protein_data$protein_acc %in% resistance.signature,]
    df.prot.sen <- protein_data[protein_data$protein_acc %in% sensitivity.signature,]

        enrichment.res_up <- enrichment.of.markers.in.list(list.of.markers = resistance.signature,
                                                        list.of.peptides = df.prot.increased$protein_acc,
                                                        background.list = background_population,
                                                        signature.name = signature)

        enrichment.res_do <- enrichment.of.markers.in.list(list.of.markers = resistance.signature,
                                                           list.of.peptides = df.prot.decreased$protein_acc,
                                                           background.list = background_population,
                                                           signature.name = signature)


        enrichment.sen_up <- enrichment.of.markers.in.list(list.of.markers = sensitivity.signature,
                                                        list.of.peptides = df.prot.increased$protein_acc,
                                                        background.list = background_population,
                                                        signature.name = signature)

        enrichment.sen_do <- enrichment.of.markers.in.list(list.of.markers = sensitivity.signature,
                                                           list.of.peptides = df.prot.decreased$protein_acc,
                                                           background.list = background_population,
                                                           signature.name = signature)


        e.resistance_in_increased[r] <- enrichment.res_up$enrichment
        e.resistance_in_decreased[r] <- enrichment.res_do$enrichment

        e.senstivity_in_decreased[r] <- enrichment.sen_do$enrichment
        e.senstivity_in_increased[r] <- enrichment.sen_up$enrichment

        p.resistance_in_increased[r] <- enrichment.res_up$pvalue
        p.resistance_in_decreased[r] <- enrichment.res_do$pvalue

        p.senstivity_in_decreased[r] <- enrichment.sen_do$pvalue
        p.senstivity_in_increased[r] <- enrichment.sen_up$pvalue


        counts.resistance_in_increased[r] <- enrichment.res_up$counts
        counts.resistance_in_decreased[r] <- enrichment.res_do$counts

        counts.senstivity_in_decreased[r] <- enrichment.sen_do$counts
        counts.senstivity_in_increased[r] <- enrichment.sen_up$counts

        resistance.markers[r] <- enrichment.res_do$proteins
        sensitivity.markers[r] <- enrichment.sen_up$proteins


    if (nrow(df.prot.res)>2 & nrow(df.prot.sen)>2){


      ks.p <- ks.test(df.prot.res[,2],df.prot.sen[,2])
        pvalue.ks[r] <- ks.p$p.value

        pvalue.bws[r] <-  BWStest::bws_test(df.prot.res[,2],df.prot.sen[,2],
                                            method = 'BWS',
                                            alternative = 'two.sided')$p.value

        zscore.resistance[r] <- (mean(df.prot.res[,2])-mean.population)/sd.population
        zscore.sensitivity[r] <- (mean(df.prot.sen[,2])-mean.population)/sd.population


        perturbagens[r] <- df.signatures[r,]$perturbagen
        perturbagen.signature[r] <- signature

        aa <-  median(df.prot.sen$protein_rank)
        bb <- median(df.prot.res$protein_rank)


        #aa <- ReMEA::gm_mean(df.prot.sen$protein_rank)
        #bb <- ReMEA::gm_mean(df.prot.res$protein_rank)

        av.rank.resistance[r] <- bb
        av.rank.sensitivity[r] <- aa



    }


      }
    }


  df <- data.frame(perturbagens,
                   perturbagen.signature,
                   e.resistance_in_decreased,
                   e.resistance_in_increased,
                   e.senstivity_in_decreased,
                   e.senstivity_in_increased,
                   counts.resistance_in_decreased,
                   counts.resistance_in_increased,
                   counts.senstivity_in_decreased,
                   counts.senstivity_in_increased,
                   p.resistance_in_decreased,
                   p.resistance_in_increased,
                   p.senstivity_in_decreased,
                   p.senstivity_in_increased,
                   zscore.resistance,
                   zscore.sensitivity,
                   pvalue.ks,
                   pvalue.bws,
                   resistance.markers,
                   sensitivity.markers,
                   av.rank.resistance,
                   av.rank.sensitivity,
                   signature.type)

df <- subset(df,df$perturbagen.signature !="")


return(df)
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

enrichment.of.markers.in.list <- function(list.of.markers, #signature
                                          list.of.peptides,
                                          background.list,
                                          split.lists=F,
                                          signature.name
                                          ){


  if (split.lists==T){
    background.list <- trimws(unlist(strsplit(background.list,";")))
    list.of.peptides <- trimws(unlist(strsplit(list.of.peptides,";")))
    list.of.markers <- trimws(unlist(strsplit(list.of.markers,";")))
  }
  q <- length(intersect(list.of.markers,list.of.peptides))
  j <- length(intersect(list.of.markers,background.list))
  m <- length(background.list)
  k <- length(list.of.peptides)

  if (q>0){
    n <- m-j
    prots <- intersect(list.of.markers,list.of.peptides)
    pvalue <- 1-phyper(q-1,j,n,k,lower.tail = TRUE, log.p = F)
    enrichment <-(q/k)/(j/m)
    counts <- q
    data.size <- k
    counts.bg<- j
    bg.size <- m
    proteins <- paste(prots,collapse = ";")
  }else{

    pvalue <- 1
    enrichment <-1
    counts <- 0
    data.size <- k
    counts.bg<- j
    bg.size <- m
    proteins <- ""

  }

  results <- data.frame(signature.name, pvalue,enrichment,counts, data.size,counts.bg,bg.size, proteins)
  return(results)
}
