




pathway_enrichment_remea <- function(genes,
                               background.genes,
                               prot_dbs=c("kegg","hallmark.genes","nci","process"))
  {

  library(foreach)
  library(doParallel)

  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  t1 <- Sys.time()
  enrich.up <- foreach(db = prot_dbs, .combine = "rbind")%dopar%{

    e <- protools2::enrichment.from.list(list.of.peptides=genes,
                                         background.genes,
                                         prot_db = db)
    if (nrow(e)>0){
      return(e)
    }
  }
  stopCluster(cl)
  t2 <- Sys.time()

  return(enrich.up)


}
