#' Build network between identifications based on kegg network database.
#'
#' @param kegg_id a vector of strings indicating KEGG ID of putative ID.
#' @return a binary matrix of network of KEGG IDs.

get_kegg_network <- function(kegg_id){

  ## initialize Wk
  c <- length(kegg_id)
  Wk <- Matrix::Matrix(0,ncol=c,nrow=c,sparse = TRUE)

  ## deal with muiltiple KEGG IDs and missing KEGG IDs
<<<<<<< HEAD
  kegg_id <- strsplit(kegg_id, ' ')
  kegg_id <- lapply(kegg_id, unique)
  ids <- unlist(kegg_id)
  indexes <- rep(seq_along(kegg_id), lengths(kegg_id))

  ## subset useful kegg_network and replace them with indexes
  kegg_network<-get("kegg_network")
=======
  ids <- c()
  indexes <- c()
  for (i in 1:c){
    if (kegg_id[i]!=''){
      ## multiple IDs
      if(grepl(' ',kegg_id[i])){
        add <- unique(strsplit(kegg_id[i],' ')[[1]])
        ids <- c(ids,add)
        indexes <- c(indexes,rep(i,length(add)))
      } else {
        ids <- c(ids,kegg_id[i])
        indexes <- c(indexes,i)
      }
    }
  }

  ## subset useful kegg_network and replace them with indexes
>>>>>>> eb66f9edd10e572e7e0edd24e2e01afcc3cdfdf5
  sub_netdb <- kegg_network[(kegg_network$r1 %in% ids)&(kegg_network$r2 %in% ids),]
  if (nrow(sub_netdb)==0){
    return(Wk)
  }
  for (i in 1: dim(sub_netdb)[1]){
    sub_netdb$r1[i] <- indexes[which(sub_netdb$r1[i]==ids)][1]
    sub_netdb$r2[i] <- indexes[which(sub_netdb$r2[i]==ids)][1]
  }

  ## turn index pairs into matrix
  obj <- igraph::graph_from_data_frame(sub_netdb,directed = FALSE)
  w <- igraph::as_adjacency_matrix(obj,type = 'both')
  tmp <- as.integer(colnames(w))
  Wk[tmp,tmp] <- w

  return(Wk)
}





#' Build network between identifications based on tanimoto score.
#'
#' @param pubchem_cid a vector of strings indicating PubChem CID of putative ID.
#' @return a binary matrix of network of tanimoto scores.

get_tani_network <- function(pubchem_cid){

  ### initialize W
  t <- pubchem_cid
  c <- length(t)
  Wt <- Matrix::Matrix(0,ncol=c,nrow=c,sparse = TRUE)

  ##### Turn t(cid) to fpset
  fpset <- list()
  size <- 100
  groups <- as.integer(c/size)
  remain <- c%%size
  if (groups != 0){
    for (g in 1:groups){
      fpset[[g]] <- get_sdf(t[(1+(g-1)*size):(g*size)])
    }
  }
  if (remain!=0){
    fpset[[groups+1]] <- get_sdf(t[(c-remain+1):c])
    groups <- groups+1
  }

  ##### get tanimoto score
  for (g1 in 1:groups){
    for (g2 in g1:groups){
      x <- fpset[[g1]]
      s1 <- length(x)
      y <- fpset[[g2]]
      s2 <- length(y)
      ind1 <- (g1-1)*100+1
      ind2 <- (g2-1)*100+1

      if (g1==g2){
        for(i in 1:s1){
          Wt[(ind1+i-1),(ind1+i-1):(ind1+s1-1)] <- fpSim(x=x[i],y=x[i:s1],method='Tanimoto')
        }
      }
      else{
        ## g1!=g2
        for(i in 1:s1){
          Wt[(ind1+i-1),ind2:(ind2+s2-1)] <- fpSim(x=x[i],y=y,method='Tanimoto')
        }
      }
    }
  }

  return(Wt+t(Wt))
}



################################ helper function ###########################

get_sdf <- function(ids){
  compounds <- getIds(ids)
  cid(compounds) <- sdfid(compounds)
  return(fp2bit(compounds))
}
<<<<<<< HEAD
=======


>>>>>>> eb66f9edd10e572e7e0edd24e2e01afcc3cdfdf5
