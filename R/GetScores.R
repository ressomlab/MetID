#' Get scores for metabolite putative IDs by LC-MS .
#'
#' @param filename the name of the file which the data are to be read from. Its type should be chosen
#'             in 'type' parameter. Also, it should have columns named exactly 'metid' (IDs for peaks),
#'             'query_m.z' (query mass of peaks), 'exact_m.z' (exact mass of putative IDs),
#'             'kegg_id' (IDs of putative IDs from KEGG Database), 'pubchem_cid' (CIDs of putative IDs
#'             from PubChem Database). Otherwise, this function would not work.
#' @param type string indicating the type of the file. It can be a 'data.frame' which is already loaded
#'             into R, or some other specified types like a csv file.
#' @param na a character vector of strings which are to be interpreted as NA values.
#' @param sep a character value which seperates multiple IDs in kegg_id or pubchem_cid field, if there
#'            are multiple IDs.
#' @param mode string indicating the mode of metabolites. It can be positive mode (POS) or negative mode
#'             (NEG).
#' @param Size an integer which indicates sample size in Gibbs sampling.
#' @param delta a hyper-parameter representing the mean value of mass ratio.
#' @param gamma_mass a hyper-parameter representing the accuracy of mass measurement.
#' @param iterations ask user to input number of interations,default 500
#' @return A dataframe which contains input data together with a
#'         column of scores in the end. In the
#'         score column, if the row contains NA values or does not has a PubChem cid, the score would be
#'         '-', which stands for missing value. Otherwise, each score would be from 0 to 1.
#' @examples
#' ## check if colnames of dataset meet requirement
#' names(demo1)
#' ## change colnames
#' colnames(demo1) <- c('query_m.z','name','formula','exact_m.z','pubchem_cid','kegg_id')
#' ## get scores
#' out <- get_scores_for_LC_MS(demo1, type = 'data.frame', na='-', mode='POS')
#'
#' @export
#' @import ChemmineR
#' @importFrom stringr str_replace_all str_replace str_trim
#' @importFrom igraph graph_from_data_frame as_adjacency_matrix
#' @importFrom Matrix Matrix t
#' @importFrom stats dnorm rmultinom
#' @importFrom utils data read.csv read.table write.csv
#' @import devtools

get_scores_for_LC_MS <- function(filename, type = c('data.frame','csv','txt'), na = 'NA', sep = ';',
                                 mode = c('POS','NEG'), Size=2000, delta=1, gamma_mass=10,iterations=500){

  ## Preprocess data
  list_from_get_cleaned <- get_cleaned(filename, type = type, na = na, sep = sep)
  df <- list_from_get_cleaned$df
  mass <- list_from_get_cleaned$mass
  ID <- list_from_get_cleaned$ID


  ## Build identification network
  message("Start building network: it may take several minutes......")
  Wk <- get_kegg_network(ID$kid)
  Wt <- get_tani_network(ID$cid)
  Wt[Wt>=0.7] <- 1
  Wt[Wt<0.7] <- 0
  W <- pmax(Wt,Wk)


  ## Gibbs sampling
  # adjust mass according to mode
  proto_mass = 1.00727646677
  if (mode=='POS'){
    qmass<-mass$qmass-proto_mass
  } else if (mode=='NEG'){
    qmass<-mass$qmass+proto_mass
  }

  # initialize Z
  m <- dim(mass)[1]
  c <- dim(ID)[1]
  Z <- Matrix::Matrix(1,nrow=c,ncol=m,sparse=TRUE)

  # I: binary matrix of identification~mass
  I <- Matrix::Matrix(0,nrow=c,ncol=m)
  for (i in 1:m){
    I[which(ID$metid==mass$metid[i]),i] <- 1
  }

  # load m.z data
  X <- matrix(qmass,nrow=1)
  Y <- as.matrix(as.numeric(ID$emass))

  # Gibbs Samplings -- burn-in
  message('Start getting random samples: it may take several minutes......')
  for (s in 1:iterations){
    beta_temp <- as.vector(W%*%Z%*%Matrix(1,ncol=1,nrow=m,sparse=TRUE))
    beta <- Matrix::Matrix(beta_temp,ncol=m,nrow=c,sparse=TRUE)-W%*%Z
    beta_sum <- Matrix::Matrix(apply(beta,1,sum),ncol=m,nrow=c,sparse=TRUE)
    prior <- (delta+beta)/(c*delta+beta_sum)
    post <- (dnorm((1/Y)%*%X,delta,gamma_mass/3*10^(-6))*prior)*I
    post <- t(t(post)/apply(post,2,sum))
    Z <- Matrix::Matrix(apply(post,2,function(x){rmultinom(1,1,x)}),sparse=TRUE)
  }

  # Gibbs Samplings
  prob <- Matrix::Matrix(0,nrow=c,ncol=m,sparse=TRUE)
  for (s in 1:(Size-iterations)){
    beta_temp <- as.vector(W%*%Z%*%Matrix(1,ncol=1,nrow=m,sparse=TRUE))
    beta <- Matrix(beta_temp,ncol=m,nrow=c,sparse=TRUE)-W%*%Z
    beta_sum <- Matrix::Matrix(apply(beta,1,sum),ncol=m,nrow=c,sparse=TRUE)
    prior <- (delta+beta)/(c*delta+beta_sum)
    post <- (dnorm((1/Y)%*%X,delta,gamma_mass/3*10^(-6))*prior)*I
    post <- t(t(post)/apply(post,2,sum))
    Z <- Matrix::Matrix(apply(post,2,function(x){rmultinom(1,1,x)}),sparse=TRUE)
    prob <- Z+prob
  }
  prob <- prob/(Size-iterations)


  ## add score column to the original file
  message('Start writing scores......')
  index_empty <- list_from_get_cleaned$index_na
  df_dup <- list_from_get_cleaned$clean_data
  df$score <- rep(0,dim(df)[1])
  df$score[index_empty] <- '-'
  for (i in 1:m){
    subdf <- df_dup[df_dup$metid==mass$metid[i],]
    inchikeys <- subdf[!duplicated(subdf$inchikey),]$inchikey
    p <- prob[,which(mass$metid==mass$metid[i])]
    p <- p[p!=0]
    for (j in 1:length(inchikeys)){
      ind <- as.numeric(rownames(subdf[subdf$inchikey==inchikeys[j],]))
      df$score[ind] <- p[j]
    }
  }
  df$score[is.na(df$score)] <- 0
  message('Completed!')


  write.csv(df,file = 'scores.csv')
  return(df)
}
