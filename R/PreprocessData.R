#' Preprocess input file.
#'
#' @param filename the name of the file which the data are to be read from. Its type should be chosen
#'             in 'type' parameter. Also, it should have columns named exactly as 'metid' (IDs for peaks),
#'             'query_m.z' (query mass of peaks), 'exact_m.z' (exact mass of putitative IDs),
#'             'kegg_id' (IDs of putitative IDs from KEGG Database), 'pubchem_cid' (CIDs of putitative IDs
#'             from PubChem Database). Otherwise, this function would not work.
#' @param type string indicating the type of the file. It can be a 'data.frame' which is already loaded
#'             into R, or some other types like a csv file.
#' @param na a character vector of strings which are to be interpreted as NA values.
#' @param sep a character value which seperates multiple IDs in kegg_id or pubchem_cid field, if there
#'            are multiple IDs.
#' @return get_cleaned returns a list containing the following components:
#'        \item{df}{a data frame which is the original input data.}
#'        \item{clean_data}{a data frame with unuseful observations and features removed.}
#'        \item{mass}{a data frame with unique query peak, along with query mass.}
#'        \item{ID}{a data frame with unique putitative IDs, along with PubChem ID, KEGG ID, exact mass.}
#'        \item{index_na}{a vector of row indexes which contains NA values.}

get_cleaned <- function(filename, type = c('data.frame','csv','txt'), na, sep){


  ## Accept different types of files
  if (type == 'data.frame'){
    df <- filename
  } else if (type == 'csv'){
      df <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE, na.strings = na)
  } else if (type == 'txt'){
        df <- read.table(filename, header = TRUE, stringsAsFactors = FALSE, na.strings = na)
  } else {
          stop('Please provide a data frame in R or a csv/text file!')
    }
  # Check column names of data
  if (!all(c("query_m.z", "exact_m.z", "kegg_id", "pubchem_cid") %in% colnames(df))){
    stop('Please provide features needed, also with right column names!')
  }


  ## Load data
  qmass <- df$query_m.z
  emass <- df$exact_m.z
  kid <- df$kegg_id
  cid <- df$pubchem_cid
  # Add metid column
  metid <- rep(0,length(qmass))
  mass_set <- qmass[!duplicated(qmass)]
  for (i in 1:length(qmass)){
      metid[qmass==mass_set[i]] <- i
  }
  # Format id columns
  kid <- format_id(kid, na = na, sep = sep)
  cid <- format_id(cid, na = na, sep = sep)
  cid <- as.numeric(get_first(cid))
  # Get inchikey
  inchikey <- get_inchikey(cid)
  df$inchikey <- inchikey
  df$metid <- metid


  ## Format data
  df_formated <- data.frame(metid,qmass,emass,kid,cid,inchikey,stringsAsFactors = FALSE)
  # remove NAs
  index_empty <- which(is.na(cid))
  df_formated <- df_formated[!is.na(cid),]
  # Combine rows with same inchikey
  metid <- unique(df_formated$metid)
  compound <- df_formated[df_formated$metid==metid[1],]
  df_combined <- combine_inchikey(compound)
  if (length(metid)>1){
    for (i in 2:length(metid)){
      compound <- df_formated[df_formated$metid==metid[i],]
      df_combined <- rbind(df_combined,combine_inchikey(compound))
    }
  }


  ## get mass data -- with m query_m.z
  mass <- df_combined[!duplicated(df_combined$metid),]
  mass <- subset(mass,select = c(metid,qmass))


  ## get ID data -- with c identifications
  ID <- df_combined[!duplicated(df_combined$cid),]
  ID <- subset(ID,select = c(metid,kid,cid,emass))

  return(list(df=df, clean_data = df_formated, mass = mass, ID = ID, index_na = index_empty))
}






################################## helper functions ##################################

format_id <- function(id, na, sep){
  # This function formats id columns.
  id[is.na(id)] <- ''
  for (i in 1:length(id)){
    id[i] <- stringr::str_replace_all(id[i], pattern = na,'')
    if (startsWith(id[i], sep)){id[i] <- stringr::str_replace(id[i], pattern = sep, '')}
    id[i] <- stringr::str_replace_all(id[i], pattern = ' ','')
    id[i] <- stringr::str_replace_all(id[i], pattern = paste0(sep,sep,'|',sep), ' ')
    id[i] <- stringr::str_trim(id[i])
  }
  return(id)
}


get_first <- function(id){
  # This function gets the first id if there are multiple ids.
  for (i in 1:length(id)){
    if (id[i]==''){
      id[i] <- ''
    } else {
      temp <- strsplit(id[i],' ')[[1]]
      id[i] <- temp[1]
    }
  }
  return(id)
}


get_inchikey <- function(id){
  # This function gets the fist 14 charactors of InchiKey from PubChem database.
  n <- length(id)
  response <- rep(NA,n)
  for (i in 1:n){
    if (!is.na(id[i])){
      if (!id[i] %in% InchiKey$CID){
        response[i] <- as.character(id[i])
      }
      else{
        response[i] <- substr(InchiKey$InchiKey[which(InchiKey$CID %in% id[i])],1,14)
      }
    }
  }
  return(response)
}


combine_inchikey <- function(compound){
  # This function combines identifications with same InchiKey.
  group <- compound$inchikey[!duplicated(compound$inchikey)]
  combined <- compound[1:length(group),]
  for (j in 1:length(group)){
    sub <- compound[compound$inchikey==group[j],]
    combined[j,] <- sub[1,]
    combined$kid[j] <- stringr::str_trim(paste(sub$kid, collapse = ' '))
  }
  return(combined)
}

