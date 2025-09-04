#Here are some modified functions for UpSetR
#from: https://github.com/hms-dbmi/UpSetR/issues/85
#by: docmanny
fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}
get_intersect_members <- function (x, ...){
  require(dplyr)
  require(tibble)
  require(rlang)
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n))
  #ar[[1]] <- x
  i=1
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  filter_str<-paste(unlist(ar),collapse=" & ")
  filter_expr<-parse_expr(filter_str)
  x %>% filter(!!filter_expr) %>% column_to_rownames() -> x
  return(x)
}

#This function takes a gene name, a matrix of counts from edgeR and the table of sample data.
#the output is a table of count data across all samples for that gene.
#Useful to get plotting data for a given gene.
gene_data<-function(x,gene,sample_dat){
  counts<-x[gene,]
  sample_dat$counts<-as.vector(counts[match(sample_dat$sample,names(counts))])
  return(sample_dat)
}
