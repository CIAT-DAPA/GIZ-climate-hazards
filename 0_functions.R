# function to clean names
cleanVC <- function(x){
  x <- unlist(strsplit(x, ','))
  return(tolower(trimws(x)))
}  
