# function to clean names
cleanVC <- function(x){
  x <- tolower(trimws(unlist(strsplit(x, ','))))
  return(sort(unique(x)))
}  
