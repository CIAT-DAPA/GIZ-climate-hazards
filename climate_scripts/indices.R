# Agro-climatic indices main functions
# H. Achicanoy and A. Esquivel
# CIAT, 2020

rsum.lapply <- function(x, n=3L) # Calculate rollin sum
{
  lapply(1:(length(x)-n+1), function(i)
  {
    # Sum for n consecutive days
    z <- sum(x[i:(i+n-1)])
    # Indices used to calculate the sum
    seq.sum <- as.numeric(i:(i+n-1))
    # List with SUM and INDICES
    results <- list(z, seq.sum)
    return(results)
  })
}
cumulative.r.sum <- function(results){ unlist(lapply(results, function(x){z <- x[[1]]; return(z)})) } # Extract the SUM
is.leapyear <- function(year){ return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)) } # Function to identify leap years


## CDD. Maximum number of consecutive dry days
calc_cdd <- function(PREC, p_thresh=1){
  runs <- rle(PREC < p_thresh)
  cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
  return(cons_days)
}
calc_cddCMP <- cmpfun(calc_cdd)

## P5D. Maximum 5-day running average precipitation
calc_p5d <- function(PREC){
  runAvg <- caTools::runmean(PREC, k=5, endrule='NA')
  runAvg <- max(runAvg, na.rm=TRUE)
  return(runAvg)
}
calc_p5dCMP <- cmpfun(calc_p5d)

## NT35. Number of days with max. temperature above 35?C
calc_hts <- function(TMAX, t_thresh=35) {
  hts <- length(which(TMAX >= t_thresh))
  return(hts)
}
calc_htsCMP <- cmpfun(calc_hts)

## P95. 95th percentile of daily precipitation
calc_p95 <- function(PREC){
  quantile(PREC, probs = .95, na.rm = T)
}
calc_p95CMP <- cmpfun(calc_p95)
