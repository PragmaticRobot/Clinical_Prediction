reorder_LS <- function(LS_pred)
{
  outback <- data.frame(matrix(vector(),dim(LS_pred[[2]])[1],100))
  for (j in 1:100) # loop over cross-validations
  {
    for (i in 1:dim(LS_pred[[2]])[1])
    {
      if (is.na(LS_pred[[2]][i,j])) {
        next 
      }
      else {
        outback[LS_pred[[2]][i,j],j] <- LS_pred[[3]][i,j]
      }
    }
  }
  
  return(outback)
}

reorder_RF <- function(RF_pred)
{
  outback <- data.frame(matrix(vector(),dim(RF_pred[[3]])[1],100))
  for (j in 1:100) # loop over cross-validations
  {
    for (i in 1:dim(RF_pred[[3]])[1])
    {
      outback[RF_pred[[2]][i,j],j] <- RF_pred[[3]][i,j]
    }
  }
  return(outback)
}

