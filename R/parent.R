#' @param temp

parent =  function(temp)
{
  p= rep(1,length(temp$pos))
  p[1] = 0
  for(i  in 2:length(temp$pos))
  {
    ind = which(temp$pos[1:(i-1)] == (temp$pos[i]-1))
    p[i] =  ind[length(ind)]
  }
  return(p)
}