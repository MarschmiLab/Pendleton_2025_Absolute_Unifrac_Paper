subset_dist_matrix <- function(dist_mat, keeps){
  
  mat <- as.matrix(dist_mat)
  
  return(as.dist(mat[keeps, keeps]))
  
}
