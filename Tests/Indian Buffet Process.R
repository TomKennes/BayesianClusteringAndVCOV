#Indian Buffet Process

rbernouilli <- function(n, p){
  u <- runif(n, 0, 1)
  x <- rep(0, n)
  x[u > p] = 1
  return(x)
}

alpha = 10
K <- 10
N <- 10
random_binary_matrix <- function(alpha, N, K){
  pi <- rbeta(N, alpha/K, 1)
  result <- matrix(0, nrow = N, ncol = K)
  for(i in 1:N){
    result[i,] <- rbernouilli(K, pi[i])
  }
  return(result)
}

random_binary_triangular_matrix <- function(alpha, N, K){
  pi <- rbeta(N, alpha/K, 1)
  result <- matrix(0, nrow = N, ncol = K)
  for(i in 1:N){
    result[i,] <- rbernouilli(K, pi[i])
  }
  for(i in 1:(N- 1)){
    for(j in (i + 1):K){
      result[i,j] = 0
    }
  }
  
  return(t(result))
}
