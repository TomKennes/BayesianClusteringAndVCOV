#Following the Generalization on Wikipedia
alpha = 0 #Discount
theta = 0.5 #Strength

B <- c(1)
N <- 20000

for(i in 2:N){
  new = (theta + length(B)*alpha)/((i - 1) + theta)
  old = (B - alpha)/((i - 1) + theta)
  p <- c(old, new)
  draw <- rmultinom(1,1, p)
  if(draw[length(draw)] == 1){
    B <- c(B, 1)
  } else {
    B[which(draw == 1)] = B[which(draw == 1)] +  1
  }
  print(B)
}






#First we generate data distributed over 3 clusters
x <- sort(rep(c(1,8), 100))
y <- -x + rnorm(length(x), 0, 1)
x <- x + rnorm(length(x), 0, 1)
x <- sample(x)
plot(x)
# Now the algorithm has to figure out how many clusters there are
posterior <- function(theta, sigma, y){
  N <- length(y)
  #-2log(y_i - theta_1) probabilty that y is from theta_1
  mat <- vapply(y, theta)
  like <- -2*N*sum(log(mat))
  mat = rowsum(mat^(-2))
  like <- like - length(theta)*log(mat)

}




vapply <- function(m, x){
  m_ <- matrix(0, nrow = length(m), ncol = length(x))
  for(i in 1:length(x)){
    m_[,i] <- m - x[i]
  }
  return(m_)
}


normal_loglik <- function(y, mu, sigma){
  ret <- -1/2*log(2*pi) - log(sigma) - 1/(2 * sigma^2)*sum((y - mu)^2)
  return(ret)
}


assign_new_cluster <- function(B, draws, theta, alpha){
  new = (theta + length(B)*alpha)/((draws - 1) + theta)
  old = (B - alpha)/((draws - 1) + theta)
  p <- c(old, new)
  return(p)
}


p <- assign_new_cluster(B, draws, theta, alpha)
weights_assign(p)



