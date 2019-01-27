#Stick Breaking Simulations
levels <- seq(1, 50, by = 0.01)
result <- rep(0, length(levels))
for(i in 1:length(levels)){
  beta = levels[i]
  N = 10000
  W <- rgamma(N, 1, beta)
  # beta^1 / gamma(1) x^(1 - 1) exp(-beta*x)
  mean(W)
  result[i] <- max(W)
}
plot(seq(1, 50, by = 0.01), result)
abline(h = 1, col = "red")
min(which(result <= 1))



setwd("C:/Users/Tom Kennes/Desktop/Research Master/Thesis/Exploratory Coding")
pdf("stick breaking with various beta levels.pdf")
levels <- seq(1, 25, by = 1)
for(i in 1:length(levels)){
  beta = levels[i]
  N = 100
  print(beta)
  W <- rgamma(N, 1, beta)
  # beta^1 / gamma(1) x^(1 - 1) exp(-beta*x)
  V <- rep(0, (N + 1))
  V[1] <- W[1]
  for(i in 2:N){
    V[i] <- W[i]*V[i - 1]/W[i - 1]*(1 - W[i - 1])
  }
  V[N + 1] <- 1 - sum(V[1:N])
  plot(V, main = paste("beta = ", beta))
}
dev.off()





beta = 15
N = 100
hist(stick_breaking(N, beta = beta)$V, breaks = 50)
stick_breaking <- function(N, alpha = 1, beta){
  W <- rgamma(N, alpha, beta)
  # beta^1 / gamma(1) x^(1 - 1) exp(-beta*x)
  V <- rep(0, (N + 1))
  V[1] <- W[1]
  for(i in 2:N){
    V[i] <- W[i]*V[i - 1]/W[i - 1]*(1 - W[i - 1])
  }
  V[N + 1] <- 1 - sum(V[1:N])
  ret <- list(W, V)
  names(ret) <- c("W", "V")
  return(ret)
}

M <- 10000
F_ <- rep(0, M)
p = 0.5
for(j in 1:M){
  for(i in 1:N){
    U <- runif(1,0,1)
    if(U >= p){
      F_[j] <- F_[j] + V[i]
    }
  }
  print(j/M)
}
hist(F_, breaks = 100, probability = T)

