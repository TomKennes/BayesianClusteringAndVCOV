#Functions 
rnormal_inverse_gamma <- function(s, S, y, m, tau){
  if(is.null(ncol(y))){
    y <- matrix(y, ncol = length(y))
  }
  mu <- matrix(rep(0, ncol(y)*nrow(y)), ncol = ncol(y), nrow = nrow(y))
  V <- array(0, dim = c(ncol(y), ncol(y), nrow(y)))
  for(i in 1:nrow(y)){
    S <- S + ((y[i,] - m)^2)/(1 + tau)
    tmp <- 1/(rgamma(ncol(y), (1 + s)/2, S))
    V[, , i] <- diag(tmp)
    x <- (m + tau*y[i,])/(1 + tau)
    mu[i, ] <-  mvrnorm(1, mu = x, Sigma = tau/(1 + tau)*V[, ,i])
  }
  ret <- list(mu, V)
  names(ret) <- c("mu", "V")
  return(ret)
}


mvc_ <- function(s, p){
  ret <- gamma((p + s)/2)*(gamma(s/2)^(-1))*s^(-p/2)
  return(ret)
}


# Student T distributions!!
q_0 <- function(s, y, m, M){
  inv_M <- solve(M)
  det_M <- det(M)
  if(is.matrix(y)){
    p <- ncol(y)
    n <- nrow(y)
    q <- rep(mvc_(s, p)*det_M^(-1/2), nrow(y))
    for(i in 1:nrow(y)){
      q[i] <- q[i] * (1 + (1/s)*(t(y[i,] - m)) %*% inv_M %*% (y[i,] - m))^(-(p+s)/2)
    }
    return(q)
  } else {
    p <- length(y)
    q <- rep(mvc_(s, p)*det_M^(-1/2), 1)
    q <-q * (1 + (1/s)*(t(y - m)) %*% inv_M %*% (y - m))^(-(p+s)/2)
    return(q)
  }
  
}


q_j <- function(y, mu, V){
  det_V <- det(V)
  inv_V <- solve(V)
  ret <- exp(-(1/2)*(t(y - mu) %*% inv_V %*%(y - mu)))*(2 * det_V)^(-1/2)
  return(ret)
}


clustered <- function(theta1, theta2, epsilon){
  epsilon <- rep(epsilon, length(theta1))
  if(sum((theta1 - theta2)^2 > epsilon) == 0){
    return(T)
  } else {
    return(F)
  }
}


n_levels <- function(x){
  levs <- levels(as.factor(x))
  levs <- cbind(levs, rep(0, length(levs)))
  for(i in 1:nrow(levs)){
    count = 0
    for(j in 1:length(x)){
      if(as.character(x[j]) == levs[i,1]){
        count = count + 1
      }
    }
    levs[i,2] <- count
  }
  return(levs)
}





mat_mean <- function(x){
  d <- dim(x)
  mat <- matrix(0, nrow = d[1], ncol = d[2])
  for(i in 1:d[1]){
    for(j in 1:d[2]){
      mat[i,j] <- mean(x[i,j,])
    }
  }
  return(mat)
}

colMins <- function(x){
  container <- rep(0, ncol(x))
  for(i in 1:ncol(x)){
    container[i] <- min(x[,i])
  }
  return(container)
}


# Clustering a la escobar and west 1995
y <- sort(rep(seq(-3,3,by = 2), 30))
x <- sort(rep(seq(-3,3,by = 2), 30))
x <- sample(x)
y <- sample(y)
y <- y + rnorm(length(y), 0, 0.1)
x <- x + rnorm(length(x), 0, 0.1)
y <- cbind(x,y)
for(i in 1:8){
  x <- sort(rep(seq(-3,3,by = 1), 100))
  x <- x + rnorm(length(x), 0, 0.1)
  y <- cbind(x,y)
}
y <- y[sample(1:nrow(y)),]
plot(y)
abline(h = 0)
abline(v = 0)


if(!("MASS" %in% installed.packages())){
  install.packages("MASS")
  library("MASS")
} else {
  library("MASS")
}

s = 12
S = 0.1
tau = 0.5
m <- colMins(y)
M <- diag(rep(1, ncol(y)))
N <- 100    # Number of Iterations in gibbs
y <- as.matrix(y)
M_y <- ncol(y)
N_y <- nrow(y)
iterations_m <- array(0, dim = c(N_y, M_y, N))
iterations_V <- array(0, dim = c(M_y, M_y, N_y, N))
tmp <- rnormal_inverse_gamma(s, S, y, m, tau)
iterations_m[, , 1] <- tmp$mu
iterations_V[, , , 1] <- tmp$V
clusters <- matrix(1, ncol = N_y, nrow = N)
clusters[1,] <- 1:N_y
nclusters <- rep(0, N)
nclusters[1] <- length(levels(as.factor(clusters[1,])))
alpha = 1e-7
colnames(clusters) <- rep(paste("y", as.character(seq(1, 700, by = 1)), sep = ""))

for(k in 2:N){
  iterations_m[, ,k] <- iterations_m[, , k - 1]
  iterations_V[, , , k] <- iterations_V[, , ,k - 1]
  clusters[k,] <- clusters[k-1,]
  for(i in 1:N_y){
     #Set Up Probalities
     qj <- rep(0, nrow(y))
     for(j in 1:nrow(y)){
       qj[j] <- q_j(y[i,], iterations_m[j, , k], iterations_V[, ,j, k])
     }
     qj[i] <- alpha*q_0(s, y[i,], iterations_m[i, ,k], M)
     qj <- qj/sum(qj)
     
     #Cluster Probabilities
     cluster_names <- sort(as.numeric(levels(as.factor(clusters[k,]))))
     markov <- rep(0, length(cluster_names))
     for(j in 1:length(cluster_names)){
       markov[j] <- sum(qj[clusters[k, ] == cluster_names[j]]) - qj[i]*(cluster_names[j] == clusters[k,i])
     }
     if(sum(clusters[k, ] == clusters[k,i]) > 1){
       markov <- c(markov, qj[i])
       cluster_names <- c(cluster_names, max(as.numeric(clusters[k,])) + 1)
     }
     markov <- markov/sum(markov)
     
     transition <- rmultinom(1, 1, markov)
     transition <- cluster_names[which(transition == 1)]
     clusters[k,i] <- transition
     if(transition == max(clusters[k,])){
       # Meaning we have created the start of a new cluster
       
       # How are we defining the new cluster?
       tmp <- rnormal_inverse_gamma(s, S, y[i,], y[i,], tau)
       iterations_m[i, , k] <- tmp$mu
       iterations_V[, , i, k] <- tmp$V
     } else {
       # Meaning that we are adding it to another cluster
       y_clust <- y[clusters[k, ] == clusters[k,i],]
       tmp <- rnormal_inverse_gamma(s, S, y[i,], colMeans(y_clust), tau)
       iterations_m[i, ,k] <- tmp$mu
       iterations_V[, ,i, k] <- tmp$V
     }
  }
  print(clusters[k,])
}

type <- levels(as.factor(clusters[k,]))
plot(y[clusters[k,] == type[1],1], y[clusters[k,] == type[1],2],
     ylim <- c(min(y)*1.1, max(y)*1.1),
     xlim <- c(min(x)*1.1, max(x)*1.1),
     type = "l", 
     ylab = "y",
     xlab = "observations",
     col = "black")
for(i in 2:length(type)){
  rgb_ <- runif(n = 3, min= 0, max = 255)
  points(y[clusters[k,] == type[i],1], y[clusters[k,] == type[i],2], type = "l", col = 
           rgb(rgb_[1], rgb_[2], rgb_[3], maxColorValue = 255))
}
write.csv(clusters, file = "clusters.csv")

nclusters <- rep(0, N)
for(i in 1:N){
  nclusters[i] <- length(levels(as.factor(clusters[i,])))
}
plot(nclusters)

setwd("C:/Users/Tom Kennes/Desktop/Research Master/Thesis/Functions and Codes")
dat <- read.csv("clusters.csv")
max <- max(dat)
min = 0
plot(dat[,2], ylim = c(min, max), type = "l")
for(i in 3:ncol(dat)){
  rgb_ <- runif(n = 3, min= 0, max = 255)
  points(dat[,i], type = "l", col = 
           rgb(rgb_[1], rgb_[2], rgb_[3], maxColorValue = 255))
}


web <- cbind(iterations_m[,,n], clusters[n,])
web <- web[order(web[,ncol(web)]),]
n_levels(web[,ncol(web)])

gen_boolean <- function(x){
  tmp <- matrix(0, nrow = length(x), ncol = length(x))
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      if(x[i] == x[j]){
        tmp[i,j] = 1
      }
    }
  }
  return(tmp)
}

probs <- matrix(0, ncol = nrow(y), nrow = nrow(y))
start = 50
end = N
for(i in start:end){
  probs <- probs + gen_boolean(clusters[i,])/(end - start + 1)
}
probs
setwd("C:/Users/Tom Kennes/Desktop/Research Master/Thesis/Functions and Codes")
write.csv(probs, "probs.csv")





