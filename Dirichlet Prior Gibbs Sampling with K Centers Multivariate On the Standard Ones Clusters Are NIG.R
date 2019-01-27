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


plot.lot <- function(y){
  m = max(y)
  n = min(y)
  colours <- matrix(runif(3*ncol(y), 0, 250), nrow = ncol(y), ncol = 3)
  plot(y[,1], 
       col = rgb(red = colours[1,1], green = colours[1,2], blue = colours[1,3], maxColorValue = 250), 
       type ="l",
       ylim = c(n,m))
  for(i in 2:ncol(y)){
    points(y[,i], col = rgb(red = colours[i,1], green = colours[i,2], blue = colours[i,3], maxColorValue = 250), type ="l")
  }
}

plot.timewindow <- function(y, pdf, period = 6){
  n <- floor(nrow(y)/period)
  pdf(pdf)
  for(i in 0:(n - 1)){
    x = y[(0 + period*i):(period*(i + 1)),]
    plot.lot(x)
  }
  plot.lot(y[n:nrow(y),])
  dev.off()
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

if(!("MASS" %in% installed.packages())){
  install.packages("MASS")
  library("MASS")
} else {
  library("MASS")
}
if(!("newtorkD3" %in% installed.packages())){
  install.packages("networkD3")
  library("networkD3")
} else {
  library("NetworkD3")
}
if(!("matrixcalc" %in% installed.packages())){
  install.packages("matrixcalc")
  library("matrixcalc")
} else {
  library("matrixcalc")
}
if(!("riverplot" %in% installed.packages())){
  install.packages("riverplot")
  library("riverplot")
} else {
  library("riverplot")
}
install.packages("fBasics")
require('fBasics')


sankey <- function(clusters, notepadding = 50, notewidth = 10, fontsize = 10, var_names = as.character(seq(1, ncol(clusters), by = 1))){
  # Set up from-to-node
  counter = 1
  from_to_node <- c(0,0,0,1,0,0)
  for(i in 1:ncol(clusters)){
    from_to_node <- rbind(from_to_node, c(0,clusters[1,i], counter, 1,i,1))
    counter = counter + 1
  }
  colnames(from_to_node) <- c("from", "to", "node", "weight", "t", "series")
  
  for(i in 2:nrow(clusters)){
    for(j in 1:ncol(clusters)){
      from_to_node <- rbind(from_to_node, c(clusters[i-1,j], clusters[i,j],counter,1,j,i))
      counter = counter + 1
    }
  }
  from_to_node <- from_to_node[2:nrow(from_to_node),]
  #from_to_node
  # Make corresponding adjustments
  # Combine nodes 
  n = 1
  while(n < nrow(from_to_node)){
    for(i in n:(n + ncol(clusters) - 2)){
      for(j in (i+1):(n + ncol(clusters) - 1)){
        #print(paste("i: ", i, " j: ", j))
        if(from_to_node[i,2] == from_to_node[j,2]){
          from_to_node[j,3] = from_to_node[i,3]
          if(from_to_node[i,1] == from_to_node[j,1]){
            from_to_node[i,4] <- from_to_node[i, 4] + 1
            from_to_node[j,4] <- from_to_node[j, 4] + 1
          }
        }
      }
    }
    n = n + ncol(clusters)
  }
  
  from_to_node <- from_to_node[,c(1,3:ncol(from_to_node))]
  #from_to_node
  
  for(i in (ncol(clusters)+1):nrow(from_to_node)){
    from_to_node[i,1] <- from_to_node[i - ncol(clusters),2]
  }
  #from_to_node
  # Scale down the to-nodes
  current = 1
  counter = 1
  for(i in seq(1,nrow(from_to_node),by = ncol(clusters))){
    subset <- from_to_node[i:(i+ncol(clusters) - 1),2]
    m <- max(subset)
    while(current <= m){
      if(sum(subset == current) == 0){
        current = current+ 1
      } else {
        subset[subset == current] <- counter
        counter = counter + 1
        current = current + 1
      }
    }
    from_to_node[i:(i+ncol(clusters) - 1),2] <- subset
  }
  #from_to_node
  
  
  for(i in 1:nrow(from_to_node)){
    from_to_node[i,3] = (from_to_node[,1] == from_to_node[i,1]) %*% (from_to_node[,2] == from_to_node[i,2])
  }
  #from_to_node
  
  for(i in (ncol(clusters) + 1):nrow(from_to_node)){
    from_to_node[i,1] <- from_to_node[i - ncol(clusters),2]
  }
  #from_to_node
  
  p = max(from_to_node[,2])
  
  nodes = data.frame("name" = 
                       paste("Node", 0:p, sep = " "))
  group <- as.character(from_to_node[,4])
  tmp <- from_to_node[,1:3]
  tmp <- vec(t(tmp))
  links = as.data.frame(matrix(tmp,# The third number is the value of the node
                               byrow = TRUE, ncol = 3))
  names(links) = c("source", "target", "value")
  occs <- rep("", nrow(links))
  
  for(i in 1:nrow(links)){
    occs[i] = var_names[from_to_node[i,4]]
  }
  links = cbind(links, occs)
  names(links) = c("source", "target", "value", "group")
  
  
  return(sankeyNetwork(Links = links, 
                       Nodes = nodes,
                       Source = "source", 
                       Target = "target",
                       Value = "value", 
                       NodeID = "name",
                       LinkGroup = "group",
                       fontSize= fontsize, 
                       nodeWidth = notewidth, 
                       nodePadding = notepadding,
                       #,NodeGroup = "group",
                       sinksRight = T
  ))
}



# Clustering a la escobar and west 1995
setwd("C:/Users/Tom Kennes/Desktop/Research Master/Thesis/Core/Functions and Codes")
y = read.csv("Monthly 49 Value Weighted.csv", sep = ";")
y = y[,2:ncol(y)]
min = 0
for(i in 1:ncol(y)){
  x = length(which(y[,i] == -99.99))
  if(max(x) > min){
    min = max(x)
    print(min)
  }
}

#Kick out the first min values.
y = y[(min + 1):nrow(y),]


t = 3
plot.timewindow(y[,2:ncol(y)], pdf = "Different Time Intervals.pdf", period = t)
y = y[20:(20+t),2:ncol(y)]
plot.lot(y)
y = t(y)
y <- as.matrix(y)


s = 4
S = 0.1
tau = 2
m <- colMeans(y)
M <- diag(colStdevs(y))*0.05
iterations = 200


#Gibbs Fashion Clustering Attempt
N <- iterations   # Number of Iterations in gibbs
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
alpha = 10E-6
fixed_variance = mean(diag(M))
use_fixed_variance = T
if(use_fixed_variance){
  iterations_V[, , , 1] <- diag(nrow(tmp$V))*fixed_variance
}
#colnames(clusters) <- rep(paste("y", as.character(seq(1, 700, by = 1)), sep = ""))

for(k in 2:N){
  if(length(levels(as.factor(as.character(clusters[k-1,])))) == 1){break}
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
      tmp <- rnormal_inverse_gamma(s, S, y[i,], y[i,], tau)
      iterations_m[i, , k] <- tmp$mu
      if(use_fixed_variance){
        iterations_V[, , i, k] <- diag(nrow(tmp$V))*fixed_variance
      } else {
        iterations_V[, , i, k] <- tmp$V
      }
    } else {
      # Meaning that we are adding it to another cluster
      y_clust <- y[clusters[k, ] == clusters[k,i],]
      tmp <- rnormal_inverse_gamma(s, S, y[i,], colMeans(y_clust), tau)
      iterations_m[i, ,k] <- tmp$mu
      if(use_fixed_variance){
        iterations_V[, , i, k] <- diag(nrow(tmp$V))*fixed_variance
      } else {
        iterations_V[, , i, k] <- tmp$V
      }
    }
  }
  nclusters[k] <- as.numeric(length(levels(as.factor(as.character(clusters[k,])))))
  print(clusters[k,])
}

write.csv(clusters, file = "clusters.csv")
dat <- read.csv("clusters.csv")
max <- max(clusters)
min = 0
plot(dat[,2], ylim = c(min, max), type = "l")
for(i in 3:ncol(dat)){
  rgb_ <- runif(n = 3, min= 0, max = 255)
  points(dat[,i], type = "l", col = 
           rgb(rgb_[1], rgb_[2], rgb_[3], maxColorValue = 255))
}


plot(nclusters, main = "Number of Clusters during the Simulation")
#p = sankey(clusters)
#p



#type <- levels(as.factor(clusters[k,]))
#plot(y[clusters[k,] == type[1],1], y[clusters[k,] == type[1],2],
#     ylim <- c(min(y)*1.1, max(y)*1.1),
#     xlim <- c(min(x)*1.1, max(x)*1.1),
#     type = "l", 
#     ylab = "y",
#     xlab = "observations",
#     col = "black")
#for(i in 2:length(type)){
#  rgb_ <- runif(n = 3, min= 0, max = 255)
#  points(y[clusters[k,] == type[i],1], y[clusters[k,] == type[i],2], type = "l", col = 
#           rgb(rgb_[1], rgb_[2], rgb_[3], maxColorValue = 255))
#}



#


#web <- cbind(iterations_m[,,n], clusters[n,])
#web <- web[order(web[,ncol(web)]),]
#n_levels(web[,ncol(web)])


#probs <- matrix(0, ncol = nrow(y), nrow = nrow(y))
#start = 50
#end = N
#for(i in start:end){
#  probs <- probs + gen_boolean(clusters[i,])/(end - start + 1)
#}
#probs
#setwd("C:/Users/Tom Kennes/Desktop/Research Master/Thesis/Functions and Codes")
#write.csv(probs, "probs.csv")





