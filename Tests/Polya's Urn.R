#Polya Urn Simulation
draws_ <- 100000
poss <- rep(0, draws_)
for(j in 1:draws_){
  draws <- 1000
  p_ <- rep(0, draws)
  balls <- c(80,10)
  p_[1] <- balls[1]/sum(balls)
  added_balls <- 1
  ratio <- rep(0, draws)
  ratio[1] <- balls[1]/balls[2]
  for(i in 2:draws){
    U <- runif(1,0,1)
    if(U <= p_[i - 1]){
      balls <- balls + c(added_balls,0)
    } else {
      balls <- balls + c(0, added_balls)
    }
    p_[i] <- balls[1]/sum(balls)
    ratio[i] <- balls[1]/balls[2]
  }
  poss[j] <- p_[draws]
  #par(mfrow = c(2,1))
  #plot(p_, type = "l", ylim = c(0,1))
  #plot(ratio, type = "l", col = "black")
  #abline(h = ratio[1], col = "red")
  print(j/draws_)
}
hist(poss, breaks = 100)
poss <- sort(poss, decreasing = T)
abline(v = poss[round(length(poss)*0.05)], col = "red")
abline(v = poss[round(length(poss)*0.95)], col = "red")








polyas_urn <- function(init = c(1,1), new = 1, cycles = 1000, init_p = c(0.5, 0.5)){
  if(sum(init_p) != 1){
    stop("Initial Probability measure does not sum to one")
  }
  if(length(init) != length(init_p)){
    stop("Initial dimensions of urns and probabilities are not matching")
  }
  container_p <- matrix(0, ncol = length(init_p), nrow = cycles)
  container_p[1,] <- init_p
  balls <- matrix(0, ncol = length(init), nrow = cycles)
  balls[1,] <- init
  for(i in 2:cycles){
    U <- runif(1,0,1)
    p_U <- cumsum(container_p[i-1,])
    sample <- min(which(p_U > U))
    draw <- rep(0, length(p_U))
    draw[sample] = 1
    balls[i,] <- balls[i - 1, ] + draw
    for(p in 1:length(init_p)){
      container_p[i, p] <- balls[i,p]/sum(balls[i,])
    }
    container_p[i, ] <- container_p[i,]
  }
  ret <- list(container_p, balls)
  names(ret) <- c("P", "draws")
  return(ret)
}
x <- polyas_urn(init = rep(1, 10), new = 1, cycles = 100000, init_p = rep(1/10, 10))


