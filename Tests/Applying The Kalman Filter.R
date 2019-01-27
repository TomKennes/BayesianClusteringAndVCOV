C  = diag(1)
T_ = diag(1)
R <- diag(1)
H <- diag(1)*10
V <- diag(1)
N = 100
test <- simulate_kalman(Z = C, T_, R, H, V, N)
plot(test$y)
points(test$a, col = "red", type = "l")
#y[t]   = Za[t] + e[t],  e[t] ~ N(0, H)
#a[t+1] = Ta[t] + Rw[t], w[t] ~ N(0, V)
#e.g. General multivariate Guassian Kalman Filter

setwd("C:/Users/tajmk/Desktop/Research Master/Thesis/Core/Functions and Codes")
y = read.csv("Monthly 49 Value Weighted.csv", sep = ";")
y = y[,1:ncol(y)]
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

#Set up Kalman Filter elements
Z  = diag(1)
T_ = diag(1)
R <- diag(1)
H <- diag(1)*10
V <- diag(1)
a0 = 0
P0 = 1
#y[t]   = Za[t] + e[t],  e[t] ~ N(0, H)
#a[t+1] = Ta[t] + Rw[t], w[t] ~ N(0, V)
p = c()
for(i in 2:ncol(y)){
  if(i == 1){
    p = Kalman_Filter(a0, P0, T_, H, Z, R, V, y[,i])$a_
  } else {
    p = cbind(p, Kalman_Filter(a0, P0, T_, H, Z, R, V, y[,i])$a_)
  }
}
colnames(p) <- colnames(y)[2:ncol(y)]
#write.csv(p, "Kalman Filtered Series.csv")

pdf("Kalman Filtered Variables.pdf")
for(i in 1:ncol(p)){
  plot(y[,i+1], main = colnames(y)[i+1], xlab = "Observations over time",ylab = "values")
  points(p[,i], col = "red", type = "l")
}
dev.off()




setwd("C:/Users/tajmk/Desktop/Research Master/Thesis/Core/Functions and Codes")
y = read.csv("Monthly 49 Value Weighted.csv", sep = ";")
y = y[,1:ncol(y)]
store = matrix("", ncol = 2, nrow = (ncol(y) - 1))
for(i in 2:ncol(y)){
  if(y[1,i] == -99.99){
    j = which.max(y[which(y[,i] == -99.99),1])
    store[i-1,] = c(y[j+1,1], colnames(y)[i])
  } else {
    store[i-1,] = c(y[1,1], colnames(y)[i])
  }
}
write.table(store, "portfolios and starting dates.txt")






Z  = diag(1)
T_ = diag(1)
R <- diag(1)
H <- diag(1)
V <- diag(1)*10
a0 = 0
P0 = 1
setwd("C:/Users/tajmk/Desktop/Research Master/Thesis/Core/Functions and Codes")
y = read.csv("Monthly 49 Value Weighted.csv", sep = ";")
y = y[,1:ncol(y)]
min = 0
for(i in 1:ncol(y)){
  x = length(which(y[,i] == -99.99))
  if(max(x) > min){
    min = max(x)
    print(min)
  }
}
y <- y[min:nrow(y),"Agric"]
plot(y)
p = Kalman_Filter(a0, P0, T_, H, Z, R, V, y)
pa = p$a_
pP = p$P
plot(y, col = "red", main = "Agriculture Portfolio. Blue: States, Red: Observations", 
     ylab = "values", xlab = "Observations over time",
     type = "p")
points(seq(1,length(y), 1), pa, col = "blue", type = "l")
plot(pP, type = "l", col = "blue",
     xlab = "Observations over Time",
     ylab = "Values P",
     main = "Agriculture Portfolio. Unobserved Variance (P)")
