
c <- matrix(0, ncol = 1, nrow = 10000)
for(i in 1:10000){
x <- rnorm(1000, 0, 5)
y <- 1 + x + rnorm(1000, 0, 1)
c[i] <- cov(x, y)
}

hist(c, breaks = 100)


install.packages("plot3D")
library("plot3D")
data("iris")
x <- sep.l <- iris$Sepal.Length
y <- pet.l <- iris$Petal.Length
z <- sep.w <- iris$Sepal.Width
# type ="h" for vertical lines
types <- c("b", "b2", "f", "g", "bl", "bl2", "u", "n")
for(i in 1:length(types)){
  scatter3D(x, y, z, phi = 25, bty = types[i],  type = "h", 
            ticktype = "detailed", pch = 20, cex = 0.5)
}
#b2 is the clearest
scatter3D(x, y, z,  bty = "f", colkey = FALSE,
          pch = 19, cex = 0.5)
scatter3D(x, y, z, bty = "f", colkey = FALSE, main ="bty= 'f'")

scatter3D(x, y, z, phi = 25, bty = "g",
          pch = 20, cex = 2, ticktype = "detailed", col.panel = "steelblue")
scatter3D(x, y, z, bty = "u", colkey = FALSE, 
          main ="bty= g", col.panel ="grey", expand =0.4, type = "j",
          col.grid = "darkblue", ticktype = "detailed", pch = 20, cex= 0.5, phi = 0)
scatter3D(x, y, z, phi = 0, bty = "g", pch = 20, cex = 0.5)
# Add text
text3D(x, y, z,  labels = rownames(iris),
       add = TRUE, colkey = FALSE, cex = 0.5)
data(VADeaths)
hist3D(z = VADeaths, scale = FALSE, expand = 0.01, bty = "b2", phi = 20,
      border = "black", shade = 0.2, theta = -45,
      space = 0.3, ticktype = "detailed", d = 2)


















